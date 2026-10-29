"""Round 6l, priorities 2 and 4: black-hole production at hadron colliders with real parton distributions, event rates above the daughter
thresholds, interior formation against evaporation, and the observables that separate low-scale gravity, the interior classes, and what no
exterior measurement can see.

Parton distributions: LHAPDF6 ``lhagrid1`` files as shipped with Pythia 8.3 (default NNPDF2.3 QCD+QED LO, alpha_s = 0.130; cross-checks with
NNPDF3.1 LO and CT14qed).  ``read_lhagrid`` is a self-contained reader (bicubic in ln x, ln Q^2); ``lumi_table`` tabulates the summed parton
luminosity dL/dtau = int dx/x F(x) F(tau/x), F = sum over quarks, antiquarks and gluons (black-hole production is flavour-blind), and the
table is stored in DATA/round6l_lumi.json so the rest needs no PDF file.

Cross-section: sigma(pp -> BH) = int dtau dL/dtau(tau, Q) sigma_hat(sqrt(tau s)), sigma_hat = F_n pi R_S(sqrt(s_hat))^2 with the
Yoshino-Nambu form factors F_n (n = 1..7: 1.084 ... 1.883) and the Dimopoulos-Landsberg radius (G_D = 1/M_D^{n+2}).  Mass trapped:
optimistic M_BH = sqrt(s_hat); conservative M_BH = y_n sqrt(s_hat) with y_n the head-on (b = 0) Eardley-Giddings/Yoshino-Nambu area bound
(computed here: y_0 = 1/sqrt2).  Scale Q = M_BH or Q = 1/R_S.

Thresholds: model U horizon (daughter) threshold M_ext = kappa_n c_l^{n+1} M_D (unified6l.model_u_thresholds) and the FMM threshold
R_S > c_l / M_D of Round 6k.

Evaporation: brane-dominated Stefan-Boltzmann emission dM/dt = -(pi^2/120) g_* pi r_c^2 T^4, r_c = ((n+3)/2)^{1/(n+1)} sqrt((n+3)/(n+1)) r_+
(geometric-optics capture radius, Emparan-Horowitz-Myers), g_* = 106.75; GR (Schwarzschild-Tangherlini) or model U (T -> 0 at M_ext:
a remnant).

Results (ROADMAP_ROUND6L.md, sections 2 and 4):
* [Proposition-N 6l-6] The cross-section depends on the model essentially through the threshold M_ext.  With conservative trapping a
  100 TeV collider (30 ab^-1) gives >= 10 events up to M_ext ~ 47-50 TeV; the LHC (13 TeV, 138 fb^-1, none seen) requires M_ext >~ 9-10 TeV.
  The window M_ext ~ 9-50 TeV is c_l ~ 0.55-0.77 (Hayward tower) or 0.64-0.92 (Dymnikova tower), depending on n and M_D.
* Interior bounce ~1e-28 s after horizon formation; evaporation to the remnant takes ~1e-27-1e-24 s (ratio 20-7e3); the remnant keeps
  7-36 % of the mass; model U emits 1.5-2.5x more, softer quanta than GR; T_U is within 10 % of GR only once r_+ >~ 2.5-3.6 l.
"""


from __future__ import annotations

import json
import os
from math import gamma, pi

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator
from scipy.optimize import brentq

from . import unified6l as U

GEV2_TO_PB = 3.8938e8
HBAR_GEV_S = 6.582119569e-25
YN_FORM_FACTOR = {0: 0.647, 1: 1.084, 2: 1.341, 3: 1.515, 4: 1.642, 5: 1.741, 6: 1.819, 7: 1.883}
PDF_SEARCH = ["/usr/local/share/Pythia8/pdfdata/NNPDF23_lo_as_0130_qed_0000.dat"]
MACHINES = {"LHC 13.6": (13.6, 0.3), "HL-LHC 14": (14.0, 3.0), "HE-LHC 27": (27.0, 15.0), "FCC-hh 85": (85.0, 20.0),
            "SPPC/FCC 100": (100.0, 30.0)}                                     # sqrt(s) [TeV], integrated luminosity [ab^-1]


# --------------------------------------------------------------------------- parton distributions


def read_lhagrid(path):
    """Blocks of an LHAPDF6 lhagrid1 file: list of dicts (x, Q, flavours, xf[ix, iq, f])."""
    txt = open(path).read().split("\n---")
    blocks = []
    for chunk in txt[1:]:
        lines = [ln for ln in chunk.strip().split("\n") if ln.strip()]
        if len(lines) < 4:
            continue
        x = np.array(lines[0].split(), float)
        Q = np.array(lines[1].split(), float)
        fl = [int(v) for v in lines[2].split()]
        data = np.array(" ".join(lines[3:3 + len(x) * len(Q)]).split(), float).reshape(len(x), len(Q), len(fl))
        blocks.append(dict(x=x, Q=Q, flavours=fl, xf=data))
    return blocks


class PDF:
    """Summed parton number density F(x, Q) = sum_i f_i(x, Q) over quarks, antiquarks and gluons (photon excluded), plus single flavours."""

    def __init__(self, path):
        self.path = path
        self.blocks = read_lhagrid(path)
        self._spl = []
        for b in self.blocks:
            keep = [k for k, f in enumerate(b["flavours"]) if f != 22]
            tot = b["xf"][:, :, keep].sum(axis=2)
            self._spl.append(RectBivariateSpline(np.log(b["x"]), np.log(b["Q"] ** 2), tot, kx=3, ky=3))

    def _block(self, Q):
        for k, b in enumerate(self.blocks):
            if b["Q"][0] <= Q <= b["Q"][-1]:
                return k
        return len(self.blocks) - 1

    def xF(self, x, Q):
        k = self._block(Q)
        return np.maximum(self._spl[k].ev(np.log(x), np.log(Q**2) * np.ones_like(np.log(x))), 0.0)

    def xf_flavour(self, x, Q, pid):
        k = self._block(Q)
        b = self.blocks[k]
        j = b["flavours"].index(pid)
        s = RectBivariateSpline(np.log(b["x"]), np.log(b["Q"] ** 2), b["xf"][:, :, j], kx=3, ky=3)
        return s.ev(np.log(x), np.log(Q**2) * np.ones_like(np.log(x)))

    def sum_rules(self, Q=100.0):
        """Momentum sum int sum_i x f_i dx and valence numbers int (u - ubar), int (d - dbar)."""
        xs = np.exp(np.linspace(np.log(1e-7), np.log(0.999), 6000))
        mom = np.trapezoid(self.xF(xs, Q), xs)
        uv = np.trapezoid((self.xf_flavour(xs, Q, 2) - self.xf_flavour(xs, Q, -2)) / xs, xs)
        dv = np.trapezoid((self.xf_flavour(xs, Q, 1) - self.xf_flavour(xs, Q, -1)) / xs, xs)
        return dict(Q=Q, momentum=float(mom), u_valence=float(uv), d_valence=float(dv))

    def lumi(self, tau, Q, n=600):
        """dL/dtau = int_tau^1 dx/x F(x) F(tau/x), F = xF/x."""
        y = np.linspace(np.log(tau), 0.0, n)
        x1 = np.exp(y)
        x2 = tau / x1
        F1 = self.xF(x1, Q) / x1
        F2 = self.xF(np.minimum(x2, 1.0), Q) / x2
        return float(np.trapezoid(F1 * F2, y))


def lumi_table(pdf: PDF, taus=None, Qs=None) -> dict:
    taus = np.exp(np.linspace(np.log(1e-4), np.log(0.95), 160)) if taus is None else taus
    Qs = np.array([300., 600., 1e3, 2e3, 4e3, 8e3, 1.5e4, 3e4, 6e4, 1e5]) if Qs is None else Qs
    L = np.array([[pdf.lumi(t, q) for q in Qs] for t in taus])
    return dict(pdf=os.path.basename(pdf.path), tau=taus.tolist(), Q=Qs.tolist(), lnL=np.log(np.maximum(L, 1e-300)).tolist())


class Lumi:
    def __init__(self, table: dict):
        self.t = table
        self._i = RegularGridInterpolator((np.log(table["tau"]), np.log(table["Q"])), np.array(table["lnL"]), bounds_error=False,
                                          fill_value=None)
        self.tau_max = table["tau"][-1]

    def __call__(self, tau, Q):
        Q = np.clip(Q, self.t["Q"][0], self.t["Q"][-1])
        return np.exp(self._i(np.column_stack([np.log(np.atleast_1d(tau)), np.log(np.atleast_1d(Q))])))


def load_lumi(path) -> Lumi:
    return Lumi(json.load(open(path)))


# --------------------------------------------------------------------------- formation


def rs_GeV(M_GeV, n, M_D_GeV):
    """Dimopoulos-Landsberg horizon radius [GeV^-1] of a D = 4 + n black hole of mass M."""
    x = 8 * gamma((n + 3) / 2) / (n + 2) * (M_GeV / M_D_GeV)
    return x ** (1 / (n + 1)) / (np.sqrt(pi) * M_D_GeV)


def trapped_fraction_headon(n):
    """Head-on (b = 0) apparent-horizon mass bound: two Aichelburg-Sexl shocks of energy sqrt(s)/2 each in D = n + 4; trapped surface of
    radius r_c^{D-3} = 8 pi G mu / Omega_{D-3}, area 2 Omega_{D-3} r_c^{D-2}/(D-2); M_AH from the Tangherlini area-mass relation."""
    D = n + 4
    G, mu = 1.0, 0.5
    om = U.omega
    rc = (8 * pi * G * mu / om(D - 3)) ** (1 / (D - 3))
    A = 2 * om(D - 3) * rc ** (D - 2) / (D - 2)
    rh = (A / om(D - 2)) ** (1 / (D - 2))
    M = (D - 2) * om(D - 2) * rh ** (D - 3) / (16 * pi * G)
    return float(M)


def sigma_pp(lumi: Lumi, sqrt_s_TeV, n, M_D_TeV, M_th_TeV, trapped="conservative", scale="M", F=None, n_tau=400) -> float:
    """sigma(pp -> BH with M_BH > M_th) in pb."""
    s = (sqrt_s_TeV * 1e3) ** 2
    MD = M_D_TeV * 1e3
    y = trapped_fraction_headon(n) if trapped == "conservative" else 1.0
    tau_min = (M_th_TeV * 1e3 / y) ** 2 / s
    if tau_min >= lumi.tau_max:
        return 0.0
    F = YN_FORM_FACTOR[n] if F is None else F
    taus = np.exp(np.linspace(np.log(tau_min), np.log(lumi.tau_max), n_tau))
    shat = np.sqrt(taus * s)
    R = rs_GeV(shat, n, MD)
    Q = y * shat if scale == "M" else 1.0 / R
    sig_hat = F * pi * R**2
    integrand = lumi(taus, Q) * sig_hat * taus
    return float(np.trapezoid(integrand, np.log(taus)) * GEV2_TO_PB)


def events(sigma_pb, L_ab):
    return sigma_pb * L_ab * 1e6


def thresholds(n, M_D_TeV, c_l=1.0, kind="hayward") -> dict:
    t = U.model_u_thresholds(n, c_l, kind)
    return dict(model_U_TeV=t["M_ext_over_MD"] * M_D_TeV, fmm_TeV=t["M_fmm_over_MD"] * M_D_TeV, M_D_TeV=M_D_TeV, n=n, c_l=c_l, kind=kind)


def c_l_max_reachable(lumi: Lumi, n, M_D_TeV, sqrt_s_TeV, L_ab, N_min=10.0, kind="hayward", trapped="conservative") -> float:
    """Largest c_l (limiting length in units of 1/M_D) for which at least N_min model-U black holes above M_ext are produced."""
    kap = U.model_u_thresholds(n, 1.0, kind)["M_ext_over_MD"]

    def excess(c):
        M_th = kap * c ** (n + 1) * M_D_TeV
        return np.log(max(events(sigma_pp(lumi, sqrt_s_TeV, n, M_D_TeV, M_th, trapped), L_ab), 1e-30)) - np.log(N_min)

    if excess(1e-3) < 0:
        return 0.0
    if excess(1.5) > 0:
        return 1.5
    return float(brentq(excess, 1e-3, 1.5, xtol=1e-4))


def c_l_min_lhc(lumi: Lumi, n, M_D_TeV, sqrt_s_TeV=13.0, L_fb=138.0, N_max=3.0, kind="hayward") -> float:
    """Smallest c_l compatible with the LHC: fewer than N_max model-U black holes above M_ext at 13 TeV with 138 fb^-1 (optimistic trapping,
    no background: a simplified stand-in for the CMS black-hole search, which found none).  Below it M_ext falls inside the LHC's reach."""
    kap = U.model_u_thresholds(n, 1.0, kind)["M_ext_over_MD"]

    def excess(c):
        M_th = kap * c ** (n + 1) * M_D_TeV
        return np.log(max(sigma_pp(lumi, sqrt_s_TeV, n, M_D_TeV, M_th, "optimistic") * L_fb * 1e3, 1e-30)) - np.log(N_max)

    if excess(1.5) > 0:
        return 1.5
    if excess(1e-3) < 0:
        return 0.0
    return float(brentq(excess, 1e-3, 1.5, xtol=1e-4))


# --------------------------------------------------------------------------- evaporation and interior formation


def _T_rplus(mu_hat, n, kind):
    D = n + 4
    if kind == "gr":
        rp = mu_hat ** (1 / (n + 1))
        return (n + 1) / (4 * pi * rp), rp
    tp = U.temperature(mu_hat, D, 1.0, kind)
    return tp["T"], tp["r_plus"]


def evaporation(n=6, M_D_TeV=8.0, c_l=1.0, M0_over_MD=None, kind="hayward", g_star=106.75, eps_remnant=0.01) -> dict:
    """Mass loss from M0 to the remnant (model U: M_ext (1 + eps)) or to M_D (GR, end of the semiclassical regime).  Units: lengths in
    l = c_l/M_D.  Returns lifetime [s], emitted quanta (<E> ~ 3T), peak temperature, and the remnant mass."""
    D = n + 4
    MD = M_D_TeV * 1e3
    l = c_l / MD                                                                   # GeV^-1
    mu_of_M = lambda M: (rs_GeV(M, n, MD) / l) ** (n + 1)                          # mu in units of l^{n+1}
    M_of_mu = lambda mu: brentq(lambda M: mu_of_M(M) - mu, 1e-8 * MD, 1e12 * MD)
    M_ext = M_of_mu(U.extremal_mu(D, 1.0, kind)) if kind != "gr" else 0.0
    if M0_over_MD is None:
        M0_over_MD = 3.0 * max(M_ext / MD, 5.0)
    M0 = M0_over_MD * MD
    M_end = M_ext * (1 + eps_remnant) if kind != "gr" else MD
    if M0 <= M_end:
        return dict(n=n, kind=kind, formed=False)
    rc_fac = ((n + 3) / 2) ** (1 / (n + 1)) * np.sqrt((n + 3) / (n + 1))
    Ms = np.exp(np.linspace(np.log(M_end), np.log(M0), 1500))
    T_GeV, rate = [], []
    for M in Ms:
        T_hat, rp_hat = _T_rplus(mu_of_M(M), n, kind)
        T = T_hat / l
        rc = rc_fac * rp_hat * l
        T_GeV.append(T)
        rate.append((pi**2 / 120) * g_star * pi * rc**2 * T**4)                   # GeV^2  (dM/dt in GeV per GeV^-1)
    T_GeV, rate = np.array(T_GeV), np.array(rate)
    life = np.trapezoid(1.0 / np.maximum(rate, 1e-300), Ms)                      # GeV^-1
    N = np.trapezoid(1.0 / (3 * np.maximum(T_GeV, 1e-300)), Ms)
    return dict(n=n, kind=kind, M_D_TeV=M_D_TeV, c_l=c_l, M0_TeV=M0 / 1e3, M_remnant_TeV=M_ext / 1e3 if kind != "gr" else None,
                lifetime_s=float(life * HBAR_GEV_S), quanta=float(N), T_peak_GeV=float(T_GeV.max()), T_at_M0_GeV=float(T_GeV[-1]),
                remnant_energy_fraction=float(M_ext / M0) if kind != "gr" else 0.0)


def interior_time(n=6, mu_over_ext=10.0, kind="hayward", start=1.5) -> dict:
    """Proper time (units of l) for matter released from rest at start * r_+ to go from the horizon to its bounce in the core."""
    D = n + 4
    mu = mu_over_ext * U.extremal_mu(D, 1.0, kind)
    rm, rp = U.static_horizons(mu, D, 1.0, kind)
    E2 = float(U.static_f(start * rp, mu, D, 1.0, kind))
    fF = lambda r: float(U.static_f(r, mu, D, 1.0, kind))
    R_turn = brentq(lambda r: fF(r) - E2, 1e-6, rm * (1 - 1e-12))
    w = lambda t: 2 * (rp - R_turn) * np.sin(t) * np.cos(t) / np.sqrt(max(E2 - fF(R_turn + (rp - R_turn) * np.sin(t) ** 2), 1e-300))
    tau = quad(w, 1e-9, pi / 2 - 1e-9, limit=400)[0]
    return dict(n=n, mu_over_ext=mu_over_ext, r_plus=rp, r_minus=rm, R_turn=R_turn, tau_horizon_to_bounce_l=float(tau))


def formation_vs_evaporation(n=6, M_D_TeV=8.0, c_l=0.7, M_over_ext=2.0, kind="hayward") -> dict:
    MD = M_D_TeV * 1e3
    l = c_l / MD
    th = U.model_u_thresholds(n, c_l, kind)
    M0 = M_over_ext * th["M_ext_over_MD"] * MD
    mu_hat = (rs_GeV(M0, n, MD) / l) ** (n + 1)
    it = interior_time(n, mu_hat / U.extremal_mu(n + 4, 1.0, kind), kind)
    ev = evaporation(n, M_D_TeV, c_l, M0 / MD, kind)
    t_int = it["tau_horizon_to_bounce_l"] * l * HBAR_GEV_S
    return dict(n=n, M_D_TeV=M_D_TeV, c_l=c_l, M_TeV=M0 / 1e3, interior_s=float(t_int), lifetime_to_remnant_s=ev["lifetime_s"],
                ratio=float(ev["lifetime_s"] / t_int), quanta=ev["quanta"], remnant_TeV=ev["M_remnant_TeV"])


# --------------------------------------------------------------------------- what experiments can and cannot distinguish


def temperature_curves(n=6, kind="hayward", n_pts=60) -> dict:
    """T(M) for model U and GR at equal mass (units: l = 1, mu in l^{n+1}): the remnant signature and where they differ by 10%."""
    D = n + 4
    me = U.extremal_mu(D, 1.0, kind)
    rows = []
    for x in np.exp(np.linspace(np.log(1.0005), np.log(1e6), n_pts)):
        mu = x * me
        TU, _ = _T_rplus(mu, n, kind)
        TG, _ = _T_rplus(mu, n, "gr")
        rows.append(dict(mu_over_ext=float(x), T_U=float(TU), T_GR=float(TG), ratio=float(TU / TG)))
    r = np.array([q["ratio"] for q in rows])
    x10 = rows[int(np.argmax(r > 0.9))]["mu_over_ext"] if np.any(r > 0.9) else np.inf
    # the same mass ratio expressed as a radius ratio r_+/l (model U approaches GR only once r_+ >> l)
    rp10 = U.static_horizons(x10 * me, D, 1.0, kind)[1] if np.isfinite(x10) else np.inf
    return dict(n=n, kind=kind, rows=rows, mu_over_ext_where_within_10pct=float(x10), r_plus_over_l_where_within_10pct=float(rp10))


def discriminants() -> list:
    """Which observable tests what (Round 6l, section 4)."""
    return [
        dict(observable="missing-energy monojet / monophoton (KK gravitons)", tests="P1: n and M_D", interior_sensitive=False),
        dict(observable="virtual-graviton dilepton / diphoton tail", tests="P1: cut-off scale M_S", interior_sensitive=False),
        dict(observable="onset of high-multiplicity, high-sum-pT, democratic events", tests="P1 (black holes exist) and, through M_th/M_D, the"
             " interior class (GR-like: few M_D; model U: kappa_n c_l^{n+1} M_D)", interior_sensitive=True),
        dict(observable="T(M) from spectra vs reconstructed mass", tests="T ~ M^{-1/(n+1)} (GR) vs T -> 0 at M_ext (regular core)",
             interior_sensitive=True),
        dict(observable="stable remnants: missing energy, heavy stable charged tracks", tests="regular core / mimetic (remnant) vs complete"
             " evaporation", interior_sensitive=True),
        dict(observable="final burst at M ~ M_D", tests="not semiclassical: unreliable discriminator", interior_sensitive=False),
        dict(observable="anything about the region beyond the core (new universe, inflating de Sitter, crunch)", tests="nothing: the exterior"
             " is identical for interiors that agree up to the inner horizon (exterior-equivalence lemma)", interior_sensitive=False),
    ]
