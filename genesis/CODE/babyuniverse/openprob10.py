"""Round 10, open problems of the model-U programme (other than the double-null code, ``doublenull10``).

* ``causality_threshold``  P-R6l-b.  The tower h = psi + alpha psi^2 + ... contains Gauss-Bonnet with coefficient alpha = l^2 (the k = 2
  quasi-topological density in D >= 5 is Gauss-Bonnet).  Camanho, Edelstein, Maldacena & Zhiboedov (JHEP 02 (2016) 020, arXiv:1407.5597):
  a Gauss-Bonnet correction to the graviton three-point coupling violates causality unless an infinite tower of massive higher-spin states
  appears at m_HS <~ 1/l (up to O(1)).  Model U is therefore an effective theory whose limiting length also sets a new-physics threshold,
  M_HS <~ M_D / c_l: it does not fix c_l, but it ties it to a second observable.  [Assessment]
* ``scalar_from_bounce`` / ``scalar_bounce_scan``  P-R6l-e.  Closed FRW interior with dust and a scalar field in a false/true-vacuum
  potential through repeated model-U bounces (generalised Friedmann h(H^2 + 1/a^2) = rho), a fraction eps of the bounce energy given to
  the field (energy conserved).  A single bounce needs V_F > Lambda_c(C) (Proposition 6l-2), but the field's work over successive cycles
  (cosmological hysteresis, cf. Sahni & Toporensky) pumps the cycle maxima until the interior crosses into the de Sitter branch: lasting
  inflation for V_F well below Lambda_c once the field is trapped in the false vacuum.  Controls: eps = 0, or V = 0 in the field's vacuum,
  cycle for ever.  [Proposition-N, homogeneous interior]
* ``spin_down``  P-R6n-b.  Angular momentum of collider holes near threshold (from the generator) and the energy a rotating hole can
  release by superradiance even at T = 0, using the singly spinning Myers-Perry proxy E_rot/M = 1 - (1 + a*^2)^{-1/(n+2)} (no rotating
  model-U solution is known).  [Assessment]
* ``remnant_discharge``  P-R9-b.  Schwinger discharge of a remnant's electric and colour charge against hadronisation: integrating the
  pair-production rate (eE)^2/(24 pi) over the Coulomb field outside r_+ gives dQ/dt ~ -Q^2 alpha^2 / (6 r_+), so tau ~ 6 r_+/(Q alpha^2).
  Colour (alpha_s(1/r_+) ~ 0.09) discharges in ~0.1 fm/c, before hadronisation (~1 fm/c); SM triality-charge correlation then makes the
  colour-singlet remnant integer-charged, and e+e- emission neutralises it within ~1e-22 s.  Baseline: neutral remnants.  Alternative
  (colour discharge suppressed): R-hadron-like dressing, charged fraction ~0.3-0.5 by flavour counting.  [Assessment]
"""

from __future__ import annotations

from math import pi

import numpy as np
from scipy.integrate import solve_ivp

from . import collider6l as C
from . import pipeline6o as P
from . import unified6l as U

HBAR_GEV_S = 6.582119569e-25
HBARC_GEV_FM = 0.1973269804


# ----------------------------------------------------------------------------------------------------------- P-R6l-b

def causality_threshold(c_ls=(0.45, 0.6, 0.7, 0.8), M_D_TeV=8.0) -> dict:
    rows = [dict(c_l=c, l_inv_TeV=M_D_TeV / c, M_HS_bound_over_MD=1.0 / c) for c in c_ls]
    return dict(statement="Gauss-Bonnet coefficient alpha = l^2 in h(psi) = psi + alpha psi^2 + ...; CEMZ causality requires higher-spin "
                          "states at m_HS <~ 1/l = M_D/c_l (O(1) factor not fixed)",
                rows=rows, reference="Camanho, Edelstein, Maldacena, Zhiboedov, JHEP 02 (2016) 020, arXiv:1407.5597")


# ----------------------------------------------------------------------------------------------------------- P-R6l-e

def _V(phi, VF, A):
    s = np.clip(phi, 0.0, 1.0)
    step = s * s * (3.0 - 2.0 * s)
    return A * phi**2 * (phi - 1.0) ** 2 + VF * step


def _dV(phi, VF, A):
    s = np.clip(phi, 0.0, 1.0)
    dstep = 6.0 * s * (1.0 - s) if 0.0 < phi < 1.0 else 0.0
    return A * 2.0 * phi * (phi - 1.0) * (2.0 * phi - 1.0) + VF * dstep


def scalar_bounce(D=10, C_dust=None, VF_over_Lc=2.0, A=None, phi0=0.0, dphi0=0.0, a0=None, kind="hayward", t_max=400.0) -> dict:
    """Closed FRW: h(H^2 + 1/a^2) = C/a^{D-1} + phidot^2/2 + V(phi), Hdot = 1/a^2 - (D-1)(rho + p)/(2 h'(psi)), field eq. with friction.
    Starts contracting at a0; classifies the end state (lasting inflation in the false vacuum, recollapse, or cycling)."""
    C_dust = 3.0 * U.extremal_mu(D, 1.0, kind) if C_dust is None else C_dust
    Lc = U.vacuum_threshold(C_dust, D, 1.0, kind)["Lam_c"]
    VF = VF_over_Lc * Lc
    A = 16.0 * VF if A is None else A                    # barrier ~ VF
    a0 = 1.5 * C_dust ** (1.0 / (D - 3)) if a0 is None else a0

    def rho(a, phi, dphi):
        return C_dust / a ** (D - 1) + 0.5 * dphi**2 + _V(phi, VF, A)

    def rhs(t, y):
        a, H, phi, dphi = y
        r_ = rho(a, phi, dphi)
        psi = float(U.h_inv(r_, 1.0, kind))
        x = max(1.0 - psi, 1e-14)
        h1 = 1.0 / x**2 if kind == "hayward" else 1.0 / x
        rp = C_dust / a ** (D - 1) + dphi**2
        return [a * H, 1.0 / a**2 - (D - 1) * rp / (2.0 * h1), dphi, -(D - 1) * H * dphi - _dV(phi, VF, A)]

    psi0 = float(U.h_inv(rho(a0, phi0, dphi0), 1.0, kind))
    H0 = -np.sqrt(max(psi0 - 1.0 / a0**2, 0.0))
    ev_big = lambda t, y: y[0] - 1e3 * a0
    ev_big.terminal = True
    sol = solve_ivp(rhs, (0.0, t_max), [a0, H0, phi0, dphi0], method="DOP853", rtol=1e-9, atol=1e-12, events=ev_big, max_step=0.5)
    a, H, phi = sol.y[0], sol.y[1], sol.y[2]
    bounced = bool(np.any(H > 0))
    inflating = bool(sol.status == 1 or (a[-1] > 20 * a0 and H[-1] > 0))
    in_false = bool(abs(phi[-1] - 1.0) < 0.25)
    fate = "lasting inflation (false vacuum)" if inflating and in_false else ("expands" if inflating else "recollapses / cycles")
    return dict(D=D, C=C_dust, Lambda_c=Lc, VF=VF, VF_over_Lc=VF_over_Lc, bounced=bounced, a_min=float(a.min()), a_end=float(a[-1]),
                phi_end=float(phi[-1]), in_false_vacuum=in_false, fate=fate)


def scalar_from_bounce(D=10, C_dust=None, VF_over_Lc=2.0, eps=0.05, phi_b=0.0, sgn=1.0, kind="hayward", t_max=600.0, max_evals=3_000_000) -> dict:
    """Start at the dust bounce (H = 0, h(1/a^2) = rho) with a fraction eps of the bounce density given to the scalar as kinetic energy
    at field value phi_b (the bounce density, ~ 1/l^2, exceeds any sub-cutoff barrier), then evolve the re-expansion."""
    C_dust = 3.0 * U.extremal_mu(D, 1.0, kind) if C_dust is None else C_dust
    Lc = U.vacuum_threshold(C_dust, D, 1.0, kind)["Lam_c"]
    VF = VF_over_Lc * Lc
    A = 16.0 * VF
    a_b = U.closed_dust(C_dust, D, 1.0, kind)["a_b"]
    rho_b = C_dust / a_b ** (D - 1)
    K0 = eps * rho_b                                       # taken from the dust: energy is conserved at the bounce
    C_dust = (1.0 - eps) * C_dust - _V(phi_b, VF, A) * a_b ** (D - 1)
    a0 = a_b                                               # total density unchanged, so the bounce radius is unchanged
    dphi0 = sgn * np.sqrt(2.0 * K0)

    count = [0]

    class _Budget(Exception):
        pass

    def rhs(t, y):
        count[0] += 1
        if count[0] > max_evals:
            raise _Budget()
        a, H, phi, dphi = y
        r_ = C_dust / a ** (D - 1) + 0.5 * dphi**2 + _V(phi, VF, A)
        psi = float(U.h_inv(r_, 1.0, kind))
        x = max(1.0 - psi, 1e-14)
        h1 = 1.0 / x**2 if kind == "hayward" else 1.0 / x
        return [a * H, 1.0 / a**2 - (D - 1) * (C_dust / a ** (D - 1) + dphi**2) / (2.0 * h1), dphi, -(D - 1) * H * dphi - _dV(phi, VF, A)]

    ev_big = lambda t, y: y[0] - 1e3 * C_dust ** (1.0 / (D - 3))
    ev_big.terminal = True
    ev_crunch = lambda t, y: y[0] - 0.5 * a_b
    ev_crunch.terminal = True
    try:                                                                   # integration budget (stiff corners)
        sol = solve_ivp(rhs, (0.0, t_max), [a0, 1e-9, phi_b, dphi0], method="DOP853", rtol=1e-9, atol=1e-12,
                        events=[ev_big, ev_crunch], max_step=0.5)
    except _Budget:
        return dict(VF_over_Lc=VF_over_Lc, eps=eps, phi_b=phi_b, a_max=float("nan"), phi_end=float("nan"), in_false_vacuum=False,
                    fate="undetermined (integration budget)", n_cycles=0, cycle_maxima=[], t_end=float("nan"))
    a, H, phi = sol.y[0], sol.y[1], sol.y[2]
    grew = bool(len(sol.t_events[0]) > 0)
    in_false = bool(abs(phi[-1] - 1.0) < 0.3)
    fate = "lasting inflation (false vacuum)" if grew and in_false else ("expands" if grew else "recollapses")
    # cycle maxima (H: + -> -): cosmological hysteresis shows up as growing maxima
    turn = np.where((H[:-1] > 0) & (H[1:] <= 0))[0]
    maxima = [float(a[i]) for i in turn]
    return dict(VF_over_Lc=VF_over_Lc, eps=eps, phi_b=phi_b, a_max=float(a.max()), phi_end=float(phi[-1]), in_false_vacuum=in_false, fate=fate,
                n_cycles=len(maxima), cycle_maxima=maxima[:12], t_end=float(sol.t[-1]))


def _scan_job(args):
    D, C_d, q, eps, phi_b, sgn, kind, t_max = args
    r = scalar_from_bounce(D, C_d, q, eps, phi_b, sgn, kind, t_max=t_max)
    return q, eps, r["fate"].startswith("lasting"), r["in_false_vacuum"], r["n_cycles"], r["fate"].startswith("undetermined")


def scalar_bounce_scan(D=10, ratios=(0.5, 0.9, 1.5, 3.0), eps_list=(0.01, 0.1), n_phase=16, kind="hayward", t_max=3000.0,
                       workers=4) -> list:
    """Fraction of bounce-time field values (uniform in [0, 1], both velocity signs) that end in lasting inflation (parallel)."""
    import multiprocessing as mp
    C_d = 3.0 * U.extremal_mu(D, 1.0, kind)
    jobs = [(D, C_d, q, eps, float(pb), sgn, kind, t_max) for q in ratios for eps in eps_list
            for pb in np.linspace(0.0, 1.0, n_phase) for sgn in (1.0, -1.0)]
    with mp.Pool(workers) as pool:
        res = list(pool.imap_unordered(_scan_job, jobs, chunksize=1))
    rows = []
    for q in ratios:
        for eps in eps_list:
            rr = [r for r in res if r[0] == q and r[1] == eps]
            cyc = [r[4] for r in rr if r[2]]
            rows.append(dict(VF_over_Lc=q, eps=eps, fraction_lasting_inflation=float(np.mean([r[2] for r in rr])),
                             fraction_in_false_vacuum=float(np.mean([r[3] for r in rr])), n=len(rr), n_undetermined=int(sum(r[5] for r in rr)),
                             median_cycles_to_escape=float(np.median(cyc)) if cyc else None))
    return rows


def lambda_c_physical(n=6, c_l=0.6, M_D_TeV=8.0, mu_over_ext=3.0, kind="hayward") -> dict:
    """Lambda_c for a collider-sized hole in units of the curvature bound 1/l^2 and of M_D^2."""
    D = n + 4
    C_d = mu_over_ext * U.extremal_mu(D, 1.0, kind)
    Lc = U.vacuum_threshold(C_d, D, 1.0, kind)["Lam_c"]
    return dict(n=n, c_l=c_l, Lambda_c_l2=float(Lc), Lambda_c_over_MD2=float(Lc / c_l**2),
                note="the false vacuum must have rho comparable to the curvature-bound density (cutoff scale)")


# ----------------------------------------------------------------------------------------------------------- P-R6n-b

def rotational_fraction(a_star, n):
    a_star = np.asarray(a_star, dtype=float)
    return 1.0 - (1.0 + a_star**2) ** (-1.0 / (n + 2))


def spin_down(lumi: C.Lumi, sqrt_s=100.0, n=6, M_D_TeV=8.0, c_l=0.6, kind="hayward", n_events=4000, seed=3) -> dict:
    g = P.Generator(lumi, sqrt_s, n, M_D_TeV, c_l, kind, seed=seed)
    a = []; M = []
    for _ in range(n_events):
        e = g.event()
        if e["horizon"]:
            a.append(e["a_star"]); M.append(e["M"])
    a = np.array(a); M = np.array(M)
    f = rotational_fraction(a, n)
    Erot = f * M
    return dict(n=n, c_l=c_l, sqrt_s=sqrt_s, n_holes=len(a), mean_a_star=float(a.mean()), median_a_star=float(np.median(a)),
                mean_Erot_fraction=float(f.mean()), mean_Erot_GeV=float(Erot.mean()), frac_Erot_gt_1TeV=float(np.mean(Erot > 1e3)),
                a_star=a.tolist()[:2000], Erot_GeV=Erot.tolist()[:2000],
                note="Myers-Perry proxy: superradiant spin-down can release up to E_rot as visible TeV-scale quanta even at T = 0")


# ----------------------------------------------------------------------------------------------------------- P-R9-b

def remnant_discharge(n=6, M_D_TeV=8.0, c_l=0.6, kind="hayward", alpha_s=0.09, alpha_em=1 / 128.0, C_colour=4.0 / 3.0) -> dict:
    D = n + 4
    l = c_l / (M_D_TeV * 1e3)                                           # GeV^-1
    me = U.extremal_mu(D, 1.0, kind)
    rp = U.static_horizons(me * 1.001, D, 1.0, kind)[1] * l             # GeV^-1
    tau_em = 6.0 * rp / alpha_em**2 * HBAR_GEV_S                         # s, Q = 1
    tau_col = 6.0 * rp / (C_colour * alpha_s**2) * HBAR_GEV_S
    tau_had = 1.0 / 0.2 * HBAR_GEV_S                                    # ~ 1/Lambda_QCD ~ 1 fm/c
    E_discharge = alpha_s * C_colour / rp                               # GeV: potential at the horizon (energy scale of the emitted parton)
    # R-hadron-like dressing by flavour counting: P(u) = P(d) = 1/(2 + s), P(s) = s/(2 + s), s ~ 0.3
    s = 0.3
    pu = pd = 1.0 / (2.0 + s); ps = s / (2.0 + s)
    charged_meson = 2 * pu * pd + 2 * pu * ps                           # u dbar, d ubar, u sbar, s ubar (gluino-like, octet)
    triplet_charged = None                                              # depends on the remnant's own charge; quoted as a range
    return dict(r_plus_m=float(rp * HBARC_GEV_FM * 1e-15), tau_electric_s=float(tau_em), tau_colour_s=float(tau_col),
                tau_hadronisation_s=float(tau_had), colour_before_hadronisation=bool(tau_col < tau_had), ratio=float(tau_had / tau_col),
                E_discharge_GeV=float(E_discharge), dressed_charged_fraction_octet=float(0.9 * charged_meson),
                baseline="neutral remnant (colour discharges in ~%.2f fm/c, then e+e- emission removes integer charge)" % (tau_col / HBAR_GEV_S * HBARC_GEV_FM),
                alternative="if colour discharge were suppressed: R-hadron-like dressing with charged fraction ~0.3-0.5")
