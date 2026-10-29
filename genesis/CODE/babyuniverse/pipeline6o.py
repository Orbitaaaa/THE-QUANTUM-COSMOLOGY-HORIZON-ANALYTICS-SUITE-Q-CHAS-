"""Round 6o: one consistent description of a collider event in model U, from the parton collision through the interior to evaporation and the
fate of the daughter, and an event generator for the observables.

Per event (ROADMAP_ROUND6O.md):
1. Formation.  sqrt(s_hat) from the NNPDF2.3 luminosity x the geometric cross-section; impact parameter b^2 uniform below
   b_max = sqrt(F_n) r_h(sqrt(s_hat)); trapped mass M = y_n sqrt(s_hat) (head-on area bound, conservative); spin parameter recorded.
   A horizon forms iff M > M_ext.
2. Interior.  The colliding partons are an ingoing null shell along v = 0.  The centre is regular, so the shell passes through r = 0 and
   continues as an outgoing null shell in the core.
3. Evaporation.  Model-U Vaidya geometry f(v, r) = 1 - r^2 psi(mu(v)/r^{D-1}) with dmu/dv = -eps_n r_+^2 T^4 (brane Stefan-Boltzmann, in
   units of l).  ``vaidya_history`` computes r_+(v), r_-(v), the outgoing shell r(v) and its affine parameter.
   [Proposition-N 6o-1] T -> 0 at mu_ext, so mu(v) - mu_ext ~ 1/v: the two trapping horizons approach but never meet (the trapped region
   never closes; closure lemma of Round 6l); the outgoing shell stays inside the untrapped core, below r_-(v), and approaches the
   degenerate horizon with FINITE affine parameter (since int kappa_- dv diverges): it reaches the boundary of the parent's chart in
   finite affine time, so the daughter region beyond persists.  In GR-like complete evaporation the trapped region would close and the
   matter would return to the parent.
4. Emission.  Monte Carlo Hawking emission on the brane: Standard Model species with their degrees of freedom (bosons weight 1, fermions
   3/4 in number flux), energies from x^3/(e^x -+ 1) at T(M), masses respected, isotropic in the hole's frame, until M reaches the remnant
   (model U) or M_D (GR, then a two-body final decay).  Observables: multiplicity, visible and missing energy (remnant + neutrinos),
   species fractions, sphericity, temperature history.
"""

from __future__ import annotations

from math import gamma, pi

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import unified6l as U
from . import collider6l as C

# Standard Model species: (name, dof, spin-statistics, mass GeV, charge, visible)
SPECIES = [
    ("quark", 72, "F", 0.0, 1, True), ("top", 12, "F", 173.0, 1, True), ("gluon", 16, "B", 0.0, 0, True),
    ("W", 6, "B", 80.4, 1, True), ("Z", 3, "B", 91.2, 0, True), ("photon", 2, "B", 0.0, 0, True), ("higgs", 1, "B", 125.0, 0, True),
    ("charged_lepton", 12, "F", 0.0, 1, True), ("neutrino", 6, "F", 0.0, 0, False),
]


def a_n(n, c_l):
    """M = a_n M_D mu_hat with mu_hat = (R_S / l)^{n+1} (Dimopoulos-Landsberg)."""
    return (np.sqrt(pi) * c_l) ** (n + 1) * (n + 2) / (8 * gamma((n + 3) / 2))


def eps_rate(n, c_l, g_star=106.75):
    rc = ((n + 3) / 2) ** (1 / (n + 1)) * np.sqrt((n + 3) / (n + 1))
    return (pi**3 / 120) * g_star * rc**2 / (c_l * a_n(n, c_l))


def _T_rp(mu_hat, D, kind):
    if kind == "gr":
        rp = mu_hat ** (1 / (D - 3))
        return (D - 3) / (4 * pi * rp), rp, None
    tp = U.temperature(mu_hat, D, 1.0, kind)
    return tp["T"], tp["r_plus"], tp["r_minus"]


def vaidya_history(n=6, c_l=0.65, mu0_over_ext=3.0, kind="hayward", v_max=None, n_pts=3000) -> dict:
    """Evaporating model-U hole (units l = 1): mu(v), horizons, and the outgoing shell from r = 0 at v = 0 with its affine parameter."""
    D = n + 4
    me = U.extremal_mu(D, 1.0, kind)
    eps = eps_rate(n, c_l)
    mu0 = mu0_over_ext * me

    def dmu(v, y):
        mu = y[0]
        if mu <= me * (1 + 1e-12):
            return [0.0]
        T, rp, _ = _T_rp(mu, D, kind)
        return [-eps * rp**2 * T**4]

    T0, rp0, rm0 = _T_rp(mu0, D, kind)
    t_scale = (mu0 - me) / (eps * rp0**2 * T0**4)
    v_max = v_max or 400 * t_scale
    vs = np.concatenate([[0.0], np.geomspace(1e-3 * t_scale, v_max, n_pts)])
    sol = solve_ivp(dmu, (0, v_max), [mu0], t_eval=vs, rtol=1e-10, atol=1e-14 * mu0, method="LSODA")
    mus = sol.y[0]
    hz = [U.static_horizons(m, D, 1.0, kind) if m > me * (1 + 1e-9) else (None, None) for m in mus]
    rm = np.array([h[0] if h and h[0] else np.nan for h in hz])
    rp = np.array([h[1] if h and h[1] else np.nan for h in hz])
    # outgoing null shell from the regular centre: dr/dv = f/2; affine parameter: dlambda/dv = exp(int d_r f / 2 dv)
    mu_of_v = lambda v: float(np.interp(v, vs, mus))
    f = lambda v, r: float(U.static_f(max(r, 1e-9), mu_of_v(v), D, 1.0, kind))

    def rhs(v, y):
        r, lnw = y
        h = 1e-6 * max(r, 1e-3)
        dfdr = (f(v, r + h) - f(v, r - h)) / (2 * h)
        return [0.5 * f(v, r), 0.5 * dfdr]

    s2 = solve_ivp(rhs, (0, v_max), [1e-6, 0.0], t_eval=vs, rtol=1e-10, atol=1e-13, method="LSODA")
    r_shell = s2.y[0]
    w = np.exp(s2.y[1])
    lam = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(s2.t))])
    gap = rm[: len(r_shell)] - r_shell
    sep = rp - rm
    tail = slice(int(0.8 * len(vs)), None)
    # power-law fits of the late approach mu - mu_ext ~ v^-p and horizon separation ~ v^-q
    vv = vs[tail]
    p_mu = -np.polyfit(np.log(vv), np.log(np.maximum(mus[tail] - me, 1e-300)), 1)[0]
    q_sep = -np.polyfit(np.log(vv), np.log(np.maximum(sep[tail], 1e-300)), 1)[0]
    return dict(n=n, c_l=c_l, kind=kind, mu0_over_ext=mu0_over_ext, eps=eps, t_scale=float(t_scale), v=vs.tolist(),
                mu_over_ext=(mus / me).tolist(), r_plus=rp.tolist(), r_minus=rm.tolist(), r_shell=r_shell.tolist(), affine=lam.tolist(),
                horizons_never_meet=bool(np.all(sep[np.isfinite(sep)] > 0)), shell_below_inner_horizon=bool(np.nanmin(gap) > -1e-9),
                affine_final=float(lam[-1]), affine_increment_last_decade=float(lam[-1] - np.interp(v_max / 10, s2.t, lam)),
                mu_minus_ext_power=float(p_mu), separation_power=float(q_sep),
                final_mu_over_ext=float(mus[-1] / me), final_gap=float(gap[-1]))


# --------------------------------------------------------------------------- event generator


class Generator:
    def __init__(self, lumi: C.Lumi, sqrt_s_TeV=100.0, n=6, M_D_TeV=8.0, c_l=0.65, kind="hayward", trapped="conservative", seed=1):
        self.lumi, self.rs, self.n, self.MD, self.c_l, self.kind = lumi, sqrt_s_TeV, n, M_D_TeV, c_l, kind
        self.D = n + 4
        self.rng = np.random.default_rng(seed)
        self.y = C.trapped_fraction_headon(n) if trapped == "conservative" else 1.0
        self.MDg = M_D_TeV * 1e3
        self.l = c_l / self.MDg
        self.A = a_n(n, c_l)
        self.me = U.extremal_mu(self.D, 1.0, kind) if kind != "gr" else 0.0
        self.M_ext = self.A * self.MDg * self.me if kind != "gr" else 0.0
        self._build_sampler()
        # temperature table T(M) [GeV] in terms of mu_hat
        if kind != "gr":
            mus = self.me * (1 + np.geomspace(1e-6, 1e4, 600))
            self._Ttab_mu = mus
            self._Ttab_T = np.array([_T_rp(m, self.D, kind)[0] for m in mus]) / self.l

    def _build_sampler(self):
        s = (self.rs * 1e3) ** 2
        M_th = self.M_ext if self.kind != "gr" else 2.0 * self.MDg
        self.M_th = max(M_th, 1e-3)
        tau_min = (self.M_th / self.y) ** 2 / s
        self.tau_min = tau_min
        taus = np.geomspace(min(tau_min, self.lumi.tau_max * 0.999), self.lumi.tau_max, 3000)
        shat = np.sqrt(taus * s)
        R = C.rs_GeV(shat, self.n, self.MDg)
        dens = self.lumi(taus, self.y * shat) * C.YN_FORM_FACTOR[self.n] * pi * R**2 * taus        # per d ln tau
        cdf = np.concatenate([[0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(np.log(taus)))])
        self.sigma_pb = float(cdf[-1] * C.GEV2_TO_PB)
        self._cdf = cdf / cdf[-1] if cdf[-1] > 0 else cdf
        self._lt = np.log(taus)

    def T_of_M(self, M):
        if self.kind == "gr":
            R = C.rs_GeV(M, self.n, self.MDg)
            return (self.n + 1) / (4 * pi * R)
        mu = M / (self.A * self.MDg)
        if mu <= self.me:
            return 0.0
        return float(np.interp(mu, self._Ttab_mu, self._Ttab_T))

    def _emit(self, T):
        """One Hawking quantum at temperature T: species by number-flux weight, energy from the Bose/Fermi spectrum (masses respected)."""
        w = np.array([d * (1.0 if st == "B" else 0.75) * (np.exp(-m / T) if m > 0 else 1.0) for (_, d, st, m, _, _) in SPECIES])
        k = self.rng.choice(len(SPECIES), p=w / w.sum())
        name, d, st, m, q, vis = SPECIES[k]
        for _ in range(200):
            x = self.rng.gamma(4.0, 1.0)                    # proposal ~ x^3 e^-x
            acc = (1.0 / (1 - np.exp(-x))) if st == "B" else (1.0 / (1 + np.exp(-x)))
            if self.rng.random() < acc / (1.0 / (1 - np.exp(-1e-3)) if st == "B" else 1.0) or st == "F":
                if st == "F" and self.rng.random() > acc:
                    continue
                E = x * T
                if E > m:
                    return name, E, vis, q
        return name, max(3 * T, m * 1.01), vis, q

    def event(self):
        u = self.rng.random()
        tau = float(np.exp(np.interp(u, self._cdf, self._lt)))
        shat = np.sqrt(tau) * self.rs * 1e3
        M = self.y * shat
        rh = C.rs_GeV(shat, self.n, self.MDg)
        b = np.sqrt(self.rng.random()) * np.sqrt(C.YN_FORM_FACTOR[self.n]) * rh
        J = 0.5 * b * shat
        a_star = (self.n + 2) * J / (2 * M * C.rs_GeV(M, self.n, self.MDg))
        out = dict(sqrt_shat=shat, M=M, b=b, a_star=a_star, horizon=bool(M > self.M_ext) if self.kind != "gr" else True)
        if not out["horizon"]:
            out.update(N=0, E_vis=0.0, E_miss=0.0, remnant=0.0)
            return out
        M_end = self.M_ext * 1.0005 if self.kind != "gr" else self.MDg
        quanta = []
        Mc = M
        T0 = self.T_of_M(M)
        while True:
            T = self.T_of_M(Mc)
            if T <= 0:
                break
            name, E, vis, q = self._emit(T)
            if Mc - E < M_end:
                break
            Mc -= E
            ct = 2 * self.rng.random() - 1
            ph = 2 * pi * self.rng.random()
            st = np.sqrt(1 - ct**2)
            quanta.append((name, E, vis, E * st * np.cos(ph), E * st * np.sin(ph), E * ct))
        remnant = 0.0
        if self.kind == "gr":
            # Planck phase: the remaining energy in a two-body decay (visible)
            for sgn in (1, -1):
                ct = 2 * self.rng.random() - 1
                quanta.append(("final", Mc / 2, True, sgn * Mc / 2 * np.sqrt(1 - ct**2), 0.0, sgn * Mc / 2 * ct))
        else:
            remnant = Mc
        E_vis = sum(qq[1] for qq in quanta if qq[2])
        E_nu = sum(qq[1] for qq in quanta if not qq[2])
        P = np.array([[qq[3], qq[4], qq[5]] for qq in quanta if qq[2]]) if any(qq[2] for qq in quanta) else np.zeros((1, 3))
        S = (P.T @ P) / max(np.sum(P**2), 1e-300)
        sph = 1.5 * (np.sort(np.linalg.eigvalsh(S))[:2].sum())
        counts = {}
        for qq in quanta:
            counts[qq[0]] = counts.get(qq[0], 0) + 1
        out.update(N=len(quanta), N_vis=sum(1 for qq in quanta if qq[2]), E_vis=E_vis, E_nu=E_nu, remnant=remnant,
                   E_miss=remnant + E_nu, f_miss=(remnant + E_nu) / M, T0=T0, sphericity=float(sph), counts=counts,
                   mean_E=float(np.mean([qq[1] for qq in quanta])) if quanta else 0.0)
        return out

    def run(self, n_events=2000) -> dict:
        ev = [self.event() for _ in range(n_events)]
        bh = [e for e in ev if e["horizon"]]
        arr = lambda k: np.array([e[k] for e in bh])
        species = {}
        for e in bh:
            for k, v in e.get("counts", {}).items():
                species[k] = species.get(k, 0) + v
        tot = sum(species.values()) or 1
        return dict(kind=self.kind, n=self.n, M_D_TeV=self.MD, c_l=self.c_l, sqrt_s_TeV=self.rs, sigma_pb=self.sigma_pb,
                    M_ext_TeV=self.M_ext / 1e3, M_th_TeV=self.M_th / 1e3, n_events=n_events, frac_horizon=len(bh) / n_events,
                    M_TeV=(arr("M") / 1e3).tolist(), N=arr("N").tolist(), N_vis=arr("N_vis").tolist(), f_miss=arr("f_miss").tolist(),
                    E_vis_TeV=(arr("E_vis") / 1e3).tolist(), remnant_TeV=(arr("remnant") / 1e3).tolist(), T0_GeV=arr("T0").tolist(),
                    sphericity=arr("sphericity").tolist(), a_star=arr("a_star").tolist(), mean_E_GeV=arr("mean_E").tolist(),
                    species_fractions={k: v / tot for k, v in species.items()})
