"""Round 6n: model U beyond idealised collapse.  Spherically symmetric perfect-fluid collapse with pressure, inhomogeneity and shocks.

Equations (ROADMAP_ROUND6N.md, section 1).  Comoving coordinates, ds^2 = -e^{2 Phi} dt^2 + e^{2 lambda} dr^2 + R^2 dOmega_{D-2}^2, fluid at rest in
r, U = D_tau R, Gamma = R' e^{-lambda}.  In theories whose spherical sector has the Lovelock structure (field equations
h'(psi)(grad_a grad_b r - ...) + g_ab(...) = T_ab, with the generalised Misner-Sharp mass), which includes the quasi-topological towers of
model U:
  (1) generalised Misner-Sharp relation  Gamma^2 - U^2 = (grad R)^2 = 1 - R^2 psi(y),  y = mu / R^{D-1}
  (2) unified first law                    mu' = (D-1) rho R^{D-2} R',   D_tau mu = -(D-1) p R^{D-2} U        (units kappa~ = 1)
  (3) tr equation = h'(psi) x (its GR form) => lambda_dot R' = Rdot' - Phi' Rdot, hence D_tau Gamma = -U Gamma p'/((rho + p) R')
  (4) Euler                               Phi' = -p'/(rho + p)
which give the acceleration equation
      D_tau U = -Gamma^2 p'/((rho + p) R') + R psi(y) - ((D-1)/2) R psi'(y) (y + p).
GR limit (psi = y): D_tau U = -Gamma^2 p'/((rho+p)R') - M/R^2 - 4 pi p R (D = 4).  Homogeneous dust: addot = a [psi - ((D-1)/2) psi' y],
the generalised Friedmann equation of Round 6l.  Both limits are tests.

Scheme: Lagrangian (May-White type) finite volumes; rest mass per zone conserved; mu evolved by (2) and checked against the constraint;
Gamma from (1); ideal gas p = (gamma - 1) n eps; artificial viscosity q in compression; RK2 with a Courant step.

Results [Proposition-N 6n-1 .. 6n-4]:
* dust OS ball reproduces the analytic bounce radius and time of Round 6l to < 1 %;
* inhomogeneous dust: every shell bounces at R_b = sqrt(k alpha) (k = 1 - Gamma^2); shells cross after the bounce (multi-stream), handled by
  pressure;
* with pressure (p/rho up to ~0.2 at the start, shocks): horizon forms, every shell bounces inside an untrapped core, a shock forms, and the
  matter re-expands through anti-trapped (white-hole) spheres: the daughter forms; psi stays <= 1/alpha throughout;
* enough pressure prevents collapse altogether (no horizon, no daughter needed): the threshold is computed.
"""

from __future__ import annotations

import numpy as np

from . import unified6l as U


def psi_and_derivative(y, alpha=1.0, kind="hayward"):
    y = np.asarray(y, dtype=float)
    if kind == "hayward":
        return y / (1 + alpha * y), 1.0 / (1 + alpha * y) ** 2
    if kind == "dymnikova":
        e = np.exp(-alpha * y)
        return (1 - e) / alpha, e
    if kind == "gr":
        return y, np.ones_like(y)
    raise ValueError(kind)


class Star:
    """A spherically symmetric fluid ball in model U, D dimensions, units l = sqrt(alpha) = 1 and kappa~ = 16 pi G/((D-1)(D-2)) = 1."""

    def __init__(self, D=5, N=200, R_s=3.0, rho_c=1.0, profile="uniform", eps0=0.0, gamma_ad=None, alpha=1.0, kind="hayward",
                 c_q=2.0, c_lin=0.3):
        self.D, self.N, self.alpha, self.kind = D, N, alpha, kind
        self.gamma = gamma_ad if gamma_ad is not None else 1 + 1.0 / (D - 1)
        self.c_q, self.c_lin = c_q, c_lin
        R = np.linspace(0.0, R_s, N + 1)
        Rc = 0.5 * (R[1:] + R[:-1])
        x = Rc / R_s
        if profile == "uniform":
            rho = rho_c * np.ones(N)
        elif profile == "parabolic":
            rho = rho_c * (1 - 0.9 * x**2)
        elif profile == "gaussian":
            rho = rho_c * np.exp(-3 * x**2)
        elif profile == "smooth":
            rho = rho_c * (1 - x**2) ** 2 + 1e-6 * rho_c
        elif profile == "cored":
            rho = rho_c * (1 - x**4) ** 2 + 1e-6 * rho_c
        else:
            raise ValueError(profile)
        # specific internal energy tapered to zero at the surface, so p -> 0 there (a free surface)
        eps = eps0 * (1 - x**2) if np.isscalar(eps0) else np.asarray(eps0)
        n = rho / (1 + eps)
        V = (R[1:] ** (D - 1) - R[:-1] ** (D - 1)) / (D - 1)
        mu = np.concatenate([[0.0], np.cumsum(rho * (R[1:] ** (D - 1) - R[:-1] ** (D - 1)))])
        Ub = np.zeros(N + 1)
        G2 = 1 + Ub**2 - R**2 * psi_and_derivative(np.where(R > 0, mu / np.maximum(R, 1e-300) ** (D - 1), 0.0), alpha, kind)[0]
        if np.any(G2[1:] <= 0):
            raise ValueError("initial data trapped")
        Gam = np.sqrt(np.maximum(G2, 0))
        Gz = 0.5 * (Gam[1:] + Gam[:-1])
        self.db = n * V / Gz                                          # rest mass per zone (per unit solid angle), conserved
        self.R, self.Ub, self.mu, self.eps = R, Ub, mu, eps
        self.t = 0.0
        self.R_s0 = R_s
        self._mu0 = float(mu[-1])

    # ---------------------------------------------------------------- derived quantities
    def state(self, R, Ub, mu, eps):
        """Mass from the constraint mu' = (D-1) rho R^{D-2} R' (energy-consistent), Gamma from the Misner-Sharp relation, n from the conserved
        rest mass; the three are solved together by fixed-point iteration."""
        D, a, k = self.D, self.alpha, self.kind
        V = (R[1:] ** (D - 1) - R[:-1] ** (D - 1)) / (D - 1)
        dR1 = R[1:] ** (D - 1) - R[:-1] ** (D - 1)
        Gam = getattr(self, "_Gam", None)
        if Gam is None or len(Gam) != len(R):
            Gam = np.ones_like(R)
        for _ in range(6):
            Gz = 0.5 * (Gam[1:] + Gam[:-1])
            n = Gz * self.db / np.maximum(V, 1e-300)
            rho = n * (1 + eps)
            mu = np.concatenate([[0.0], np.cumsum(rho * dR1)])
            y = np.where(R > 0, mu / np.maximum(R, 1e-300) ** (D - 1), 0.0)
            psi, dpsi = psi_and_derivative(y, a, k)
            G2 = 1 + Ub**2 - R**2 * psi
            Gam = np.sqrt(np.maximum(G2, 1e-14))
        self._Gam = Gam
        p = (self.gamma - 1) * n * eps
        # artificial viscosity on non-homologous compression only (tensor-viscosity form): homologous collapse, U proportional to R, is
        # untouched, so the dust Oppenheimer-Snyder solution is exact in the scheme
        Rbar = 0.5 * (R[1:] + R[:-1])
        Ubar = 0.5 * (Ub[1:] + Ub[:-1])
        dU = (Ub[1:] - Ub[:-1]) - Ubar * (R[1:] - R[:-1]) / np.maximum(Rbar, 1e-300)
        cs = np.sqrt(np.maximum(self.gamma * p / np.maximum(rho + p, 1e-300), 0))
        q = np.where(dU < 0, rho * (self.c_q**2 * dU**2 + self.c_lin * cs * np.abs(dU)), 0.0)
        P = p + q
        H = rho + P
        # lapse: Phi = 0 at the surface; Phi' = -P'/H
        Phi_z = np.zeros(self.N)
        Phi_z[-1] = -P[-1] / H[-1] * 0.5
        for i in range(self.N - 2, -1, -1):
            Phi_z[i] = Phi_z[i + 1] + (P[i + 1] - P[i]) / (0.5 * (H[i] + H[i + 1]))
        Phi_b = np.concatenate([[Phi_z[0]], 0.5 * (Phi_z[1:] + Phi_z[:-1]), [0.0]])
        Rz = 0.5 * (R[1:] + R[:-1])
        dPdR = np.zeros(self.N + 1)
        Hb = np.ones(self.N + 1)
        Pb = np.zeros(self.N + 1)
        dPdR[1:-1] = (P[1:] - P[:-1]) / np.maximum(Rz[1:] - Rz[:-1], 1e-300)
        dPdR[-1] = (0 - P[-1]) / max(R[-1] - Rz[-1], 1e-300)
        Hb[1:-1] = 0.5 * (H[1:] + H[:-1]); Hb[-1] = H[-1]; Hb[0] = H[0]
        Pb[1:-1] = 0.5 * (P[1:] + P[:-1]); Pb[-1] = 0.0
        acc = -G2 * dPdR / Hb + R * psi - 0.5 * (D - 1) * R * dpsi * (y + Pb)
        acc[0] = 0.0
        return dict(y=y, psi=psi, G2=G2, Gam=Gam, V=V, n=n, p=p, rho=rho, q=q, P=P, H=H, Phi_b=Phi_b, Phi_z=Phi_z, acc=acc, Pb=Pb, mu=mu)

    def rhs(self, R, Ub, mu, eps):
        s = self.state(R, Ub, mu, eps)
        lap = np.exp(s["Phi_b"])
        D = self.D
        Rdot = lap * Ub
        Udot = lap * s["acc"]
        mudot = -(D - 1) * s["Pb"] * R ** (D - 2) * Rdot
        mudot[0] = 0.0
        # d(1/n)/dt for the energy equation: 1/n = V / (Gz db)
        Vdot = R[1:] ** (D - 2) * Rdot[1:] - R[:-1] ** (D - 2) * Rdot[:-1]
        psi, dpsi = psi_and_derivative(s["y"], self.alpha, self.kind)
        ydot = np.where(R > 0, mudot / np.maximum(R, 1e-300) ** (D - 1) - (D - 1) * s["y"] * Rdot / np.maximum(R, 1e-300), 0.0)
        G2dot = 2 * Ub * Udot - (2 * R * Rdot * psi + R**2 * dpsi * ydot)
        Gdot = G2dot / (2 * s["Gam"])
        Gz = 0.5 * (s["Gam"][1:] + s["Gam"][:-1])
        Gzdot = 0.5 * (Gdot[1:] + Gdot[:-1])
        inv_n_dot = (Vdot / Gz - s["V"] * Gzdot / Gz**2) / self.db
        epsdot = -s["P"] * inv_n_dot
        return Rdot, Udot, mudot, epsdot, s

    def dt_courant(self, s, cfl=0.25):
        dRp = (self.R[1:] - self.R[:-1]) / np.maximum(0.5 * (s["Gam"][1:] + s["Gam"][:-1]), 1e-6)
        cs = np.sqrt(np.maximum(self.gamma * s["p"] / np.maximum(s["rho"] + s["p"], 1e-300), 0)) + 1e-3
        lapz = np.exp(s["Phi_z"])
        dU = np.abs(self.Ub[1:] - self.Ub[:-1])
        dt1 = np.min(dRp / (lapz * (cs + dU + 1e-3)))
        dt2 = 0.02 * np.min(np.maximum(self.R[1:], 1e-3) / (np.abs(self.Ub[1:]) * np.exp(s["Phi_b"][1:]) + 1e-2))
        dt3 = 0.05 / np.max(np.abs(s["acc"]) * np.exp(s["Phi_b"]) + 1e-3)
        return cfl * min(dt1, dt2, dt3)

    def step(self):
        R, Ub, mu, eps = self.R, self.Ub, self.mu, self.eps
        k1 = self.rhs(R, Ub, mu, eps)
        dt = self.dt_courant(k1[4])
        R2, U2, m2, e2 = R + 0.5 * dt * k1[0], Ub + 0.5 * dt * k1[1], mu + 0.5 * dt * k1[2], eps + 0.5 * dt * k1[3]
        k2 = self.rhs(R2, U2, m2, np.maximum(e2, 0))
        self.R = R + dt * k2[0]
        self.Ub = Ub + dt * k2[1]
        self.eps = np.maximum(eps + dt * k2[3], 0.0)
        self.R[0] = 0.0; self.Ub[0] = 0.0
        self.t += dt
        s_new = self.state(self.R, self.Ub, self.mu, self.eps)
        self.mu = s_new["mu"]
        return dt, s_new

    def diagnostics(self, s):
        D = self.D
        gradR2 = 1 - self.R**2 * s["psi"]
        trapped = (gradR2 < 0) & (self.Ub < 0)
        anti = (gradR2 < 0) & (self.Ub > 0)
        if not hasattr(self, "_mu0"):
            self._mu0 = float(self.mu[-1])
        resid = abs(float(s["mu"][-1]) - self._mu0) / self._mu0                  # total-energy conservation (surface pressure ~ 0)
        return dict(t=self.t, R_surface=float(self.R[-1]), U_surface=float(self.Ub[-1]), rho_c=float(s["rho"][0]),
                    psi_max=float(np.max(s["psi"])), n_trapped=int(trapped.sum()), n_anti=int(anti.sum()),
                    n_bounced=int(np.sum(self.Ub[1:] > 0)), min_dR=float(np.min(np.diff(self.R))), constraint=float(resid),
                    p_over_rho_max=float(np.max(s["p"] / np.maximum(s["rho"], 1e-300))), mu_total=float(self.mu[-1]),
                    G2_min=float(np.min(s["G2"][1:])))


def run(star: Star, t_max=40.0, n_save=400, max_steps=400000, stop_when_reexpanded=True) -> dict:
    """Evolve and record the history; classify the outcome."""
    hist = []
    first_trapped = first_anti = None
    all_bounced_t = None
    save_every = max(1, max_steps // n_save)
    k = 0
    status = "ok"
    while star.t < t_max and k < max_steps:
        try:
            dt, s = star.step()
        except FloatingPointError:
            status = "floating point error"
            break
        k += 1
        if not np.all(np.isfinite(star.R)) or np.min(np.diff(star.R)) <= 0:
            status = "shell crossing or breakdown"
            break
        if k % 20 == 0 or k < 5:
            d = star.diagnostics(s)
            if d["n_trapped"] and first_trapped is None:
                first_trapped = d["t"]
            if d["n_anti"] and first_anti is None:
                first_anti = d["t"]
            if d["n_bounced"] == star.N and all_bounced_t is None:
                all_bounced_t = d["t"]
            if k % (20 * max(1, save_every // 20)) == 0 or k < 5:
                hist.append(d)
            if stop_when_reexpanded and all_bounced_t is not None and d["n_trapped"] == 0 and d["R_surface"] > star.R_s0 * 0.98:
                break
    s = star.state(star.R, star.Ub, star.mu, star.eps)
    final = star.diagnostics(s)
    hist.append(final)
    psi_max = max(h["psi_max"] for h in hist)
    outcome = ("daughter: horizon formed, all shells bounced, re-expansion through anti-trapped spheres"
               if first_trapped is not None and first_anti is not None and all_bounced_t is not None else
               "no horizon: pressure-supported or bounced without trapping" if first_trapped is None else "incomplete")
    return dict(D=star.D, N=star.N, kind=star.kind, gamma=star.gamma, steps=k, status=status, t_final=star.t, first_trapped=first_trapped,
                first_anti_trapped=first_anti, all_bounced=all_bounced_t, psi_max_alpha=psi_max * star.alpha,
                max_constraint=max(h["constraint"] for h in hist), outcome=outcome, history=hist,
                mu_over_ext=float(star.mu[-1] / U.extremal_mu(star.D, star.alpha, star.kind)))


def radius_for(mu, D, sin2chi=0.7, kind="hayward"):
    """Initial radius of a uniform ball of mass parameter mu at rest with R^2 psi(mu/R^{D-1}) = sin^2 chi (the outer root: untrapped)."""
    from scipy.optimize import brentq
    g = lambda R: R**2 * float(psi_and_derivative(mu / R ** (D - 1), 1.0, kind)[0]) - sin2chi
    _, r_peak = U._g_peak(mu, D, 1.0, kind)
    hi = r_peak * 2
    while g(hi) > 0:
        hi *= 2
    return brentq(g, r_peak, hi)


def os_check(D=5, N=120, R_s=None, rho_c=None, mu_over_ext=20.0, t_max=None, kind="hayward", sin2chi=0.7) -> dict:
    """Homogeneous dust: the comoving code against the analytic closed FRW (surface R = a(tau) sin chi_b, a_max = 1/sqrt(psi(rho0)))."""
    me = U.extremal_mu(D, 1.0, kind)
    if R_s is None:
        R_s = radius_for(mu_over_ext * me, D, sin2chi, kind)
    if rho_c is None:
        rho_c = mu_over_ext * me / R_s ** (D - 1)
    st = Star(D, N, R_s, rho_c, "uniform", 0.0, kind=kind)
    psi0 = float(psi_and_derivative(rho_c, 1.0, kind)[0])
    a_max = 1 / np.sqrt(psi0)
    sinchi = R_s / a_max
    C = rho_c * a_max ** (D - 1)
    an = U.closed_dust(C, D, 1.0, kind)
    Rmin_an = an["a_b"] * sinchi
    tmin_an = an["tau_half_cycle"]
    Rmin, tmin = np.inf, None
    while st.t < 2.2 * tmin_an:
        dt, s = st.step()
        if st.R[-1] < Rmin:
            Rmin, tmin = st.R[-1], st.t
    return dict(D=D, C=C, sin_chi=float(sinchi), R_min=float(Rmin), R_min_analytic=float(Rmin_an), t_bounce=float(tmin),
                t_bounce_analytic=float(tmin_an), rel_err_R=float(abs(Rmin / Rmin_an - 1)), rel_err_t=float(abs(tmin / tmin_an - 1)),
                mu_over_ext=float(st.mu[-1] / me))


def ltb_bounce_radii(D=5, profile="gaussian", R_s=2.0, rho_c=None, mu_over_ext=20.0, kind="hayward", N=100) -> dict:
    """Inhomogeneous dust: each shell obeys U^2 = Gamma^2 - 1 + R^2 psi(mu/R^{D-1}) with Gamma, mu fixed (LTB).  Bounce radius vs the
    small-radius law R_b = sqrt(k alpha), k = 1 - Gamma^2, and the proper time of each shell's bounce (shells crossing after the bounce)."""
    st = make_star(D, N, mu_over_ext, profile, 0.0, 0.7, kind)
    from scipy.optimize import brentq
    from scipy.integrate import quad
    rows = []
    for j in range(1, N + 1):
        R0, mu = st.R[j], st.mu[j]
        G2 = 1 - R0**2 * float(psi_and_derivative(mu / R0 ** (D - 1), 1.0, kind)[0])
        kk = 1 - G2
        F = lambda R: G2 - 1 + R**2 * float(psi_and_derivative(mu / R ** (D - 1), 1.0, kind)[0])
        try:
            Rb = brentq(F, 1e-6, R0 * (1 - 1e-9))
        except ValueError:
            continue
        w = lambda t: 2 * (R0 - Rb) * np.sin(t) * np.cos(t) / np.sqrt(max(F(Rb + (R0 - Rb) * np.sin(t) ** 2), 1e-300))
        tb = quad(w, 1e-9, np.pi / 2 - 1e-9, limit=200)[0]
        rows.append(dict(R0=float(R0), mu=float(mu), k=float(kk), R_b=float(Rb), R_b_law=float(np.sqrt(kk)), tau_bounce=float(tb)))
    tb = np.array([r["tau_bounce"] for r in rows])
    return dict(D=D, profile=profile, rows=rows, bounce_time_monotonic=bool(np.all(np.diff(tb) > -1e-9)),
                spread_of_bounce_times=float(tb.max() / tb.min()))


def make_star(D=5, N=100, mu_over_ext=20.0, profile="uniform", eps0=0.0, sin2chi=0.7, kind="hayward", **kw) -> Star:
    """A ball of total mass parameter mu_over_ext * mu_ext at rest, with the given density profile, starting untrapped."""
    mu = mu_over_ext * U.extremal_mu(D, 1.0, kind)
    R_s = radius_for(mu, D, sin2chi, kind)
    for _ in range(60):
        R = np.linspace(0.0, R_s, N + 1)
        x = 0.5 * (R[1:] + R[:-1]) / R_s
        shape = {"uniform": np.ones(N), "parabolic": 1 - 0.9 * x**2, "gaussian": np.exp(-3 * x**2),
                 "smooth": (1 - x**2) ** 2 + 1e-6, "cored": (1 - x**4) ** 2 + 1e-6}[profile]
        rho_c = mu / np.sum(shape * (R[1:] ** (D - 1) - R[:-1] ** (D - 1)))
        try:
            return Star(D, N, R_s, rho_c, profile, eps0, kind=kind, **kw)
        except ValueError:
            R_s *= 1.1
    raise ValueError("could not build untrapped initial data")


class ConservativeStar(Star):
    """Energy-conserving variant (May-White form): the boundary masses mu_j are evolved by the unified first law
    D_tau mu = -(D-1) P R^{D-2} U (P including viscosity), so the zone energy density rho_i = (mu_{i+1} - mu_i)/((D-1) V_i) satisfies the
    constraint exactly and the total mass is conserved to round-off (P = 0 at the surface).  Gamma follows algebraically from mu, the rest-mass
    density from Gamma and the conserved rest mass, and eps = rho/n - 1 (viscous heating appears automatically)."""

    def __init__(self, *a, dust=False, **kw):
        super().__init__(*a, **kw)
        self.dust = dust

    def state(self, R, Ub, mu, eps=None):
        D, a, k = self.D, self.alpha, self.kind
        V = (R[1:] ** (D - 1) - R[:-1] ** (D - 1)) / (D - 1)
        y = np.where(R > 0, mu / np.maximum(R, 1e-300) ** (D - 1), 0.0)
        psi, dpsi = psi_and_derivative(y, a, k)
        G2 = 1 + Ub**2 - R**2 * psi
        Gam = np.sqrt(np.maximum(G2, 1e-14))
        Gz = 0.5 * (Gam[1:] + Gam[:-1])
        n = Gz * self.db / np.maximum(V, 1e-300)
        rho = (mu[1:] - mu[:-1]) / np.maximum((D - 1) * V, 1e-300)
        epsz = np.maximum(rho / np.maximum(n, 1e-300) - 1.0, 0.0)
        p = np.zeros_like(rho) if self.dust else (self.gamma - 1) * n * epsz
        Rbar = 0.5 * (R[1:] + R[:-1])
        Ubar = 0.5 * (Ub[1:] + Ub[:-1])
        dU = (Ub[1:] - Ub[:-1]) - Ubar * (R[1:] - R[:-1]) / np.maximum(Rbar, 1e-300)
        cs = np.sqrt(np.maximum(self.gamma * p / np.maximum(rho + p, 1e-300), 0))
        q = np.where(dU < 0, rho * (self.c_q**2 * dU**2 + self.c_lin * cs * np.abs(dU)), 0.0)
        P = p + q
        H = rho + P
        Phi_z = np.zeros(self.N)
        Phi_z[-1] = -P[-1] / H[-1] * 0.5
        for i in range(self.N - 2, -1, -1):
            Phi_z[i] = Phi_z[i + 1] + (P[i + 1] - P[i]) / (0.5 * (H[i] + H[i + 1]))
        Phi_b = np.concatenate([[Phi_z[0]], 0.5 * (Phi_z[1:] + Phi_z[:-1]), [0.0]])
        Rz = 0.5 * (R[1:] + R[:-1])
        dPdR = np.zeros(self.N + 1); Hb = np.ones(self.N + 1); Pb = np.zeros(self.N + 1)
        dPdR[1:-1] = (P[1:] - P[:-1]) / np.maximum(Rz[1:] - Rz[:-1], 1e-300)
        dPdR[-1] = (0 - P[-1]) / max(R[-1] - Rz[-1], 1e-300)
        Hb[1:-1] = 0.5 * (H[1:] + H[:-1]); Hb[-1] = H[-1]; Hb[0] = H[0]
        Pb[1:-1] = 0.5 * (P[1:] + P[:-1])
        acc = -G2 * dPdR / Hb + R * psi - 0.5 * (D - 1) * R * dpsi * (y + Pb)
        acc[0] = 0.0
        return dict(y=y, psi=psi, G2=G2, Gam=Gam, V=V, n=n, p=p, rho=rho, q=q, P=P, H=H, Phi_b=Phi_b, Phi_z=Phi_z, acc=acc, Pb=Pb, mu=mu,
                    eps=epsz)

    def rhs(self, R, Ub, mu, eps=None):
        s = self.state(R, Ub, mu)
        lap = np.exp(s["Phi_b"])
        Rdot = lap * Ub
        Udot = lap * s["acc"]
        mudot = -(self.D - 1) * s["Pb"] * R ** (self.D - 2) * Rdot
        mudot[0] = 0.0; mudot[-1] = 0.0
        return Rdot, Udot, mudot, None, s

    def step(self):
        R, Ub, mu = self.R, self.Ub, self.mu
        k1 = self.rhs(R, Ub, mu)
        dt = self.dt_courant(k1[4])
        k2 = self.rhs(R + 0.5 * dt * k1[0], Ub + 0.5 * dt * k1[1], mu + 0.5 * dt * k1[2])
        self.R = R + dt * k2[0]
        self.Ub = Ub + dt * k2[1]
        self.mu = mu + dt * k2[2]
        self.R[0] = 0.0; self.Ub[0] = 0.0; self.mu[0] = 0.0
        self.t += dt
        s = self.state(self.R, self.Ub, self.mu)
        self.eps = s["eps"]
        return dt, s


def make_conservative(D=5, N=100, mu_over_ext=20.0, profile="smooth", eps0=0.0, sin2chi=0.7, kind="hayward", **kw) -> ConservativeStar:
    st = make_star(D, N, mu_over_ext, profile, eps0, sin2chi, kind)
    c = ConservativeStar.__new__(ConservativeStar)
    c.__dict__.update(st.__dict__)
    c.dust = (np.isscalar(eps0) and eps0 == 0.0)
    for k_, v in kw.items():
        setattr(c, k_, v)
    return c


def frw_bounce_any_fluid(w=1 / 3, C=1e4, D=5, alpha=1.0, kind="hayward") -> dict:
    """[Proposition 6n-1] Closed FRW with p = w rho, rho = C a^{-(D-1)(1+w)}: adot^2 = a^2 psi(rho) - 1.  Because psi <= 1/alpha, adot^2 < 0
    for a < sqrt(alpha): every closed universe with any w > -1 bounces at a_b >= sqrt(alpha) (and reaches a_max where GR recollapse sets in)."""
    from scipy.optimize import brentq
    F = lambda a: a**2 * float(psi_and_derivative(C * a ** (-(D - 1) * (1 + w)), alpha, kind)[0]) - 1.0
    grid = np.geomspace(1e-3, 1e8, 20000)
    vals = np.array([F(a) for a in grid])
    idx = np.where(np.diff(np.sign(vals)) != 0)[0]
    roots = [brentq(F, grid[i], grid[i + 1]) for i in idx]
    return dict(w=w, C=C, D=D, a_b=float(roots[0]) if roots else None, a_max=float(roots[1]) if len(roots) > 1 else None,
                bounces=bool(roots and roots[0] >= np.sqrt(alpha) * (1 - 1e-9)), n_roots=len(roots))


def merge_thin_zones(star: ConservativeStar, frac=0.02, min_zones=20) -> int:
    """Inelastic rezoning: a zone thinner than frac x (initial mean width) is merged with its thinner neighbour.  Rest mass and energy
    (the boundary masses mu_j) are conserved exactly; the boundary between the two zones is removed and its kinetic energy becomes heat
    (eps = rho/n - 1 of the merged zone).  This lets the evolution continue through the pile-up where the bounced core meets the envelope."""
    merged = 0
    thr = frac * star.R_s0 / max(len(star.db), 1) * (len(star.db) / star.N0 if hasattr(star, "N0") else 1)
    while len(star.db) > min_zones:
        dR = np.diff(star.R)
        i = int(np.argmin(dR))
        if dR[i] > thr:
            break
        # merge zone i with neighbour j (the thinner one); remove the boundary between them
        if i == 0:
            j = 1
        elif i == len(dR) - 1:
            j = i - 1
        else:
            j = i - 1 if dR[i - 1] < dR[i + 1] else i + 1
        lo, hi = min(i, j), max(i, j)
        b = hi                                                  # boundary index between zones lo and hi
        star.db = np.concatenate([star.db[:lo], [star.db[lo] + star.db[hi]], star.db[hi + 1:]])
        star.R = np.delete(star.R, b)
        star.Ub = np.delete(star.Ub, b)
        star.mu = np.delete(star.mu, b)
        star.N = len(star.db)
        star._Gam = None
        merged += 1
    return merged


def run_rezoned(star: ConservativeStar, t_max=60.0, max_steps=200000, frac=0.02, n_save=600, wall=1500.0) -> dict:
    """Evolution with inelastic rezoning; classifies whether the bounced matter reaches anti-trapped (white-hole) spheres and re-expands."""
    import time
    star.N0 = star.N
    t0 = time.time()
    hist, k, merges = [], 0, 0
    first_trapped = first_anti = None
    max_R_after_bounce = 0.0
    status = "ok"
    while star.t < t_max and k < max_steps and time.time() - t0 < wall:
        dt, s = star.step()
        k += 1
        if not np.all(np.isfinite(star.R)) or np.min(np.diff(star.R)) <= 0:
            status = "breakdown"
            break
        m = merge_thin_zones(star, frac)
        merges += m
        if m:
            s = star.state(star.R, star.Ub, star.mu)
        if k % 25 == 0:
            d = star.diagnostics(s)
            d["N_zones"] = star.N
            if d["n_trapped"] and first_trapped is None:
                first_trapped = d["t"]
            if d["n_anti"] and first_anti is None:
                first_anti = d["t"]
            hist.append(d)
    s = star.state(star.R, star.Ub, star.mu)
    d = star.diagnostics(s); d["N_zones"] = star.N
    hist.append(d)
    return dict(D=star.D, status=status, steps=k, merges=merges, t_final=star.t, first_trapped=first_trapped, first_anti_trapped=first_anti,
                psi_max_alpha=max(h["psi_max"] for h in hist) * star.alpha, energy_error=max(h["constraint"] for h in hist), history=hist,
                final=d)


def run_to_throat(D=5, N=50, mu_over_ext=20.0, profile="smooth", eps0=0.1, frac=0.02, max_steps=4000, kind="hayward") -> dict:
    """Rezoned evolution up to the first failure of the comoving (Gamma > 0) slicing, then a classification of the last good slice:
    'throat' if an anti-trapped shell (U > 0, (grad R)^2 < 0) sits next to a trapped one with Gamma -> 0 between them (the slice passes
    through a neck: the daughter's white-hole region attached to the parent's black-hole region), 'shell crossing' otherwise."""
    st = make_conservative(D, N, mu_over_ext, profile, eps0, kind=kind)
    st.N0 = st.N
    prev = None
    first_trapped = first_bounce = first_anti = None
    rho_c_max = 0.0
    for k in range(max_steps):
        dt, s = st.step()
        if (not np.all(np.isfinite(st.R))) or np.min(np.diff(st.R)) <= 0:
            break
        if merge_thin_zones(st, frac):
            s = st.state(st.R, st.Ub, st.mu)
        g = 1 - st.R**2 * s["psi"]
        if first_trapped is None and np.any((g < 0) & (st.Ub < 0)):
            first_trapped = st.t
        if first_bounce is None and np.any(st.Ub[1:] > 0):
            first_bounce = st.t
        if first_anti is None and np.any((g < 0) & (st.Ub > 0)):
            first_anti = st.t
        rho_c_max = max(rho_c_max, float(s["rho"][0]))
        prev = (st.R.copy(), st.Ub.copy(), s, st.t, g)
    R, Ub, s, t, g = prev
    anti = (g < 0) & (Ub > 0)
    trap = (g < 0) & (Ub < 0)
    adj = [j for j in range(len(R) - 1) if (anti[j] and trap[j + 1]) or (trap[j] and anti[j + 1])]
    G2_at = float(min(s["G2"][j] for j in adj)) if adj else None
    kind_end = "throat" if adj and G2_at < 0.05 else "shell crossing / shock"
    return dict(D=D, N=N, profile=profile, eps0=eps0, t_end=float(t), first_trapped=first_trapped, first_bounce=first_bounce,
                first_anti_trapped=first_anti, rho_c_max=rho_c_max, rho_c_end=float(s["rho"][0]), psi_max_alpha=float(np.max(s["psi"]) * st.alpha),
                energy_error=float(abs(s["mu"][-1] - st._mu0) / st._mu0), end=kind_end, G2_at_neck=G2_at,
                n_anti=int(anti.sum()), n_trap=int(trap.sum()),
                profile_end=dict(R=R.tolist(), U=Ub.tolist(), gradR2=g.tolist(), G2=s["G2"].tolist()))


def pressure_threshold(D=5, N=50, mu_over_ext=20.0, profile="smooth", eps_list=(0.1, 0.5, 1.0, 2.0, 4.0, 8.0), t_max=40.0, kind="hayward") -> list:
    """Does the star form a trapped surface, or does pressure stop the collapse (no horizon: no daughter needed)?"""
    rows = []
    for e in eps_list:
        st = make_conservative(D, N, mu_over_ext, profile, e, kind=kind)
        st.N0 = st.N
        trapped = None
        R0 = st.R[-1]
        Rmin = R0
        for k in range(20000):
            dt, s = st.step()
            if (not np.all(np.isfinite(st.R))) or np.min(np.diff(st.R)) <= 0:
                break
            merge_thin_zones(st, 0.05)
            g = 1 - st.R**2 * st.state(st.R, st.Ub, st.mu)["psi"]
            Rmin = min(Rmin, st.R[-1])
            if np.any((g < 0) & (st.Ub < 0)):
                trapped = st.t
                break
            if st.t > t_max or (st.Ub[-1] > 0 and st.R[-1] > R0):
                break
        p_over_rho = float(np.max(st.state(st.R, st.Ub, st.mu)["p"] / np.maximum(st.state(st.R, st.Ub, st.mu)["rho"], 1e-300)))
        rows.append(dict(eps0=e, horizon=trapped is not None, t_trapped=trapped, R_min_over_R0=float(Rmin / R0)))
    return rows


def trajectories(D=5, N=50, mu_over_ext=20.0, profile="smooth", eps0=0.1, frac=0.02, every=5, max_steps=4000, kind="hayward") -> dict:
    """Shell radii and trapping status along the evolution (for spacetime diagrams); shells are labelled by enclosed rest mass."""
    st = make_conservative(D, N, mu_over_ext, profile, eps0, kind=kind)
    st.N0 = st.N
    b_edges0 = np.concatenate([[0], np.cumsum(st.db)])
    labels = np.linspace(0, b_edges0[-1], 21)[1:]
    out = dict(t=[], R=[], status=[], psi_max=[], rho_c=[], energy=[])
    for k in range(max_steps):
        dt, s = st.step()
        if (not np.all(np.isfinite(st.R))) or np.min(np.diff(st.R)) <= 0:
            break
        if merge_thin_zones(st, frac):
            s = st.state(st.R, st.Ub, st.mu)
        if k % every:
            continue
        b_edges = np.concatenate([[0], np.cumsum(st.db)])
        Rl = np.interp(labels, b_edges, st.R)
        Ul = np.interp(labels, b_edges, st.Ub)
        psil = np.interp(labels, b_edges, s["psi"])
        g = 1 - Rl**2 * psil
        status = np.where(g >= 0, 0, np.where(Ul < 0, 1, 2))          # 0 untrapped, 1 trapped, 2 anti-trapped
        out["t"].append(st.t); out["R"].append(Rl.tolist()); out["status"].append(status.tolist())
        out["psi_max"].append(float(np.max(s["psi"]))); out["rho_c"].append(float(s["rho"][0]))
        out["energy"].append(float(abs(s["mu"][-1] - st._mu0) / st._mu0))
    return out
