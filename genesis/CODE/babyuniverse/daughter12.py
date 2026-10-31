"""Round 12: a finite daughter that forms, expands and recollapses or ends, inside the region fixed by the field equations.

A. Cauchy-development daughter (``cdd``).  A homogeneous dust ball (closed FRW patch, comoving radius chi_s) collapses in model U
   (any D >= 4, Hayward or Dymnikova tower, alpha = 1, l = 1).  Its interior obeys adot^2 = a^2 psi(a) - 1 with psi = h^{-1}(C/a^{D-1}),
   C = mu / sin^{D-1} chi_s, and bounces at a_b >= l.  The exterior is the static model-U hole of mass mu.  After the bounce the ball's
   surface leaves the inner core through r = r_- at the event E.  E lies on the null hypersurface v = infinity (the Cauchy horizon of the
   exterior), and inside the dust the Cauchy horizon continues as the ingoing null cone from E, chi_cone(eta) = chi_s - (eta - eta_E).
   The daughter region is
       D = { post-bounce events of the dust interior with chi < chi_cone(eta) }.
   D lies in the maximal globally hyperbolic development: no continuation through a Cauchy horizon is needed.  It lies inside the event
   horizon, and it contains a de Sitter-like (accelerating) phase.  The code measures its expansion, volume, duration and acceleration.

B. Lemma (no recollapse in the core).  For these towers f decreases on (0, r_-), so the ball's surface (adot^2 = E^2 - f) has no turning
   point between the bounce and r_-: a single ball cannot recollapse before leaving the core (``core_monotone``).

C. Loophole 1, compact trapped region (``compact_trapped_region``).  Ingoing Vaidya with mu(v) that forms the hole and then evaporates
   below mu_ext in finite advanced time.  The code computes the trapped region, the inner-horizon blueshift exponent int kappa_- dv
   (Carballo-Rubio et al. 2026: the s-wave RSET grows like exp(int kappa_- dv) and is finite for a compact trapped region), outgoing null
   rays (is there an event horizon?) and a bounced test body.

D. Loophole 2, inner extremality as an attractor (``extremality_attractor``).  Model-U evaporation drives kappa_- -> 0.  The code finds
   the rate, and the blueshift exponent int kappa_- dv, for power-law and exponentially decaying perturbations, and for finite-time
   extremality.
"""

from __future__ import annotations

from math import gamma, pi

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar

from . import hydro6n as H
from . import unified6l as U


# --------------------------------------------------------------------------------------------------------------- A. FRW interior

def psi_y(a, C, D, kind):
    y = C / a ** (D - 1)
    psi, dpsi = H.psi_and_derivative(y, 1.0, kind)
    return float(psi), float(dpsi), y


def turning_points(C, D, kind="hayward"):
    """Bounce a_b and maximum a_max of the closed FRW interior: roots of a^2 psi(a) = 1."""
    g = lambda la: np.exp(2 * la) * psi_y(np.exp(la), C, D, kind)[0] - 1.0
    res = minimize_scalar(lambda la: -g(la), bounds=(np.log(0.5), np.log(10 * C ** (1.0 / max(D - 3, 1)) + 10)), method="bounded")
    lp = res.x
    if g(lp) <= 0:
        return None
    a_b = np.exp(brentq(g, np.log(1e-3), lp, xtol=1e-14))
    hi = lp + 1.0
    while g(hi) > 0:
        hi += 1.0
    a_max = np.exp(brentq(g, lp, hi, xtol=1e-14))
    return float(a_b), float(a_max)


def _I(chi, D):
    """int_0^chi sin^{D-2} x dx (comoving volume of a ball of radius chi on the unit S^{D-1}, without the Omega_{D-2} factor)."""
    return quad(lambda x: np.sin(x) ** (D - 2), 0, chi)[0]


def evolve_ball(C, D, kind="hayward", rtol=1e-11):
    """Interior from rest at a_max through the bounce and back to a_max.  Returns dense arrays in proper time tau."""
    tp = turning_points(C, D, kind)
    if tp is None:
        raise ValueError("no oscillation: C below the critical value")
    a_b, a_max = tp

    def rhs(t, y):
        a, ad, eta = y
        psi, dpsi, yy = psi_y(a, C, D, kind)
        return [ad, a * (psi - 0.5 * (D - 1) * dpsi * yy), 1.0 / a]

    def ev_bounce(t, y):
        return y[1]
    ev_bounce.direction = 1
    ev_bounce.terminal = False

    # half period tau(a_b -> a_max) with a = a_b + (a_max - a_b) sin^2(th) to remove the endpoint singularities
    def dtau(th):
        a = a_b + (a_max - a_b) * np.sin(th) ** 2
        g = a * a * psi_y(a, C, D, kind)[0] - 1.0
        return 2 * (a_max - a_b) * np.sin(th) * np.cos(th) / np.sqrt(max(g, 1e-300))
    T_half = quad(dtau, 0, np.pi / 2, limit=400)[0]
    a0 = a_max * (1 - 1e-13)
    def ev_max(t, y):
        return y[1]
    ev_max.direction = -1
    ev_max.terminal = False
    # collapse, bounce, expansion to a_max, recollapse, second bounce (3 half periods)
    sol = solve_ivp(rhs, (0, 3 * T_half * (1 + 1e-3)), [a0, -1e-12, 0.0], method="DOP853", rtol=rtol, atol=1e-13, events=(ev_bounce, ev_max),
                    dense_output=True, max_step=max(T_half / 400, 1e-3))
    tb = sol.t_events[0][0]
    tm = [t for t in sol.t_events[1] if t > tb]
    tb2 = [t for t in sol.t_events[0] if t > tb + 0.5 * T_half]
    t_max = tm[0] if tm else float("nan")
    t_b2 = tb2[0] if tb2 else float("nan")
    return dict(sol=sol, a_b=a_b, a_max=a_max, tau_b=float(tb), eta_b=float(sol.sol(tb)[2]), tau_end=float(sol.t[-1]), T_half=float(T_half),
                tau_max=float(t_max), eta_max=float(sol.sol(t_max)[2]) if np.isfinite(t_max) else float("nan"),
                tau_b2=float(t_b2), eta_b2=float(sol.sol(t_b2)[2]) if np.isfinite(t_b2) else float("nan"))


def cdd(D=5, mu_over_ext=20.0, chi_s=0.5, kind="hayward", n_grid=4000) -> dict:
    """The Cauchy-development daughter of a homogeneous dust ball (see module docstring)."""
    me = U.extremal_mu(D, 1.0, kind)
    mu = mu_over_ext * me
    C = mu / np.sin(chi_s) ** (D - 1)
    rm, rp = U.static_horizons(mu, D, 1.0, kind)
    ev = evolve_ball(C, D, kind)
    sol, a_b, a_max, tb = ev["sol"], ev["a_b"], ev["a_max"], ev["tau_b"]
    R0 = a_max * np.sin(chi_s)
    untrapped_initially = bool(R0 > rp)
    # dense post-bounce grid
    tau = np.linspace(tb, ev["tau_end"], n_grid)
    tau = np.unique(np.concatenate([tau, tb + (ev["tau_end"] - tb) * np.geomspace(1e-9, 1e-2, 400)]))
    Y = sol.sol(tau)
    a, ad, eta = Y
    R = a * np.sin(chi_s)
    iE = np.where(R >= rm)[0]
    if len(iE) == 0:
        return dict(error="surface never leaves the core")
    k = iE[0]
    # refine E by root finding on R(tau) - r_-
    tE = brentq(lambda t: sol.sol(t)[0] * np.sin(chi_s) - rm, tau[max(k - 1, 0)], tau[k], xtol=1e-14)
    aE, _, etaE = sol.sol(tE)
    eta_cone = etaE + chi_s
    # acceleration phase after the bounce
    acc = np.array([a_ * (psi_y(a_, C, D, kind)[0] - 0.5 * (D - 1) * psi_y(a_, C, D, kind)[1] * psi_y(a_, C, D, kind)[2]) for a_ in a])
    ia = np.where(acc <= 0)[0]
    a_acc = float(a[ia[0]]) if len(ia) else float("nan")
    # center reached by the cone?
    if eta_cone <= eta[-1]:
        tc = brentq(lambda t: sol.sol(t)[2] - eta_cone, tau[0], tau[-1], xtol=1e-13)
        a_cone = float(sol.sol(tc)[0])
        tau_D_center = float(tc - tb)
    else:
        tc, a_cone, tau_D_center = float("nan"), float("nan"), float("nan")
    # volume of D
    chiD = np.where(eta < etaE, chi_s, np.clip(chi_s - (eta - etaE), 0, None))
    Ivals = np.array([_I(c, D) for c in chiD[::max(1, len(chiD) // 600)]])
    aa = a[::max(1, len(chiD) // 600)]
    V = aa ** (D - 1) * Ivals
    V_b = a_b ** (D - 1) * _I(chi_s, D)
    j = int(np.argmax(V))
    out = dict(D=D, kind=kind, mu_over_ext=mu_over_ext, chi_s=chi_s, C=float(C), r_minus=float(rm), r_plus=float(rp), a_b=a_b, a_max=a_max,
               R0=float(R0), untrapped_initially=untrapped_initially, eta_b=ev["eta_b"], eta_E=float(etaE), a_E=float(aE),
               eta_cone=float(eta_cone), dS_prediction_eta_cone_minus_eta_b=float(np.pi / 2),
               eta_cone_minus_eta_b=float(eta_cone - ev["eta_b"]), a_cone=a_cone, a_cone_over_a_b=a_cone / a_b,
               a_cone_over_C_root=a_cone / C ** (1.0 / (D - 1)), a_acc=a_acc, efolds_accelerating=float(np.log(a_acc / a_b)),
               efolds_to_cone=float(np.log(a_cone / a_b)), fraction_of_max=a_cone / a_max, tau_daughter_center=tau_D_center,
               V_max_over_V_bounce=float(V[j] / V_b), log10_V_max_over_V_bounce=float(np.log10(V[j] / V_b)),
               eta_max_minus_eta_b=float(ev["eta_max"] - ev["eta_b"]), eta_b2_minus_eta_b=float(ev["eta_b2"] - ev["eta_b"]),
               recollapse_begins_inside_D=bool(eta_cone > ev["eta_max"]), full_cycle_inside_D=bool(eta_cone > ev["eta_b2"]),
               a_at_cone_phase="recollapsing" if eta_cone > ev["eta_max"] else "expanding",
               hidden=True)
    out["series"] = dict(eta=(eta - ev["eta_b"])[::8].tolist(), a=a[::8].tolist(), chiD=chiD[::8].tolist())
    return out


def core_monotone(D=5, mu_over_ext=20.0, kind="hayward", n=4000) -> dict:
    """Lemma B: f is strictly decreasing on (0, r_-]."""
    me = U.extremal_mu(D, 1.0, kind)
    mu = mu_over_ext * me
    rm, rp = U.static_horizons(mu, D, 1.0, kind)
    r = np.linspace(1e-4, rm, n)
    f = U.static_f(r, mu, D, 1.0, kind)
    return dict(D=D, kind=kind, mu_over_ext=mu_over_ext, r_minus=float(rm), f_decreasing=bool(np.all(np.diff(f) < 0)),
                f0=float(f[0]), f_at_rm=float(f[-1]))


def physical_scales(D, l_m, mu_over_ext=None, r_g_m=None, chi_s=0.5, kind="hayward") -> dict:
    """Run ``cdd`` and convert to physical units (l in metres).  Give either mu/mu_ext or the gravitational radius r_g (mu = r_g^{D-3})."""
    me = U.extremal_mu(D, 1.0, kind)
    if r_g_m is not None:
        mu_over_ext = (r_g_m / l_m) ** (D - 3) / me
    r = cdd(D, mu_over_ext, chi_s, kind)
    c = 2.99792458e8
    r.update(l_m=l_m, a_b_m=r["a_b"] * l_m, a_cone_m=r["a_cone"] * l_m, a_max_m=r["a_max"] * l_m,
             tau_daughter_center_s=r["tau_daughter_center"] * l_m / c)
    return r


# --------------------------------------------------------------------------------------------------------------- C. compact trapped region

def kappa_minus(mu, D, kind="hayward"):
    hz = U.static_horizons(mu, D, 1.0, kind)
    if hz is None:
        return 0.0, None, None
    rm, rp = hz
    h = 1e-7 * rm
    fp = (U.static_f(rm + h, mu, D, 1.0, kind) - U.static_f(rm - h, mu, D, 1.0, kind)) / (2 * h)
    return float(abs(fp) / 2), float(rm), float(rp)


def compact_trapped_region(D=5, mu0_over_ext=5.0, v_form=5.0, T_life=40.0, v_evap=20.0, mu_end_over_ext=0.5, kind="hayward",
                           n_rays=24, v_max=None) -> dict:
    """Ingoing Vaidya ds^2 = -f(v, r) dv^2 + 2 dv dr + r^2 dOmega with f = 1 - r^2 psi(mu(v)/r^{D-1}).
    mu(v): smooth rise to mu0 over v_form, plateau, then a smooth fall to mu_end < mu_ext over v_evap, starting at v_form + T_life.
    Returns horizons, trapped-region lifetime, int kappa_- dv, whether outgoing rays all escape (no event horizon), and a bounced test
    body (radial timelike geodesic from near the centre with E = 0.9)."""
    me = U.extremal_mu(D, 1.0, kind)
    mu0, mu1 = mu0_over_ext * me, mu_end_over_ext * me
    s = lambda x: 0.5 * (1 + np.tanh(x))
    v1 = v_form + T_life

    def mu_of(v):
        return mu0 * s((v - 0.5 * v_form) / (0.15 * v_form)) - (mu0 - mu1) * s((v - v1 - 0.5 * v_evap) / (0.15 * v_evap))

    def dmu(v):
        return (mu0 / (0.15 * v_form) * 0.5 / np.cosh((v - 0.5 * v_form) / (0.15 * v_form)) ** 2
                - (mu0 - mu1) / (0.15 * v_evap) * 0.5 / np.cosh((v - v1 - 0.5 * v_evap) / (0.15 * v_evap)) ** 2)

    v_max = v_max or v1 + v_evap + 60.0
    f = lambda v, r: float(U.static_f(max(r, 1e-9), max(mu_of(v), 1e-12), D, 1.0, kind))

    def fr(v, r, h=1e-6):
        return (f(v, r + h) - f(v, r - h)) / (2 * h)

    vv = np.linspace(0, v_max, 1500)
    rows = [kappa_minus(mu_of(v), D, kind) for v in vv]
    km = np.array([r[0] for r in rows]); rmv = np.array([r[1] if r[1] else np.nan for r in rows]); rpv = np.array([r[2] if r[2] else np.nan for r in rows])
    trapped = np.isfinite(rmv)
    v_open, v_close = (float(vv[trapped][0]), float(vv[trapped][-1])) if trapped.any() else (None, None)
    blueshift_exponent = float(np.trapezoid(km, vv))
    # outgoing null rays dr/dv = f/2 from points inside the core and the trapped region during the plateau
    rays = []
    for v0 in np.linspace(v_form + 2, v1, 4):
        hz = U.static_horizons(mu_of(v0), D, 1.0, kind)
        for r0 in np.linspace(0.3 * hz[0], 0.95 * hz[1], n_rays // 4):
            sol = solve_ivp(lambda v, y: [0.5 * f(v, y[0])], (v0, v_max), [r0], max_step=0.5, rtol=1e-9, atol=1e-11)
            rays.append(dict(v0=float(v0), r0=float(r0), r_end=float(sol.y[0, -1]), escaped=bool(sol.y[0, -1] > 3 * np.nanmax(rpv)),
                             v=sol.t[::10].tolist(), r=sol.y[0, ::10].tolist()))
    all_escape = all(r["escaped"] for r in rays)

    # bounced test body: radial timelike geodesic, Lagrangian L = -f vdot^2 + 2 vdot rdot.  Geodesic equations:
    #   vddot = -(1/2) f_r vdot^2 ,  rddot = -(1/2) f_v vdot^2 ... derived from E-L:  d/dtau(-f vdot + rdot) = -(1/2) f_v vdot^2,
    #   d/dtau(vdot) = -(1/2) f_r vdot^2 (from the r equation).
    def geo(tau, y):
        v, r, vd, rd = y
        fv = (f(v + 1e-6, r) - f(v - 1e-6, r)) / 2e-6
        frr = fr(v, r)
        vdd = -0.5 * frr * vd**2
        # from d/dtau(-f vdot + rdot) = -(1/2) f_v vdot^2:  -f_v vd^2 - f_r rd vd - f vdd + rdd = -(1/2) f_v vd^2
        rdd = 0.5 * fv * vd**2 + frr * rd * vd + f(v, r) * vdd
        return [vd, rd, vdd, rdd]

    v0 = v_form + 3.0
    r0 = 0.3
    fv0 = f(v0, r0)
    rd0 = 0.1                                                  # outward, bound (static energy < 1), from the bounce region
    vd0 = (rd0 + np.sqrt(rd0**2 + fv0)) / fv0                 # -f vd^2 + 2 vd rd = -1, future-directed root
    stop = lambda t, y: y[0] - v_max
    stop.terminal = True
    sol = solve_ivp(geo, (0, 2000), [v0, r0, vd0, rd0], max_step=0.05, rtol=1e-9, atol=1e-11, events=[stop])
    body = dict(v=sol.y[0, ::20].tolist(), r=sol.y[1, ::20].tolist(), r_max=float(np.max(sol.y[1])),
                v_at_rmax=float(sol.y[0, int(np.argmax(sol.y[1]))]),
                recollapsed=bool(sol.y[1, -1] < 0.9 * np.max(sol.y[1])))
    return dict(D=D, kind=kind, mu0_over_ext=mu0_over_ext, T_life=T_life, v_open=v_open, v_close=v_close,
                lifetime=(v_close - v_open) if v_open is not None else None, blueshift_exponent=blueshift_exponent,
                rset_growth_factor=float(np.exp(min(blueshift_exponent, 700))), no_event_horizon=all_escape, rays=rays, body=body,
                v=vv.tolist(), r_minus=rmv.tolist(), r_plus=rpv.tolist(), kappa_minus=km.tolist())


# --------------------------------------------------------------------------------------------------------------- D. extremality attractor

def extremality_attractor(n=6, c_l=0.70, mu0_over_ext=3.0, kind="hayward", alpha_decay=None, p_tail=3.0) -> dict:
    """Model-U evaporation (Round 6o Vaidya, brane Stefan-Boltzmann; greybody correction ~ 0.9 near extremality is a constant factor).
    kappa_-(v) -> 0 is an attractor; fit kappa_- ~ v^-q and Lambda(v) = int kappa_- dv.  Ori-type inner mass function
    m(v) ~ tail(v) exp(Lambda(v)) for power-law tails v^-p or exponential decay exp(-alpha v) (a positive cosmological constant)."""
    from . import pipeline6o as P
    vh = P.vaidya_history(n, c_l, mu0_over_ext, kind, n_pts=1500)
    D = n + 4
    me = U.extremal_mu(D, 1.0, kind)
    v = np.array(vh["v"]); mu = np.array(vh["mu_over_ext"]) * me
    km = np.array([kappa_minus(m, D, kind)[0] for m in mu])
    tail = slice(int(0.7 * len(v)), None)
    ok = (v[tail] > 0) & (km[tail] > 0)
    q = -np.polyfit(np.log(v[tail][ok]), np.log(km[tail][ok]), 1)[0]
    Lam = np.concatenate([[0.0], np.cumsum(0.5 * (km[1:] + km[:-1]) * np.diff(v))])
    c_fit = float(np.polyfit(v[tail] ** (1 - q), Lam[tail], 1)[0])
    lnm_power = -p_tail * np.log(np.maximum(v, 1.0)) + Lam
    out = dict(n=n, c_l=c_l, kappa_minus_power_q=float(q), Lambda_coefficient=c_fit, Lambda_final=float(Lam[-1]), v_final=float(v[-1]),
               lnm_power_tail_final=float(lnm_power[-1]), lnm_power_tail_increasing=bool(lnm_power[-1] > lnm_power[len(v) // 2]),
               t_scale=vh["t_scale"], v=v[::10].tolist(), kappa_minus=km[::10].tolist(), Lambda=Lam[::10].tolist())
    if alpha_decay is not None:
        # m ~ exp(Lambda(v) - alpha v); with Lambda ~ c v^{1-q}, the maximum is at v* = ((1-q) c / alpha)^{1/q}
        v_star = ((1 - q) * c_fit / alpha_decay) ** (1.0 / q)
        out.update(alpha_decay=alpha_decay, v_star=float(v_star), bounded_with_exponential_decay=True,
                   lnm_peak=float(c_fit * v_star ** (1 - q) - alpha_decay * v_star))
    return out


def finite_time_extremality(D=10, kind="hayward", q_law=0.5, mu0_over_ext=3.0, rate=0.05, v_max=400.0) -> dict:
    """If the luminosity near extremality scaled as (mu - mu_ext)^{q_law} with q_law < 1, extremality would be reached at finite v_e and
    kappa_- = 0 afterwards, so int kappa_- dv is finite: the blueshift is bounded (inner-extremal end state).  Model U's T^4 law has
    effective q_law = 2 (kappa ~ sqrt(mu - mu_ext)), so this needs a non-thermal emission channel (assumption)."""
    me = U.extremal_mu(D, 1.0, kind)

    def rhs(v, y):
        x = max(y[0] - me, 0.0)
        return [-rate * me * (x / me) ** q_law]

    sol = solve_ivp(rhs, (0, v_max), [mu0_over_ext * me], max_step=0.2, rtol=1e-10, atol=1e-14)
    mu = sol.y[0]
    km = np.array([kappa_minus(m, D, kind)[0] if m > me * (1 + 1e-10) else 0.0 for m in mu])
    reach = np.where(mu <= me * (1 + 1e-9))[0]
    Lam = float(np.trapezoid(km, sol.t))
    return dict(D=D, q_law=q_law, v_extremal=float(sol.t[reach[0]]) if len(reach) else None, Lambda_total=Lam,
                blueshift_bounded=bool(len(reach) > 0))


# --------------------------------------------------------------------------------------------------------------- A'. quadrature version

def _eta_integral(a1, a2, a_b, a_max, C, D, kind):
    """int_{a1}^{a2} da / (a sqrt(a^2 psi - 1)) on [a_b, a_max], with a = a_b + (a_max - a_b) sin^2(th) near the turning points."""
    def integrand(th):
        a = a_b + (a_max - a_b) * np.sin(th) ** 2
        g = a * a * psi_y(a, C, D, kind)[0] - 1.0
        return 2 * (a_max - a_b) * np.sin(th) * np.cos(th) / (a * np.sqrt(max(g, 1e-300)))
    th = lambda a: np.arcsin(np.sqrt(np.clip((a - a_b) / (a_max - a_b), 0, 1)))
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return quad(integrand, th(a1), th(a2), limit=800, epsabs=1e-12, epsrel=1e-10)[0]


def cdd_fast(D=10, mu_over_ext=1.05, chi_s=1.4, kind="hayward") -> dict:
    """Independent quadrature check of the timing: cone budget = (eta_E - eta_b) + chi_s, half period, full period."""
    me = U.extremal_mu(D, 1.0, kind)
    mu = mu_over_ext * me
    C = mu / np.sin(chi_s) ** (D - 1)
    hz = U.static_horizons(mu, D, 1.0, kind)
    tp = turning_points(C, D, kind)
    if hz is None or tp is None:
        return dict(valid=False)
    rm, rp = hz
    a_b, a_max = tp
    aE = rm / np.sin(chi_s)
    budget = _eta_integral(a_b, min(aE, a_max), a_b, a_max, C, D, kind) + chi_s
    half = _eta_integral(a_b, a_max, a_b, a_max, C, D, kind)
    if not (np.isfinite(budget) and np.isfinite(half) and 0 < budget < 3 * np.pi and 0 < half < 3 * np.pi):
        r = cdd(D, mu_over_ext, chi_s, kind)                    # quadrature failed near chi_s -> pi/2: fall back to the ODE
        budget, half = r["eta_cone_minus_eta_b"], r["eta_max_minus_eta_b"]
    return dict(valid=True, D=D, kind=kind, mu_over_ext=mu_over_ext, chi_s=chi_s, untrapped_initially=bool(a_max * np.sin(chi_s) > rp),
                R0_over_rplus=float(a_max * np.sin(chi_s) / rp), a_b=a_b, a_max=a_max, cone_budget=float(budget), half_period=float(half),
                chi_reaches_max=float(max(budget - half, 0.0)), chi_full_cycle=float(max(budget - 2 * half, 0.0)),
                max_inside_D=bool(budget > half), full_cycle_inside_D=bool(budget > 2 * half))


def cycle_map(D=10, kind="dymnikova", mus=None, chis=None) -> dict:
    mus = np.array(mus if mus is not None else np.concatenate([np.linspace(1.005, 1.5, 34), np.geomspace(1.6, 30, 16)]))
    chis = np.array(chis if chis is not None else np.linspace(0.3, 1.565, 40))
    budget = np.full((len(mus), len(chis)), np.nan); half = budget.copy(); untr = np.zeros_like(budget, dtype=bool)
    for i, m in enumerate(mus):
        for j, c in enumerate(chis):
            r = cdd_fast(D, m, c, kind)
            if r["valid"]:
                budget[i, j], half[i, j], untr[i, j] = r["cone_budget"], r["half_period"], r["untrapped_initially"]
    return dict(D=D, kind=kind, mu_over_ext=mus.tolist(), chi_s=chis.tolist(), cone_budget=budget.tolist(), half_period=half.tolist(),
                untrapped=untr.tolist(), frac_max_inside=float(np.mean((budget > half) & untr)),
                frac_full_inside=float(np.mean((budget > 2 * half) & untr)))
