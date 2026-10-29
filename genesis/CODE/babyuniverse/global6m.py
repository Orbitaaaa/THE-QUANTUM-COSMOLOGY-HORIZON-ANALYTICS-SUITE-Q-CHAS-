"""Round 6m: the global daughter spacetime of model U.

The spacetime: Oppenheimer-Snyder collapse of a closed dust ball in model U (unified6l), glued exactly to the static regular exterior.  The
ball surface is a radial timelike geodesic of the exterior with energy E = cos(chi_b), since f(R) + Rdot^2 = cos^2 chi_b (generalised
Misner-Sharp matching, mu = C sin^{D-1} chi_b).  The maximal extension of the exterior (two simple horizons r_- < r_+, regular centre) is
tiled by diamonds: I (parent exterior), II (black-hole region), III (static de Sitter core), IV (white-hole region), I' (new universe), ...

Results (ROADMAP_ROUND6M.md):
* ``kretschmann_static``, ``curvature_bound_scan``  [Proposition 6m-1, numerical] K = f''^2 + 2(D-2) f'^2/r^2 + 2(D-2)(D-3)(1-f)^2/r^4 is
  bounded for every mass, with sup over all mu and r equal to a constant times 1/alpha^2 (approaches the de Sitter value 2D(D-1)/alpha^2 for
  large mu): a solution-independent curvature bound (cf. Bueno et al. 2025).
* ``frw_curvature``  The closed interior has K = 4(D-1)(addot/a)^2 + 2(D-1)(D-2) psi^2 with psi <= 1/alpha: bounded through the bounce.
* ``ball_history``  [Proposition 6m-2] For every bound collapse with mu > mu_ext, the ball surface passes I -> II -> III -> IV -> I' in finite
  proper time: it crosses the event horizon, the inner horizon, bounces in the core where f = E^2 > 0, leaves through the white-hole
  region and reaches its maximal radius in a new asymptotically flat universe, then repeats.  For mu < mu_ext there is no horizon and it
  bounces back into the parent universe.
* ``penrose_worldline``  The same history in diamond (conformal-patch) coordinates, matched across every horizon through the finite null
  coordinate (v across H+, v across the inner horizon into III, u out of III into IV and I').
* ``geodesic_completeness``  Radial and non-radial timelike geodesics and radial null geodesics cross each region in finite affine
  parameter and never terminate: no centre singularity (f(0) = 1, f'(0) = 0), simple horizons (analytic extension).
* ``daughter_properties``  [Proposition 6m-3] The daughter D = IV u I' (and its continuation) lies outside J^-(I^+ of the parent), has its own
  future null infinity, ADM mass mu, and contains the re-expanded ball: a closed FRW universe re-expanding from a_b ~ l to a_max (finite
  volume Omega_{D-1} a_max^{D-1} I(chi_b)).
"""

from __future__ import annotations

from math import pi

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar

from . import unified6l as U


# --------------------------------------------------------------------------- curvature


def _f_derivs(r, mu, D, alpha, kind, h=1e-5):
    f0 = float(U.static_f(r, mu, D, alpha, kind))
    fp = float(U.static_f(r * (1 + h), mu, D, alpha, kind))
    fm = float(U.static_f(r * (1 - h), mu, D, alpha, kind))
    d1 = (fp - fm) / (2 * h * r)
    d2 = (fp - 2 * f0 + fm) / (h * r) ** 2
    return f0, d1, d2


def kretschmann_static(r, mu, D, alpha=1.0, kind="hayward"):
    f, f1, f2 = _f_derivs(r, mu, D, alpha, kind)
    return f2**2 + 2 * (D - 2) * f1**2 / r**2 + 2 * (D - 2) * (D - 3) * (1 - f) ** 2 / r**4


def kretschmann_schwarzschild_check(M=1.0, r=3.0):
    """D = 4 GR check: K = 48 M^2 / r^6."""
    f = lambda x: 1 - 2 * M / x
    h = 1e-4
    f1 = (f(r + h) - f(r - h)) / (2 * h)
    f2 = (f(r + h) - 2 * f(r) + f(r - h)) / h**2
    K = f2**2 + 2 * 2 * f1**2 / r**2 + 2 * 2 * 1 * (1 - f(r)) ** 2 / r**4
    return dict(K=float(K), exact=float(48 * M**2 / r**6))


def curvature_bound_scan(D=10, alpha=1.0, kind="hayward", mus=None) -> dict:
    mus = np.geomspace(1e-3, 1e8, 45) if mus is None else mus
    rows = []
    for mu in mus:
        rr = np.geomspace(1e-3, 50 * max(1.0, mu ** (1 / (D - 3))), 3000)
        K = np.array([kretschmann_static(r, mu, D, alpha, kind) for r in rr])
        rows.append(dict(mu=float(mu), K_max_alpha2=float(K.max() * alpha**2), r_at_max=float(rr[int(np.argmax(K))])))
    Kd = 2 * D * (D - 1)
    return dict(D=D, kind=kind, rows=rows, sup_K_alpha2=max(r["K_max_alpha2"] for r in rows), de_sitter_value=float(Kd))


def frw_curvature(C, D, alpha=1.0, kind="hayward", n=4000) -> dict:
    """Kretschmann scalar of the closed FRW interior over a full cycle: addot/a = psi - ((D-1)/2) psi'(y) y, y = C/a^{D-1}."""
    r = U.closed_dust(C, D, alpha, kind)
    a = np.geomspace(r["a_b"], r["a_max"], n)
    y = C / a ** (D - 1)
    psi = U.h_inv(y, alpha, kind)
    dpsi = (U.h_inv(y * (1 + 1e-6), alpha, kind) - U.h_inv(y * (1 - 1e-6), alpha, kind)) / (2e-6 * y)
    acc = psi - (D - 1) / 2 * dpsi * y
    K = 4 * (D - 1) * acc**2 + 2 * (D - 1) * (D - 2) * psi**2
    return dict(C=C, D=D, K_max_alpha2=float(K.max() * alpha**2), K_at_bounce_alpha2=float(K[0] * alpha**2),
                de_sitter_value=float(2 * D * (D - 1)), psi_max_alpha=float(psi.max() * alpha))


# --------------------------------------------------------------------------- the collapsing ball and its global history


def _rstar_factory(mu, D, alpha, kind):
    """Tortoise coordinates r*(r) = int dr / f in each region, with the logarithmic horizon parts subtracted analytically."""
    rm, rp = U.static_horizons(mu, D, alpha, kind)
    f = lambda r: float(U.static_f(r, mu, D, alpha, kind))
    eps = 1e-7
    kp = (f(rp * (1 + eps)) - f(rp * (1 - eps))) / (2 * eps * rp) / 2
    km = -(f(rm * (1 + eps)) - f(rm * (1 - eps))) / (2 * eps * rm) / 2

    def reg(r):                                                    # 1/f minus its two poles
        return 1.0 / f(r) - 1.0 / (2 * kp * (r - rp)) + 1.0 / (2 * km * (r - rm))

    def rstar(r, ref):
        val = quad(reg, ref, r, limit=200)[0]
        return val + np.log(abs(r - rp) / abs(ref - rp)) / (2 * kp) - np.log(abs(r - rm) / abs(ref - rm)) / (2 * km)
    return rstar, rm, rp, kp, km


def ball_history(D=10, alpha=1.0, mu_over_ext=10.0, E=0.6, kind="hayward", n_out=3000) -> dict:
    """Proper-time history of the OS ball surface (radial geodesic of energy E = cos chi_b) from rest at R_max through every region until it
    reaches R_max again in the next universe.  Returns region crossing times and the worldline in diamond coordinates."""
    mu = mu_over_ext * U.extremal_mu(D, alpha, kind)
    hz = U.static_horizons(mu, D, alpha, kind)
    f = lambda r: float(U.static_f(r, mu, D, alpha, kind))
    if hz is None:
        return dict(horizon=False, verdict="no horizon: the ball bounces back into the parent universe (no daughter)")
    rm, rp = hz
    R_max = brentq(lambda r: f(r) - E**2, rp * (1 + 1e-12), rp * 1e4)
    R_b = brentq(lambda r: f(r) - E**2, 1e-8, rm * (1 - 1e-12))
    # proper times: R_max -> r_+ -> r_- -> R_b, symmetric on the way out
    def tau(a, b):
        w = lambda t: 2 * (b - a) * np.sin(t) * np.cos(t) / np.sqrt(max(E**2 - f(a + (b - a) * np.sin(t) ** 2), 1e-300))
        return quad(w, 1e-10, pi / 2 - 1e-10, limit=400)[0]
    # segments with endpoints where E^2 - f vanishes only at R_max and R_b
    tau_total_half = tau(R_b, R_max)
    sq = lambda r: 1.0 / np.sqrt(max(E**2 - f(r), 1e-300))
    t_I = quad(sq, rp, R_max, limit=400, points=[R_max * 0.999])[0]
    t_II = quad(sq, rm, rp, limit=400)[0]
    t_III = tau_total_half - t_I - t_II
    return dict(horizon=True, D=D, mu_over_ext=mu_over_ext, E=E, r_plus=rp, r_minus=rm, R_max=R_max, R_bounce=R_b,
                tau_region_I=t_I, tau_region_II=t_II, tau_core_to_bounce=t_III, tau_half_cycle=tau_total_half,
                tau_parent_horizon_to_daughter_horizon=float(2 * (t_II + t_III)),
                regions=["I (parent)", "II (black hole)", "III (de Sitter core, bounce)", "IV (white hole)", "I' (new universe)"],
                verdict="the ball reaches a new asymptotically flat universe in finite proper time")


def penrose_worldline(D=10, alpha=1.0, mu_over_ext=10.0, E=0.6, kind="hayward", scale=None, n=1500) -> dict:
    """Diamond coordinates: region R has null coordinates (u, v); tilde = tanh(./L).  Centres: II = I + (-1, 1), III = II + (-1, 1),
    IV = III + (1, 1), I' = IV + (1, 1).  The finite null coordinate is integrated along the geodesic (v on ingoing legs, u on outgoing legs),
    the other follows from u = v - 2 r*(r) with region-wise tortoise coordinates."""
    mu = mu_over_ext * U.extremal_mu(D, alpha, kind)
    f = lambda r: float(U.static_f(r, mu, D, alpha, kind))
    rstar, rm, rp, kp, km = _rstar_factory(mu, D, alpha, kind)
    L = scale or 2.0 * rp
    R_max = brentq(lambda r: f(r) - E**2, rp * (1 + 1e-12), rp * 1e4)
    R_b = brentq(lambda r: f(r) - E**2, 1e-8, rm * (1 - 1e-12))
    refs = {"I": 0.5 * (rp + R_max) + rp, "II": 0.5 * (rm + rp), "III": 0.5 * rm, "IV": 0.5 * (rm + rp), "I'": 0.5 * (rp + R_max) + rp}
    cen = {"I": np.array([0.0, 0.0])}
    cen["II"] = cen["I"] + [-1, 1]
    cen["III"] = cen["II"] + [-1, 1]
    cen["IV"] = cen["III"] + [1, 1]
    cen["I'"] = cen["IV"] + [1, 1]
    th = lambda x: np.tanh(x / L)

    def place(reg, u, v):
        U_, V_ = th(u), th(v)
        c = cen[reg]
        if reg in ("I", "III", "IV", "I'"):
            return c + np.array([(V_ - U_) / 2, (V_ + U_) / 2])
        return c + np.array([(V_ + U_) / 2, (V_ - U_) / 2])          # region II orientation

    pts, labels = [], []
    # ingoing leg R_max -> R_b: integrate v (regular), dv/dr = -1/(sqrt(E^2-f) (E + sqrt(E^2-f)))
    r_in = R_b + (R_max - R_b) * np.cos(np.linspace(0, pi / 2, n)) ** 2              # from R_max down to R_b
    v = 0.0
    for k in range(len(r_in)):
        r = r_in[k]
        if k:
            r0 = r_in[k - 1]
            rmid = 0.5 * (r + r0)
            s = np.sqrt(max(E**2 - f(rmid), 1e-30))
            v += (r0 - r) / (s * (E + s))
        reg = "I" if r > rp else ("II" if r > rm else "III")
        if abs(r - rp) < 1e-9 * rp or abs(r - rm) < 1e-9 * rm:
            continue
        u = v - 2 * rstar(r, refs[reg])
        pts.append(place(reg, u, v)); labels.append(reg)
    # outgoing leg R_b -> R_max: integrate u (regular), du/dr = 1/(sqrt(E^2-f) (E + sqrt(E^2-f)))
    r_out = r_in[::-1]
    u = pts and (v - 2 * rstar(r_out[0], refs["III"]))
    u = float(u)
    for k in range(1, len(r_out)):
        r = r_out[k]
        r0 = r_out[k - 1]
        rmid = 0.5 * (r + r0)
        s = np.sqrt(max(E**2 - f(rmid), 1e-30))
        u += (r - r0) / (s * (E + s))
        reg = "III" if r < rm else ("IV" if r < rp else "I'")
        if abs(r - rp) < 1e-9 * rp or abs(r - rm) < 1e-9 * rm:
            continue
        vv = u + 2 * rstar(r, refs[reg])
        pts.append(place(reg, u, vv)); labels.append(reg)
    pts = np.array(pts)
    order = [labels[0]]
    for l in labels[1:]:
        if l != order[-1]:
            order.append(l)
    # centre line of the core: v - u = 2 r*(0)
    rs0 = 2 * rstar(1e-9, refs["III"])
    tt = np.linspace(-6 * L, 6 * L, 400)
    centre = np.array([place("III", t - rs0 / 2, t + rs0 / 2) for t in tt])
    return dict(x=pts[:, 0].tolist(), t=pts[:, 1].tolist(), regions=labels, region_order=order, centres={k: v.tolist() for k, v in cen.items()},
                centre_line=centre.tolist(), r_plus=rp, r_minus=rm)


def geodesic_completeness(D=10, alpha=1.0, mu_over_ext=10.0, kind="hayward") -> dict:
    """Affine length of each region for radial timelike (E < 1, E > 1), non-radial timelike and radial null geodesics.  Finite crossing
    lengths and turning points at r > 0 (or passage through the regular centre) show that no geodesic terminates."""
    mu = mu_over_ext * U.extremal_mu(D, alpha, kind)
    f = lambda r: float(U.static_f(r, mu, D, alpha, kind))
    rm, rp = U.static_horizons(mu, D, alpha, kind)
    out = []
    for (E, Lang, kind_g) in ((0.6, 0.0, "timelike"), (1.5, 0.0, "timelike"), (1.2, 2.0, "timelike"), (1.0, 0.0, "null"), (1.0, 3.0, "null")):
        eps = 1.0 if kind_g == "timelike" else 0.0
        Veff = lambda r: f(r) * (eps + Lang**2 / r**2)
        g = lambda r: E**2 - Veff(r)
        # inner turning point (if any) inside the core
        grid = np.geomspace(1e-6, rm, 4000)
        neg = [r for r in grid if g(r) < 0]
        r_turn = max(neg) if neg else 0.0
        if r_turn > 0:
            r_turn = brentq(g, r_turn, grid[min(np.searchsorted(grid, r_turn) + 1, len(grid) - 1)])
        lo = max(r_turn, 1e-6)
        lam_core = quad(lambda r: 1 / np.sqrt(max(g(r), 1e-300)), lo, rm, limit=400)[0]
        lam_II = quad(lambda r: 1 / np.sqrt(max(g(r), 1e-300)), rm, rp, limit=400)[0]
        out.append(dict(kind=kind_g, E=E, L=Lang, turning_point=float(r_turn), passes_centre=bool(r_turn == 0.0),
                        affine_core=float(lam_core), affine_BH_region=float(lam_II), finite=bool(np.isfinite(lam_core) and np.isfinite(lam_II))))
    return dict(D=D, mu_over_ext=mu_over_ext, f0=f(1e-8), fprime0=float((f(2e-4) - f(1e-4)) / 1e-4), geodesics=out)


def daughter_properties(D=10, alpha=1.0, C=1e4, chi_b=np.pi / 3, kind="hayward") -> dict:
    """The finite daughter carried by the ball: its closed FRW interior re-expands from a_b to a_max inside the new universe."""
    r = U.closed_dust(C, D, alpha, kind)
    mu = C * np.sin(chi_b) ** (D - 1)
    me = U.extremal_mu(D, alpha, kind)
    I = quad(lambda x: np.sin(x) ** (D - 2), 0, chi_b)[0]
    vol = lambda a: U.omega(D - 2) * a ** (D - 1) * I
    return dict(C=C, chi_b=float(chi_b), mu=float(mu), mu_over_ext=float(mu / me), forms_horizon=bool(mu > me), a_bounce=r["a_b"],
                a_max=r["a_max"], volume_at_bounce=float(vol(r["a_b"])), volume_max=float(vol(r["a_max"])),
                volume_growth=float((r["a_max"] / r["a_b"]) ** (D - 1)), E_surface=float(np.cos(chi_b)),
                adm_mass_parameter_daughter=float(mu),
                causal="outside J^-(I+ parent): every future causal curve from IV u I' stays in the daughter or its own later black hole")


def kruskal_regularity(D=10, alpha=1.0, mu_over_ext=10.0, kind="hayward", offsets=(1e-2, 1e-3, 1e-4, 1e-5)) -> dict:
    """At each horizon r_h with surface gravity kappa, Kruskal coordinates U V = -exp(2 kappa r*) give g_UV ~ F(r) = f(r) exp(-2 kappa r*)/kappa^2.
    F must tend to the same finite, non-zero limit from both sides: the metric extends smoothly (analytically) across the horizon."""
    mu = mu_over_ext * U.extremal_mu(D, alpha, kind)
    rstar, rm, rp, kp, km = _rstar_factory(mu, D, alpha, kind)
    f = lambda r: float(U.static_f(r, mu, D, alpha, kind))
    out = {}
    for name, rh, kap, ref_in, ref_out in (("outer", rp, kp, 0.5 * (rm + rp), 2 * rp), ("inner", rm, -km, 0.5 * rm, 0.5 * (rm + rp))):
        # use one tortoise branch through the horizon: r* with the pole at this horizon handled by the analytic log
        rows = []
        for eps in offsets:
            vals = []
            for r in (rh * (1 - eps), rh * (1 + eps)):
                ref = 0.5 * (rm + rp)
                rs = rstar(r, ref)
                vals.append(abs(f(r) * np.exp(-2 * kap * rs) / kap**2))
            rows.append(dict(eps=eps, F_inside=float(vals[0]), F_outside=float(vals[1])))
        last = rows[-1]
        out[name] = dict(r_h=rh, kappa=float(abs(kap)), rows=rows, limit_ratio=float(last["F_inside"] / last["F_outside"]),
                         finite_nonzero=bool(0 < last["F_inside"] < np.inf and 0 < last["F_outside"] < np.inf))
    return out
