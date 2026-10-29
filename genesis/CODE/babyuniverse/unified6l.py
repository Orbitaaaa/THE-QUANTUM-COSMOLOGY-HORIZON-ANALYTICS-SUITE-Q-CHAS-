"""Round 6l, priorities 1 and 3: one consistent higher-dimensional model that contains both doors, whether its black-hole interiors inflate,
and whether one black hole can hold more than one daughter.

Model U (ROADMAP_ROUND6L.md, section 1).  D = 4 + n dimensional Einstein gravity (P1: large extra dimensions, fundamental scale M_D) plus an
infinite tower of quasi-topological curvature corrections (P2: singularity resolution), the pure-gravity theories of Bueno, Cano & Hennigar
(arXiv:2403.04827) whose spherically symmetric sector has second-order equations, a Birkhoff theorem and Israel-type junction conditions
(Bueno, Cano, Hennigar & Murcia, PRD 111, 104009 (2025); PRL 134, 181401 (2025)).  On spherical symmetry everything is fixed by one function
h(psi), psi = (1 - (grad r)^2)/r^2, through the generalised Misner-Sharp relation

    h(psi) = mu / r^{D-1},      mu = 16 pi G_D M / ((D-2) Omega_{D-2})     (GR: h = psi, f = 1 - mu/r^{D-3})

and, for a closed FRW region (psi = H^2 + 1/a^2), the generalised Friedmann equation h(H^2 + 1/a^2) = 16 pi G_D rho/((D-1)(D-2)) (the form
used by Bueno, Cano, Hennigar, Murcia & Vicente-Cano, arXiv:2505.09680).  Two towers with a pole at psi_max = 1/alpha are implemented:
"hayward"   h = psi/(1 - alpha psi),        psi = y/(1 + alpha y)                   (D-dimensional Hayward black hole)
"dymnikova" h = -ln(1 - alpha psi)/alpha,  psi = (1 - exp(-alpha y))/alpha          (D-dimensional Dymnikova-type core)
The limiting length is l = sqrt(alpha) = c_l / M_D.

Results [labels as in RESULTS.md]:
* ``static_horizons``, ``extremal_mu``, ``temperature``   Regular black holes: de Sitter core f ~ 1 - r^2/alpha, two horizons above an extremal
  mass, T -> 0 there (a remnant).  Hayward: mu_ext = alpha^{(D-3)/2} ((D-1)/(D-3))^{(D-1)/2} (D-3)/2, T = ((D-3) r^2 - (D-1) alpha)/(4 pi r^3).
* ``closed_dust``  [Proposition 6l-1] A closed dust interior obeys adot^2 = a^2 psi(C/a^{D-1}) - 1.  Its turning points are the horizon radii
  of the static hole of mass C (same polynomial); it bounces at a_b -> sqrt(alpha) (universal: any tower with a pole at 1/alpha), goes
  through a transient accelerating phase of N_acc = ln(a_infl/a_b) e-folds, a_infl^{D-1} = 2 alpha C/(D-3) (Hayward, exact), i.e.
  N_acc -> ((D-3)/(D-1)) ln(a_max/l) + ln(2/(D-3))/(D-1), re-expands to its original size and cycles.  No lasting inflation.
* ``vacuum_threshold`` [Proposition 6l-2] With interior vacuum energy (Lambda term in the density) the bounced region inflates for ever iff
  Lambda > Lambda_c(C); GR limit Lambda_c = [2 ((D-3)/(D-1))^{(D-1)/2} / ((D-3) C)]^{2/(D-3)} (D = 4: l_Lambda < (3 sqrt3/2) C), the same
  function as Hayward extremality with alpha -> 1/Lambda.  Because model U has no singularities, the Farhi-Guth past-singularity obstruction
  is absent: a trapped false-vacuum region need not come from a singularity.  The vacuum must still be as dense as the hole (Nariai-partner
  scale).
* ``vacuum_collapse_interior``  For vacuum collapse (colliding shock waves, thin shells) Birkhoff fixes the interior to the static regular hole:
  trapped region -> inner horizon -> static de Sitter core of radius < l.  Not inflating; matter bounces in the core and emerges through a
  white-hole region into a new asymptotic region (Bueno et al.).
* ``bounce_untrapped_lemma``  [Proposition 6l-3, kinematic, any theory] At a bounce of a radial timelike worldline (dR/dtau = 0),
  (grad R)^2 = (n.grad R)^2 >= 0: the sphere is untrapped.  A bounce inside a black hole therefore needs an untrapped core bounded by an inner
  apparent horizon.
* ``ecsk_os_mass_deficit``  [Proposition 6l-4] ECSK homogeneous Oppenheimer-Snyder: the ball's Misner-Sharp mass falls from (sin^3 chi/2) A to
  (sin^3 chi/2) a_b at the bounce (spin pressure does work), so it cannot be matched to a vacuum exterior of fixed mass: the bounce requires an
  inhomogeneous, untrapped core inside still-falling layers (corrects the homogeneous picture of Round 6k, F61).
* ``two_epoch_infall``  [Proposition 6l-5] In ingoing Eddington-Finkelstein coordinates the bounced surface of the first collapse reaches the
  inner horizon only at v = +infinity, so any later infall at finite advanced time v meets it inside the same core: it enters and modifies the
  existing daughter.  Computed for test shells in model U.
* ``fmm_patch_separation``  For a limiting-curvature (FMM) spacelike transition joined to the expanding branch of de Sitter, two infall events
  separated by Delta v enter causally disjoint regions iff Delta v > l ln((r0 + l)/(r0 - l)) ~ 2 l (l/r_s)^{1/3}; on the contracting branch
  (extrinsic-curvature continuous, the FMM choice) the daughter is one closed de Sitter world and their futures meet.
* ``add_scales``  ADD radius (2 pi R)^n = M_P^2 / M_D^{n+2}: collider holes (R_S ~ 1/M_D) are 10^4-10^5 times smaller, so the D-dimensional
  asymptotically flat description applies.
"""

from __future__ import annotations

from math import gamma, pi

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar

M_PL_GEV = 1.22089e19
HBARC_GEV_M = 1.97327e-16


# --------------------------------------------------------------------------- the tower


def h_inv(y, alpha=1.0, kind="hayward"):
    """psi(y) solving h(psi) = y."""
    y = np.asarray(y, dtype=float)
    if kind == "hayward":
        return y / (1.0 + alpha * y)
    if kind == "dymnikova":
        return -np.expm1(-alpha * y) / alpha
    if kind == "gr":
        return y
    raise ValueError(kind)


def h_fun(psi, alpha=1.0, kind="hayward"):
    psi = np.asarray(psi, dtype=float)
    if kind == "hayward":
        return psi / (1.0 - alpha * psi)
    if kind == "dymnikova":
        return -np.log1p(-alpha * psi) / alpha
    if kind == "gr":
        return psi
    raise ValueError(kind)


def omega(k):
    """Area of the unit k-sphere."""
    return 2 * pi ** ((k + 1) / 2) / gamma((k + 1) / 2)


# --------------------------------------------------------------------------- static regular black holes


def static_f(r, mu, D, alpha=1.0, kind="hayward"):
    r = np.asarray(r, dtype=float)
    return 1.0 - r**2 * h_inv(mu / r ** (D - 1), alpha, kind)


def _g_peak(mu, D, alpha, kind):
    """Maximum of r^2 psi(r) (the horizon function) and where it sits."""
    res = minimize_scalar(lambda s: -float(np.exp(2 * s) * h_inv(mu / np.exp(s * (D - 1)), alpha, kind)),
                          bounds=(np.log(1e-4 * np.sqrt(alpha)), np.log(1e4 * max(np.sqrt(alpha), mu ** (1 / (D - 3))))),
                          method="bounded", options=dict(xatol=1e-12))
    return float(-res.fun), float(np.exp(res.x))


def static_horizons(mu, D, alpha=1.0, kind="hayward"):
    """(r_minus, r_plus) or None below the extremal mass."""
    gmax, rp = _g_peak(mu, D, alpha, kind)
    if gmax <= 1.0:
        return None
    g = lambda r: float(r**2 * h_inv(mu / r ** (D - 1), alpha, kind)) - 1.0
    r_minus = brentq(g, 1e-6 * np.sqrt(alpha), rp, xtol=1e-14)
    hi = rp
    while g(hi) > 0:
        hi *= 2.0
    r_plus = brentq(g, rp, hi, xtol=1e-14)
    return r_minus, r_plus


def extremal_mu(D, alpha=1.0, kind="hayward"):
    if kind == "hayward":
        return alpha ** ((D - 3) / 2) * ((D - 1) / (D - 3)) ** ((D - 1) / 2) * (D - 3) / 2
    f = lambda lmu: _g_peak(np.exp(lmu), D, alpha, kind)[0] - 1.0
    return float(np.exp(brentq(f, np.log(1e-3 * alpha ** ((D - 3) / 2)), np.log(1e4 * alpha ** ((D - 3) / 2)), xtol=1e-13)))


def temperature(mu, D, alpha=1.0, kind="hayward"):
    """Hawking temperature f'(r_+)/4 pi (0 at and below extremality) and the inner-horizon surface gravity."""
    hz = static_horizons(mu, D, alpha, kind)
    if hz is None:
        return dict(T=0.0, kappa_minus=0.0, r_plus=None, r_minus=None)
    rm, rp = hz
    d = lambda r: (static_f(r * (1 + 1e-6), mu, D, alpha, kind) - static_f(r * (1 - 1e-6), mu, D, alpha, kind)) / (2e-6 * r)
    return dict(T=float(d(rp) / (4 * pi)), kappa_minus=float(abs(d(rm)) / 2), r_plus=rp, r_minus=rm)


def hayward_T_analytic(r_plus, D, alpha=1.0):
    return ((D - 3) * r_plus**2 - (D - 1) * alpha) / (4 * pi * r_plus**3)


# --------------------------------------------------------------------------- closed FRW interiors (Oppenheimer-Snyder)


def _adot2(a, C, D, alpha, kind, Lam):
    a = np.asarray(a, dtype=float)
    return a**2 * h_inv(C / a ** (D - 1) + Lam, alpha, kind) - 1.0


def closed_dust(C, D, alpha=1.0, kind="hayward", Lam=0.0, a_hi_factor=1e3) -> dict:
    """Bounce, accelerating phase, turnaround and proper times of a closed dust (+ vacuum) interior."""
    sa = np.sqrt(alpha)
    a_scale = max(sa, C ** (1 / (D - 3)) if C > 0 else sa)
    grid = np.exp(np.linspace(np.log(1e-3 * sa), np.log(a_hi_factor * a_scale), 40000))
    F = _adot2(grid, C, D, alpha, kind, Lam)
    s = np.sign(F)
    idx = np.where(np.diff(s) != 0)[0]
    roots = [brentq(lambda a: float(_adot2(a, C, D, alpha, kind, Lam)), grid[i], grid[i + 1], xtol=1e-15 * grid[i] + 1e-300) for i in idx]
    if not roots or F[0] > 0:
        return dict(C=C, D=D, kind=kind, Lam=Lam, a_b=None, fate="no bounce")
    a_b = roots[0]
    a_max = roots[1] if len(roots) > 1 else None
    G = np.gradient(grid**2 * h_inv(C / grid ** (D - 1) + Lam, alpha, kind), grid)
    sel = np.where((grid > a_b) & (np.diff(np.sign(G), append=np.sign(G[-1])) != 0))[0]
    a_infl = float(grid[sel[0]]) if len(sel) and (a_max is None or grid[sel[0]] < a_max) else None
    N_acc = float(np.log(a_infl / a_b)) if a_infl else None

    def tau_between(a1, a2):
        w = lambda t: 2 * (a2 - a1) * np.sin(t) * np.cos(t) / np.sqrt(max(float(_adot2(a1 + (a2 - a1) * np.sin(t) ** 2, C, D, alpha, kind, Lam)), 1e-300))
        return quad(w, 1e-9, pi / 2 - 1e-9, limit=400)[0]

    tau_half = tau_between(a_b, a_max) if a_max else None
    fate = "bounces, accelerates for N_acc e-folds, re-expands to a_max and cycles" if a_max else "bounces and inflates for ever"
    out = dict(C=C, D=D, kind=kind, Lam=Lam, a_b=float(a_b), a_max=float(a_max) if a_max else None, a_infl=a_infl, N_acc=N_acc,
               tau_half_cycle=tau_half, fate=fate, H2_max_bound=1.0 / alpha)
    if kind == "hayward" and Lam == 0.0:
        out["a_infl_exact"] = float((2 * alpha * C / (D - 3)) ** (1 / (D - 1)))
        if a_max:
            out["N_acc_asymptotic"] = float((D - 3) / (D - 1) * np.log(a_max / sa) + np.log(2 / (D - 3)) / (D - 1))
    return out


def vacuum_threshold(C, D, alpha=1.0, kind="hayward") -> dict:
    """Smallest interior vacuum term Lambda (in the density units of h) for which the bounced region inflates for ever."""
    def recollapses(Lam):
        return closed_dust(C, D, alpha, kind, Lam)["a_max"] is not None
    lo, hi = 0.0, 1.0 / alpha
    while recollapses(hi):
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if recollapses(mid) else (lo, mid)
    Lam_gr = (2 * ((D - 3) / (D - 1)) ** ((D - 1) / 2) / ((D - 3) * C)) ** (2 / (D - 3))
    return dict(C=C, D=D, Lam_c=float(hi), Lam_c_gr_limit=float(Lam_gr), ratio=float(hi / Lam_gr), l_Lam_c=float(1 / np.sqrt(hi)))


def vacuum_collapse_interior(mu, D, alpha=1.0, kind="hayward") -> dict:
    """What a vacuum (shock-wave or thin-shell) collapse leaves inside, by Birkhoff: the static regular hole."""
    hz = static_horizons(mu, D, alpha, kind)
    if hz is None:
        return dict(mu=mu, horizon=False, interior="no horizon: horizonless regular object (matter bounces and re-emerges)")
    rm, rp = hz
    return dict(mu=mu, horizon=True, r_plus=rp, r_minus=rm, core_radius_over_l=float(rm / np.sqrt(alpha)),
                interior="trapped shell r_- < r < r_+, then a static de Sitter core of radius r_- < l: not inflating",
                core_static=True, inflating=False)


# --------------------------------------------------------------------------- daughters: kinematics


def bounce_untrapped_lemma(n_samples=2000, seed=1) -> dict:
    """(grad R)^2 = -(u.grad R)^2 + (n.grad R)^2 in an orthonormal radial frame; at a bounce u.grad R = 0 so (grad R)^2 >= 0.  Random checks
    of the frame identity for 2D metrics, then the lemma."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(n_samples):
        A = rng.normal(size=(2, 2))
        g = np.array([[-np.exp(rng.normal()), 0.0], [0.0, np.exp(rng.normal())]])
        L = A if abs(np.linalg.det(A)) > 1e-3 else np.eye(2)
        gm = L.T @ g @ L                                         # a generic Lorentzian 2-metric
        dR = rng.normal(size=2)
        ginv = np.linalg.inv(gm)
        w, V = np.linalg.eig(gm)
        # orthonormal frame from the eigen-decomposition
        e = [V[:, i] / np.sqrt(abs(w[i])) for i in range(2)]
        t_i = int(np.argmin(w))
        u, nvec = e[t_i], e[1 - t_i]
        lhs = dR @ ginv @ dR
        rhs = -(u @ dR) ** 2 + (nvec @ dR) ** 2
        worst = max(worst, abs(lhs - rhs) / (1 + abs(lhs)))
    return dict(frame_identity_max_rel_error=float(worst),
                lemma="at dR/dtau = 0 along a radial timelike worldline, (grad R)^2 = (n.grad R)^2 >= 0: the bounce sphere is untrapped")


def ecsk_os_mass_deficit(A=1.0, S=1e-6, chi_b=np.pi / 4) -> dict:
    """ECSK closed dust with spin, adot^2 = A/a - S/a^4 - 1 (Round 6k units, G = c = 1).  Misner-Sharp mass of the ball:
    M_MS = (R/2)(1 - (grad R)^2) = (sin^3 chi_b / 2)(A - S/a^3).  At the bounce A - S/a_b^3 = a_b."""
    a_b = brentq(lambda a: A / a - S / a**4 - 1.0, (S / A) ** (1 / 3) * 0.5, (S / A) ** (1 / 3) * 4)
    M0 = np.sin(chi_b) ** 3 / 2 * A
    Mb = np.sin(chi_b) ** 3 / 2 * (A - S / a_b**3)
    R_b = a_b * np.sin(chi_b)
    return dict(A=A, S=S, chi_b=chi_b, a_b=float(a_b), M_initial=float(M0), M_at_bounce=float(Mb), deficit_fraction=float(1 - Mb / M0),
                compactness_at_bounce=float(2 * Mb / R_b), exterior_compactness_if_vacuum=float(2 * M0 / R_b),
                verdict="the homogeneous ball cannot sit in a vacuum exterior of fixed mass: an inhomogeneous untrapped core is required")


def _geodesic_vr(mu, D, alpha, kind, E, R0, v0, tau_max, n_out=4000):
    """Radial timelike geodesic of energy E in ingoing EF coordinates (v, r), from rest at R0 (f(R0) = E^2): falls, bounces where f = E^2
    in the core, then moves out toward the inner horizon, which it reaches only as v -> +infinity.  Returns the ingoing and outgoing legs."""
    fF = lambda r: float(static_f(r, mu, D, alpha, kind))
    rm, rp = static_horizons(mu, D, alpha, kind)
    R_turn = brentq(lambda r: fF(r) - E**2, 1e-6 * np.sqrt(alpha), rm * (1 - 1e-12))

    def rhs_in(t, y):
        v, r = y
        k = np.sqrt(max(E**2 - fF(r), 0.0))
        return [1.0 / (E + k), -k]

    r_start = R0 * (1 - 1e-9)
    ev = lambda t, y: y[1] - R_turn * (1 + 1e-7)
    ev.terminal = True
    s_in = solve_ivp(rhs_in, (0, tau_max), [v0, r_start], events=ev, max_step=tau_max / n_out, rtol=1e-10, atol=1e-12)

    def rhs_out(t, y):
        v, r = y
        f = fF(r)
        k = np.sqrt(max(E**2 - f, 0.0))
        return [(E + k) / f, k]

    ev2 = lambda t, y: (rm - y[1]) / rm - 1e-9
    ev2.terminal = True
    y0 = [s_in.y[0, -1], R_turn * (1 + 2e-7)]
    s_out = solve_ivp(rhs_out, (0, tau_max), y0, events=ev2, max_step=tau_max / n_out, rtol=1e-10, atol=1e-12)
    return dict(v_in=s_in.y[0], r_in=s_in.y[1], tau_in=s_in.t, v_out=s_out.y[0], r_out=s_out.y[1], R_turn=R_turn, r_minus=rm, r_plus=rp)


def two_epoch_infall(D=10, alpha=1.0, mu_over_ext=20.0, E1=0.7, delay=50.0, kind="hayward") -> dict:
    """First collapse (surface = geodesic E1 from rest) and a later test shell (same orbit shape, released a coordinate time `delay` later,
    in units of l).  Find where the later shell meets the re-expanding first surface."""
    mu = mu_over_ext * extremal_mu(D, alpha, kind)
    fF = lambda r: float(static_f(r, mu, D, alpha, kind))
    rm, rp = static_horizons(mu, D, alpha, kind)
    R0 = brentq(lambda r: fF(r) - E1**2, rp * (1 + 1e-9), rp * 1e3)
    tau_max = 200 * rp
    g1 = _geodesic_vr(mu, D, alpha, kind, E1, R0, 0.0, tau_max)
    g2 = _geodesic_vr(mu, D, alpha, kind, E1, R0, delay, tau_max)
    # entry of shell 2 into the core (crossing r_-) on its ingoing leg
    i_in = np.where(g2["r_in"] < rm)[0]
    v2_entry = float(g2["v_in"][i_in[0]]) if len(i_in) else None
    # meeting: shell 2 ingoing r2(v) vs surface 1 outgoing r1(v).  Near r_- the outgoing leg obeys r_- - r1 ~ exp(-kappa_- v) (v diverges
    # only logarithmically), so it is continued analytically beyond the last integrated point.
    kap = temperature(mu, D, alpha, kind)["kappa_minus"]
    v1, r1 = g1["v_out"], g1["r_out"]
    v2, r2 = g2["v_in"], g2["r_in"]

    def r1_of_v(v):
        v = np.asarray(v, dtype=float)
        tail = rm - (rm - r1[-1]) * np.exp(-kap * (v - v1[-1]))
        return np.where(v <= v1[-1], np.interp(v, v1, r1), tail)

    meet = None
    if v2_entry is not None:
        v_start = float(g2["v_in"][max(i_in[0] - 1, 0)])                  # last point outside the core, where r2 >= r_- > r1
        vv = np.concatenate([np.linspace(v_start, v2_entry, 2000), np.linspace(v2_entry, v2[-1], 200000)[1:]])
        d = np.interp(vv, v2, r2) - r1_of_v(vv)
        k = np.where(np.diff(np.sign(d)) != 0)[0]
        if len(k):
            vm = float(vv[k[0]])
            log10_gap = (np.log10(rm - np.interp(vm, v1, r1)) if vm <= v1[-1]
                         else float(np.log10(rm - r1[-1]) - kap * (vm - v1[-1]) / np.log(10)))
            meet = dict(v=vm, r=float(r1_of_v(vm)), log10_rminus_minus_r=float(log10_gap))
    return dict(D=D, mu_over_ext=mu_over_ext, r_plus=rp, r_minus=rm, kappa_minus=kap, R_turn=g1["R_turn"], v1_last=float(v1[-1]),
                r1_last=float(r1[-1]), v2_entry_core=v2_entry, meeting=meet, meets_inside_core=bool(meet is not None and np.isfinite(meet["log10_rminus_minus_r"])),
                verdict="later infall meets the bounced matter inside the same core: it joins and modifies the existing daughter")


def fmm_patch_separation(r_s, l=1.0) -> dict:
    """Spacelike FMM junction at r0 = (r_s l^2)^{1/3} (continuous induced metric => Delta t_dS = Delta t_S = Delta v).  On the expanding de
    Sitter branch, points of the junction separated by Delta v have disjoint futures iff Delta v > l ln((r0 + l)/(r0 - l))."""
    r0 = (r_s * l**2) ** (1 / 3)
    if r0 <= l:
        return dict(r_s=r_s, valid=False)
    dv = l * np.log((r0 + l) / (r0 - l))
    return dict(r_s=r_s, r0=float(r0), valid=True, delta_v_disjoint=float(dv), approx_2l_l_over_rs_13=float(2 * l * (l / r_s) ** (1 / 3)),
                induced_metric_match=float(abs((r_s / r0 - 1) - (r0**2 / l**2 - 1))))


# --------------------------------------------------------------------------- P1 in numbers


def add_scales(n=6, M_D_TeV=8.0) -> dict:
    """ADD (Dimopoulos-Landsberg convention G_D = 1/M_D^{n+2}): (2 pi R)^n = M_P^2/M_D^{n+2}."""
    MD = M_D_TeV * 1e3
    R = (M_PL_GEV**2 / MD ** (n + 2)) ** (1 / n) / (2 * pi)
    return dict(n=n, M_D_TeV=M_D_TeV, R_extra_m=float(R * HBARC_GEV_M), inv_MD_m=float(HBARC_GEV_M / MD), ratio=float(R * MD))


def rs_over_MD_to_mass(rs_MD, n):
    """Inverse of the Dimopoulos-Landsberg radius: M/M_D for R_S M_D = rs_MD."""
    return (np.sqrt(pi) * rs_MD) ** (n + 1) * (n + 2) / (8 * gamma((n + 3) / 2))


def model_u_thresholds(n=6, c_l=1.0, kind="hayward") -> dict:
    """Horizon (daughter) threshold of model U in units of M_D: R_S,ext = mu_ext^{1/(D-3)} with alpha = (c_l/M_D)^2."""
    D = n + 4
    rhat = extremal_mu(D, 1.0, kind) ** (1 / (D - 3))
    return dict(n=n, c_l=c_l, kind=kind, rs_ext_over_l=float(rhat), M_ext_over_MD=float(rs_over_MD_to_mass(c_l * rhat, n)),
                M_fmm_over_MD=float(rs_over_MD_to_mass(c_l, n)))


# --------------------------------------------------------------------------- late infall meets the bounced matter: DTR and the curvature bound


def dtr_crossing(mu1, mu2, r_X, D, alpha=1.0, kind="hayward", log10_gap=None) -> dict:
    """Dray-'t Hooft-Redmount relation f_P f_F = f_L f_R at the crossing of the bounced (outgoing) first collapse, a thin shell of mass mu1
    with flat interior (f_L = 1), and a later ingoing shell mu2 (region R outside both: mu1 + mu2; past region P between them: mu1).  The
    relation is kinematic (continuity of r and a C^0 metric across both shells), so it holds in model U.  Model U bounds f from below,
    f_F > 1 - r_X^2/alpha (psi < 1/alpha); GR does not.  If log10_gap is given, r_X = r_-(mu1) - 10^log10_gap and f_P is taken from the
    near-horizon expansion f_P = 2 kappa_- (r_- - r_X)."""
    if log10_gap is not None:
        tp = temperature(mu1, D, alpha, kind)
        r_X = tp["r_minus"]
        f_P = 2 * tp["kappa_minus"] * 10.0**log10_gap
    else:
        f_P = float(static_f(r_X, mu1, D, alpha, kind))
    f_R = float(static_f(r_X, mu1 + mu2, D, alpha, kind))
    f_F = f_R / f_P if f_P != 0 else -np.inf
    f_min = 1.0 - r_X**2 / alpha
    allowed = bool(f_F > f_min)
    mu_F_U = float(r_X ** (D - 1) * h_fun((1 - f_F) / r_X**2, alpha, kind)) if allowed else np.inf
    mu_F_GR = float(r_X ** (D - 3) * (1 - f_F))
    return dict(r_X=float(r_X), f_P=float(f_P), f_R=f_R, f_F_required=float(f_F), f_F_min_model_U=float(f_min), allowed_in_model_U=allowed,
                mu_F_model_U=mu_F_U, mu_F_GR=mu_F_GR, mu_F_GR_over_mu1=float(mu_F_GR / mu1))


def dtr_critical_delay(D=10, alpha=1.0, mu_over_ext=20.0, mu2_over_mu1=1e-3, E1=0.7, kind="hayward") -> dict:
    """Delay (advanced time, units of l) after which the late shell's crossing would need curvature beyond the bound (mu_F -> infinity at the
    critical delay: region F saturates at psi = 1/alpha, i.e. becomes locally de Sitter; beyond it no thin-shell solution exists)."""
    mu1 = mu_over_ext * extremal_mu(D, alpha, kind)
    mu2 = mu2_over_mu1 * mu1

    def margin(delay):
        m = two_epoch_infall(D, alpha, mu_over_ext, E1, delay, kind)["meeting"]
        c = dtr_crossing(mu1, mu2, m["r"], D, alpha, kind, log10_gap=m["log10_rminus_minus_r"]) if m["log10_rminus_minus_r"] < -6 \
            else dtr_crossing(mu1, mu2, m["r"], D, alpha, kind)
        return c["f_F_required"] - c["f_F_min_model_U"], c

    lo, hi = 0.5, 1.0
    while margin(hi)[0] > 0:
        lo, hi = hi, hi * 2
        if hi > 1e4:
            return dict(found=False)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if margin(mid)[0] > 0 else (lo, mid)
    kap = temperature(mu1, D, alpha, kind)["kappa_minus"]
    before = margin(0.9 * lo)[1]
    return dict(found=True, D=D, mu_over_ext=mu_over_ext, mu2_over_mu1=mu2_over_mu1, delay_c=float(hi), kappa_minus=kap,
                mu_F_U_just_before=before["mu_F_model_U"], mu_F_GR_just_before=before["mu_F_GR"],
                verdict="beyond delay_c the crossing region would need psi > 1/alpha: it saturates at the limiting curvature (locally de Sitter)")
