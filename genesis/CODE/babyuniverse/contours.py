"""Phase III: the reduced Hamilton-Jacobi function on the complex R-plane.

On shell (fixed M) the canonical momentum is

    P(R) = (R/2G) [ log((1+x_out)/(1-x_out)) - log((1+x_in)/(1-x_in)) ],
    x = Rdot/A(R),  Rdot = sqrt(-V(R)),

an analytic function on a Riemann surface whose finite singular points are

* the turning points (zeros of V): square-root branch points;
* the horizons (f_in = 0 or f_out = 0, where Rdot^2 = A^2 and x = +-1):
  logarithmic branch points with monodromy log -> log + 2 pi i, i.e.
  P -> P +- i pi R / G;
* R = 0.

The flip radii (A = 0) are regular points: x -> infinity and
(1+x)/(1-x) -> -1 stays away from 0 and infinity.

`continue_P` integrates P dR along an arbitrary complex path with continuous
branch tracking (square root by nearest-value continuation, logarithms by
unwrapping their phase), so that horizon windings can be measured directly.
"""

from __future__ import annotations

import numpy as np

from .model import ThinWallModel


def V_complex(model: ThinWallModel, R):
    c2, e2, k2 = model.chi**2, model.eta**2, model.kappa**2
    GM = model.G * model.M
    return 1 - model.L4 / (4 * k2) * R**2 + GM * (c2 - e2 - k2) / (k2 * R) - GM**2 / (k2 * R**4)


def A_complex(model: ThinWallModel, R, side: str):
    fin = 1 - model.chi**2 * R**2
    fout = 1 - 2 * model.G * model.M / R - model.eta**2 * R**2
    s = 1 if side == "in" else -1
    return (fin - fout + s * (model.kappa * R) ** 2) / (2 * model.kappa * R)


def continue_P(model: ThinWallModel, path: np.ndarray, sqrt_sign: complex = None):
    """Analytically continue P along `path` (complex array).
    Returns (integral of P dR, P along the path)."""
    R = np.asarray(path, dtype=complex)
    s = np.sqrt(-V_complex(model, R))
    if sqrt_sign is not None and abs(s[0]) > 0:
        if abs(s[0] - sqrt_sign) > abs(-s[0] - sqrt_sign):
            s[0] = -s[0]
    for i in range(1, len(s)):
        if abs(s[i] - s[i - 1]) > abs(-s[i] - s[i - 1]):
            s[i] = -s[i]
    logs = []
    for side in ("out", "in"):
        x = s / A_complex(model, R, side)
        z = (1 + x) / (1 - x)
        logs.append(np.log(np.abs(z)) + 1j * np.unwrap(np.angle(z)))
    P = R / (2 * model.G) * (logs[0] - logs[1])
    integral = np.sum(0.5 * (P[1:] + P[:-1]) * np.diff(R))
    return integral, P


def _seg(a, b, n):
    t = np.linspace(0.0, 1.0, n)
    return a + (b - a) * t


def horizon_loop_path(model: ThinWallModel, r_h: float, turns: int = 1, n: int = 20000):
    """Path from just above R1 to just below R2 along the real barrier, with a
    loop of winding number `turns` around the horizon point r_h attached at the
    barrier midpoint through a tail in the upper half-plane (the tail passes
    over, never around, the turning points)."""
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]
    eps = 1e-7 * (R2 - R1)
    a, b, c = R1 + eps, R2 - eps, 0.5 * (R1 + R2)
    dist = min(abs(r_h - R1), abs(r_h - R2))
    rad = 0.45 * dist
    top = 1j * 0.3 * (R2 - R1)
    tail = np.concatenate([_seg(c, c + top, n // 5), _seg(c + top, r_h + top, n // 2),
                           _seg(r_h + top, r_h + 1j * rad, n // 5)])
    t = np.linspace(0, 2 * np.pi * turns, 2 * n)
    loop = r_h + rad * np.exp(1j * (np.pi / 2 + t))
    return np.concatenate([_seg(a + 0j, c + 0j, n // 2), tail, loop, tail[::-1], _seg(c + 0j, b + 0j, n // 2)])


def real_barrier_path(model: ThinWallModel, n: int = 40000):
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]
    eps = 1e-7 * (R2 - R1)
    return _seg(R1 + eps + 0j, R2 - eps + 0j, n)


def B_of_path(model: ThinWallModel, path) -> float:
    """WKB exponent 2 |Im int P dR| for a path between the turning points,
    oriented so that the real barrier contour gives B > 0."""
    I0, _ = continue_P(model, real_barrier_path(model))
    I, _ = continue_P(model, path)
    sgn = -np.sign(np.imag(I0))
    return float(-2 * sgn * np.imag(I) * -1)


def monodromy_shift(model: ThinWallModel, r_h: float, turns: int = 1) -> dict:
    """Change of B and of Re int P dR when a horizon loop is inserted, together
    with the prediction (A(R2) - A_h)/4G per winding."""
    I0, _ = continue_P(model, real_barrier_path(model))
    I, _ = continue_P(model, horizon_loop_path(model, r_h, turns))
    sgn = np.sign(np.imag(I0))
    dB = 2 * sgn * (np.imag(I) - np.imag(I0))
    R2 = model.turning_points()[-1]
    pred = np.pi * (R2**2 - r_h**2) / model.G
    return dict(r_h=r_h, turns=turns, delta_B=float(dB), prediction_per_turn=float(pred),
                delta_Re=float(np.real(I - I0)), B_real_contour=float(2 * abs(np.imag(I0))))


def principal_branch_constraint_violation(model: ThinWallModel, n: int = 400) -> dict:
    """Theorem III-2 check.  The analytic-branch (complex-angle) P_E satisfies the
    sign-free constraint cos(G P_E / R) = C on the whole barrier.  The
    principal-branch (real-angle, artanh) P_E satisfies instead the
    sign-dressed constraint cos(G P_E / R) = eps_in eps_out C, whose right-hand
    side jumps at the flip radius."""
    from . import tunnelling as T
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]
    Rs = np.linspace(R1, R2, n)[1:-1]
    fi, fo = model.f_in(Rs), model.f_out(Rs)
    C = (fi + fo - model.kappa**2 * Rs**2) / (2 * np.sqrt(fi * fo))
    ee = np.sign(model.A_in(Rs)) * np.sign(model.A_out(Rs))
    pa = T.euclidean_momentum(model, Rs, "analytic")
    pp = T.euclidean_momentum(model, Rs, "principal")
    ca = np.cos(model.G * pa / Rs)
    cp = np.cos(model.G * pp / Rs)
    return dict(analytic_max_violation=float(np.max(np.abs(ca - C))),
                principal_max_violation=float(np.max(np.abs(cp - C))),
                principal_dressed_violation=float(np.max(np.abs(cp - ee * C))),
                analytic_dressed_violation=float(np.max(np.abs(ca - ee * C))))


# ------------------------------------------------- zero-net horizon periods
def two_horizon_loop_path(model: ThinWallModel, rh1: float, t1: int, rh2: float, t2: int,
                          n: int = 20000):
    """Real barrier contour with two horizon loops (windings t1 around rh1 and
    t2 around rh2), both attached at the barrier midpoint."""
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]
    eps = 1e-7 * (R2 - R1)
    a, b, c = R1 + eps, R2 - eps, 0.5 * (R1 + R2)
    top = 1j * 0.3 * (R2 - R1)

    def block(rh, t):
        rad = 0.45 * min(abs(rh - R1), abs(rh - R2))
        tail = np.concatenate([_seg(c, c + top, n // 5), _seg(c + top, rh + top, n // 2),
                               _seg(rh + top, rh + 1j * rad, n // 5)])
        ang = np.linspace(0, 2 * np.pi * t, 2 * n)
        return np.concatenate([tail, rh + rad * np.exp(1j * (np.pi / 2 + ang)), tail[::-1]])

    return np.concatenate([_seg(a + 0j, c + 0j, n // 2), block(rh1, t1), block(rh2, t2),
                           _seg(c + 0j, b + 0j, n // 2)])


def zero_net_shift(model: ThinWallModel, rh1: float, rh2: float) -> dict:
    """Windings around two horizons whose log-monodromies cancel: the momentum
    returns to the analytic branch at R2 (no change of sector), but B changes
    by a period +-(A_h2 - A_h1)/4G of the differential P dR.  The relative
    orientation that cancels is model dependent; both are tried."""
    I0, P0 = continue_P(model, real_barrier_path(model))
    best = None
    for t2 in (-1, 1):
        I, P = continue_P(model, two_horizon_loop_path(model, rh1, 1, rh2, t2))
        dP = abs(P[-1] - P0[-1])
        if best is None or dP < best["endpoint_dP"]:
            best = dict(t1=1, t2=t2, endpoint_dP=float(dP),
                        delta_B=float(2 * np.sign(np.imag(I0)) * (np.imag(I) - np.imag(I0))))
    best["period"] = float(np.pi * abs(rh2**2 - rh1**2) / model.G)
    return best


# ---------------------------------------- true Hamiltonian and its branch locus
def true_hamiltonian(model: ThinWallModel, R, Pu, root: int):
    """Solve the sign-free constraint 2 sqrt(f_in) s cosh(G Pu) = f_in + s^2 - kappa^2 R^2
    for s = sqrt(f_out): s_pm = F c +- sqrt(F^2 (c^2-1) + kappa^2 R^2), and return
    H_pm = (R/2G)(1 - eta^2 R^2 - s_pm^2), the mass conjugate to the asymptotic clock."""
    R = np.asarray(R, dtype=complex)
    F = np.sqrt(1 - model.chi**2 * R**2)
    c = np.cosh(model.G * np.asarray(Pu, dtype=complex))
    s = F * c + root * np.sqrt(F**2 * (c**2 - 1) + (model.kappa * R) ** 2)
    return R / (2 * model.G) * (1 - model.eta**2 * R**2 - s**2)


def branch_locus_check(model: ThinWallModel) -> dict:
    """Theorem III-4: on the barrier the classical trajectory (analytic branch,
    G Pu = i theta) meets the branch locus kappa^2 R^2 = f_in sin^2 theta of the
    true Hamiltonian exactly at the exterior flip radius."""
    from . import tunnelling as T
    Rf = model.flip_radius_out()
    th = float(T.euclidean_momentum(model, Rf) * model.G / Rf)
    rad = (model.kappa * Rf) ** 2 - model.f_in(Rf) * np.sin(th) ** 2
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]
    Rs = np.linspace(R1, R2, 2001)[1:-1]
    ths = T.euclidean_momentum(model, Rs) * model.G / Rs
    rads = (model.kappa * Rs) ** 2 - model.f_in(Rs) * np.sin(ths) ** 2
    sign_changes = int(np.sum(np.diff(np.sign(rads)) != 0))
    Hm = true_hamiltonian(model, Rf, 1j * th / model.G, -1)
    Hp = true_hamiltonian(model, Rf, 1j * th / model.G, +1)
    return dict(R_flip=float(Rf), relative_radicand=float(rad / (model.kappa * Rf) ** 2),
                radicand_sign_changes_on_barrier=sign_changes,
                H_minus_minus_M=float(abs(Hm - model.M)), H_plus_minus_M=float(abs(Hp - model.M)))


def sector_hamiltonians(model: ThinWallModel, frac: float = 1e-3) -> dict:
    """Which sheet H_pm carries the Lorentzian states next to each turning point,
    in the real-angle (artanh) and complex-angle (Sorkin) conventions."""
    from .quantum import canonical_momentum
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]
    out = {}
    for lab, R in (("R1", R1 * (1 - frac)), ("R2", R2 * (1 + frac))):
        rdot = np.sqrt(max(-model.V(R), 0.0))
        P_art = float(np.real(canonical_momentum(model, R, rdot)))
        n = (1 if model.A_out(R) < 0 else 0) - (1 if model.A_in(R) < 0 else 0)
        Pu_art = P_art / R
        Pu_sor = Pu_art + 1j * np.pi * n / model.G
        res = {}
        for conv, Pu in (("artanh", Pu_art), ("sorkin", Pu_sor)):
            res[conv] = {r: float(abs(true_hamiltonian(model, R, Pu, r) - model.M)) for r in (-1, 1)}
        out[lab] = dict(R=float(R), n=n, eps=(int(np.sign(model.A_in(R))), int(np.sign(model.A_out(R)))),
                        **res)
    return out


# ----------------------------------------- sector switching (squared constraint)
def flip_class(model: ThinWallModel) -> str:
    """'single': eps_out flips on the barrier, eps_in does not (sector (+,-) at R2);
    'double': both flip ((-,-) at R2); 'none' otherwise."""
    R2 = model.turning_points()[-1]
    eo, ei = model.A_out(R2) < 0, model.A_in(R2) < 0
    return "double" if (eo and ei) else ("single" if eo else "none")


def B_sector_switch(model: ThinWallModel) -> float:
    """Exponent of the sector-switching connection admitted by the sector-blind
    (squared) constraint cos^2 theta = C^2: theta follows arccos C up to the
    coalescence point C = 0 (theta = pi/2, where the branches theta and pi - theta
    meet) and pi - arccos C beyond.  Defined for the single-flip class."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    if flip_class(model) != "single":
        return float("nan")
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]

    def C(R):
        fi, fo = model.f_in(R), model.f_out(R)
        return (fi + fo - model.kappa**2 * R**2) / (2 * np.sqrt(fi * fo))

    Rq = brentq(C, R1 * (1 + 1e-12), R2 * (1 - 1e-12), xtol=1e-15)
    th = lambda R: np.arccos(np.clip(C(R), -1.0, 1.0))
    a = quad(lambda R: R * th(R), R1, Rq, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    b = quad(lambda R: R * (np.pi - th(R)), Rq, R2, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    return 2 * (a + b) / model.G


# --------------------------------------- sector-blind quantization (Prop. III-8)
def coalescence_numerator(model: ThinWallModel, R):
    """N(R) = f_in + f_out - kappa^2 R^2 = 2 - (chi^2 + eta^2 + kappa^2) R^2 - 2GM/R,
    strictly concave on R > 0 (Lemma III-9); C = N / (2 sqrt(f_in f_out))."""
    R = np.asarray(R, dtype=float)
    return 2 - (model.chi**2 + model.eta**2 + model.kappa**2) * R**2 - 2 * model.G * model.M / R


def coalescence_points(model: ThinWallModel, n: int = 4001) -> list:
    """Zeros of C (theta = pi/2) inside the barrier: points where the squared
    constraint's branches theta and pi - theta coalesce."""
    from scipy.optimize import brentq
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]
    Rs = np.linspace(R1, R2, n)[1:-1]
    N = coalescence_numerator(model, Rs)
    idx = np.nonzero(np.diff(np.sign(N)) != 0)[0]
    return [brentq(lambda r: float(coalescence_numerator(model, r)), Rs[i], Rs[i + 1], xtol=1e-15) for i in idx]


def B_sector_blind(model: ThinWallModel) -> float:
    """Exponent of the sector-blind (squared-constraint) quantization:
    B_sq = (2/G) int R arccos|C| dR  (min(theta, pi - theta) on the barrier).
    Equals B_sector_switch in the single-flip class and B_an in the double-flip class."""
    from scipy.integrate import quad
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]

    def C(R):
        fi, fo = model.f_in(R), model.f_out(R)
        return (fi + fo - model.kappa**2 * R**2) / (2 * np.sqrt(fi * fo))

    pts = coalescence_points(model) or None
    v = quad(lambda R: R * np.arccos(min(1.0, abs(C(R)))), R1, R2, points=pts, limit=400,
             epsabs=1e-13, epsrel=1e-12)[0]
    return 2 * v / model.G


def weber_multiplier(nu: float, z: float = 30.0) -> dict:
    """Stokes multiplier of the local normal form at a coalescence point.
    The solution recessive on the left, D_nu(-z), contains the dominant branch
    e^{z^2/4} z^{-nu-1} on the right with coefficient sqrt(2 pi)/Gamma(-nu)
    (up to the phase e^{i pi nu}); it vanishes iff nu = 0, 1, 2, ...  Checked
    against mpmath's parabolic-cylinder function at large z."""
    import mpmath as mp
    mp.mp.dps = 50
    nu_m, zz = mp.mpf(nu), mp.mpf(z)
    measured = mp.pcfd(nu_m, -zz) / (mp.e ** (zz**2 / 4) * zz ** (-nu_m - 1))
    predicted = mp.mpf(0) if (mp.isint(nu_m) and nu_m >= 0) else mp.sqrt(2 * mp.pi) / mp.gamma(-nu_m)
    return dict(nu=float(nu), measured_abs=float(abs(measured)), predicted_abs=float(abs(predicted)))
