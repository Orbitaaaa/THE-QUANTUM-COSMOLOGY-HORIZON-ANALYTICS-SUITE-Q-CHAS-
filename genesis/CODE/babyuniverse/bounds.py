"""Roadmap Stage 4 (task 30) and Stage 1 (P-T5): analytic suppression bounds.

Everything here is for laboratory transitions: Schwarzschild exterior (eta = 0), regime A.  Notation as in
model.py; gamma, E are the dimensionless parameters, E_AB = 2 - 2 gamma is the A/B boundary.

* ``sine_bound``      Lemma S1:  B_an >= B_sq >= B_mem = (2 kappa/G) int R^2 sqrt(V/(f_in f_out)) dR.
* ``elementary_bound`` Theorem S2: B >= 2 pi kappa^2 R1^2 f_out(R1) / (G (L^2 + 4GM/R1^3) sqrt(f_in(R1) f_out(R2))).
* ``Phi``             the universal thin-barrier function of Theorem S3, B = B(0) Phi(E) (1 + O(kappa^2/chi^2)).
* ``M_sd``            Theorem S4: the transition is double-flip (FGG-unambiguous) iff kappa < chi and M < M_sd.
* ``omega_condition`` Theorem S5: a pointwise sufficient condition for Delta t_E^W < 0 (hence dB/dM > 0 in the
                      double-flip class, B >= B_CdL).
* ``boundary_limit``  Proposition S5-b: lim_{E -> E_AB+} Delta t_E^W < 0 for every gamma (the A/B end of the strip).
* ``laboratory_law``  Corollary S6: B_lab = 4 pi^2 sigma^2 R_M^2 / rho  (R_M >> R_CdL, kappa << chi).

See ROADMAP_STAGE1_2_4.md, Part A.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from .model import ThinWallModel
from . import tunnelling as T


# ------------------------------------------------------------------ Lemma S1 / Theorem S2
def sine_bound(model: ThinWallModel) -> float:
    """B_mem = (2/G) int R sin(theta) dR with sin(theta) = kappa R sqrt(V/(f_in f_out)).
    theta >= sin(theta) and arccos|C| >= sqrt(1 - C^2) = sin(theta) give B_an >= B_sq >= B_mem."""
    Ra, Rb = T.barrier(model)

    def f(R):
        fi, fo = float(model.f_in(R)), float(model.f_out(R))
        V = max(float(model.V(R)), 0.0)
        return R * R * np.sqrt(V / (fi * fo))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        v = quad(f, Ra, Rb, limit=600, epsabs=0.0, epsrel=1e-11)[0]
    return 2 * model.kappa * v / model.G


def elementary_bound(model: ThinWallModel) -> float:
    """Theorem S2 (eta = 0, regime A).  On [R1, R'] with A_out(R') = -sqrt(f_out(R1)):
    R^2 >= R1^2, V >= f_out(R1) - A_out^2, f_in f_out <= f_in(R1) f_out(R2), |A_out'| <= (L^2 + 4GM/R1^3)/(2 kappa);
    the w = A_out integral of sqrt(f1 - w^2) is pi f1 / 2."""
    if model.eta != 0:
        raise ValueError("elementary_bound needs a Schwarzschild exterior")
    R1, R2 = T.barrier(model)
    GM = model.G * model.M
    return float(2 * np.pi * model.kappa**2 * R1**2 * model.f_out(R1)
                 / (model.G * (model.L2 + 4 * GM / R1**3) * np.sqrt(model.f_in(R1) * model.f_out(R2))))


# ------------------------------------------------------------------ Theorem S3: thin barrier
def _g1(z):
    return (z**3 - 1) ** 2 / z**4


def Phi(E: float) -> float:
    """Phi(E) = (16 / (pi E^2)) int z^2 sqrt(E - (z^3-1)^2/z^4) dz over {g_1 < E}.
    Computed in the variable s = (z^3-1)/z^2 (monotone, g_1 = s^2), z^2 dz = z^5/(z^3+2) ds,
    with s = sqrt(E) sin(phi): no endpoint singularities."""
    rE = np.sqrt(E)

    def z_of_s(s):
        # z^3 - s z^2 - 1 = 0 has exactly one positive root
        lo, hi = 1e-12, max(2.0, abs(s) + 2.0)
        return brentq(lambda z: z**3 - s * z**2 - 1, lo, hi, xtol=1e-15, rtol=1e-15, maxiter=200)

    def integrand(phi):
        s = rE * np.sin(phi)
        z = z_of_s(s)
        return z**5 / (z**3 + 2) * E * np.cos(phi) ** 2

    v = quad(integrand, -np.pi / 2, np.pi / 2, limit=400, epsabs=0.0, epsrel=1e-12)[0]
    return 16 / (np.pi * E**2) * v


def Phi_small_E(E):
    """Phi(E) = 8/(3E) + 2/3 + O(E)   (h(s) = z^5/(z^3+2) = (1 + 4s/3 + s^2)/3 + O(s^3))."""
    return 8 / (3 * np.asarray(E)) + 2 / 3


def Phi_large_E(E):
    """Phi(E) = 1 + (16/pi) E^(-3/2) + o(E^(-3/2))   (first-order effect of the -2/z term of g)."""
    return 1 + 16 / np.pi * np.asarray(E) ** -1.5


def thin_barrier_B(model: ThinWallModel) -> float:
    """Leading thin-barrier prediction B(0) Phi(E), with B(0) the exact closed form T1."""
    B0 = T.B_canonical(ThinWallModel(model.chi, model.kappa, 0.0, model.eta, model.G))
    return B0 * Phi(model.E)


def R_flip(model: ThinWallModel) -> float:
    return (2 * model.G * model.M / model.L2) ** (1 / 3)


def thin_barrier_asymptote(model: ThinWallModel) -> float:
    """E -> 0 end of Theorem S3: B -> (2 pi / 3) (kappa/chi)^2 R_f^2 / G  = (2/3)(kappa/chi)^2 A_f / 4G."""
    return 2 * np.pi / 3 * (model.kappa / model.chi) ** 2 * R_flip(model) ** 2 / model.G


# ------------------------------------------------------------------ Theorem S4: flip class
def M_sd(model: ThinWallModel) -> float:
    """Single/double-flip boundary: R2 = 1/chi, i.e. 2 G M chi = 1 - kappa^2/chi^2 (eta = 0).
    Double flip iff kappa < chi and M < M_sd; returns 0 if kappa >= chi."""
    if model.chi <= model.kappa:
        return 0.0
    return (model.chi**2 - model.kappa**2) / (2 * model.G * model.chi**3)


def E_sd(gamma: float) -> float:
    """Dimensionless form of M_sd: E_sd = (1 - gamma^2) / gamma^(2/3) (gamma > 0)."""
    return (1 - gamma**2) / gamma ** (2 / 3) if gamma > 0 else np.inf


def flip_class_predicted(model: ThinWallModel) -> str:
    return "double" if model.M < M_sd(model) else "single"


# ------------------------------------------------------------------ Theorem S5: sign of Delta t_E
def _dimless_barrier(gamma: float, E: float):
    g = lambda z: z**2 - 2 * gamma / z + 1 / z**4
    zs = ThinWallModel.z_star(gamma)
    z1 = brentq(lambda z: g(z) - E, 1e-10, zs, xtol=1e-15)
    z2 = brentq(lambda z: g(z) - E, zs, 10 * np.sqrt(E) + 10, xtol=1e-15)
    return z1, z2


def omega(z):
    """omega(z) = z^2/12 + 7/(6 z): Delta t_E^W < 0 whenever E / E_AB >= max(omega(z1), omega(z2))."""
    return z**2 / 12 + 7 / (6 * z)


def omega_condition(gamma: float, E: float) -> bool:
    z1, z2 = _dimless_barrier(gamma, E)
    return E / (2 - 2 * gamma) >= max(omega(z1), omega(z2))


def pointwise_D_min(gamma: float, E: float, n: int = 4001) -> float:
    """min over the barrier of D = 6z(W - s^2)W - b(z^3+2)(W/2 - s^2) (b = E_AB, W = E - b/z, s = (z^3-1)/z^2).
    D >= 0 on the barrier implies Delta t_E^W < 0 (proof of Theorem S5)."""
    z1, z2 = _dimless_barrier(gamma, E)
    b = 2 - 2 * gamma
    z = np.linspace(z1, z2, n)[1:-1]
    W = E - b / z
    s = (z**3 - 1) / z**2
    return float((6 * z * (W - s**2) * W - b * (z**3 + 2) * (W / 2 - s**2)).min())


def delta_tE_dimless(gamma: float, E: float) -> float:
    """J(gamma, E) = int (1 - z^3) / (z^2 f_out sqrt(E - g)) dz, proportional (positive factor) to Delta t_E^W.

    Near E_AB the outer horizon pinches the turning point z1 (W(z1) = E f_out(z1) = s(z1)^2 -> 0) and the integrand has a
    spike of width w = W(z1) z1^2 / b at z1 (it carries +pi in regime A).  The integrand is therefore written in y = z - z1
    without cancellations, with E taken as g(z1) (an O(1e-16) shift):
        E - g = -y [2 z1 + y + 2 gamma/(z1 z) - (z1 + z)(z1^2 + z^2)/(z1^4 z^4)],   W = s(z1)^2 + b y/(z1 z),
    and the range is broken geometrically from z1 on the scale w, and at s = 0; cos-substitution on each piece."""
    z1, z2 = _dimless_barrier(gamma, E)
    b = 2 - 2 * gamma
    W1 = ((z1**3 - 1) / z1**2) ** 2

    def integrand(y):
        z = z1 + y
        r = -y * (2 * z1 + y + 2 * gamma / (z1 * z) - (z1 + z) * (z1**2 + z**2) / (z1**4 * z**4))
        if r <= 0:           # rounding at the far turning point; the true contribution there is O(sqrt(r)) -> 0
            return 0.0
        W = W1 + b * y / (z1 * z)
        return (1 - z**3) * E / (z**2 * W * np.sqrt(r))

    def piece(a, c):
        return quad(lambda t: integrand(a + (c - a) * (1 - np.cos(t)) / 2) * (c - a) * np.sin(t) / 2, 0, np.pi, limit=800)[0]

    L = z2 - z1
    w = max(W1 * z1**2 / b, 1e-300)
    pts = {0.0, L, *((1.0 - z1,) if z1 < 1 < z2 else ())}
    p = w
    while p < 0.5 * L:
        pts.add(p)
        p *= 10
    pts = sorted(pts)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sum(piece(a, c) for a, c in zip(pts[:-1], pts[1:]))


# ------------------------------------------------------------------ Proposition S5-b: the A/B end of the strip
def _h_ab(z):
    """h(z) = (z^3 - 1)(z^2 + z + 1)/z^3, so that E_AB - g(z) = (z - 1)(E_AB - h(z))/z."""
    return z**2 + z + 1 - 1 / z - 1 / z**2 - 1 / z**3


def _F_ab(z):
    """F(z) = 2 (z^2 + z + 1)^2 / (z^2 h'(z)); F(1) = 2 and F > 2 for z > 1 (boundary_certificate)."""
    hp = 2 * z + 1 + z**-2 + 2 * z**-3 + 3 * z**-4
    return 2 * (z * z + z + 1) ** 2 / (z * z * hp)


def boundary_J0(gamma: float) -> float:
    """J0(gamma) = J evaluated at E = E_AB (barrier [1, z2]):
        J0 = -int_1^{z2} (z^2 + z + 1) dz / (z sqrt(E_AB - g)) = -int_0^{pi/2} F(z(phi)) dphi,   h(z(phi)) = E_AB sin^2 phi."""
    b = 2 - 2 * gamma

    def z_of(ph):
        t = b * np.sin(ph) ** 2
        return 1.0 if t <= 0 else brentq(lambda z: _h_ab(z) - t, 1.0, 10 + 2 * np.sqrt(b), xtol=1e-15)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return -quad(lambda ph: _F_ab(z_of(ph)), 0, np.pi / 2, limit=200)[0]


def boundary_limit(gamma: float) -> float:
    """lim_{E -> E_AB+} J(gamma, E) = J0(gamma) + pi  (the pinch of the outer horizon with s = 0 at z = 1 contributes +pi).
    Negative for every gamma in (-1, 1): |J0| = int F dphi > (pi/2) * 2 because F > 2 on z > 1."""
    return boundary_J0(gamma) + np.pi


def boundary_certificate() -> dict:
    """(z^2 + z + 1)^2 - z^2 h'(z) = (z - 1)(z + 1)(z^2 - z + 3)(z^2 + z + 1)/z^2 > 0 for z > 1  =>  F(z) > 2."""
    import sympy as sp
    z = sp.symbols("z", positive=True)
    h = z**2 + z + 1 - 1 / z - 1 / z**2 - 1 / z**3
    lhs = (z**2 + z + 1) ** 2 - z**2 * sp.diff(h, z)
    rhs = (z - 1) * (z + 1) * (z**2 - z + 3) * (z**2 + z + 1) / z**2
    h_ok = sp.simplify(h - (z**3 - 1) * (z**2 + z + 1) / z**3) == 0
    g = z**2 - 2 * sp.Symbol("gamma") / z + 1 / z**4
    split_ok = sp.simplify((2 - 2 * sp.Symbol("gamma")) - g - (z - 1) * ((2 - 2 * sp.Symbol("gamma")) - h) / z) == 0
    return dict(factorization=bool(sp.simplify(lhs - rhs) == 0), h_definition=bool(h_ok), barrier_split=bool(split_ok))


def strip_lp_certificate(gamma: float, E: float, deg: int = 5, n: int = 800, verify_n: int = 20000) -> dict:
    """Generalised integration by parts for -J/E = int s k (W - s^2)^(-1/2) ds, k = q/W: for any smooth phi,
    int [s k + phi'(W - s^2) + phi (W'/2 - s)] (W - s^2)^(-1/2) ds = int s k (W - s^2)^(-1/2) ds (it adds d(phi sqrt(W - s^2))/ds).
    phi = sum_j (a_j k + b_j q + c_j q sqrt(W)) x^j, x the rescaled s.  A linear programme maximises the minimum of the
    modified integrand on n Chebyshev-clustered points (relative to max|s k|); the solution is then re-checked on verify_n
    points.  Positive margins certify J < 0 on the samples (sampling-based, not interval-verified).
    Theorem S5 is the single choice phi = k."""
    from scipy.optimize import linprog
    b = 2 - 2 * gamma
    z1, z2 = _dimless_barrier(gamma, E)
    s1, s2 = (z1**3 - 1) / z1**2, (z2**3 - 1) / z2**2
    sm, hw = (s1 + s2) / 2, (s2 - s1) / 2

    def columns(t):
        S = s1 + (s2 - s1) * (1 - np.cos(t)) / 2
        Z = np.array([brentq(lambda z: z**3 - si * z**2 - 1, 1e-12, abs(si) + 3, xtol=1e-15) for si in S])
        W = E - b / Z
        dzds = Z**3 / (Z**3 + 2)
        q, dq = Z**3 / (Z**3 + 2), 6 * Z**2 / (Z**3 + 2) ** 2 * dzds
        dW = b / Z**2 * dzds
        k, dk = q / W, dq / W - q * dW / W**2
        r, dr = q * np.sqrt(W), dq * np.sqrt(W) + q * dW / (2 * np.sqrt(W))
        u = W - S**2
        x = (S - sm) / hw
        cols = []
        for base, dbase in ((k, dk), (q, dq), (r, dr)):
            for j in range(deg + 1):
                m = x**j
                dm = j * x ** (j - 1) / hw if j > 0 else 0 * x
                cols.append((dbase * m + base * dm) * u + base * m * (dW / 2 - S))
        return S * k, np.array(cols).T

    sk, A = columns(np.linspace(0, np.pi, n)[1:-1])
    scale = np.abs(sk).max()
    nv = A.shape[1]
    c = np.zeros(nv + 1)
    c[-1] = -1
    res = linprog(c, A_ub=np.hstack([-A, scale * np.ones((len(sk), 1))]), b_ub=sk,
                  bounds=[(-50, 50)] * nv + [(None, 1)], method="highs")
    if res.status != 0:
        return dict(margin=float("nan"), verified_margin=float("nan"))
    sk2, A2 = columns(np.linspace(0, np.pi, verify_n)[1:-1])
    return dict(margin=float(-res.fun), verified_margin=float((sk2 + A2 @ res.x[:-1]).min() / scale))


# ------------------------------------------------------------------ Corollary S6: laboratory law
def laboratory_law(E_vac_GeV: float, E_wall_GeV: float, M_planck: float) -> dict:
    """B_lab = 4 pi^2 sigma^2 R_M^2 / rho with R_M = (3M / 4 pi rho)^(1/3)  (Planck units),
    and the equivalent form 4 pi sqrt(2/27) sqrt(B_flat) (3M/4pi)^(2/3) rho^(-1/6) with
    B_flat = 27 pi^2 sigma^4 / (2 rho^3)."""
    from .observability import planck_params
    p = planck_params(E_vac_GeV, E_wall_GeV)
    rho, sigma = p["rho"], p["sigma"]
    RM = (3 * M_planck / (4 * np.pi * rho)) ** (1 / 3)
    R0 = 3 * sigma / rho
    B_flat = 27 * np.pi**2 * sigma**4 / (2 * rho**3)
    B_law = 4 * np.pi**2 * sigma**2 * RM**2 / rho
    B_phi = B_flat * Phi((R0 / RM) ** 2)
    return dict(R_M_lP=float(RM), R_CdL_lP=float(R0), x=float(RM / R0), B_flat=float(B_flat),
                B_law=float(B_law), B_Phi=float(B_phi), kappa_over_chi=float(p["k"]),
                B_law_identity=float(4 * np.pi * np.sqrt(2 / 27) * np.sqrt(B_flat)
                                     * (3 * M_planck / (4 * np.pi)) ** (2 / 3) * rho ** (-1 / 6)))


def minimum_exponent(B_survive: float, M_planck: float, rho_planck: float) -> float:
    """Smallest thin-barrier laboratory exponent compatible with metastability B_CdL >= B_survive,
    a black hole of at least M_planck Planck masses and vacuum energy density at most rho_planck:
    B_min = 4 pi sqrt(2 B_survive / 27) (3 M / 4 pi)^(2/3) rho^(-1/6)."""
    return float(4 * np.pi * np.sqrt(2 * B_survive / 27) * (3 * M_planck / (4 * np.pi)) ** (2 / 3)
                 * rho_planck ** (-1 / 6))


# ------------------------------------------------------------------ Corollary S7: the semiclassical frontier
def h_even_defect(s: float) -> float:
    """(h(s) + h(-s))/2 - 1/3 - s^2/3 with h(s) = z^5/(z^3+2), z^3 - s z^2 - 1 = 0 (numerically non-negative; ~0.48 s^4 at 0).
    The theorem Phi(E) >= 8/(3E) + 2/3 needs only the integrated form, proved in ``phi_bound_certificate``."""
    def z(x):
        return brentq(lambda zz: zz**3 - x * zz**2 - 1, 1e-12, abs(x) + 3, xtol=1e-15)

    def h(x):
        zz = z(x)
        return zz**5 / (zz**3 + 2)
    return float(0.5 * (h(s) + h(-s)) - 1 / 3 - s * s / 3)


FRONTIER_CONST = 4 * np.pi * np.sqrt(2 / 27) * np.sqrt(3 / (4 * np.pi))     # = (4/3) sqrt(pi/2) = 1.6711


def frontier_bound(B_CdL: float, M_planck: float, R_lP: float) -> float:
    """B_lab >= (4/3) sqrt(pi/2) sqrt(B_CdL * M/m_P * R_M/l_P) for every thin-barrier laboratory transition
    (from B_lab >= 4 pi^2 sigma^2 R_M^2 / rho, B_CdL = 27 pi^2 sigma^4 / (2 rho^3), M = (4 pi/3) rho R_M^3)."""
    return float(FRONTIER_CONST * np.sqrt(B_CdL * M_planck * R_lP))


def frontier_optimum(B_survive: float, M_min: float, R_min: float) -> dict:
    """Vacuum scale at which the frontier bound is attained: rho* = 3 M_min / (4 pi R_min^3)."""
    rho = 3 * M_min / (4 * np.pi * R_min**3)
    B = frontier_bound(B_survive, M_min, R_min)
    sigma = (2 * B_survive * rho**3 / (27 * np.pi**2)) ** 0.25
    return dict(B_min=B, log10_P_max=float(-B / np.log(10)), E_vac_GeV=float(rho**0.25 * 1.22e19),
                E_wall_GeV=float(sigma ** (1 / 3) * 1.22e19), R_CdL_lP=float(3 * sigma / rho),
                kappa_over_chi=float(np.sqrt(6 * np.pi * sigma**2 / rho)))


def phi_bound_certificate() -> dict:
    """Proof of Phi(E) >= max(1, 8/(3E) + 2/3) (ROADMAP A3).  With x = z(s) > 1, y = z(-s) < 1, F = z^3 and
    G(s) = F(s) - F(-s) = s (x^2 + y^2), integration by parts gives
        int (h - 1/3 - s^2/3) sqrt(E - s^2) ds = (1/3) int_0^sqrt(E) K(s) s / sqrt(E - s^2) ds,  K = G - 2s - 2s^3/3,
    so it suffices that Q = x^2 + y^2 >= 2 + 2 s^2/3.  From x - 1/x^2 = s = 1/y^2 - y one gets, with P = xy,
    Q^2 = P^4 Q + 2 P^5 (P >= 1 forced by Q >= 2P) and s^2 = (Q - 2P)(Q + P)^2 / Q^2, so the claim is
    f = Q^3 - 6 Q^2 + 6 P^2 Q + 4 P^3 = c1 Q + c0 >= 0.  The squared form reduces to
    2 c1^2 P^5 - c0^2 - c0 c1 P^4 = 8 P^6 (P-1)^2 (2P^7 + 7P^6 + 12P^5 + 20P^4 + 22P^3 + 6P^2 - 4P - 2) >= 0 for P >= 1.
    Phi > 1: the same argument with K2 = G - s^3 = s(x^2 + y^2 - s^2) > 0 since 0 < s < x.  Returns the sympy checks."""
    import sympy as sp
    P = sp.symbols("P", positive=True)
    c1 = P**8 + 2 * P**5 - 6 * P**4 + 6 * P**2
    c0 = 2 * P**9 - 12 * P**5 + 4 * P**3
    key = sp.expand(2 * c1**2 * P**5 - c0**2 - c0 * c1 * P**4)
    cof = 2 * P**7 + 7 * P**6 + 12 * P**5 + 20 * P**4 + 22 * P**3 + 6 * P**2 - 4 * P - 2
    identity = sp.simplify(key - 8 * P**6 * (P - 1) ** 2 * cof) == 0
    cof_roots = [float(r) for r in sp.Poly(cof, P).real_roots()]
    c1_roots = [float(r) for r in sp.Poly(c1, P).real_roots() if r >= 1]
    # f = c1 Q + c0 with c1 Q from the root of Q^2 - P^4 Q - 2P^5 = 0: check f(Q) and the relation symbolically
    Q = sp.symbols("Q", positive=True)
    f_reduced = sp.expand(sp.rem(sp.expand(Q**3 - 6 * Q**2 + 6 * P**2 * Q + 4 * P**3), Q**2 - P**4 * Q - 2 * P**5, Q))
    return dict(identity=bool(identity), cofactor_real_roots=cof_roots, cofactor_positive_on_P_ge_1=bool(max(cof_roots) < 1),
                c1_roots_ge_1=c1_roots, f_reduced_matches=bool(sp.simplify(f_reduced - (c1 * Q + c0)) == 0))
