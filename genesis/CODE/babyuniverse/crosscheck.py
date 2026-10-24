"""Stage 1 task 3: an independent check of the universal suppression law (Theorem S3) and of its parameter range.

Nothing here imports model.py, tunnelling.py or bounds.py.  Everything is rebuilt from the junction data in mpmath:

    f_in = 1 - chi^2 R^2,  f_out = 1 - 2GM/R,  A_out = (f_in - f_out - kappa^2 R^2) / (2 kappa R),  V = f_out - A_out^2,
    cos(theta) = (f_in + f_out - kappa^2 R^2) / (2 sqrt(f_in f_out)),  sin(theta) = kappa R sqrt(V / (f_in f_out)),
    B_an = (2/G) int_{V>0} R theta dR   (theta = atan2(sin, cos), the analytic determination on the barrier),
    B_sq = (2/G) int_{V>0} R arccos|cos theta| dR,

and Phi(E) = (16 / (pi E^2)) int_{g1<E} z^2 sqrt(E - g1) dz with g1 = (z^3 - 1)^2 / z^4, directly in z.  The mass is fixed by
E = 4 kappa^2 / (L^(8/3) (2GM)^(2/3)), L^2 = chi^2 + kappa^2 (eta = 0).  B(0) is the same quadrature at M = 0.

Round 6c (P-S3, Theorem 6c-S3): the correction coefficient in B_an = B(0) Phi(E) [1 + c(E) kappa^2/chi^2 + O(kappa^4/chi^4)] in closed
form.  In units R = kappa rho, 2GM = kappa^3 m (m = 4 E^(-3/2)) the barrier is V = 1 - X^2 + kappa^2 V_1 with X = m/rho^2 - rho/2 and
V_1 = -P, P = rho^2/2 + m/rho; expanding the junction angle and B(0) = pi kappa^4 (1 - 2 kappa^2) / (8G) to O(kappa^2) gives

    c(E) = 2 - 2 m I_0'(m)/I_0 + I_1/I_0,   I_0 = int rho^2 sqrt(1 - X^2) d rho,   I_1 = int rho^2 P (sqrt(1 - X^2) - 1/(2 sqrt(1 - X^2))) d rho,

over the barrier |X| < 1, with I_0' = -int X / sqrt(1 - X^2) d rho.  As E -> 0 (m -> infinity) the barrier is the thin shell
rho = rho_*(1 + y), (1 + y)^-2 - (1 + y) = 2 X / rho_*, and the exact asymptotics I_0 = pi rho_*^2/3 + pi/3, I_1 = -25 pi rho_*^2/54 - 2 pi/3
give c(0+) = 2/3 - 25/18 = -13/18 (the Round 3 conjecture).  As E -> infinity (m -> 0): I_0 = pi/2 + 2m, I_0' -> 2, I_1 + 2 I_0 = 2m
(the terms 20m/3, -2m/3 and -4m from the shift of X, the m/rho part of P and the finite part of the 1/sqrt(1 - X^2) integral), so
c = -8m/pi + 4m/pi = -(16/pi) E^(-3/2) + O(E^(-9/4)); the Round 3 fit -5 E^(-3/2) was this law with 16/pi = 5.093.
"""

from __future__ import annotations

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def _junction(chi, kappa, M, G=1):
    f_in = lambda R: 1 - chi**2 * R**2
    f_out = lambda R: 1 - 2 * G * M / R
    A_out = lambda R: (f_in(R) - f_out(R) - kappa**2 * R**2) / (2 * kappa * R)
    V = lambda R: f_out(R) - A_out(R) ** 2
    cos_t = lambda R: (f_in(R) + f_out(R) - kappa**2 * R**2) / (2 * mp.sqrt(f_in(R) * f_out(R)))
    sin_t = lambda R: kappa * R * mp.sqrt(max(V(R), 0) / (f_in(R) * f_out(R)))
    return f_in, f_out, V, cos_t, sin_t


def _barrier(V, lo, hi, n=1500):
    """Endpoints of {V > 0} inside (lo, hi): logarithmic scan, then bisect each sign change."""
    xs = [lo * (hi / lo) ** (mp.mpf(k) / n) for k in range(1, n)]
    vs = [V(x) for x in xs]
    roots = [mp.findroot(V, (xs[k], xs[k + 1]), solver="bisect") for k in range(len(xs) - 1) if vs[k] * vs[k + 1] < 0]
    return roots


def exponents(chi, kappa, E, dps=30, G=1) -> dict:
    """B_an, B_sq for the laboratory transition with dimensionless energy E, and B(0)."""
    mp.mp.dps = dps
    chi, kappa, E = mp.mpf(chi), mp.mpf(kappa), mp.mpf(E)
    L2 = chi**2 + kappa**2
    M = (4 * kappa**2 / (E * L2 ** (mp.mpf(4) / 3))) ** (mp.mpf(3) / 2) / (2 * G)
    out = {}
    for label, m in (("M", M), ("zero", mp.mpf(0))):
        f_in, f_out, V, c, s = _junction(chi, kappa, m, G)
        lo = 2 * G * m * (1 + mp.mpf(10) ** (-dps // 2)) if m > 0 else mp.mpf(10) ** (-dps // 2)
        r = _barrier(V, lo, (1 - mp.mpf(10) ** (-dps // 2)) / chi)
        if m > 0:
            if len(r) != 2:
                raise ValueError("expected a two-sided barrier (regime A, M < M_cr)")
            Ra, Rb = r
        else:
            Ra, Rb = mp.mpf(0), r[-1]
        an = 2 / G * mp.quad(lambda R: R * mp.atan2(s(R), c(R)), [Ra, Rb])
        sq = 2 / G * mp.quad(lambda R: R * mp.acos(min(abs(c(R)), 1)), [Ra, Rb])
        out[label] = (an, sq, Ra, Rb)
    return dict(M=M, B_an=out["M"][0], B_sq=out["M"][1], B0=out["zero"][0], R1=out["M"][2], R2=out["M"][3])


def Phi_direct(E, dps=30):
    mp.mp.dps = dps
    E = mp.mpf(E)
    g1 = lambda z: (z**3 - 1) ** 2 / z**4
    za = mp.findroot(lambda z: g1(z) - E, (mp.mpf("1e-6"), mp.mpf(1)), solver="bisect")
    zb = mp.findroot(lambda z: g1(z) - E, (mp.mpf(1), 10 + 2 * mp.sqrt(E)), solver="bisect")
    return 16 / (mp.pi * E**2) * mp.quad(lambda z: z**2 * mp.sqrt(max(E - g1(z), 0)), [za, 1, zb])


def s3_check(eps, E, dps=30) -> dict:
    """Relative deviation of B_an from B(0) Phi(E) at kappa/chi = eps, and the same deviation divided by eps^2."""
    ex = exponents(1, eps, E, dps)
    ph = Phi_direct(E, dps)
    dev = ex["B_an"] / (ex["B0"] * ph) - 1
    gamma = (1 - eps**2) / (1 + eps**2)
    return dict(eps=float(eps), E=float(E), E_over_EAB=float(E / (2 - 2 * gamma)), B_an=float(ex["B_an"]), B_sq=float(ex["B_sq"]),
                B0=float(ex["B0"]), Phi=float(ph), rel_dev=float(dev), rel_dev_over_eps2=float(dev / eps**2),
                double_flip=bool(abs(ex["B_an"] - ex["B_sq"]) < mp.mpf(10) ** (-dps // 2) * ex["B_an"]))


def s3_coefficient(E: float) -> float:
    """Theorem 6c-S3: c(E) from the closed form.  1 - X^2 = (r - lo)(hi - r) g(r) with the divided differences
    (1 - X)/(r - lo) = m (r + lo)/(lo^2 r^2) + 1/2 and (1 + X)/(hi - r) = m (r + hi)/(r^2 hi^2) + 1/2, so the square-root endpoints of the
    barrier are carried exactly by algebraic quadrature weights."""
    m = 4.0 * E ** -1.5
    X = lambda r: m / r**2 - r / 2
    r_star = (2 * m) ** (1 / 3)                                   # X(r_star) = 0, X decreasing
    lo = brentq(lambda r: X(r) - 1, 1e-14 * r_star, r_star, xtol=1e-16 * r_star, rtol=1e-15)
    hi = brentq(lambda r: X(r) + 1, r_star, r_star + 4, xtol=1e-15, rtol=1e-15)
    g = lambda r: (m * (r + lo) / (lo**2 * r**2) + 0.5) * (m * (r + hi) / (r**2 * hi**2) + 0.5)
    P = lambda r: r**2 / 2 + m / r
    kw = dict(a=lo, b=hi, weight="alg", limit=400, epsabs=0, epsrel=1e-12)
    I0 = quad(lambda r: r**2 * np.sqrt(g(r)), wvar=(0.5, 0.5), **kw)[0]
    I1 = quad(lambda r: r**2 * P(r) * np.sqrt(g(r)), wvar=(0.5, 0.5), **kw)[0] \
        - quad(lambda r: r**2 * P(r) / (2 * np.sqrt(g(r))), wvar=(-0.5, -0.5), **kw)[0]
    dI0 = -quad(lambda r: X(r) / np.sqrt(g(r)), wvar=(-0.5, -0.5), **kw)[0]
    return float(2 - 2 * m * dI0 / I0 + I1 / I0)


def s3_small_E_limit():
    """Exact c(0+): series inversion of the thin-shell barrier and symbolic Chebyshev integrals.  Returns a sympy Rational."""
    import sympy as sp
    Xs, eps, u = sp.symbols("X epsilon u", real=True)
    N = 6
    a = sp.symbols(f"a1:{N + 1}")
    ser = sum(a[k] * u ** (k + 1) for k in range(N))
    expr = sp.expand(sp.series((1 + ser) ** -2 - (1 + ser) - u, u, 0, N + 1).removeO())
    y = ser.subs(sp.solve([expr.coeff(u, k) for k in range(1, N + 1)], a, dict=True)[0])
    rho = ((1 + y) / eps).subs(u, 2 * eps * Xs)              # rho_* = 1/eps, m = rho_*^3/2, u = 2 X/rho_*
    jac = -sp.diff(rho, Xs)
    P = rho**2 / 2 + (1 / eps) ** 3 / (2 * rho)
    sq = sp.sqrt(1 - Xs**2)
    integ = lambda e, w: sp.integrate(sp.expand(sp.series(e, eps, 0, 3).removeO()) * w, (Xs, -1, 1))
    I0 = integ(rho**2 * jac, sq)
    I1 = integ(rho**2 * P * jac, sq) - integ(rho**2 * P * jac, 1 / (2 * sq))
    return sp.Rational(2, 3) + sp.limit(sp.simplify(I1 / I0), eps, 0)
