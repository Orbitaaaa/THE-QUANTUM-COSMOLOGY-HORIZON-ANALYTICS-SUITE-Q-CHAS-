"""Stage 1 task 3: an independent check of the universal suppression law (Theorem S3) and of its parameter range.

Nothing here imports model.py, tunnelling.py or bounds.py.  Everything is rebuilt from the junction data in mpmath:

    f_in = 1 - chi^2 R^2,  f_out = 1 - 2GM/R,  A_out = (f_in - f_out - kappa^2 R^2) / (2 kappa R),  V = f_out - A_out^2,
    cos(theta) = (f_in + f_out - kappa^2 R^2) / (2 sqrt(f_in f_out)),  sin(theta) = kappa R sqrt(V / (f_in f_out)),
    B_an = (2/G) int_{V>0} R theta dR   (theta = atan2(sin, cos), the analytic determination on the barrier),
    B_sq = (2/G) int_{V>0} R arccos|cos theta| dR,

and Phi(E) = (16 / (pi E^2)) int_{g1<E} z^2 sqrt(E - g1) dz with g1 = (z^3 - 1)^2 / z^4, directly in z.  The mass is fixed by
E = 4 kappa^2 / (L^(8/3) (2GM)^(2/3)), L^2 = chi^2 + kappa^2 (eta = 0).  B(0) is the same quadrature at M = 0.
"""

from __future__ import annotations

import mpmath as mp


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
