"""Round 6c (P-R6f): semiclassical spherically symmetric Einstein-Maxwell-scalar-CFT evolution in double-null coordinates.

Metric ds^2 = -e^{2s} du dv + r^2 dOmega^2 (u, v increasing to the future, v outward; normal r_u < 0 < r_v).  Sources:

* a magnetic charge Q (geometric units), G T_uv = e^{2s} Q^2/(16 pi r^4), G T_thth = Q^2/(8 pi r^2);
* a minimally coupled scalar with the tilted quartic U(phi) of Round 5 (false vacuum -v, true vacuum +v);
* c two-dimensional conformal fields in the (u, v) plane (the s-wave / lowest-Landau-level sector), with the Polyakov stress tensor
      T2_uu = -(c/12 pi)(s_u^2 - s_uu) + t_u(u),  T2_vv = -(c/12 pi)(s_v^2 - s_vv) + t_v(v),  T2_uv = -(c/12 pi) s_uv,
  and T_ab = T2_ab / (4 pi r^2), zero angular stress (conserved in 4D);
* optionally a classical outgoing null fluid N(u) of either sign (the stand-in for classical NEC violation, Route D), T2_uu += N(u).

The general equations (sympy, `scratchpad/r6c/einstein.py`, also `symbolic_check`):
    r_uv = -r_u r_v / r - e^{2s}/(4 r) + 4 pi G r T_uv,
    s_uv = r_u r_v / r^2 + e^{2s}/(4 r^2) - 4 pi G T_uv - 2 pi G e^{2s} T_thth / r^2,
    r_uu - 2 s_u r_u = -4 pi G r T_uu,   r_vv - 2 s_v r_v = -4 pi G r T_vv   (constraints).
With k = G c / 12 pi, Lambda(phi) = 8 pi G U(phi) and W = Lambda r^2 + Q^2/r^2 - 1:
    s_uv (1 - k/r^2) = r_u r_v/r^2 + e^{2s}/(4 r^2) - 4 pi G phi_u phi_v - e^{2s} Q^2/(2 r^4),
    r_uv = -r_u r_v / r + e^{2s} W/(4 r) - (k/r) s_uv,
    phi_uv = -(r_u phi_v + r_v phi_u)/r - (e^{2s}/4) U'(phi).
The state functions t_u, t_v and the null fluid enter only the constraints, i.e. the characteristic data.

Data: gauge s = s0 on both initial rays.  On v = v0 (ingoing, index i) E_uu gives r_uu = -4 pi G r phi_u^2 - (G/r)(t_u + N); on u = u0
(outgoing, index j) r_vv = -4 pi G r phi_v^2 - (G/r) t_v.  Scheme: second-order diamond X_N = X_E + X_W - X_S + h^2 F(centre) with a
three-pass predictor-corrector, as in Phase IV (`doublenull`).

The 1 - k/r^2 factor is the s-wave form of the condition S_BH(r) > c/12 (Proposition GSL-5c): the model is used only where it holds.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp

try:
    from numba import njit
except ImportError:          # pragma: no cover
    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda f: f


@njit(cache=True)
def _U(p, lam, vv, rho_f):
    x = p / vv
    return lam / 4.0 * (p * p - vv * vv) ** 2 + rho_f * (0.5 - 0.75 * x + 0.25 * x ** 3)


@njit(cache=True)
def _dU(p, lam, vv, rho_f):
    x = p / vv
    return lam * p * (p * p - vv * vv) - rho_f * (0.75 - 0.75 * x * x) / vv


@njit(cache=True)
def _rhs(r, ru, rv, s, su, sv, p, pu, pv, G, Q2, k, lam, vv, rho_f):
    a2 = math.exp(2.0 * s)
    Lam = 8.0 * math.pi * G * _U(p, lam, vv, rho_f)
    W = Lam * r * r + Q2 / (r * r) - 1.0
    As = ru * rv / (r * r) + a2 / (4.0 * r * r) - 4.0 * math.pi * G * pu * pv - a2 * Q2 / (2.0 * r ** 4)
    Fs = As / (1.0 - k / (r * r))
    Fr = -ru * rv / r + a2 * W / (4.0 * r) - (k / r) * Fs
    Fp = -(ru * pv + rv * pu) / r - 0.25 * a2 * _dU(p, lam, vv, rho_f)
    return Fr, Fs, Fp


@njit(cache=True)
def _diamond(rE, rW, rS, sE, sW, sS, pE, pW, pS, h, G, Q2, k, lam, vv, rho_f, n_iter):
    rN = rE + rW - rS
    sN = sE + sW - sS
    pN = pE + pW - pS
    for it in range(n_iter):
        rc = 0.25 * (rN + rE + rW + rS)
        sc = 0.25 * (sN + sE + sW + sS)
        pc = 0.25 * (pN + pE + pW + pS)
        ru = (rN + rW - rE - rS) / (2.0 * h)
        rv = (rN + rE - rW - rS) / (2.0 * h)
        su = (sN + sW - sE - sS) / (2.0 * h)
        sv = (sN + sE - sW - sS) / (2.0 * h)
        pu = (pN + pW - pE - pS) / (2.0 * h)
        pv = (pN + pE - pW - pS) / (2.0 * h)
        Fr, Fs, Fp = _rhs(rc, ru, rv, sc, su, sv, pc, pu, pv, G, Q2, k, lam, vv, rho_f)
        rN = rE + rW - rS + h * h * Fr
        sN = sE + sW - sS + h * h * Fs
        pN = pE + pW - pS + h * h * Fp
    return rN, sN, pN


@njit(cache=True)
def _evolve(ra, sa, pa, rb, sb, pb, h, G, Q2, k, lam, vv, rho_f, n_iter, r_exc, r_inf, phi_max):
    """ra.. on u = u0 (j = 0..N), rb.. on v = v0 (i = 0..M).  Returns full arrays (NaN beyond breakdown) and an end code per ray."""
    N = ra.shape[0] - 1
    M = rb.shape[0] - 1
    R = np.full((M + 1, N + 1), np.nan)
    S = np.full((M + 1, N + 1), np.nan)
    P = np.full((M + 1, N + 1), np.nan)
    end = np.zeros(M + 1, dtype=np.int64)
    R[0, :] = ra
    S[0, :] = sa
    P[0, :] = pa
    for i in range(1, M + 1):
        R[i, 0] = rb[i]
        S[i, 0] = sb[i]
        P[i, 0] = pb[i]
        for j in range(1, N + 1):
            if not (math.isfinite(R[i - 1, j]) and math.isfinite(R[i - 1, j - 1]) and math.isfinite(R[i, j - 1])):
                break
            rN, sN, pN = _diamond(R[i - 1, j], R[i, j - 1], R[i - 1, j - 1], S[i - 1, j], S[i, j - 1], S[i - 1, j - 1],
                                  P[i - 1, j], P[i, j - 1], P[i - 1, j - 1], h, G, Q2, k, lam, vv, rho_f, n_iter)
            if not (math.isfinite(rN) and math.isfinite(sN) and math.isfinite(pN)) or abs(sN) > 500.0 or abs(pN) > phi_max:
                end[i] = 3                       # numerical breakdown (in practice next to the singularity)
                break
            if rN < r_exc or rN * rN <= 1.05 * k:
                end[i] = 1                       # singularity, or the s-wave model's validity edge S_BH = c/12
                break
            if rN > r_inf:
                end[i] = 2                       # grown beyond the tracked range
                break
            R[i, j] = rN
            S[i, j] = sN
            P[i, j] = pN
    return R, S, P, end


class Potential:
    """Tilted quartic (Round 5): false vacuum phi = -v with U = rho_f, true vacuum +v with U = 0.  rho_f = 0, lam = 0: no scalar."""

    def __init__(self, lam: float = 0.0, v: float = 1.0, rho_f: float = 0.0):
        self.lam, self.v, self.rho_f = lam, v, rho_f

    def U(self, p):
        x = p / self.v
        return self.lam / 4 * (p**2 - self.v**2) ** 2 + self.rho_f * (0.5 - 0.75 * x + 0.25 * x**3)


def characteristic_data(Q: float, c: float, r0: float, ru0: float, rv0: float, h: float, Nv: int, Nu: int,
                        phi_v=None, phi_u=None, t_u=None, N_u=None, G: float = 1.0, s0: float = 0.0, N_v=None):
    """Constraint-satisfying data in the gauge s = s0 on both rays.  phi_v(v), phi_u(u): scalar profiles on u = u0 and v = v0 (must
    agree at the corner); t_u(u): CFT state flux; N_u(u): classical null-fluid flux.  Returns (ra, sa, pa), (rb, sb, pb)."""
    v = np.arange(Nv + 1) * h
    u = np.arange(Nu + 1) * h
    phi_v = phi_v or (lambda x: np.zeros_like(np.asarray(x, float)))
    phi_u = phi_u or (lambda x: np.zeros_like(np.asarray(x, float)))
    t_u = t_u or (lambda x: np.zeros_like(np.asarray(x, float)))
    N_u = N_u or (lambda x: np.zeros_like(np.asarray(x, float)))
    N_v = N_v or (lambda x: np.zeros_like(np.asarray(x, float)))
    eps = 1e-6 * h

    def dphi(f, x):
        return (f(x + eps) - f(x - eps)) / (2 * eps)

    def rhs_v(x, y):
        return [y[1], -4 * np.pi * G * y[0] * dphi(phi_v, x) ** 2 - G / y[0] * N_v(x)]

    def rhs_u(x, y):
        return [y[1], -4 * np.pi * G * y[0] * dphi(phi_u, x) ** 2 - G / y[0] * (t_u(x) + N_u(x))]

    sa = solve_ivp(rhs_v, (0, v[-1]), [r0, rv0], t_eval=v, method="DOP853", rtol=1e-12, atol=1e-14, max_step=h / 2)
    sb = solve_ivp(rhs_u, (0, u[-1]), [r0, ru0], t_eval=u, method="DOP853", rtol=1e-12, atol=1e-14, max_step=h / 2)
    if sa.y.shape[1] < Nv + 1 or sb.y.shape[1] < Nu + 1:
        raise ValueError("an initial ray reached r = 0")
    ra, rb = sa.y[0], sb.y[0]
    return (ra, np.full(Nv + 1, s0), phi_v(v)), (rb, np.full(Nu + 1, s0), phi_u(u)), dict(rv_a=sa.y[1], ru_b=sb.y[1])


def evolve(data_a, data_b, Q: float, c: float, pot: Potential, h: float, G: float = 1.0, n_iter: int = 3,
           r_exc: float = 1e-3, r_inf: float = 1e6, phi_max: float = None) -> dict:
    (ra, sa, pa), (rb, sb, pb) = data_a, data_b
    k = G * c / (12 * np.pi)
    R, S, P, end = _evolve(np.asarray(ra, float), np.asarray(sa, float), np.asarray(pa, float), np.asarray(rb, float),
                           np.asarray(sb, float), np.asarray(pb, float), h, G, Q * Q, k, pot.lam, pot.v, pot.rho_f, n_iter,
                           r_exc, r_inf, phi_max if phi_max is not None else 1e300)
    return dict(r=R, s=S, phi=P, end=end, h=h, Q=Q, c=c, k=k, G=G)


# ------------------------------------------------------------------ diagnostics
def derivatives(res: dict) -> dict:
    """Centred r_u, r_v, s_u, s_v (NaN where neighbours are missing)."""
    R, S, h = res["r"], res["s"], res["h"]
    ru = np.full_like(R, np.nan)
    rv = np.full_like(R, np.nan)
    ru[1:-1, :] = (R[2:, :] - R[:-2, :]) / (2 * h)
    rv[:, 1:-1] = (R[:, 2:] - R[:, :-2]) / (2 * h)
    su = np.full_like(S, np.nan)
    sv = np.full_like(S, np.nan)
    su[1:-1, :] = (S[2:, :] - S[:-2, :]) / (2 * h)
    sv[:, 1:-1] = (S[:, 2:] - S[:, :-2]) / (2 * h)
    return dict(ru=ru, rv=rv, su=su, sv=sv)


def misner_sharp(res: dict, d: dict = None) -> np.ndarray:
    """Mass parameter M = (r/2)(1 + Q^2/r^2 + 4 r_u r_v e^{-2s}) (Reissner-Nordstrom M in electrovacuum)."""
    d = d or derivatives(res)
    R = res["r"]
    return 0.5 * R * (1 + res["Q"] ** 2 / R**2 + 4 * d["ru"] * d["rv"] * np.exp(-2 * res["s"]))


def constraint_residual(res: dict, t_u=None, N_u=None, t_v=None) -> dict:
    """Relative residuals of E_vv (along outgoing rays) and E_uu (along ingoing rays), including the CFT terms, in the interior."""
    R, S, P, h, G, c = res["r"], res["s"], res["phi"], res["h"], res["G"], res["c"]
    M, N = R.shape
    u = np.arange(M) * h
    v = np.arange(N) * h
    tu = (t_u(u) if t_u else 0 * u) + (N_u(u) if N_u else 0 * u)
    tv = t_v(v) if t_v else 0 * v
    rvv = (R[:, 2:] - 2 * R[:, 1:-1] + R[:, :-2]) / h**2
    rv = (R[:, 2:] - R[:, :-2]) / (2 * h)
    sv = (S[:, 2:] - S[:, :-2]) / (2 * h)
    svv = (S[:, 2:] - 2 * S[:, 1:-1] + S[:, :-2]) / h**2
    pv = (P[:, 2:] - P[:, :-2]) / (2 * h)
    r = R[:, 1:-1]
    T2vv = -(c / (12 * np.pi)) * (sv**2 - svv) + tv[None, 1:-1]
    Evv = rvv - 2 * sv * rv + 4 * np.pi * G * r * pv**2 + G / r * T2vv
    scale_v = np.abs(rvv) + np.abs(2 * sv * rv) + np.abs(4 * np.pi * G * r * pv**2) + np.abs(G / r * T2vv) + 1e-300
    ruu = (R[2:, :] - 2 * R[1:-1, :] + R[:-2, :]) / h**2
    ru = (R[2:, :] - R[:-2, :]) / (2 * h)
    su = (S[2:, :] - S[:-2, :]) / (2 * h)
    suu = (S[2:, :] - 2 * S[1:-1, :] + S[:-2, :]) / h**2
    pu = (P[2:, :] - P[:-2, :]) / (2 * h)
    r2 = R[1:-1, :]
    T2uu = -(c / (12 * np.pi)) * (su**2 - suu) + tu[1:-1, None]
    Euu = ruu - 2 * su * ru + 4 * np.pi * G * r2 * pu**2 + G / r2 * T2uu
    scale_u = np.abs(ruu) + np.abs(2 * su * ru) + np.abs(4 * np.pi * G * r2 * pu**2) + np.abs(G / r2 * T2uu) + 1e-300
    return dict(Evv=np.nanmax(np.abs(Evv)), Euu=np.nanmax(np.abs(Euu)),
                Evv_rel=np.nanmedian(np.abs(Evv) / scale_v), Euu_rel=np.nanmedian(np.abs(Euu) / scale_u))


def inflated(res: dict, cl: dict, r_ref: float, factor: float = 5.0, min_area: float = 1.0) -> dict:
    """Daughter growth detector: a coherent set (coordinate area >= min_area, |s| < 30, so not the breakdown edge next to a singularity) of
    anti-trapped points in the false vacuum (phi < 0) with r > factor * r_ref, and the largest such r."""
    P, R, S = res["phi"], res["r"], res["s"]
    good = (cl["cls"] == 1) & (P < 0) & np.isfinite(R) & (np.abs(S) < 30)
    big = good & (R > factor * r_ref)
    return dict(grew=bool(big.sum() * res["h"] ** 2 >= min_area), r_max_antifalse=float(np.nanmax(np.where(good, R, np.nan))) if good.any() else None,
                n_big=int(big.sum()), n_good=int(good.sum()))


def classify(res: dict, d: dict = None) -> dict:
    """Region map: +1 anti-trapped (r_u, r_v > 0), -1 trapped, 0 normal; summary of the anti-trapped set and of growth."""
    d = d or derivatives(res)
    ru, rv, R = d["ru"], d["rv"], res["r"]
    cls = np.where((ru > 0) & (rv > 0), 1, np.where((ru < 0) & (rv < 0), -1, 0)).astype(float)
    cls[~np.isfinite(ru) | ~np.isfinite(rv)] = np.nan
    anti = cls == 1
    trapped = cls == -1
    return dict(cls=cls, n_anti=int(np.sum(anti)), n_trapped=int(np.sum(trapped)),
                r_max_anti=float(np.nanmax(np.where(anti, R, np.nan))) if anti.any() else None,
                r_max=float(np.nanmax(R)), r_min=float(np.nanmin(R)))


# ------------------------------------------------------------------ Round 6c: the Nariai throat (positive result)
def conformal_pulse(c: float, A: float, up: float, w: float, L: float, W2: float):
    """A CFT state obtained from the vacuum by the conformal map X(u) with y = ln X', y' = A s1(u) s2(u): a rise of width w at up and a
    fall of width W2 at up + L.  Its flux t_u = (c/24 pi)(y'^2/2 - y'') saturates the strong QNEC (Proposition 6c-7) and obeys ANEC
    (int t_u = (c/48 pi) int y'^2 >= 0) and the Fewster-Hollands null QEI.  The negative lobe -(c/24 pi) A (rise) comes first; the
    compensating positive energy is spread over L.  Returns t_u(u) and y'(u)."""
    def yp(u):
        u = np.asarray(u, float)
        return A * 0.25 * (1 + np.tanh((u - up) / w)) * (1 - np.tanh((u - up - L) / W2))

    def ypp(u):
        u = np.asarray(u, float)
        s1, s2 = 0.5 * (1 + np.tanh((u - up) / w)), 0.5 * (1 - np.tanh((u - up - L) / W2))
        return A * (0.5 / w / np.cosh((u - up) / w) ** 2 * s2 - s1 * 0.5 / W2 / np.cosh((u - up - L) / W2) ** 2)

    return (lambda u: c / (24 * np.pi) * (0.5 * yp(u) ** 2 - ypp(u))), yp


def nariai_false_vacuum(l_H: float = 30.0, delta: float = 1e-8, vphi: float = 0.005, eps: float = 2.0, wall_frac: float = 0.95,
                        h: float = 0.08, Lu: float = 200.0, Lv: float = 200.0, Nflux: float = 0.0, up: float = 2.0, width: float = 1.0,
                        c: float = 0.0, Y: float = 0.0, wp: float = 1.0, tu_func=None, u_inner: float = None,
                        w_inner: float = 2.0, r0_factor: float = 1.0, Q: float = 0.0, Nv_func=None, r_inf: float = 1e4) -> dict:
    """Near-Nariai configuration (Proposition 6c-9): a Schwarzschild-de Sitter black hole of mass M = M_N (1 - delta), M_N = 1/(3 sqrt Lambda),
    Lambda = 3/l_H^2, inside a false-vacuum region (tilted quartic, light thin walls) whose outer wall sits at wall_frac of the way from
    r_N = 1/sqrt Lambda to the cosmological horizon r_c (so the data are normal).  Corner at r_N with symmetric expansions; false vacuum
    on the whole ingoing ray (the hole is inside the false vacuum).  Optional outgoing negative-energy pulse: classical null fluid of
    amplitude Nflux or a QNEC-saturating CFT step (c, Y).  "delta_m" is the pulse's local energy W = int|N| du = Nflux sqrt(pi) width
    in the corner's geodesic frame; "dm_misner_sharp" = sqrt(f0) W is its invariant (Misner-Sharp) size.  Returns the evolution and its
    diagnostics."""
    from scipy.optimize import brentq
    Lam = 3 / l_H**2
    rho_f = Lam / (8 * np.pi)
    lam = 4 * rho_f / (eps * vphi**4)
    w = np.sqrt(2) / (np.sqrt(lam) * vphi)
    pot = Potential(lam, vphi, rho_f)
    if Q == 0:
        rN = 1 / np.sqrt(Lam)
        MN = 1 / (3 * np.sqrt(Lam))
    else:                                   # Nariai branch of Reissner-Nordstrom-de Sitter (Lambda Q^2 < 1/12)
        rN = np.sqrt((1 + np.sqrt(1 - 12 * Lam * Q * Q)) / (2 * Lam))
        MN = 0.5 * (rN + Q * Q / rN - Lam * rN**3 / 3)
    M = MN * (1 - delta)
    f = lambda r: 1 - 2 * M / r + Q * Q / r**2 - Lam * r**2 / 3
    f0 = f(rN)
    r0 = rN * r0_factor
    if r0_factor != 1.0:                    # corner moved (semiclassical re-tuning): same corner expansions, mass fixed by the corner
        M = 0.5 * r0 * (1 - f0 + Q * Q / r0**2 - Lam * r0**2 / 3)
        f = lambda r: 1 - 2 * M / r + Q * Q / r**2 - Lam * r**2 / 3
    rc = brentq(f, r0, 3 * rN)
    rv0 = np.sqrt(f0) / 2
    rgrid = np.linspace(0.3 * rN, r0, 4001)
    sgn = np.where(np.diff(np.sign(f(rgrid))) != 0)[0]
    rb = brentq(f, rgrid[sgn[-1]], rgrid[sgn[-1] + 1]) if len(sgn) else np.nan
    v2 = wall_frac * (rc - r0) / rv0
    phi_v = lambda vv: vphi * np.tanh((np.asarray(vv, float) - v2) / w)
    if u_inner is None:
        phi_u = lambda uu: -vphi + 0 * np.asarray(uu, float)
    else:                                   # inner wall: true vacuum deeper than u_inner along the ingoing ray (a hole formed by collapse)
        phi_u = lambda uu: -vphi * np.tanh((u_inner - np.asarray(uu, float)) / w_inner) / np.tanh(u_inner / w_inner)
    tu = tu_func
    if tu is None and c > 0 and Y != 0:
        tu = lambda u: c / (24 * np.pi) * (0.5 * (0.5 * Y * (1 + np.tanh((u - up) / wp))) ** 2 - 0.5 * Y / wp / np.cosh((u - up) / wp) ** 2)
    Nu = (lambda u: -Nflux * np.exp(-((u - up) / width) ** 2)) if Nflux else None
    a, b, info = characteristic_data(Q, c, r0, -rv0, rv0, h, int(Lv / h), int(Lu / h), phi_v=phi_v, phi_u=phi_u, t_u=tu, N_u=Nu,
                                     N_v=Nv_func)
    res = evolve(a, b, Q, c, pot, h, r_inf=r_inf, phi_max=10 * vphi)
    cl = classify(res)
    g = inflated(res, cl, rN)
    P, R, C, S = res["phi"], res["r"], cl["cls"], res["s"]
    good = (C == 1) & (P < 0) & (np.abs(S) < 30)
    rmax_u = np.nanmax(np.where(good, R, np.nan), axis=1)
    k = np.where(rmax_u > 5 * rN)[0]
    tv = (P > 0) & np.isfinite(R)
    uu = np.arange(int(Lu / h) + 1) * h
    flux = (tu(uu) if tu is not None else 0 * uu) + (Nu(uu) if Nu is not None else 0 * uu)
    return dict(res=res, cls=C, grew=g["grew"], n_big=g["n_big"], r_max=g["r_max_antifalse"], r_N=float(rN), M_N=float(MN),
                neg_energy=float(-np.trapezoid(np.minimum(flux, 0), uu)), net_energy=float(np.trapezoid(flux, uu)),
                dm_misner_sharp=float(np.sqrt(f0) * -np.trapezoid(np.minimum(flux, 0), uu)), r_c=float(rc), M=float(M),
                r_b=float(rb), u_hole=float((r0 - rb) / rv0), rv_data=info["rv_a"], v_wall=float(v2),
                f0=float(f0), delta_m=float(Nflux * np.sqrt(np.pi) * width) if Nflux else (float(c * Y / (24 * np.pi)) if c else 0.0),
                u_growth=float(k[0] * h) if len(k) else None, true_vacuum_trapped=int(np.sum(tv & (C == -1))),
                true_vacuum_untrapped_large_r=int(np.sum(tv & (C >= 0) & (R > 5 * rN))), n_anti=cl["n_anti"])
