"""Round 6d: an independent solver for the Round 6c throat problem (verification of the Nariai-throat daughter).

Same field equations as ``throat_evolution`` (double-null, gauge-free form, Polyakov CFT, magnetic charge, tilted quartic), but a different
numerical method and an independently written data integrator:

* **Method of lines along the characteristics.**  On each ingoing line u = const the first derivatives Z = (r_u, s_u, phi_u) obey the linear
  ODE system  dZ/dv = F(Y, Y_v, Z)  (the u-v equations read as ODEs in v), integrated outward from the ingoing data ray v = 0 with Heun's
  method; Y_v comes from second-order finite differences.  Y = (r, s, phi) is then advanced in u with Heun's predictor-corrector,
  Y(u + du) = Y + du/2 [Z(Y) + Z(Y*)], du = h / n_sub.  No diamond stencil, no cell-centre averaging.  The u-update applies a Volterra
  operator whose norm is ~ m^2 L_v / 4 for a field of mass m on a line of length L_v (m^2 = 21, L_v = 200 in the Nariai runs), so the
  explicit step needs du * m^2 L_v / 4 <~ 1: n_sub = 64 at h = 0.08.
* **Data** from a fixed-step RK4 integration of the null constraints on a 4x finer grid (``throat_evolution`` uses adaptive DOP853).

Breakdown handling: once r leaves (r_exc, r_inf), |phi| > phi_max or s is not finite at some v on a line, that point and every later v on the
line are discarded (the region beyond an inflating region's future null infinity or a singularity), exactly as in the diamond code.
"""

from __future__ import annotations

import math

import numpy as np

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
def _F(r, s, p, rv, sv, pv, ru, su, pu, G, Q2, k, lam, vv, rho_f):
    a2 = math.exp(2.0 * s)
    Lam = 8.0 * math.pi * G * _U(p, lam, vv, rho_f)
    Fs = (ru * rv / (r * r) + a2 / (4.0 * r * r) - 4.0 * math.pi * G * pu * pv - a2 * Q2 / (2.0 * r ** 4)) / (1.0 - k / (r * r))
    Fr = -ru * rv / r + a2 * (Lam * r * r + Q2 / (r * r) - 1.0) / (4.0 * r) - (k / r) * Fs
    Fp = -(ru * pv + rv * pu) / r - 0.25 * a2 * _dU(p, lam, vv, rho_f)
    return Fr, Fs, Fp


@njit(cache=True)
def _vderiv(Y, h, n):
    d = np.empty(n)
    d[0] = (-3.0 * Y[0] + 4.0 * Y[1] - Y[2]) / (2.0 * h)
    for j in range(1, n - 1):
        d[j] = (Y[j + 1] - Y[j - 1]) / (2.0 * h)
    d[n - 1] = (3.0 * Y[n - 1] - 4.0 * Y[n - 2] + Y[n - 3]) / (2.0 * h)
    return d


@njit(cache=True)
def _line(R, S, P, n, ru0, su0, pu0, h, G, Q2, k, lam, vv, rho_f):
    """Integrate dZ/dv = F along one line for the valid prefix j < n; returns Z arrays."""
    RV = _vderiv(R, h, n); SV = _vderiv(S, h, n); PV = _vderiv(P, h, n)
    RU = np.empty(n); SU = np.empty(n); PU = np.empty(n)
    RU[0] = ru0; SU[0] = su0; PU[0] = pu0
    for j in range(n - 1):
        f1 = _F(R[j], S[j], P[j], RV[j], SV[j], PV[j], RU[j], SU[j], PU[j], G, Q2, k, lam, vv, rho_f)
        zr = RU[j] + h * f1[0]; zs = SU[j] + h * f1[1]; zp = PU[j] + h * f1[2]
        f2 = _F(R[j + 1], S[j + 1], P[j + 1], RV[j + 1], SV[j + 1], PV[j + 1], zr, zs, zp, G, Q2, k, lam, vv, rho_f)
        RU[j + 1] = RU[j] + 0.5 * h * (f1[0] + f2[0])
        SU[j + 1] = SU[j] + 0.5 * h * (f1[1] + f2[1])
        PU[j + 1] = PU[j] + 0.5 * h * (f1[2] + f2[2])
    return RU, SU, PU


@njit(cache=True)
def _valid_prefix(R, S, P, n, r_exc, r_inf, phi_max):
    for j in range(n):
        if not (R[j] > r_exc and R[j] < r_inf and abs(P[j]) < phi_max and abs(S[j]) < 1e3) or R[j] != R[j]:
            return j
    return n


@njit(cache=True)
def _evolve(ra, sa, pa, rb, sb, pb, rub, sub, pub, h, n_sub, G, Q2, k, lam, vv, rho_f, r_exc, r_inf, phi_max):
    """rb, sb, pb, rub, sub, pub are given on the fine u-grid of step h/n_sub; every n_sub-th line is stored."""
    Nf = rb.shape[0] - 1
    Nu = Nf // n_sub
    Nv = ra.shape[0] - 1
    hu = h / n_sub
    R = np.full((Nu + 1, Nv + 1), np.nan); S = np.full((Nu + 1, Nv + 1), np.nan); P = np.full((Nu + 1, Nv + 1), np.nan)
    R[0, :] = ra; S[0, :] = sa; P[0, :] = pa
    Rc = ra.copy(); Sc = sa.copy(); Pc = pa.copy()
    n = Nv + 1
    for i in range(Nf):
        n = _valid_prefix(Rc, Sc, Pc, n, r_exc, r_inf, phi_max)
        if n < 4:
            break
        RU, SU, PU = _line(Rc, Sc, Pc, n, rub[i], sub[i], pub[i], h, G, Q2, k, lam, vv, rho_f)
        Rs = np.empty(n); Ss = np.empty(n); Ps = np.empty(n)
        for j in range(n):
            Rs[j] = Rc[j] + hu * RU[j]; Ss[j] = Sc[j] + hu * SU[j]; Ps[j] = Pc[j] + hu * PU[j]
        Rs[0] = rb[i + 1]; Ss[0] = sb[i + 1]; Ps[0] = pb[i + 1]
        m = _valid_prefix(Rs, Ss, Ps, n, r_exc, r_inf, phi_max)
        if m < 4:
            break
        RU2, SU2, PU2 = _line(Rs, Ss, Ps, m, rub[i + 1], sub[i + 1], pub[i + 1], h, G, Q2, k, lam, vv, rho_f)
        for j in range(m):
            Rc[j] = Rc[j] + 0.5 * hu * (RU[j] + RU2[j])
            Sc[j] = Sc[j] + 0.5 * hu * (SU[j] + SU2[j])
            Pc[j] = Pc[j] + 0.5 * hu * (PU[j] + PU2[j])
        for j in range(m, n):
            Rc[j] = np.nan; Sc[j] = np.nan; Pc[j] = np.nan
        Rc[0] = rb[i + 1]; Sc[0] = sb[i + 1]; Pc[0] = pb[i + 1]
        n = m
        if (i + 1) % n_sub == 0:
            ii = (i + 1) // n_sub
            for j in range(n):
                R[ii, j] = Rc[j]; S[ii, j] = Sc[j]; P[ii, j] = Pc[j]
    return R, S, P


def _rk4_ray(rhs, y0, h_out, n_out, sub=4):
    """Fixed-step RK4 on a grid sub times finer than the output grid; returns y at the output points."""
    hh = h_out / sub
    y = np.array(y0, float)
    out = np.empty((n_out + 1, len(y0)))
    out[0] = y
    x = 0.0
    for i in range(n_out):
        for _ in range(sub):
            k1 = rhs(x, y); k2 = rhs(x + hh / 2, y + hh / 2 * k1); k3 = rhs(x + hh / 2, y + hh / 2 * k2); k4 = rhs(x + hh, y + hh * k3)
            y = y + hh / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            x += hh
        out[i + 1] = y
    return out


def characteristic_data(r0, ru0, rv0, h, Nv, Nu, phi_v, phi_u, flux_u=None, flux_v=None, G=1.0, eps=1e-6, n_sub=1):
    """Gauge s = 0 on both rays.  Outgoing ray u = 0: r_vv = -4 pi G r phi_v'^2 - (G/r) flux_v(v); ingoing ray v = 0:
    r_uu = -4 pi G r phi_u'^2 - (G/r) flux_u(u).  flux_* are the 2D null energy fluxes T2 (CFT state plus classical null fluid)."""
    fu = flux_u or (lambda x: 0.0)
    fv = flux_v or (lambda x: 0.0)
    dpv = lambda x: (phi_v(x + eps) - phi_v(x - eps)) / (2 * eps)
    dpu = lambda x: (phi_u(x + eps) - phi_u(x - eps)) / (2 * eps)
    ya = _rk4_ray(lambda x, y: np.array([y[1], -4 * np.pi * G * y[0] * dpv(x) ** 2 - G / y[0] * fv(x)]), [r0, rv0], h, Nv)
    yb = _rk4_ray(lambda x, y: np.array([y[1], -4 * np.pi * G * y[0] * dpu(x) ** 2 - G / y[0] * fu(x)]), [r0, ru0], h / n_sub,
                  Nu * n_sub, sub=max(1, 4 // n_sub))
    v = np.arange(Nv + 1) * h
    u = np.arange(Nu * n_sub + 1) * h / n_sub
    pa = np.array([phi_v(x) for x in v]); pb = np.array([phi_u(x) for x in u])
    pub = np.array([dpu(x) for x in u])
    return dict(ra=ya[:, 0], sa=np.zeros(Nv + 1), pa=pa, rb=yb[:, 0], sb=np.zeros(len(u)), pb=pb, rub=yb[:, 1], sub=np.zeros(len(u)),
                pub=pub, rva=ya[:, 1], n_sub=n_sub)


def evolve(data, h, Q=0.0, c=0.0, pot=(0.0, 1.0, 0.0), G=1.0, r_exc=1e-3, r_inf=1e4, phi_max=1e300):
    lam, vv, rho_f = pot
    k = G * c / (12 * np.pi)
    R, S, P = _evolve(data["ra"], data["sa"], data["pa"], data["rb"], data["sb"], data["pb"], data["rub"], data["sub"], data["pub"],
                      h, data.get("n_sub", 1), G, Q * Q, k, lam, vv, rho_f, r_exc, r_inf, phi_max)
    return dict(r=R, s=S, phi=P, h=h, Q=Q, c=c, k=k, G=G)


def nariai(l_H=30.0, delta=1e-8, vphi=0.005, eps=2.0, wall_frac=0.95, h=0.08, Lu=200.0, Lv=200.0, flux_u=None, flux_v=None, c=0.0,
           u_inner=None, w_inner=2.0, Q=0.0, n_sub=64):
    """The Round 6c near-Nariai configuration, rebuilt independently (Q = 0: Schwarzschild-de Sitter; Q > 0: the RNdS Nariai branch with
    the same construction).  Returns the solution and the growth / hiding diagnostics of ``throat_evolution`` (shared region map only)."""
    from scipy.optimize import brentq
    from . import throat_evolution as TE
    Lam = 3 / l_H**2
    rho_f = Lam / (8 * np.pi)
    lam = 4 * rho_f / (eps * vphi**4)
    w = np.sqrt(2) / (np.sqrt(lam) * vphi)
    if Q == 0:
        rN = 1 / np.sqrt(Lam)
        M = (1 / (3 * np.sqrt(Lam))) * (1 - delta)
    else:
        rN = np.sqrt((1 + np.sqrt(1 - 12 * Lam * Q * Q)) / (2 * Lam))
        MN = 0.5 * (rN + Q * Q / rN - Lam * rN**3 / 3)
        M = MN * (1 - delta)
    f = lambda r: 1 - 2 * M / r + Q * Q / r**2 - Lam * r**2 / 3
    rc = brentq(f, rN, 3 * rN)
    f0 = f(rN)
    rv0 = np.sqrt(f0) / 2
    v2 = wall_frac * (rc - rN) / rv0
    phi_v = lambda x: vphi * np.tanh((x - v2) / w)
    if u_inner is None:
        phi_u = lambda x: -vphi
    else:
        phi_u = lambda x: -vphi * np.tanh((u_inner - x) / w_inner) / np.tanh(u_inner / w_inner)
    data = characteristic_data(rN, -rv0, rv0, h, int(Lv / h), int(Lu / h), phi_v, phi_u, flux_u, flux_v, n_sub=n_sub)
    res = evolve(data, h, Q=Q, c=c, pot=(lam, vphi, rho_f), r_inf=1e4, phi_max=10 * vphi)
    cl = TE.classify(res)
    g = TE.inflated(res, cl, rN)
    C, P, R = cl["cls"], res["phi"], res["r"]
    tv = (P > 0) & np.isfinite(R)
    return dict(res=res, cls=C, grew=g["grew"], n_big=g["n_big"], r_max=g["r_max_antifalse"], r_N=float(rN), f0=float(f0),
                true_vacuum_trapped=int(np.sum(tv & (C == -1))), true_vacuum_untrapped_large_r=int(np.sum(tv & (C >= 0) & (R > 5 * rN))))
