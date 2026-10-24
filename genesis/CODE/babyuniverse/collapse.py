"""Round 6d (P-R6j): preparing the near-Nariai configuration from ordinary initial data with a regular centre.

Einstein + phi (tilted quartic: false vacuum -v, true vacuum +v) + chi (massless, the collapsing ordinary matter), spherical symmetry, double null
ds^2 = -e^{2s} du dv + r^2 dOmega^2, regular centre (Phase IV problem (a), gauge u = v on the centre and r = v/2 on the initial cone u = 0).
Equations as in ``doublenull`` with the second field added (psi = r phi, xi = r chi):
    r_uv = -r_u r_v/r - e^{2s}/(4r) + 2 pi G e^{2s} r U(phi),    s_uv = r_u r_v/r^2 + e^{2s}/(4r^2) - 4 pi G (phi_u phi_v + chi_u chi_v),
    psi_uv = phi r_uv - (e^{2s}/4) r U'(phi),                  xi_uv = chi r_uv,
and on the cone (ln alpha)_v = 4 pi G r (phi_v^2 + chi_v^2).  Same second-order diamond scheme, near-axis regular limit, singularity / null
infinity cuts.  Initial data: false vacuum inside R_w (light tanh wall), true vacuum outside, and a chi shell A exp(-((r - r_s)/w_s)^2) inside the
false vacuum: no trapped or anti-trapped sphere on the cone, null energy condition everywhere, asymptotically flat, regular centre.
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
def _rhs(r, ru, rv, s, su, sv, q, qu, qv, x, xu, xv, G, lam, vv, rho_f, near_axis):
    a2 = math.exp(2.0 * s)
    p = q / r; pu = (qu - p * ru) / r; pv = (qv - p * rv) / r
    c = x / r; cu = (xu - c * ru) / r; cv = (xv - c * rv) / r
    Uc = _U(p, lam, vv, rho_f)
    if near_axis:
        rho = Uc + (pu * pu + pv * pv + cu * cu + cv * cv) / a2
        Fr = -2.0 * math.pi * G / 3.0 * a2 * rho * r + 2.0 * math.pi * G * a2 * r * Uc
        Fs = 2.0 * math.pi * G / 3.0 * a2 * rho - 4.0 * math.pi * G * (pu * pv + cu * cv)
    else:
        Fr = -ru * rv / r - a2 / (4.0 * r) + 2.0 * math.pi * G * a2 * r * Uc
        Fs = ru * rv / (r * r) + a2 / (4.0 * r * r) - 4.0 * math.pi * G * (pu * pv + cu * cv)
    Fq = p * Fr - 0.25 * a2 * r * _dU(p, lam, vv, rho_f)
    Fx = c * Fr
    return Fr, Fs, Fq, Fx


@njit(cache=True)
def _diamond(rE, rW, rS, sE, sW, sS, qE, qW, qS, xE, xW, xS, h, G, lam, vv, rho_f, n_iter, near_axis):
    rN = rE + rW - rS; sN = sE + sW - sS; qN = qE + qW - qS; xN = xE + xW - xS
    Fs = 0.0
    for it in range(n_iter):
        rc = 0.25 * (rN + rE + rW + rS); sc = 0.25 * (sN + sE + sW + sS)
        qc = 0.25 * (qN + qE + qW + qS); xc = 0.25 * (xN + xE + xW + xS)
        ru = (rN + rW - rE - rS) / (2.0 * h); rv = (rN + rE - rW - rS) / (2.0 * h)
        su = (sN + sW - sE - sS) / (2.0 * h); sv = (sN + sE - sW - sS) / (2.0 * h)
        qu = (qN + qW - qE - qS) / (2.0 * h); qv = (qN + qE - qW - qS) / (2.0 * h)
        xu = (xN + xW - xE - xS) / (2.0 * h); xv = (xN + xE - xW - xS) / (2.0 * h)
        Fr, Fs, Fq, Fx = _rhs(rc, ru, rv, sc, su, sv, qc, qu, qv, xc, xu, xv, G, lam, vv, rho_f, near_axis)
        rN = rE + rW - rS + h * h * Fr; sN = sE + sW - sS + h * h * Fs
        qN = qE + qW - qS + h * h * Fq; xN = xE + xW - xS + h * h * Fx
    return rN, sN, qN, xN, Fs


@njit(cache=True)
def _evolve(r0, s0, q0, x0, h, n_rays, G, lam, vv, rho_f, r_exc, r_inf, n_iter, k_axis):
    N = r0.shape[0] - 1
    R = np.full((n_rays + 1, N + 1), np.nan); S = np.full((n_rays + 1, N + 1), np.nan)
    Q = np.full((n_rays + 1, N + 1), np.nan); X = np.full((n_rays + 1, N + 1), np.nan)
    rp = r0.copy(); sp = s0.copy(); qp = q0.copy(); xp = x0.copy()
    vp = np.ones(N + 1, dtype=np.bool_)
    rc = np.full(N + 1, np.nan); sc = np.full(N + 1, np.nan); qc = np.full(N + 1, np.nan); xc = np.full(N + 1, np.nan)
    vc = np.zeros(N + 1, dtype=np.bool_)
    Fs_prev = np.zeros(N + 1); Fs_cur = np.zeros(N + 1)
    for j in range(N + 1):
        R[0, j] = r0[j]; S[0, j] = s0[j]; Q[0, j] = q0[j]; X[0, j] = x0[j]
    for i in range(1, n_rays + 1):
        for j in range(N + 1):
            rc[j] = np.nan; sc[j] = np.nan; qc[j] = np.nan; xc[j] = np.nan; vc[j] = False
        if not (vp[i] and vp[i - 1] and i < N):
            break
        rc[i] = 0.0; qc[i] = 0.0; xc[i] = 0.0
        sc[i] = 2.0 * sp[i] - sp[i - 1] + h * h * Fs_prev[i]
        vc[i] = True
        if vp[i + 1]:
            for corr in range(2):
                rN, sN, qN, xN, Fs = _diamond(rp[i + 1], rc[i], rp[i], sp[i + 1], sc[i], sp[i], qp[i + 1], qc[i], qp[i],
                                              xp[i + 1], xc[i], xp[i], h, G, lam, vv, rho_f, n_iter, True)
                rv_axis = (rN + rp[i]) / (2.0 * h)
                if rv_axis > 0.0:
                    sc[i] = math.log(2.0 * rv_axis)
        for j in range(i + 1, N + 1):
            if not (vc[j - 1] and vp[j] and vp[j - 1]):
                break
            na = j - i <= k_axis
            rN, sN, qN, xN, Fs = _diamond(rp[j], rc[j - 1], rp[j - 1], sp[j], sc[j - 1], sp[j - 1], qp[j], qc[j - 1], qp[j - 1],
                                          xp[j], xc[j - 1], xp[j - 1], h, G, lam, vv, rho_f, n_iter, na)
            ok = math.isfinite(rN) and math.isfinite(sN) and math.isfinite(qN) and math.isfinite(xN) and abs(sN) < 600.0
            if ok and rN < r_exc and rN < rc[j - 1]:
                ok = False
            if ok and rN > r_inf and rN > rc[j - 1] and rN > rp[j]:
                ok = False
            if ok and j - i >= 2 and rc[j - 1] < rc[j - 2] and rc[j - 1] < rp[j - 1] and rN > rc[j - 1]:
                ok = False
            if not ok:
                break
            rc[j] = rN; sc[j] = sN; qc[j] = qN; xc[j] = xN; vc[j] = True; Fs_cur[j] = Fs
        for j in range(N + 1):
            if vc[j]:
                R[i, j] = rc[j]; S[i, j] = sc[j]; Q[i, j] = qc[j]; X[i, j] = xc[j]
            rp[j] = rc[j]; sp[j] = sc[j]; qp[j] = qc[j]; xp[j] = xc[j]; vp[j] = vc[j]; Fs_prev[j] = Fs_cur[j]
    return R, S, Q, X


def cone_data(phi_r, chi_r, h, v_max, G=1.0):
    """r = v/2 on the cone; (ln alpha)_v = 4 pi G r (phi_v^2 + chi_v^2) by the trapezoid rule on a 4x finer grid."""
    N = int(round(v_max / h))
    vf = np.linspace(0.0, N * h, 4 * N + 1)
    rf = vf / 2
    pf, cf = phi_r(rf), chi_r(rf)
    dp = np.gradient(pf, vf, edge_order=2); dc = np.gradient(cf, vf, edge_order=2)
    integ = 4 * np.pi * G * rf * (dp**2 + dc**2)
    sf = np.concatenate([[0.0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(vf))])
    return rf[::4].copy(), sf[::4].copy(), pf[::4].copy(), cf[::4].copy()


def cone_data_alpha1(phi_r, chi_r, h, v_max, G=1.0, dr=1e-6):
    """Gauge alpha = 1 on the cone (r(0) = 0, r_v(0) = 1/2): r_vv = -4 pi G r (phi'(r)^2 + chi'(r)^2) r_v^2, profiles given as functions of the
    areal radius (gauge-invariant data).  Avoids the large lapse that the r = v/2 gauge produces beyond a massive shell."""
    from scipy.integrate import solve_ivp
    N = int(round(v_max / h))
    v = np.arange(N + 1) * h
    dphi = lambda r: (phi_r(r + dr) - phi_r(r - dr)) / (2 * dr)
    dchi = lambda r: (chi_r(r + dr) - chi_r(r - dr)) / (2 * dr)

    def rhs(x, y):
        r, rv = y
        return [rv, -4 * np.pi * G * r * (dphi(r) ** 2 + dchi(r) ** 2) * rv * rv]
    sol = solve_ivp(rhs, (0.0, N * h), [0.0, 0.5], t_eval=v, method="DOP853", rtol=1e-12, atol=1e-14, max_step=h / 4)
    r = sol.y[0]
    if len(r) < N + 1:
        raise ValueError("cone integration failed")
    return r, np.zeros(N + 1), phi_r(r), chi_r(r)


def evolve(phi_r, chi_r, pot, h, v_max, u_max=None, G=1.0, r_inf=1e4, n_iter=3, k_axis=2, gauge="alpha1"):
    lam, vv, rho_f = pot
    r0, s0, p0, c0 = (cone_data_alpha1 if gauge == "alpha1" else cone_data)(phi_r, chi_r, h, v_max, G)
    N = len(r0) - 1
    n_rays = N - 1 if u_max is None else min(N - 1, int(round(u_max / h)))
    R, S, Q, X = _evolve(r0, s0, r0 * p0, r0 * c0, h, n_rays, G, lam, vv, rho_f, 2 * h, r_inf, n_iter, k_axis)
    with np.errstate(invalid="ignore", divide="ignore"):
        P = np.where(R > 0, Q / R, np.nan); C = np.where(R > 0, X / R, np.nan)
    return dict(r=R, s=S, phi=P, chi=C, h=h, G=G)


def misner_sharp(res):
    """m = (r/2)(1 + 4 r_u r_v e^{-2s}) with centred differences (NaN at edges)."""
    R, S, h = res["r"], res["s"], res["h"]
    ru = np.full_like(R, np.nan); rv = np.full_like(R, np.nan)
    ru[1:-1, :] = (R[2:, :] - R[:-2, :]) / (2 * h)
    rv[:, 1:-1] = (R[:, 2:] - R[:, :-2]) / (2 * h)
    return 0.5 * R * (1 + 4 * ru * rv * np.exp(-2 * S)), ru, rv


def potential(l_H=30.0, eps=4.0, vphi=0.0015):
    Lam = 3 / l_H**2
    rho_f = Lam / (8 * np.pi)
    lam = 4 * rho_f / (eps * vphi**4)
    return (lam, vphi, rho_f), dict(Lam=Lam, rho_f=rho_f, w=vphi * np.sqrt(eps / (2 * rho_f)), r_N=1 / np.sqrt(Lam), M_N=1 / (3 * np.sqrt(Lam)))
