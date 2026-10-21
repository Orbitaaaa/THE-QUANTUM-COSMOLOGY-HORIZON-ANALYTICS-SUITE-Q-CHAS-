"""Spherically symmetric Einstein-scalar evolution in double-null coordinates (Phase IV).

Metric and matter:

    ds^2 = -alpha(u,v)^2 du dv + r(u,v)^2 dOmega^2,
    S_m  = int sqrt(-g) [ -(1/2)(grad phi)^2 - U(phi) ].

Field equations (each is a non-zero multiple of a component of G_mn - 8 pi G T_mn or of
the wave equation; verified symbolically by `symbolic_check`):

    r_uv      = -r_u r_v / r - alpha^2/(4 r) + 2 pi G alpha^2 r U(phi)          (E_uv)
    (ln a)_uv = r_u r_v / r^2 + alpha^2/(4 r^2) - 4 pi G phi_u phi_v           (E_thth)
    phi_uv    = -(r_u phi_v + r_v phi_u)/r - (alpha^2/4) U'(phi)               (wave eq.)
    r_uu - 2 (alpha_u/alpha) r_u + 4 pi G r phi_u^2 = 0                        (E_uu, constraint)
    r_vv - 2 (alpha_v/alpha) r_v + 4 pi G r phi_v^2 = 0                        (E_vv, constraint)

Misner-Sharp mass m = (r/2G)(1 + 4 r_u r_v / alpha^2).  With u and v increasing to the
future and v outward: normal region r_u < 0 < r_v, trapped r_u, r_v < 0, anti-trapped
r_u, r_v > 0.

Two characteristic problems are solved.

(a) `evolve`: regular centre.  Gauge u = v on r = 0, r = v/2 on the initial cone u = 0;
    the constraint E_vv gives (ln alpha)_v = 4 pi G r phi_v^2 there, alpha(0,0) = 1.
(b) `evolve_rect`: no centre; data on an outgoing ray u = 0 and an ingoing ray v = 0
    (used for anti-trapped, inflating data, which cannot contain a regular centre).

Scheme.  Second-order diamond scheme for X_uv = F(X, X_u, X_v) on a uniform grid
u_i = i h, v_j = j h: X_N = X_E + X_W - X_S + h^2 F(centre), N = (i,j), W = (i,j-1),
E = (i-1,j), S = (i-1,j-1), with F at the cell centre from a predictor-corrector
(three passes).  The scalar is evolved as psi = r phi,

    psi_uv = phi r_uv - (alpha^2/4) r U'(phi),

because the phi equation's -(r_u phi_v + r_v phi_u)/r term is unstable at the centre.
Centre: r = psi = 0 and alpha = 2 r_v (regularity, r_u r_v = -alpha^2/4 at r = 0).
Within `k_axis` points of the centre the 1/r^2 combination r_u r_v + alpha^2/4 =
(alpha^2/4)(2Gm/r) is replaced by its regular limit m = (4 pi/3) rho r^3 (exact in
de Sitter).  Points with r < r_exc in a trapped region (singularity) or r > r_inf in an
anti-trapped region (beyond future null infinity of an inflating region) are cut off,
together with the rest of their outgoing ray.

Accuracy (tests and `scripts/phase4_dynamics.py`): exact flat-space scalar evolution to
machine precision; second-order convergence of r, phi, ln alpha away from the centre;
r second order and ln alpha first order at the centre against exact de Sitter.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

try:
    from numba import njit
except ImportError:          # pragma: no cover - slow but correct fallback
    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda f: f


# per-ray diagnostic columns
DIAG = ("v_wall", "r_wall", "s_wall",
        "r_out", "ru_out", "rv_out", "s_out", "m_out",
        "r_in", "ru_in", "rv_in", "s_in", "m_in",
        "n_antitrapped", "n_trapped", "r_ah", "m_boundary", "j_end", "end_kind")
_ND = len(DIAG)


# ----------------------------------------------------------------- potential
@njit(cache=True)
def _U(phi, lam, vv, rho_f):
    x = phi / vv
    return lam / 4.0 * (phi * phi - vv * vv) ** 2 + rho_f * (0.5 - 0.75 * x + 0.25 * x ** 3)


@njit(cache=True)
def _dU(phi, lam, vv, rho_f):
    x = phi / vv
    return lam * phi * (phi * phi - vv * vv) - rho_f * (0.75 - 0.75 * x * x) / vv


@njit(cache=True)
def _rhs(r, ru, rv, s, su, sv, q, qu, qv, G, lam, vv, rho_f, near_axis):
    """Right-hand sides for (r, ln alpha, psi = r phi)."""
    a2 = math.exp(2.0 * s)
    p = q / r
    pu = (qu - p * ru) / r
    pv = (qv - p * rv) / r
    Uc = _U(p, lam, vv, rho_f)
    if near_axis:
        rho = Uc + (pu * pu + pv * pv) / a2
        Fr = -2.0 * math.pi * G / 3.0 * a2 * rho * r + 2.0 * math.pi * G * a2 * r * Uc
        Fs = 2.0 * math.pi * G / 3.0 * a2 * rho - 4.0 * math.pi * G * pu * pv
    else:
        Fr = -ru * rv / r - a2 / (4.0 * r) + 2.0 * math.pi * G * a2 * r * Uc
        Fs = ru * rv / (r * r) + a2 / (4.0 * r * r) - 4.0 * math.pi * G * pu * pv
    Fq = p * Fr - 0.25 * a2 * r * _dU(p, lam, vv, rho_f)
    return Fr, Fs, Fq


@njit(cache=True)
def _diamond(rE, rW, rS, sE, sW, sS, qE, qW, qS, h, G, lam, vv, rho_f, n_iter, near_axis):
    rN = rE + rW - rS; sN = sE + sW - sS; qN = qE + qW - qS
    Fs = 0.0
    for it in range(n_iter):
        rcen = 0.25 * (rN + rE + rW + rS)
        scen = 0.25 * (sN + sE + sW + sS)
        qcen = 0.25 * (qN + qE + qW + qS)
        ru = (rN + rW - rE - rS) / (2.0 * h); rv = (rN + rE - rW - rS) / (2.0 * h)
        su = (sN + sW - sE - sS) / (2.0 * h); sv = (sN + sE - sW - sS) / (2.0 * h)
        qu = (qN + qW - qE - qS) / (2.0 * h); qv = (qN + qE - qW - qS) / (2.0 * h)
        Fr, Fs, Fq = _rhs(rcen, ru, rv, scen, su, sv, qcen, qu, qv, G, lam, vv, rho_f, near_axis)
        rN = rE + rW - rS + h * h * Fr
        sN = sE + sW - sS + h * h * Fs
        qN = qE + qW - qS + h * h * Fq
    return rN, sN, qN, Fs


@njit(cache=True)
def _diag_ray(i, j0, jlast, rc, sc, qc, rp, valid_prev, h, G, vv, d):
    """Fill d[i, :] (columns DIAG) from the current ray (rc, sc, qc on j0..jlast) and the
    previous ray rp (for r_u, backward difference)."""
    jw = -1
    for j in range(jlast, j0 + 1, -1):           # outermost crossing of phi = 0
        if qc[j - 1] < 0.0 <= qc[j]:
            t = -qc[j - 1] / (qc[j] - qc[j - 1])
            d[i, 0] = (j - 1 + t) * h
            d[i, 1] = rc[j - 1] + t * (rc[j] - rc[j - 1])
            d[i, 2] = sc[j - 1] + t * (sc[j] - sc[j - 1])
            jw = j
            break
    if i > 0 and jw > 0:
        for j in range(jw, jlast):                # just outside: |phi - v| < 0.02 v
            if valid_prev[j] and valid_prev[j - 1] and abs(qc[j] / rc[j] - vv) < 0.02 * vv:
                ru = (rc[j] - rp[j]) / h
                rv = (rc[j + 1] - rc[j - 1]) / (2.0 * h)
                a2 = math.exp(2.0 * sc[j])
                d[i, 3] = rc[j]; d[i, 4] = ru; d[i, 5] = rv; d[i, 6] = sc[j]
                d[i, 7] = rc[j] / (2.0 * G) * (1.0 + 4.0 * ru * rv / a2)
                break
        for j in range(jw - 1, j0 + 1, -1):       # just inside: |phi + v| < 0.02 v
            if valid_prev[j] and valid_prev[j + 1] and abs(qc[j] / rc[j] + vv) < 0.02 * vv:
                ru = (rc[j] - rp[j]) / h
                rv = (rc[j + 1] - rc[j - 1]) / (2.0 * h)
                a2 = math.exp(2.0 * sc[j])
                d[i, 8] = rc[j]; d[i, 9] = ru; d[i, 10] = rv; d[i, 11] = sc[j]
                d[i, 12] = rc[j] / (2.0 * G) * (1.0 + 4.0 * ru * rv / a2)
                break
    if i > 0:
        nat = 0; ntr = 0
        for j in range(j0 + 2, jlast):
            if not valid_prev[j]:
                break
            ru = rc[j] - rp[j]
            rv = rc[j + 1] - rc[j - 1]
            if ru > 0.0 and rv > 0.0:
                nat += 1
            elif ru < 0.0 and rv < 0.0:
                ntr += 1
        d[i, 13] = nat; d[i, 14] = ntr
    for j in range(jlast - 1, j0 + 2, -1):        # outermost sign change of r_v
        if (rc[j + 1] - rc[j - 1]) * (rc[j] - rc[j - 2]) <= 0.0:
            d[i, 15] = rc[j]
            break


@njit(cache=True)
def _evolve_core(r0, s0, q0, rb, sb, qb, centre, h, n_rays, G, lam, vv, rho_f,
                 r_exc, r_inf, ks, kv, n_iter, k_axis):
    """Core loop for both characteristic problems.  centre=True: ray i covers j = i..N
    with a regular centre at j = i.  centre=False: ray i covers j = 0..N with (i, 0) given
    by (rb, sb, qb)."""
    N = r0.shape[0] - 1
    nr_out = n_rays // ks + 1
    nv_out = N // kv + 1
    R_out = np.full((nr_out, nv_out), np.nan)
    S_out = np.full((nr_out, nv_out), np.nan)
    P_out = np.full((nr_out, nv_out), np.nan)
    M_out = np.full((nr_out, nv_out), np.nan)
    d = np.full((n_rays + 1, _ND), np.nan)
    rp = r0.copy(); sp = s0.copy(); qp = q0.copy()
    valid_prev = np.ones(N + 1, dtype=np.bool_)
    rc = np.full(N + 1, np.nan); sc = np.full(N + 1, np.nan); qc = np.full(N + 1, np.nan)
    valid = np.zeros(N + 1, dtype=np.bool_)
    Fs_prev = np.zeros(N + 1); Fs_cur = np.zeros(N + 1)
    n_done = 0
    for i in range(0, n_rays + 1):
        j0 = i if centre else 0
        kind = 0
        if i == 0:
            for j in range(N + 1):
                rc[j] = r0[j]; sc[j] = s0[j]; qc[j] = q0[j]; valid[j] = True
        else:
            for j in range(N + 1):
                rc[j] = np.nan; sc[j] = np.nan; qc[j] = np.nan; valid[j] = False
            if centre:
                if not (valid_prev[i] and valid_prev[i - 1] and i < N):
                    break
                rc[i] = 0.0; qc[i] = 0.0
                sc[i] = 2.0 * sp[i] - sp[i - 1] + h * h * Fs_prev[i]
                valid[i] = True
                if valid_prev[i + 1]:
                    for corr in range(2):
                        rN, sN, qN, Fs = _diamond(rp[i + 1], rc[i], rp[i], sp[i + 1], sc[i], sp[i],
                                                  qp[i + 1], qc[i], qp[i], h, G, lam, vv, rho_f, n_iter, True)
                        rv_axis = (rN + rp[i]) / (2.0 * h)          # r(i, i-1) = -r(i-1, i)
                        if rv_axis > 0.0:
                            sc[i] = math.log(2.0 * rv_axis)
            else:
                rc[0] = rb[i]; sc[0] = sb[i]; qc[0] = qb[i]; valid[0] = True
            for j in range(j0 + 1, N + 1):
                if not (valid[j - 1] and valid_prev[j] and valid_prev[j - 1]):
                    break
                na = centre and (j - i <= k_axis)
                rN, sN, qN, Fs = _diamond(rp[j], rc[j - 1], rp[j - 1], sp[j], sc[j - 1], sp[j - 1],
                                          qp[j], qc[j - 1], qp[j - 1], h, G, lam, vv, rho_f, n_iter, na)
                ok = math.isfinite(rN) and math.isfinite(sN) and math.isfinite(qN) and abs(sN) < 600.0
                if not ok:
                    kind = 3
                elif rN < r_exc and rN < rc[j - 1]:
                    ok = False; kind = 1                      # singularity (trapped, r -> 0)
                elif rN > r_inf and rN > rc[j - 1] and rN > rp[j]:
                    ok = False; kind = 2                      # beyond future null infinity
                elif (j - j0 >= 2 and rc[j - 1] < rc[j - 2] and rc[j - 1] < rp[j - 1] and rN > rc[j - 1]):
                    # W is trapped (r_u, r_v < 0); with the null energy condition r_v/alpha^2 cannot
                    # increase along an outgoing ray (Raychaudhuri), so r_v > 0 here signals numerical
                    # breakdown near the singularity
                    ok = False; kind = 4
                if not ok:
                    break
                rc[j] = rN; sc[j] = sN; qc[j] = qN; valid[j] = True
                Fs_cur[j] = Fs
        jlast = j0
        for j in range(j0, N + 1):
            if valid[j]:
                jlast = j
            else:
                break
        d[i, 17] = jlast; d[i, 18] = kind
        _diag_ray(i, j0, jlast, rc, sc, qc, rp, valid_prev, h, G, vv, d)
        if i > 0 and jlast == N and valid_prev[N] and valid_prev[N - 1]:
            ru = (rc[N] - rp[N]) / h
            rv = (rc[N] - rc[N - 1]) / h
            d[i, 16] = rc[N] / (2.0 * G) * (1.0 + 4.0 * ru * rv / math.exp(2.0 * sc[N]))
        if i % ks == 0:
            io = i // ks
            for jo in range(nv_out):
                j = jo * kv
                if j >= j0 and j <= N and valid[j]:
                    R_out[io, jo] = rc[j]; S_out[io, jo] = sc[j]
                    if rc[j] > 0.0:
                        P_out[io, jo] = qc[j] / rc[j]
                    elif j + 1 <= N and valid[j + 1]:
                        P_out[io, jo] = qc[j + 1] / rc[j + 1]  # centre: limit of psi/r
                    if i > 0 and j > j0 and valid_prev[j]:
                        ru = (rc[j] - rp[j]) / h
                        rv = (rc[j] - rc[j - 1]) / h
                        M_out[io, jo] = rc[j] / (2.0 * G) * (1.0 + 4.0 * ru * rv / math.exp(2.0 * sc[j]))
        for j in range(N + 1):
            rp[j] = rc[j]; sp[j] = sc[j]; qp[j] = qc[j]; valid_prev[j] = valid[j]
            Fs_prev[j] = Fs_cur[j]
        n_done = i
        if centre and i > 0 and kind == 1 and jlast - i <= 3:
            break                                  # the singularity has reached the centre
        if not centre and not valid[0]:
            break
    return R_out, S_out, P_out, M_out, d, n_done


@dataclass
class Potential:
    """Quartic potential of thickwall.QuarticPotential: false vacuum -v (U = rho_f),
    true vacuum +v (U = 0)."""
    lam: float
    v: float
    rho_f: float

    @classmethod
    def from_quartic(cls, q):
        return cls(q.lam, q.v, q.rho_f)

    def U(self, phi):
        x = phi / self.v
        return self.lam / 4 * (phi**2 - self.v**2) ** 2 + self.rho_f * (0.5 - 0.75 * x + 0.25 * x**3)


def _pack(out, h, N, n_rays, ks, kv, G):
    R, S, P, M, d, n_done = out
    res = dict(r=R, s=S, phi=P, m=M, h=h, N=N, n_rays=n_rays, n_done=int(n_done), ks=ks, kv=kv, G=G,
               u=np.arange(R.shape[0]) * h * ks, v=np.arange(R.shape[1]) * h * kv,
               u_ray=np.arange(n_rays + 1) * h)
    for k, name in enumerate(DIAG):
        res[name] = d[:, k]
    return res


def initial_cone(profile, h: float, v_max: float, G: float = 1.0):
    """Data on u = 0 for problem (a): r = v/2, phi = profile(r), ln alpha from E_vv
    ((ln alpha)_v = 4 pi G r phi_v^2), trapezoid rule on a 4x finer grid."""
    N = int(round(v_max / h))
    vf = np.linspace(0.0, N * h, 4 * N + 1)
    rf = vf / 2
    pf = profile(rf)
    dpf = np.gradient(pf, vf, edge_order=2)
    integrand = 4 * np.pi * G * rf * dpf**2
    sf = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(vf))])
    return rf[::4].copy(), sf[::4].copy(), pf[::4].copy()


def initial_cone_alpha1(profile_v, h: float, v_max: float, G: float = 1.0):
    """Data on u = 0 for problem (a) in the gauge alpha = 1 on the cone (alpha(0,0) = 2 r_v = 1):
    phi = profile_v(v), r from the E_vv constraint r_vv = -4 pi G r phi_v^2, r(0) = 0,
    r_v(0) = 1/2.  Unlike r = v/2, this allows r_v to change sign on the cone (a wall with
    eps_out = -1, the bounce branch beyond the Einstein-Rosen bridge)."""
    from scipy.integrate import solve_ivp
    N = int(round(v_max / h))
    v = np.arange(N + 1) * h
    vf = np.linspace(0.0, N * h, 8 * N + 1)
    dp = np.gradient(profile_v(vf), vf, edge_order=2)

    def rhs(x, y):
        return [y[1], -4 * np.pi * G * y[0] * np.interp(x, vf, dp) ** 2]
    sol = solve_ivp(rhs, (0.0, N * h), [0.0, 0.5], t_eval=v, method="DOP853", rtol=1e-12, atol=1e-14,
                    max_step=h)
    r = sol.y[0]
    if len(r) < N + 1 or np.min(r[np.argmax(r):]) <= 2 * h:
        raise ValueError("the initial cone reaches r = 0 before v_max: reduce v_max")
    return r, np.zeros(N + 1), profile_v(v)


def evolve_cone(r0, s0, p0, pot: Potential, h: float, u_max: float = None, G: float = 1.0,
                r_exc: float = None, ks: int = 1, kv: int = 1, n_iter: int = 3, k_axis: int = 2) -> dict:
    """Problem (a) with explicitly supplied cone data (r, ln alpha, phi on u = 0, r(0) = 0)."""
    N = len(r0) - 1
    n_rays = N - 1 if u_max is None else min(N - 1, int(round(u_max / h)))
    r_exc = r_exc if r_exc is not None else 2 * h
    dummy = np.zeros(1)
    out = _evolve_core(r0, s0, r0 * p0, dummy, dummy, dummy, True, h, n_rays, G, pot.lam, pot.v, pot.rho_f,
                       r_exc, 1e300, ks, kv, n_iter, k_axis)
    return _pack(out, h, N, n_rays, ks, kv, G)


def evolve(profile, pot: Potential, h: float, v_max: float, u_max: float = None, G: float = 1.0,
           r_exc: float = None, ks: int = 1, kv: int = 1, n_iter: int = 3, k_axis: int = 2) -> dict:
    """Problem (a): regular centre, data phi = profile(r) on the initial cone."""
    r0, s0, p0 = initial_cone(profile, h, v_max, G)
    N = len(r0) - 1
    n_rays = N - 1 if u_max is None else min(N - 1, int(round(u_max / h)))
    r_exc = r_exc if r_exc is not None else 2 * h
    dummy = np.zeros(1)
    out = _evolve_core(r0, s0, r0 * p0, dummy, dummy, dummy, True, h, n_rays, G, pot.lam, pot.v, pot.rho_f,
                       r_exc, 1e300, ks, kv, n_iter, k_axis)
    return _pack(out, h, N, n_rays, ks, kv, G)


def anti_trapped_data(profile_v, chi: float, r0: float, h: float, v_max: float, u_max: float,
                      G: float = 1.0, rv0: float = 0.5):
    """Data for problem (b) in the gauge alpha = 1 on both initial rays.

    * u = 0 (outgoing): phi = profile_v(v); r from the E_vv constraint r_vv = -4 pi G r phi_v^2
      with r(0) = r0, r_v(0) = rv0.  The wall focuses the ray, so r_v drops across it and can
      change sign (an exterior with eps_out = -1, beyond the Einstein-Rosen bridge).
    * v = 0 (ingoing): pure false vacuum phi = profile_v(0), exact de Sitter beyond the
      cosmological horizon (r0 > 1/chi).  E_uu with alpha = 1 gives r_uu = 0, so r = r0 + r_u u,
      and the Misner-Sharp mass at the corner fixes r_u = (chi^2 r0^2 - 1)/(4 rv0) > 0."""
    from scipy.integrate import solve_ivp
    N = int(round(v_max / h))
    n_rays = int(round(u_max / h))
    v = np.arange(N + 1) * h
    vf = np.linspace(0.0, N * h, 8 * N + 1)
    dp = np.gradient(profile_v(vf), vf, edge_order=2)

    def rhs(x, y):
        return [y[1], -4 * np.pi * G * y[0] * np.interp(x, vf, dp) ** 2]

    def hit_zero(x, y):
        return y[0] - 2 * h
    hit_zero.terminal = True
    sol = solve_ivp(rhs, (0.0, N * h), [r0, rv0], t_eval=v, method="DOP853", rtol=1e-11, atol=1e-13,
                    max_step=h, events=hit_zero)
    if sol.y.shape[1] < N + 1:
        raise ValueError("the initial outgoing ray reaches r = 0 before v_max: reduce v_max")
    r_a = sol.y[0]
    p_a = profile_v(v)
    ru0 = (chi**2 * r0**2 - 1) / (4 * rv0)
    ui = np.arange(n_rays + 1) * h
    r_b = r0 + ru0 * ui
    return (r_a, np.zeros(N + 1), r_a * p_a), (r_b, np.zeros_like(ui), r_b * p_a[0])


def evolve_rect(data_u0, data_v0, pot: Potential, h: float, G: float = 1.0, r_exc: float = None,
                r_inf: float = 1e4, ks: int = 1, kv: int = 1, n_iter: int = 3) -> dict:
    """Problem (b): data on u = 0 (data_u0 = (r, ln alpha, psi) for j = 0..N) and on
    v = 0 (data_v0 for i = 0..n_rays)."""
    (r0, s0, q0), (rb, sb, qb) = data_u0, data_v0
    N, n_rays = len(r0) - 1, len(rb) - 1
    r_exc = r_exc if r_exc is not None else 2 * h
    out = _evolve_core(r0, s0, q0, rb, sb, qb, False, h, n_rays, G, pot.lam, pot.v, pot.rho_f,
                       r_exc, r_inf, ks, kv, n_iter, 0)
    return _pack(out, h, N, n_rays, ks, kv, G)


# ------------------------------------------------------------ exact solutions
def de_sitter_exact(chi: float, u, v):
    """Exact de Sitter in the gauge of problem (a) (centre at u = v, r = v/2 on u = 0):
    V(v) = (2/chi) artanh(chi v / 2), r = tanh(chi (V - U)/2)/chi,
    alpha^2 = (1 - chi^2 r^2) U'(u) V'(v), U' = 1/(1 - chi^2 u^2/4).  Valid for v < 2/chi."""
    U = 2 / chi * np.arctanh(chi * u / 2)
    V = 2 / chi * np.arctanh(chi * v / 2)
    r = np.tanh(chi * (V - U) / 2) / chi
    a2 = (1 - chi**2 * r**2) / ((1 - chi**2 * u**2 / 4) * (1 - chi**2 * v**2 / 4))
    return r, 0.5 * np.log(a2)


# ------------------------------------------------------------ wall diagnostics
def wall_trajectory(res: dict) -> dict:
    """Wall worldline (phi = 0 crossing) with proper time d tau = alpha sqrt(dv_w/du) du and
    the sign functions A = n^mu d_mu r = -udot r_u + vdot r_v just outside (A_out) and inside
    (A_in) the wall, with the outward unit normal n = -udot d_u + vdot d_v, udot =
    1/(alpha sqrt(dv_w/du)), vdot = sqrt(dv_w/du)/alpha, alpha taken at the evaluation point.
    Also f = 1 - 2Gm/r = -4 r_u r_v/alpha^2 there.  Thin-wall limit: A_in - A_out = kappa R,
    A^2 = f + Rdot^2, eps = sgn A."""
    u_all = res["u_ray"]
    vw, rw, sw = res["v_wall"], res["r_wall"], res["s_wall"]
    idx = np.nonzero(np.isfinite(vw) & np.isfinite(rw))[0]
    if len(idx) < 5:
        return dict(u=u_all[idx], r=rw[idx], tau=np.zeros(len(idx)))
    br = np.nonzero(np.diff(idx) > 1)[0]
    if len(br):
        idx = idx[:br[0] + 1]                       # leading contiguous stretch
    u, v_, r_, s_ = u_all[idx], vw[idx], rw[idx], sw[idx]
    dvdu = np.gradient(v_, u)
    a = np.exp(s_)
    dtau = a * np.sqrt(np.clip(dvdu, 0, None))
    tau = np.concatenate([[0.0], np.cumsum(0.5 * (dtau[1:] + dtau[:-1]) * np.diff(u))])
    with np.errstate(divide="ignore", invalid="ignore"):
        sq = np.sqrt(dvdu)
        # normals at the evaluation points, normalised with the local lapse (alpha varies
        # across a thick wall), coordinate velocity of the wall
        a_out, a_in = np.exp(res["s_out"][idx]), np.exp(res["s_in"][idx])
        A_out = (-res["ru_out"][idx] / sq + sq * res["rv_out"][idx]) / a_out
        A_in = (-res["ru_in"][idx] / sq + sq * res["rv_in"][idx]) / a_in
        f_out = -4 * res["ru_out"][idx] * res["rv_out"][idx] / a_out**2
        f_in = -4 * res["ru_in"][idx] * res["rv_in"][idx] / a_in**2
        Rdot = np.gradient(r_, tau)
        # the wall is resolved when both evaluation points were found close to it
        resolved = (np.abs(res["r_out"][idx] - r_) < 0.25 * r_) & (np.abs(r_ - res["r_in"][idx]) < 0.25 * r_)
    return dict(u=u, v=v_, r=r_, tau=tau, dvdu=dvdu, A_out=A_out, A_in=A_in, Rdot=Rdot, f_out=f_out, f_in=f_in,
                resolved=resolved,
                m_out=res["m_out"][idx], m_in=res["m_in"][idx], r_out=res["r_out"][idx], r_in=res["r_in"][idx])


# ------------------------------------------------------------ symbolic check
def symbolic_check() -> dict:
    """Verify with sympy that the five equations above are non-zero multiples of the
    components of G_mn - 8 pi G T_mn (and of box phi - U'(phi)) for the metric
    -alpha^2 du dv + r^2 dOmega^2.  Returns the ratios (functions of r, alpha only)."""
    import sympy as sp
    u, v, th, ph, G = sp.symbols("u v theta phi_ G", real=True)
    a = sp.Function("alpha")(u, v); r = sp.Function("r")(u, v); f = sp.Function("f")(u, v)
    Uf = sp.Function("U")
    X = [u, v, th, ph]
    g = sp.zeros(4)
    g[0, 1] = g[1, 0] = -a**2 / 2
    g[2, 2] = r**2; g[3, 3] = r**2 * sp.sin(th) ** 2
    gi = g.inv()
    Ga = [[[sp.simplify(sum(gi[l, k] * (sp.diff(g[k, m], X[n]) + sp.diff(g[k, n], X[m]) - sp.diff(g[m, n], X[k]))
                             for k in range(4)) / 2) for n in range(4)] for m in range(4)] for l in range(4)]

    def Ric(m, n):
        return sp.simplify(sum(sp.diff(Ga[l][m][n], X[l]) - sp.diff(Ga[l][m][l], X[n])
                               + sum(Ga[l][l][k] * Ga[k][m][n] - Ga[l][n][k] * Ga[k][m][l] for k in range(4))
                               for l in range(4)))
    Rm = sp.Matrix(4, 4, lambda m, n: Ric(m, n))
    Rs = sp.simplify(sum(gi[m, n] * Rm[m, n] for m in range(4) for n in range(4)))
    dphi = [sp.diff(f, x) for x in X]
    kin = sum(gi[m, n] * dphi[m] * dphi[n] for m in range(4) for n in range(4))
    T = sp.Matrix(4, 4, lambda m, n: dphi[m] * dphi[n] - g[m, n] * (kin / 2 + Uf(f)))
    E = sp.simplify(Rm - Rs * g / 2 - 8 * sp.pi * G * T)
    ru, rv, ruv = sp.diff(r, u), sp.diff(r, v), sp.diff(r, u, v)
    fu, fv = sp.diff(f, u), sp.diff(f, v)
    eq_uu = sp.diff(r, u, 2) - 2 * sp.diff(a, u) / a * ru + 4 * sp.pi * G * r * fu**2
    eq_vv = sp.diff(r, v, 2) - 2 * sp.diff(a, v) / a * rv + 4 * sp.pi * G * r * fv**2
    eq_r = ruv + ru * rv / r + a**2 / (4 * r) - 2 * sp.pi * G * a**2 * r * Uf(f)
    eq_a = sp.diff(sp.log(a), u, v) - (ru * rv / r**2 + a**2 / (4 * r**2) - 4 * sp.pi * G * fu * fv)
    sol = sp.solve(eq_r, ruv)[0]
    sqrtg = sp.sqrt(-g.det())
    box = sp.simplify(sum(sp.diff(sqrtg * gi[m, n] * dphi[n], X[m]) for m in range(4) for n in range(4)) / sqrtg)
    dU = sp.Function("dU")(f)
    eq_p = sp.diff(f, u, v) + (ru * fv + rv * fu) / r + a**2 / 4 * dU
    ratios = dict(
        uu=sp.simplify(E[0, 0] / eq_uu), vv=sp.simplify(E[1, 1] / eq_vv), uv=sp.simplify(E[0, 1] / eq_r),
        thth=sp.simplify(E[2, 2].subs(ruv, sol) / eq_a.subs(ruv, sol)),
        wave=sp.simplify((box - dU) / eq_p))
    return {k: str(v_) for k, v_ in ratios.items()}
