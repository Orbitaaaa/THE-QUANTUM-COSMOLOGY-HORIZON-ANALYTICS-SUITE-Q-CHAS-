"""Round 10 (open problems P-R6n-a, P-R6l-a, P-R6l-c): a characteristic (double-null) code for model U with a massless scalar field.

Metric ds^2 = -2 e^{2 sigma(u,v)} du dv + r(u,v)^2 dOmega_{D-2}^2, scalar phi(u,v), units alpha = l^2 = 1 and 8 pi G_D/(D-2) = 1 (so
T_uu = phi_u^2 enters with unit coupling).  With psi = (1 + 2 e^{-2 sigma} r_u r_v)/r^2 and the tower function h(psi) (generalised
Misner-Sharp mass mu = r^{D-1} h(psi)) the field equations are

    constraints   h'(psi) (r_uu - 2 sigma_u r_u) = -r phi_u^2,      h'(psi) (r_vv - 2 sigma_v r_v) = -r phi_v^2
    evolution     r_uv     = (r e^{2 sigma}/2) (2 psi - (D-1) h/h')                       (T_uv = 0 for a massless scalar)
                  phi_uv   = -((D-2)/(2r)) (r_u phi_v + r_v phi_u)
                  sigma_uv = e^{-2 sigma}/(4 h'^3 r^2) [ -(D-1)^2 h^2 h'' r^2 e^{4 sigma} + 2 h'^3 (2 r_u r_v + e^{2 sigma}) e^{2 sigma}
                             + h'^2 r^2 e^{2 sigma} ((D-1)(D-4) h e^{2 sigma} - 2 (D-2) phi_u phi_v) + 4 h'' r^2 phi_u^2 phi_v^2 ]

The uu constraint is the null projection of the first law mu_u = -2 r^{D-2} e^{-2 sigma} r_v phi_u^2; r_uv follows from it; sigma_uv is the
condition that the evolution preserves the constraint (derived with sympy, ``derive_sigma_uv``).  Checks: GR (h = psi) gives the
D-dimensional Einstein-scalar system; the static model-U hole satisfies it exactly (sigma_uv = -f f''/8).

Curvature: K = R_2^2 + 4 (D-2) (nabla nabla r)^2 / r^2 + 2 (D-2)(D-3) psi^2, with R_2 = 4 e^{-2 sigma} sigma_uv and
(nabla nabla r)^2 = 2 e^{-4 sigma} [(r_uu - 2 sigma_u r_u)(r_vv - 2 sigma_v r_v) + r_uv^2], all evaluated algebraically on shell.

Numerics: second-order predictor-corrector "diamond" scheme on a uniform (u, v) grid (Hamade-Stewart / Burko-Ori type), compiled with
numba.  Initial data: ingoing ray v = 0 in the r-gauge r = r0 - u/2 with phi = 0 (sigma constant, matched to the static hole of mass mu0);
outgoing ray u = 0 with r = r0 + v/2 and a scalar pulse; sigma and r_u along it from the vv constraint and the r_uv equation (RK4).
"""

from __future__ import annotations

import numpy as np
from numba import njit

from . import unified6l as U

KIND = {"gr": 0, "hayward": 1, "dymnikova": 2}


@njit(cache=True, error_model="numpy")
def _h012(psi, kind):
    if kind == 0:
        return psi, 1.0, 0.0
    x = 1.0 - psi
    if x < 1e-14:
        x = 1e-14
    if kind == 1:
        return psi / x, 1.0 / x**2, 2.0 / x**3
    return -np.log(x), 1.0 / x, 1.0 / x**2


@njit(cache=True, error_model="numpy")
def _rhs(r, s, ru, rv, pu, pv, D, kind):
    """r_uv, sigma_uv, phi_uv and psi from the local state."""
    E = np.exp(2.0 * s)
    psi = (1.0 + 2.0 * ru * rv / E) / r**2
    h0, h1, h2 = _h012(psi, kind)
    F = 0.5 * r * E * (2.0 * psi - (D - 1.0) * h0 / h1)
    W = -(D - 2.0) / (2.0 * r) * (ru * pv + rv * pu)
    S = (-(D - 1.0) ** 2 * h0**2 * h2 * r**2 * E**2 + 2.0 * h1**3 * (2.0 * ru * rv + E) * E
         + h1**2 * r**2 * E * ((D - 1.0) * (D - 4.0) * h0 * E - 2.0 * (D - 2.0) * pu * pv) + 4.0 * h2 * r**2 * pu**2 * pv**2) / (E * 4.0 * h1**3 * r**2)
    return F, S, W, psi


@njit(cache=True, error_model="numpy")
def _evolve(r, s, p, uu, dv, D, kind, n_iter):
    nu, nv = r.shape
    for i in range(1, nu):
        du = uu[i] - uu[i - 1]
        for j in range(1, nv):
            rA, rB, rC = r[i - 1, j - 1], r[i - 1, j], r[i, j - 1]
            sA, sB, sC = s[i - 1, j - 1], s[i - 1, j], s[i, j - 1]
            pA, pB, pC = p[i - 1, j - 1], p[i - 1, j], p[i, j - 1]
            if rB <= 0.0 or rC <= 0.0 or not np.isfinite(rB + rC + sB + sC):
                r[i, j] = np.nan; s[i, j] = np.nan; p[i, j] = np.nan
                continue
            # predictor: N from the flat-cell rule
            rN = rB + rC - rA; sN = sB + sC - sA; pN = pB + pC - pA
            for _ in range(n_iter):
                rm = 0.25 * (rA + rB + rC + rN); sm = 0.25 * (sA + sB + sC + sN)
                ru = 0.5 * ((rC - rA) + (rN - rB)) / du; rv = 0.5 * ((rB - rA) + (rN - rC)) / dv
                pu = 0.5 * ((pC - pA) + (pN - pB)) / du; pv = 0.5 * ((pB - pA) + (pN - pC)) / dv
                F, S, W, psi = _rhs(rm, sm, ru, rv, pu, pv, D, kind)
                rN = rB + rC - rA + du * dv * F
                sN = sB + sC - sA + du * dv * S
                pN = pB + pC - pA + du * dv * W
            r[i, j] = rN; s[i, j] = sN; p[i, j] = pN
    return r, s, p


@njit(cache=True, error_model="numpy")
def _diag_point(r, s, ru, rv, pu, pv, D, kind):
    """mu, psi, K and the trapped flag from the local state (on shell)."""
    E = np.exp(2.0 * s)
    psi = (1.0 + 2.0 * ru * rv / E) / r**2
    h0, h1, h2 = _h012(psi, kind)
    F, S, W, _ = _rhs(r, s, ru, rv, pu, pv, D, kind)
    Auu = -r * pu**2 / h1; Avv = -r * pv**2 / h1
    nnr2 = 2.0 * np.exp(-4.0 * s) * (Auu * Avv + F**2)
    R2 = 4.0 * np.exp(-2.0 * s) * S
    K = R2**2 + 4.0 * (D - 2.0) * nnr2 / r**2 + 2.0 * (D - 2.0) * (D - 3.0) * psi**2
    return r ** (D - 1.0) * h0, psi, K


@njit(cache=True, error_model="numpy")
def _evolve_rows(r_row0, s_row0, p_row0, r_col0, s_col0, p_col0, uu, dv, D, kind, n_iter, su, sv):
    """Row-by-row march keeping two rows; returns r, sigma, phi, mu, psi, K on the sub-grid (every su-th row, sv-th column)."""
    nu = len(uu); nv = len(r_row0)
    nsu = (nu - 1) // su + 1; nsv = (nv - 1) // sv + 1
    out = np.full((6, nsu, nsv), np.nan)
    rp_ = r_row0.copy(); sp_ = s_row0.copy(); pp_ = p_row0.copy()
    rn = np.empty(nv); sn = np.empty(nv); pn = np.empty(nv)
    for jj in range(nsv):
        j = jj * sv
        out[0, 0, jj] = rp_[j]; out[1, 0, jj] = sp_[j]; out[2, 0, jj] = pp_[j]
    for i in range(1, nu):
        du = uu[i] - uu[i - 1]
        rn[0] = r_col0[i]; sn[0] = s_col0[i]; pn[0] = p_col0[i]
        for j in range(1, nv):
            rA, rB, rC = rp_[j - 1], rp_[j], rn[j - 1]
            sA, sB, sC = sp_[j - 1], sp_[j], sn[j - 1]
            pA, pB, pC = pp_[j - 1], pp_[j], pn[j - 1]
            if not (rB > 0.0 and rC > 0.0 and np.isfinite(rB + rC + sB + sC + pB + pC)):
                rn[j] = np.nan; sn[j] = np.nan; pn[j] = np.nan
                continue
            rN = rB + rC - rA; sN = sB + sC - sA; pN = pB + pC - pA
            for _ in range(n_iter):
                rm = 0.25 * (rA + rB + rC + rN); sm = 0.25 * (sA + sB + sC + sN)
                ru = 0.5 * ((rC - rA) + (rN - rB)) / du; rv = 0.5 * ((rB - rA) + (rN - rC)) / dv
                pu = 0.5 * ((pC - pA) + (pN - pB)) / du; pv = 0.5 * ((pB - pA) + (pN - pC)) / dv
                F, S, W, psi = _rhs(rm, sm, ru, rv, pu, pv, D, kind)
                rN = rB + rC - rA + du * dv * F
                sN = sB + sC - sA + du * dv * S
                pN = pB + pC - pA + du * dv * W
            rn[j] = rN; sn[j] = sN; pn[j] = pN
        if i % su == 0:
            ii = i // su
            for jj in range(nsv):
                j = jj * sv
                out[0, ii, jj] = rn[j]; out[1, ii, jj] = sn[j]; out[2, ii, jj] = pn[j]
                if j >= 1:
                    rA, rB, rC, rN = rp_[j - 1], rp_[j], rn[j - 1], rn[j]
                    sm = 0.25 * (sp_[j - 1] + sp_[j] + sn[j - 1] + sn[j])
                    ru = 0.5 * ((rC - rA) + (rN - rB)) / du; rv = 0.5 * ((rB - rA) + (rN - rC)) / dv
                    pu = 0.5 * ((pn[j - 1] - pp_[j - 1]) + (pn[j] - pp_[j])) / du
                    pv = 0.5 * ((pp_[j] - pp_[j - 1]) + (pn[j] - pn[j - 1])) / dv
                    m, ps, K = _diag_point(0.25 * (rA + rB + rC + rN), sm, ru, rv, pu, pv, D, kind)
                    out[3, ii, jj] = m; out[4, ii, jj] = ps; out[5, ii, jj] = K
        for j in range(nv):
            rp_[j] = rn[j]; sp_[j] = sn[j]; pp_[j] = pn[j]
    return out


@njit(cache=True, error_model="numpy")
def _hinv(y, kind):
    if kind == 0:
        return y
    if kind == 1:
        return y / (1.0 + y)
    return -np.expm1(-y)


@njit(cache=True, error_model="numpy")
def _rhs_mu(r, mu, ru, rv, pu, pv, D, kind):
    """(r, mu, phi) formulation: psi = h^{-1}(mu/r^{D-1}) < 1/alpha automatically; e^{2 sigma} = -2 r_u r_v / f, f = 1 - r^2 psi.
    r_uv = -(r r_u r_v / f)(2 psi - (D-1) h/h'),  mu_u = r^{D-2} f phi_u^2 / r_u,  mu_v = r^{D-2} f phi_v^2 / r_v."""
    psi = _hinv(mu / r ** (D - 1.0), kind)
    h0, h1, h2 = _h012(psi, kind)
    f = 1.0 - r**2 * psi
    F = -(r * ru * rv / f) * (2.0 * psi - (D - 1.0) * h0 / h1)
    W = -(D - 2.0) / (2.0 * r) * (ru * pv + rv * pu)
    rD2 = r ** (D - 2.0)
    mu_u = rD2 * f * pu**2 / ru if ru != 0.0 else 0.0
    mu_v = rD2 * f * pv**2 / rv if rv != 0.0 else 0.0
    return F, W, mu_u, mu_v, psi, f


@njit(cache=True, error_model="numpy")
def _evolve_rows_mu(r_row0, m_row0, p_row0, r_col0, m_col0, p_col0, uu, dv, D, kind, n_iter, su, sv):
    """Row march in (r, mu, phi); returns r, mu, phi, psi, K, f on the sub-grid."""
    nu = len(uu); nv = len(r_row0)
    nsu = (nu - 1) // su + 1; nsv = (nv - 1) // sv + 1
    out = np.full((6, nsu, nsv), np.nan)
    rp_ = r_row0.copy(); mp_ = m_row0.copy(); pp_ = p_row0.copy()
    rn = np.empty(nv); mn = np.empty(nv); pn = np.empty(nv)
    for jj in range(nsv):
        j = jj * sv
        out[0, 0, jj] = rp_[j]; out[1, 0, jj] = mp_[j]; out[2, 0, jj] = pp_[j]
    for i in range(1, nu):
        du = uu[i] - uu[i - 1]
        rn[0] = r_col0[i]; mn[0] = m_col0[i]; pn[0] = p_col0[i]
        for j in range(1, nv):
            rA, rB, rC = rp_[j - 1], rp_[j], rn[j - 1]
            mA, mB, mC = mp_[j - 1], mp_[j], mn[j - 1]
            pA, pB, pC = pp_[j - 1], pp_[j], pn[j - 1]
            if not (rB > 0.0 and rC > 0.0 and np.isfinite(rB + rC + mB + mC + pB + pC)):
                rn[j] = np.nan; mn[j] = np.nan; pn[j] = np.nan
                continue
            rN = rB + rC - rA; pN = pB + pC - pA; mN = 0.5 * (mB + mC)
            for _ in range(n_iter):
                rm = 0.25 * (rA + rB + rC + rN); mm = 0.25 * (mA + mB + mC + mN)
                ru = 0.5 * ((rC - rA) + (rN - rB)) / du; rv = 0.5 * ((rB - rA) + (rN - rC)) / dv
                pu = 0.5 * ((pC - pA) + (pN - pB)) / du; pv = 0.5 * ((pB - pA) + (pN - pC)) / dv
                F, W, mu_u, mu_v, psi, f = _rhs_mu(rm, mm, ru, rv, pu, pv, D, kind)
                rN = rB + rC - rA + du * dv * F
                pN = pB + pC - pA + du * dv * W
                wu = ru * ru; wv = rv * rv
                mN = (wu * (mB + du * mu_u) + wv * (mC + dv * mu_v)) / (wu + wv) if (wu + wv) > 0.0 else 0.5 * (mB + mC)
            rn[j] = rN; mn[j] = mN; pn[j] = pN
        if i % su == 0:
            ii = i // su
            for jj in range(nsv):
                j = jj * sv
                out[0, ii, jj] = rn[j]; out[1, ii, jj] = mn[j]; out[2, ii, jj] = pn[j]
                if j >= 1:
                    rA, rB, rC, rN = rp_[j - 1], rp_[j], rn[j - 1], rn[j]
                    rm = 0.25 * (rA + rB + rC + rN); mm = 0.25 * (mp_[j - 1] + mp_[j] + mn[j - 1] + mn[j])
                    ru = 0.5 * ((rC - rA) + (rN - rB)) / du; rv = 0.5 * ((rB - rA) + (rN - rC)) / dv
                    pu = 0.5 * ((pn[j - 1] - pp_[j - 1]) + (pn[j] - pp_[j])) / du
                    pv = 0.5 * ((pp_[j] - pp_[j - 1]) + (pn[j] - pn[j - 1])) / dv
                    F, W, mu_u, mu_v, psi, f = _rhs_mu(rm, mm, ru, rv, pu, pv, D, kind)
                    E = -2.0 * ru * rv / f
                    if E > 0.0:
                        _, _, K = _diag_point(rm, 0.5 * np.log(E), ru, rv, pu, pv, D, kind)
                    else:
                        K = np.nan
                    out[3, ii, jj] = psi; out[4, ii, jj] = K; out[5, ii, jj] = f
        for j in range(nv):
            rp_[j] = rn[j]; mp_[j] = mn[j]; pp_[j] = pn[j]
    return out


@njit(cache=True, error_model="numpy")
def _rhs_hyb(r, s, mu, ru, rv, pu, pv, D, kind):
    """Hybrid: psi from the evolved mass (bounded by 1/alpha), metric (r, sigma) evolved with it; first-law fluxes without divisions."""
    E = np.exp(2.0 * s)
    psi = _hinv(mu / r ** (D - 1.0), kind)
    h0, h1, h2 = _h012(psi, kind)
    F = 0.5 * r * E * (2.0 * psi - (D - 1.0) * h0 / h1)
    W = -(D - 2.0) / (2.0 * r) * (ru * pv + rv * pu)
    S = (-(D - 1.0) ** 2 * h0**2 * h2 * r**2 * E**2 + 2.0 * h1**3 * r**2 * psi * E**2
         + h1**2 * r**2 * E * ((D - 1.0) * (D - 4.0) * h0 * E - 2.0 * (D - 2.0) * pu * pv) + 4.0 * h2 * r**2 * pu**2 * pv**2) / (E * 4.0 * h1**3 * r**2)
    rD2 = r ** (D - 2.0)
    mu_u = -2.0 * rD2 * rv * pu**2 / E
    mu_v = -2.0 * rD2 * ru * pv**2 / E
    return F, S, W, mu_u, mu_v, psi


@njit(cache=True, error_model="numpy")
def _evolve_rows_hyb(r_row0, s_row0, m_row0, p_row0, r_col0, s_col0, m_col0, p_col0, uu, dv, D, kind, n_iter, su, sv):
    """Row march in (r, sigma, mu, phi); sub-grid output: r, sigma, mu, phi, psi(mu), psi(metric), K."""
    nu = len(uu); nv = len(r_row0)
    nsu = (nu - 1) // su + 1; nsv = (nv - 1) // sv + 1
    out = np.full((9, nsu, nsv), np.nan)
    rp_ = r_row0.copy(); sp_ = s_row0.copy(); mp_ = m_row0.copy(); pp_ = p_row0.copy()
    rn = np.empty(nv); sn = np.empty(nv); mn = np.empty(nv); pn = np.empty(nv)
    for jj in range(nsv):
        j = jj * sv
        out[0, 0, jj] = rp_[j]; out[1, 0, jj] = sp_[j]; out[2, 0, jj] = mp_[j]; out[3, 0, jj] = pp_[j]
    for i in range(1, nu):
        du = uu[i] - uu[i - 1]
        rn[0] = r_col0[i]; sn[0] = s_col0[i]; mn[0] = m_col0[i]; pn[0] = p_col0[i]
        for j in range(1, nv):
            rA, rB, rC = rp_[j - 1], rp_[j], rn[j - 1]
            sA, sB, sC = sp_[j - 1], sp_[j], sn[j - 1]
            mA, mB, mC = mp_[j - 1], mp_[j], mn[j - 1]
            pA, pB, pC = pp_[j - 1], pp_[j], pn[j - 1]
            if not (rB > 0.0 and rC > 0.0 and np.isfinite(rB + rC + sB + sC + mB + mC + pB + pC)):
                rn[j] = np.nan; sn[j] = np.nan; mn[j] = np.nan; pn[j] = np.nan
                continue
            rN = rB + rC - rA; sN = sB + sC - sA; pN = pB + pC - pA; mN = mB + mC - mA
            for _ in range(n_iter):
                rm = 0.25 * (rA + rB + rC + rN); sm = 0.25 * (sA + sB + sC + sN); mm = 0.25 * (mA + mB + mC + mN)
                ru = 0.5 * ((rC - rA) + (rN - rB)) / du; rv = 0.5 * ((rB - rA) + (rN - rC)) / dv
                pu = 0.5 * ((pC - pA) + (pN - pB)) / du; pv = 0.5 * ((pB - pA) + (pN - pC)) / dv
                F, S, W, mu_u, mu_v, psi = _rhs_hyb(rm, sm, mm, ru, rv, pu, pv, D, kind)
                rN = rB + rC - rA + du * dv * F
                sN = sB + sC - sA + du * dv * S
                pN = pB + pC - pA + du * dv * W
                mN = 0.5 * ((mB + du * mu_u) + (mC + dv * mu_v))
            rn[j] = rN; sn[j] = sN; mn[j] = mN; pn[j] = pN
        if i % su == 0:
            ii = i // su
            for jj in range(nsv):
                j = jj * sv
                out[0, ii, jj] = rn[j]; out[1, ii, jj] = sn[j]; out[2, ii, jj] = mn[j]; out[3, ii, jj] = pn[j]
                if j >= 1:
                    rA, rB, rC, rN = rp_[j - 1], rp_[j], rn[j - 1], rn[j]
                    rm = 0.25 * (rA + rB + rC + rN); sm = 0.25 * (sp_[j - 1] + sp_[j] + sn[j - 1] + sn[j])
                    mm = 0.25 * (mp_[j - 1] + mp_[j] + mn[j - 1] + mn[j])
                    ru = 0.5 * ((rC - rA) + (rN - rB)) / du; rv = 0.5 * ((rB - rA) + (rN - rC)) / dv
                    pu = 0.5 * ((pn[j - 1] - pp_[j - 1]) + (pn[j] - pp_[j])) / du
                    pv = 0.5 * ((pp_[j] - pp_[j - 1]) + (pn[j] - pn[j - 1])) / dv
                    F, S, W, mu_u, mu_v, psi = _rhs_hyb(rm, sm, mm, ru, rv, pu, pv, D, kind)
                    h0, h1, h2 = _h012(psi, kind)
                    E = np.exp(2.0 * sm)
                    Auu = -rm * pu**2 / h1; Avv = -rm * pv**2 / h1
                    T1 = (4.0 * S / E) ** 2
                    T2 = 4.0 * (D - 2.0) * 2.0 * (Auu * Avv + F**2) / E**2 / rm**2
                    K = T1 + T2 + 2.0 * (D - 2.0) * (D - 3.0) * psi**2
                    out[4, ii, jj] = psi; out[5, ii, jj] = (1.0 + 2.0 * ru * rv / E) / rm**2; out[6, ii, jj] = K
                    out[7, ii, jj] = T1; out[8, ii, jj] = T2
        for j in range(nv):
            rp_[j] = rn[j]; sp_[j] = sn[j]; mp_[j] = mn[j]; pp_[j] = pn[j]
    return out


def run_hyb(D=5, mu0=10.0, kind="hayward", r0=None, u_max=None, v_max=60.0, nu=4000, nv=12000, amp=0.0, v1=1.0, v2=5.0,
            amp2=0.0, v3=None, v4=None, n_iter=3, u_eh=None, n_cluster=0, n_out_u=400, n_out_v=600, tail_amp=0.0, tail_p=3.0) -> dict:
    """Hybrid (r, sigma, mu, phi) evolution (the production solver).  Pulses as in ``run``; optional Price-tail influx
    phi ~ tail_amp (v/v2)^(-tail_p) for v > v2 on the initial outgoing ray."""
    k = KIND[kind]
    rm_, rp = U.static_horizons(mu0, D, 1.0, kind) if kind != "gr" else (0.0, mu0 ** (1.0 / (D - 3)))
    r0 = 1.6 * rp if r0 is None else r0
    u_max = 2.0 * (r0 - max(rm_, 0.25 * rp)) - 0.15 if u_max is None else u_max
    u = u_grid(u_max, nu, u_eh, n_cluster); v = np.linspace(0.0, v_max, nv)
    dv = v[1] - v[0]
    psi0 = U.h_inv(mu0 / r0 ** (D - 1), 1.0, kind) if kind != "gr" else mu0 / r0 ** (D - 1)
    s0 = -0.5 * np.log(2.0 * (1.0 - r0**2 * psi0))

    def phi_row(vv):
        val = 0.0
        if amp and v1 < vv < v2:
            val += amp * np.sin(np.pi * (vv - v1) / (v2 - v1)) ** 2
        if amp2 and v3 < vv < v4:
            val += amp2 * np.sin(np.pi * (vv - v3) / (v4 - v3)) ** 2
        if tail_amp and vv > v2:
            val += tail_amp * (1.0 - np.exp(-(vv - v2))) * (vv / v2) ** (-tail_p)
        return val

    p_row = np.array([phi_row(x) for x in v])
    pvv = np.gradient(p_row, dv)
    ode = _ode_outgoing(r0, s0, D, kind, v, lambda vv: float(np.interp(vv, v, pvv)))
    s_row = ode[:, 0].copy(); ru_row = ode[:, 1]
    # mu on u = 0 from the first law mu_v = -2 r^{D-2} e^{-2 sigma} r_u phi_v^2 (trapezoid)
    r_row = r0 + 0.5 * v
    dmu = -2.0 * r_row ** (D - 2) * np.exp(-2.0 * s_row) * ru_row * pvv**2
    m_row = mu0 + np.concatenate([[0.0], np.cumsum(0.5 * (dmu[1:] + dmu[:-1]) * dv)])
    su = max(1, (len(u) - 1) // n_out_u); sv = max(1, (nv - 1) // n_out_v)
    out = _evolve_rows_hyb(r_row, s_row, m_row, p_row, r0 - 0.5 * u, np.full(len(u), s0), np.full(len(u), mu0), np.zeros(len(u)), u, dv,
                           float(D), k, n_iter, su, sv)
    return dict(D=D, mu0=mu0, kind=kind, r0=r0, r_plus=rp, r_minus=rm_, u=u[::su], v=v[::sv], r=out[0], s=out[1], mu=out[2], p=out[3],
                psi=out[4], psi_metric=out[5], K=out[6], K_2d=out[7], K_mixed=out[8], nu=len(u), nv=nv)


def run_mu(D=5, mu0=10.0, kind="hayward", r0=None, u_max=None, v_max=60.0, nu=4000, nv=12000, amp=0.0, v1=1.0, v2=5.0,
           amp2=0.0, v3=None, v4=None, n_iter=3, u_eh=None, n_cluster=0, n_out_u=400, n_out_v=600, tail_p=0.0, tail_amp=0.0) -> dict:
    """(r, mu, phi) evolution.  Optional power-law influx tail phi ~ tail_amp (v/v2)^(-tail_p) for v > v2 on the initial outgoing ray
    (a Price tail; with tail_amp = 0 only the backscattered tail of the compact pulse is present)."""
    k = KIND[kind]
    rm_, rp = U.static_horizons(mu0, D, 1.0, kind) if kind != "gr" else (0.0, mu0 ** (1.0 / (D - 3)))
    r0 = 1.6 * rp if r0 is None else r0
    u_max = 2.0 * (r0 - max(rm_, 0.25 * rp)) - 0.15 if u_max is None else u_max
    u = u_grid(u_max, nu, u_eh, n_cluster); v = np.linspace(0.0, v_max, nv)
    dv = v[1] - v[0]

    def phi_row(vv):
        val = 0.0
        if amp and v1 < vv < v2:
            val += amp * np.sin(np.pi * (vv - v1) / (v2 - v1)) ** 2
        if amp2 and v3 < vv < v4:
            val += amp2 * np.sin(np.pi * (vv - v3) / (v4 - v3)) ** 2
        if tail_amp and vv > v2:
            val += tail_amp * (1.0 - np.exp(-(vv - v2))) * (vv / v2) ** (-tail_p)
        return val

    p_row = np.array([phi_row(x) for x in v])
    pv = np.gradient(p_row, dv)
    # mu along u = 0 (r = r0 + v/2): mu_v = 2 r^{D-2} f phi_v^2
    m_row = np.empty(nv); m_row[0] = mu0
    for j in range(nv - 1):
        def g(mu, jj):
            rr = r0 + 0.5 * v[jj]
            ps = float(_hinv(mu / rr ** (D - 1), k))
            return 2.0 * rr ** (D - 2) * (1.0 - rr**2 * ps) * pv[jj] ** 2
        k1 = g(m_row[j], j); k2 = g(m_row[j] + dv * k1, j + 1)
        m_row[j + 1] = m_row[j] + 0.5 * dv * (k1 + k2)
    su = max(1, (len(u) - 1) // n_out_u); sv = max(1, (nv - 1) // n_out_v)
    out = _evolve_rows_mu(r0 + 0.5 * v, m_row, p_row, r0 - 0.5 * u, np.full(len(u), mu0), np.zeros(len(u)), u, dv, float(D), k,
                          n_iter, su, sv)
    return dict(D=D, mu0=mu0, kind=kind, r0=r0, r_plus=rp, r_minus=rm_, u=u[::su], v=v[::sv], r=out[0], mu=out[1], p=out[2],
                psi=out[3], K=out[4], f=out[5], nu=len(u), nv=nv)


def _ode_outgoing(r0, s0, D, kind, v, pv_fun, ru0=-0.5):
    """Along u = 0 with r = r0 + v/2 (r_vv = 0): sigma_v = r phi_v^2 / h'(psi) (vv constraint with r_v = 1/2) and d r_u / dv = r_uv."""
    k = KIND[kind]

    def f(vv, y):
        sig, ru = y
        rr = r0 + 0.5 * vv
        pv = pv_fun(vv)
        F, S, W, psi = _rhs(rr, sig, ru, 0.5, 0.0, pv, float(D), k)
        h0, h1, h2 = _h012(psi, k)
        return np.array([rr * pv**2 / h1, F])

    y = np.array([s0, ru0]); out = np.zeros((len(v), 2)); out[0] = y
    for n in range(len(v) - 1):
        h = v[n + 1] - v[n]; t = v[n]
        k1 = f(t, y); k2 = f(t + h / 2, y + h / 2 * k1); k3 = f(t + h / 2, y + h / 2 * k2); k4 = f(t + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[n + 1] = y
    return out


def u_grid(u_max, nu, u_eh=None, n_cluster=0, gap_min=1e-11, gap_max=0.3):
    """Uniform grid plus geometric clustering on both sides of the event horizon u_eh (outgoing rays that hug the horizon fall from
    r_+ to r_- over a u-width ~ exp(-kappa_+ v); the clustering resolves them up to v ~ ln(1/gap_min)/kappa_+)."""
    base = np.linspace(0.0, u_max, nu)
    if u_eh is None or n_cluster == 0:
        return base
    g = np.geomspace(gap_min, gap_max, n_cluster // 2)
    pts = np.concatenate([base, u_eh - g, u_eh + g, [u_eh]])
    pts = np.unique(np.clip(pts, 0.0, u_max))
    keep = np.concatenate([[True], np.diff(pts) > 0.2 * gap_min])
    return pts[keep]


def run(D=5, mu0=10.0, kind="hayward", r0=None, u_max=None, v_max=60.0, nu=800, nv=2400, amp=0.0, v1=1.0, v2=5.0,
        amp2=0.0, v3=None, v4=None, n_iter=2, u_eh=None, n_cluster=0) -> dict:
    """Scalar pulse(s) on the outgoing ray u = 0 falling into a static model-U (or GR) hole of mass mu0.
    A second ("late infall") pulse with amplitude amp2 on [v3, v4] tests P-R6l-a.  u_eh: event-horizon guess for grid clustering."""
    k = KIND[kind]
    rm, rp = U.static_horizons(mu0, D, 1.0, kind) if kind != "gr" else (0.0, mu0 ** (1.0 / (D - 3)))
    r0 = 1.6 * rp if r0 is None else r0
    u_max = 2.0 * (r0 - max(rm, 0.25 * rp)) - 0.15 if u_max is None else u_max
    u = u_grid(u_max, nu, u_eh, n_cluster); v = np.linspace(0.0, v_max, nv)
    nu = len(u)
    du, dv = np.diff(u), v[1] - v[0]
    psi0 = U.h_inv(mu0 / r0 ** (D - 1), 1.0, kind) if kind != "gr" else mu0 / r0 ** (D - 1)
    f0 = 1.0 - r0**2 * psi0
    s0 = -0.5 * np.log(2.0 * f0)

    def pulse(vv, a, va, vb):
        return a * np.sin(np.pi * (vv - va) / (vb - va)) ** 2 if (a and va < vv < vb) else 0.0

    def dpulse(vv, a, va, vb):
        return a * np.pi / (vb - va) * np.sin(2 * np.pi * (vv - va) / (vb - va)) if (a and va < vv < vb) else 0.0

    pv_fun = lambda vv: dpulse(vv, amp, v1, v2) + (dpulse(vv, amp2, v3, v4) if amp2 else 0.0)
    r = np.zeros((nu, nv)); s = np.zeros((nu, nv)); p = np.zeros((nu, nv))
    r[:, 0] = r0 - 0.5 * u; s[:, 0] = s0
    ode = _ode_outgoing(r0, s0, D, kind, v, pv_fun)
    r[0, :] = r0 + 0.5 * v; s[0, :] = ode[:, 0]
    p[0, :] = [pulse(x, amp, v1, v2) + (pulse(x, amp2, v3, v4) if amp2 else 0.0) for x in v]
    _evolve(r, s, p, u, dv, float(D), k, n_iter)
    return dict(D=D, mu0=mu0, kind=kind, r0=r0, r_plus=rp, r_minus=rm, u=u, v=v, r=r, s=s, p=p, du=du, dv=dv)


def run_fast(D=5, mu0=10.0, kind="hayward", r0=None, u_max=None, v_max=60.0, nu=4000, nv=12000, amp=0.0, v1=1.0, v2=5.0,
             amp2=0.0, v3=None, v4=None, n_iter=2, u_eh=None, n_cluster=0, n_out_u=400, n_out_v=600) -> dict:
    """Same physics as ``run`` with the low-memory row march; returns r, sigma, phi, mu, psi, K on a sub-grid."""
    k = KIND[kind]
    rm, rp = U.static_horizons(mu0, D, 1.0, kind) if kind != "gr" else (0.0, mu0 ** (1.0 / (D - 3)))
    r0 = 1.6 * rp if r0 is None else r0
    u_max = 2.0 * (r0 - max(rm, 0.25 * rp)) - 0.15 if u_max is None else u_max
    u = u_grid(u_max, nu, u_eh, n_cluster); v = np.linspace(0.0, v_max, nv)
    dv = v[1] - v[0]
    psi0 = U.h_inv(mu0 / r0 ** (D - 1), 1.0, kind) if kind != "gr" else mu0 / r0 ** (D - 1)
    s0 = -0.5 * np.log(2.0 * (1.0 - r0**2 * psi0))

    def pulse(vv, a, va, vb):
        return a * np.sin(np.pi * (vv - va) / (vb - va)) ** 2 if (a and va < vv < vb) else 0.0

    def dpulse(vv, a, va, vb):
        return a * np.pi / (vb - va) * np.sin(2 * np.pi * (vv - va) / (vb - va)) if (a and va < vv < vb) else 0.0

    pv_fun = lambda vv: dpulse(vv, amp, v1, v2) + (dpulse(vv, amp2, v3, v4) if amp2 else 0.0)
    ode = _ode_outgoing(r0, s0, D, kind, v, pv_fun)
    p_row = np.array([pulse(x, amp, v1, v2) + (pulse(x, amp2, v3, v4) if amp2 else 0.0) for x in v])
    su = max(1, (len(u) - 1) // n_out_u); sv = max(1, (nv - 1) // n_out_v)
    out = _evolve_rows(r0 + 0.5 * v, ode[:, 0].copy(), p_row, r0 - 0.5 * u, np.full(len(u), s0), np.zeros(len(u)), u, dv, float(D), k,
                       n_iter, su, sv)
    return dict(D=D, mu0=mu0, kind=kind, r0=r0, r_plus=rp, r_minus=rm, u=u[::su], v=v[::sv], r=out[0], s=out[1], p=out[2], mu=out[3],
                psi=out[4], K=out[5], nu=len(u), nv=nv)


def diagnostics(res: dict) -> dict:
    """mu, psi, curvature and constraint residuals on the grid (cell centres)."""
    D, k = res["D"], KIND[res["kind"]]
    r, s, p, dv = res["r"], res["s"], res["p"], res["dv"]
    du = np.diff(res["u"])[:, None]
    rc = 0.25 * (r[1:, 1:] + r[:-1, 1:] + r[1:, :-1] + r[:-1, :-1]); sc = 0.25 * (s[1:, 1:] + s[:-1, 1:] + s[1:, :-1] + s[:-1, :-1])
    ru = 0.5 * (r[1:, 1:] - r[:-1, 1:] + r[1:, :-1] - r[:-1, :-1]) / du
    rv = 0.5 * (r[1:, 1:] - r[1:, :-1] + r[:-1, 1:] - r[:-1, :-1]) / dv
    pu = 0.5 * (p[1:, 1:] - p[:-1, 1:] + p[1:, :-1] - p[:-1, :-1]) / du
    pv = 0.5 * (p[1:, 1:] - p[1:, :-1] + p[:-1, 1:] - p[:-1, :-1]) / dv
    E = np.exp(2 * sc)
    psi = (1 + 2 * ru * rv / E) / rc**2
    x = np.clip(1 - psi, 1e-14, None)
    if k == 0:
        h0, h1, h2 = psi, np.ones_like(psi), np.zeros_like(psi)
    elif k == 1:
        h0, h1, h2 = psi / x, 1 / x**2, 2 / x**3
    else:
        h0, h1, h2 = -np.log(x), 1 / x, 1 / x**2
    mu = rc ** (D - 1) * h0
    F = 0.5 * rc * E * (2 * psi - (D - 1) * h0 / h1)
    S = (-(D - 1) ** 2 * h0**2 * h2 * rc**2 * E**2 + 2 * h1**3 * (2 * ru * rv + E) * E
         + h1**2 * rc**2 * E * ((D - 1) * (D - 4) * h0 * E - 2 * (D - 2) * pu * pv) + 4 * h2 * rc**2 * pu**2 * pv**2) / (E * 4 * h1**3 * rc**2)
    Auu = -rc * pu**2 / h1; Avv = -rc * pv**2 / h1
    nnr2 = 2 * np.exp(-4 * sc) * (Auu * Avv + F**2)
    R2 = 4 * np.exp(-2 * sc) * S
    K = R2**2 + 4 * (D - 2) * nnr2 / rc**2 + 2 * (D - 2) * (D - 3) * psi**2
    # vv-constraint residual from second differences along v (interior points)
    res_vv = np.full_like(rc, np.nan)
    rvv = (r[1:-1, 2:] - 2 * r[1:-1, 1:-1] + r[1:-1, :-2]) / dv**2
    sv = (s[1:-1, 2:] - s[1:-1, :-2]) / (2 * dv); rv_c = (r[1:-1, 2:] - r[1:-1, :-2]) / (2 * dv); pv_c = (p[1:-1, 2:] - p[1:-1, :-2]) / (2 * dv)
    ru_c = (r[2:, 1:-1] - r[:-2, 1:-1]) / (res["u"][2:, None] - res["u"][:-2, None])
    psi_c = (1 + 2 * ru_c * rv_c / np.exp(2 * s[1:-1, 1:-1])) / r[1:-1, 1:-1] ** 2
    xc = np.clip(1 - psi_c, 1e-14, None)
    h1c = np.ones_like(psi_c) if k == 0 else (1 / xc**2 if k == 1 else 1 / xc)
    cvv = h1c * (rvv - 2 * sv * rv_c) + r[1:-1, 1:-1] * pv_c**2
    return dict(mu=mu, psi=psi, K=K, rc=rc, ru=ru, rv=rv, theta_out=rv, theta_in=ru, constraint_vv=cvv,
                trapped=(ru < 0) & (rv < 0))


def static_check(D=5, mu0=10.0, kind="hayward", nu=400, nv=800, v_max=20.0, n_cluster=300) -> dict:
    """No scalar field: the evolution must keep mu = mu0 everywhere (region I and the trapped region II)."""
    rm, rp = U.static_horizons(mu0, D, 1.0, kind) if kind != "gr" else (0.0, mu0 ** (1.0 / (D - 3)))
    res = run(D, mu0, kind, nu=nu, nv=nv, v_max=v_max, amp=0.0, u_eh=2.0 * 0.6 * rp, n_cluster=n_cluster)
    d = diagnostics(res)
    m = d["mu"][np.isfinite(d["mu"])]
    return dict(D=D, kind=kind, nu=nu, nv=nv, max_rel_mass_error=float(np.max(np.abs(m / mu0 - 1))),
                trapped_fraction=float(np.mean(d["trapped"][np.isfinite(d["mu"])])))


def convergence(D=5, mu0=10.0, kind="hayward", amp=0.12, base=(200, 400), v_max=16.0) -> dict:
    """Three resolutions; the observed order from mu at a fixed physical point inside the hole."""
    vals = []
    for f in (1, 2, 4):
        res = run(D, mu0, kind, nu=base[0] * f, nv=base[1] * f, v_max=v_max, amp=amp)
        d = diagnostics(res)
        iu = int(0.85 * (base[0] * f - 1)); iv = int(0.9 * (base[1] * f - 1))
        vals.append(float(d["mu"][iu, iv]))
    order = float(np.log2(abs(vals[0] - vals[1]) / max(abs(vals[1] - vals[2]), 1e-300)))
    return dict(mu_values=vals, observed_order=order)


def derive_sigma_uv():
    """Symbolic derivation (sympy) of the sigma_uv equation from preservation of the uu constraint; returns the GR reduction."""
    import sympy as sp
    D = sp.Symbol("D")
    r, ru, rv, ruu, rvv, s, su, sv, p, pu, pv, puu = sp.symbols("r r_u r_v r_uu r_vv s s_u s_v p p_u p_v p_uu")
    Suv, RUUV, suu = sp.symbols("S_uv r_uuv s_uu")
    h0, h1, h2 = sp.symbols("h0 h1 h2")
    E = sp.exp(2 * s)
    psi = (1 + 2 * ru * rv / E) / r**2
    F = (r * E / 2) * (2 * psi - (D - 1) * h0 / h1)
    W = -(D - 2) / (2 * r) * (ru * pv + rv * pu)

    def Du(e, use_h=True):
        d = (sp.diff(e, r) * ru + sp.diff(e, ru) * ruu + sp.diff(e, rv) * F + sp.diff(e, s) * su + sp.diff(e, su) * suu
             + sp.diff(e, sv) * Suv + sp.diff(e, p) * pu + sp.diff(e, pu) * puu + sp.diff(e, pv) * W)
        if use_h:
            dpsi = Du(psi, False)
            d += sp.diff(e, h0) * h1 * dpsi + sp.diff(e, h1) * h2 * dpsi
        return d

    def Dv(e, use_h=True):
        d = (sp.diff(e, r) * rv + sp.diff(e, ru) * F + sp.diff(e, rv) * rvv + sp.diff(e, s) * sv + sp.diff(e, su) * Suv
             + sp.diff(e, p) * pv + sp.diff(e, pu) * W + sp.diff(e, ruu) * RUUV)
        if use_h:
            dpsi = Dv(psi, False)
            d += sp.diff(e, h0) * h1 * dpsi + sp.diff(e, h1) * h2 * dpsi
        return d

    Cu = h1 * (ruu - 2 * su * ru) + r * pu**2
    expr = sp.expand(Dv(Cu).subs(RUUV, Du(F)).subs(ruu, sp.solve(Cu, ruu)[0]))
    S = sp.solve(expr, Suv)[0]
    S = S.subs(rvv, 2 * sv * rv - r * pv**2 / h1)            # vv constraint
    S = sp.simplify(S)
    S_gr = sp.factor(sp.simplify(S.subs({h0: psi, h1: 1, h2: 0})))
    target = (D - 2) / (4 * r**2) * ((D - 3) * (E + 2 * ru * rv) - 2 * r**2 * pu * pv)
    return dict(sigma_uv=S, gr=S_gr, gr_matches=bool(sp.simplify(S_gr - target) == 0))
