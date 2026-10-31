"""Round 11 (G11-02): is there a stable inner horizon?  An inner-extremal model-U tower and its double-null evolution.

Construction (D = 5, x = r^2, alpha = 1).  For one tuned mass the static metric function is chosen as
    f(x) = (x - x_-)^3 (x - x_+) / Q(x),   Q(x) = x^4 + p3 x^3 + p2 x^2 + p1 x + q,
with q = x_-^3 x_+ (so f(0) = 1), p1 fixed by f'(0) = -1 (de Sitter core f ~ 1 - r^2/alpha), mu = 3 x_- + x_+ + p3 (GR fall-off
f ~ 1 - mu/r^2), and (p2, p3) chosen so that Q > 0 and psi(x) = (1 - f)/x decreases monotonically from 1.  Then h(psi(x)) = mu/x^2 defines
the tower function h on (0, 1) (any such h corresponds to a quasi-topological tower; Bueno, Cano & Hennigar), with a triply degenerate
inner horizon (kappa_- = 0, f'' = 0) and a non-degenerate outer horizon - possible only for the tuned mass, as Di Filippo, Kolar &
Kubiznak (PRD 111, L041505 (2025)) proved in general.

The double-null kernel is the hybrid (r, sigma, mu, phi) solver of ``doublenull10`` with h, h', h'' and h^{-1} read from the table.
"""

from __future__ import annotations

import numpy as np
from numba import njit

X_MINUS, X_PLUS, P2, P3 = 1.5, 9.0, 4.0, 2.5


def ie_metric(xm=X_MINUS, xp=X_PLUS, p2=P2, p3=P3):
    N = np.poly1d([1, -xm]) ** 3 * np.poly1d([1, -xp])
    q = xm**3 * xp
    p1 = N.deriv()(0.0) + q
    Q = np.poly1d([1, p3, p2, p1, q])
    return N, Q, 3 * xm + xp + p3                       # GR fall-off: 1 - f -> (3 x_- + x_+ + p3)/x


def ie_table(n=400000, x_max=1e5, x_min=1e-3):
    """Tower table: psi ascending, h(psi), h'(psi), h''(psi) for the inner-extremal metric (mass mu = 3 x_- + x_+)."""
    N, Q, mu = ie_metric()
    x = np.geomspace(x_min, x_max, n)
    f = N(x) / Q(x)
    psi = (1 - f) / x
    h = mu / x**2
    dpsi = np.gradient(psi, x); dh = np.gradient(h, x)
    h1 = dh / dpsi
    h2 = np.gradient(h1, x) / dpsi
    o = np.argsort(psi)
    return dict(psi=psi[o], h0=h[o], h1=h1[o], h2=h2[o], mu=mu, x=x[o], f=f[o])


def ie_checks():
    N, Q, mu = ie_metric()
    roots = sorted(float(r.real) for r in N.r)
    x = np.linspace(1e-6, 200, 200000)
    psi = (1 - N(x) / Q(x)) / x
    dfdx = (N.deriv() * Q - N * Q.deriv())
    k_minus = float(dfdx(X_MINUS) / Q(X_MINUS) ** 2)
    return dict(mu=mu, roots_x=roots, r_minus=float(np.sqrt(X_MINUS)), r_plus=float(np.sqrt(X_PLUS)), psi_monotone=bool(np.all(np.diff(psi) < 0)),
                psi_max=float(psi.max()), Q_positive=bool(np.all(Q(x) > 0)), df_dx_at_inner=k_minus)


@njit(cache=True, error_model="numpy")
def _tab(psi, tp, t0, t1, t2):
    return np.interp(psi, tp, t0), np.interp(psi, tp, t1), np.interp(psi, tp, t2)


@njit(cache=True, error_model="numpy")
def _rhs_tab(r, s, mu, ru, rv, pu, pv, D, tp, t0, t1, t2):
    E = np.exp(2.0 * s)
    y = mu / r ** (D - 1.0)
    psi = np.interp(y, t0, tp)                       # h is increasing: invert by interpolation
    h0, h1, h2 = _tab(psi, tp, t0, t1, t2)
    F = 0.5 * r * E * (2.0 * psi - (D - 1.0) * h0 / h1)
    W = -(D - 2.0) / (2.0 * r) * (ru * pv + rv * pu)
    S = (-(D - 1.0) ** 2 * h0**2 * h2 * r**2 * E**2 + 2.0 * h1**3 * r**2 * psi * E**2
         + h1**2 * r**2 * E * ((D - 1.0) * (D - 4.0) * h0 * E - 2.0 * (D - 2.0) * pu * pv) + 4.0 * h2 * r**2 * pu**2 * pv**2) / (E * 4.0 * h1**3 * r**2)
    rD2 = r ** (D - 2.0)
    return F, S, W, -2.0 * rD2 * rv * pu**2 / E, -2.0 * rD2 * ru * pv**2 / E, psi, h1


@njit(cache=True, error_model="numpy")
def _evolve_tab(r_row0, s_row0, m_row0, p_row0, r_col0, s_col0, m_col0, p_col0, uu, dv, D, tp, t0, t1, t2, n_iter, su, sv):
    nu = len(uu); nv = len(r_row0)
    nsu = (nu - 1) // su + 1; nsv = (nv - 1) // sv + 1
    out = np.full((7, nsu, nsv), np.nan)
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
                F, S, W, mu_u, mu_v, psi, h1 = _rhs_tab(rm, sm, mm, ru, rv, pu, pv, D, tp, t0, t1, t2)
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
                    F, S, W, mu_u, mu_v, psi, h1 = _rhs_tab(rm, sm, mm, ru, rv, pu, pv, D, tp, t0, t1, t2)
                    E = np.exp(2.0 * sm)
                    Auu = -rm * pu**2 / h1; Avv = -rm * pv**2 / h1
                    K = (4.0 * S / E) ** 2 + 4.0 * (D - 2.0) * 2.0 * (Auu * Avv + F**2) / E**2 / rm**2 + 2.0 * (D - 2.0) * (D - 3.0) * psi**2
                    out[4, ii, jj] = psi; out[5, ii, jj] = K; out[6, ii, jj] = h1
        for j in range(nv):
            rp_[j] = rn[j]; sp_[j] = sn[j]; mp_[j] = mn[j]; pp_[j] = pn[j]
    return out


def run_ie(amp=0.02, v_max=24.0, nu=2000, nv=8000, r0=None, u_max=None, u_eh=None, n_cluster=800, n_out_u=240, n_out_v=400, v1=1.0, v2=5.0):
    """Scalar pulse into the inner-extremal hole (D = 5).  Same gauge and initial data as ``doublenull10.run_hyb``."""
    from .doublenull10 import u_grid
    T = ie_table()
    D = 5.0
    tp, t0, t1, t2 = T["psi"], T["h0"], T["h1"], T["h2"]
    mu0 = T["mu"]
    rp, rm = np.sqrt(X_PLUS), np.sqrt(X_MINUS)
    r0 = 1.6 * rp if r0 is None else r0
    u_max = 2.0 * (r0 - rm - 0.08) if u_max is None else u_max
    u = u_grid(u_max, nu, u_eh, n_cluster if u_eh is not None else 0); v = np.linspace(0.0, v_max, nv)
    dv = v[1] - v[0]
    psi0 = float(np.interp(mu0 / r0 ** 4, t0, tp))
    s0 = -0.5 * np.log(2.0 * (1.0 - r0**2 * psi0))
    p_row = np.array([amp * np.sin(np.pi * (x - v1) / (v2 - v1)) ** 2 if v1 < x < v2 else 0.0 for x in v])
    pv = np.gradient(p_row, dv)
    # u = 0 ray (r = r0 + v/2): sigma_v = r phi_v^2 / h'(psi), d r_u / dv = r_uv; mu from the first law (RK2)
    s_row = np.empty(nv); ru_row = np.empty(nv); m_row = np.empty(nv)
    s_row[0], ru_row[0], m_row[0] = s0, -0.5, mu0
    for j in range(nv - 1):
        def g(sv, ruv, mv, jj):
            rr = r0 + 0.5 * v[jj]
            F, S, W, mu_u, mu_v, psi, h1 = _rhs_tab(rr, sv, mv, ruv, 0.5, 0.0, pv[jj], D, tp, t0, t1, t2)
            return np.array([rr * pv[jj] ** 2 / h1, F, mu_v])
        y = np.array([s_row[j], ru_row[j], m_row[j]])
        k1 = g(*y, j); k2 = g(*(y + dv * k1), j + 1)
        y = y + 0.5 * dv * (k1 + k2)
        s_row[j + 1], ru_row[j + 1], m_row[j + 1] = y
    su = max(1, (len(u) - 1) // n_out_u); sv = max(1, (nv - 1) // n_out_v)
    out = _evolve_tab(r0 + 0.5 * v, s_row, m_row, p_row, r0 - 0.5 * u, np.full(len(u), s0), np.full(len(u), mu0), np.zeros(len(u)), u, dv, D,
                      tp, t0, t1, t2, 3, su, sv)
    return dict(D=5, mu0=mu0, kind="inner-extremal", r0=r0, r_plus=rp, r_minus=rm, u=u[::su], v=v[::sv], r=out[0], s=out[1], mu=out[2],
                p=out[3], psi=out[4], K=out[5], h1=out[6])
