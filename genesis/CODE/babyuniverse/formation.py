"""Round 6f (items 1 and 4): one continuous formation history in spherical symmetry.

Matter, all evolved together with gravity (metric -e^{2s} du dv + r^2 dOmega^2):
* phi: the vacuum field, tilted quartic U(phi) (false vacuum -v, true +v), optionally with the condensate tilt (rho_f -> rho_f + lam_f T_c);
* a thermal plasma: a perfect fluid (temperature T, rapidity eta_f) whose free energy depends on phi,
      F(T, phi) = -a_L T^4 - a_H T^4 exp(-k g^2 (phi + v)^2 / T^2),
  light species (a_L, massless in both vacua) and heavy species (a_H, massless at the false vacuum, Boltzmann-like suppressed in the true
  vacuum).  p = -F, rho = F - T F_T.  At small mass it reduces to the high-T form (c_T/2) T^2 (phi + v)^2 with c_T = 2 a_H k g^2; in the
  true vacuum the heavy species decouple (bag), and energy leaves through the light species;
* friction between plasma and field (Ignatius-Kajantie-Kurki-Suonio-Laine): box(phi) = U' + F_phi + Gamma u_f . d phi, with the
  matching source in the fluid so that the total stress tensor is conserved;
* the ghost condensate (IR fluid of ``condensate``: rapidity eta_c, charge D, clock T_c; tilt source in the false vacuum);
* an ingoing classical null shell N_v(v) (ordinary radiation; enters only the constraints).

Fluid in double-null form (A = r^2 T_uv, X = r^2 e^{2s} (rho + p)/4, B_u = X e^{2 eta}, B_v = X e^{-2 eta}):
    A_u + (B_u)_v = 2 s_u A - r r_u p e^{2s} - (e^{2s} r^2/2) Q_u,     A_v + (B_v)_u = 2 s_v A - r r_v p e^{2s} - (e^{2s} r^2/2) Q_v,
    Q_b = -(F_phi + Gamma u^a d_a phi) d_b phi,   r_uv += 4 pi G A / r,   s_uv += -4 pi G A / r^2 - 2 pi G e^{2s} p.
Box scheme for both conservation laws (2x2 Newton for ln T and eta_f at the new point); the scalar potential force (U + F) is implicit.
"""

from __future__ import annotations

import math

import numpy as np

from .throat_evolution import njit
from .condensate import _chi, _Uc, _dUc, _d2Uc


@njit(cache=True)
def _eos(T, ph, aL, aH, kg2, vv):
    """(rho, p, F_phi, F_phiphi) of the plasma at temperature T and field phi."""
    m2 = (ph + vv) * (ph + vv)
    x = kg2 * m2 / (T * T)
    E = math.exp(-x)
    T4 = T * T * T * T
    p = aL * T4 + aH * T4 * E
    rho = 3.0 * aL * T4 + aH * E * (3.0 * T4 + 2.0 * kg2 * m2 * T * T)
    Fp = 2.0 * aH * kg2 * T * T * E * (ph + vv)
    Fpp = 2.0 * aH * kg2 * T * T * E * (1.0 - 2.0 * x)
    return rho, p, Fp, Fpp


@njit(cache=True)
def _fluid_vars(r, s, T, eta, ph, aL, aH, kg2, vv):
    rho, p, Fp, Fpp = _eos(T, ph, aL, aH, kg2, vv)
    e2s = math.exp(2.0 * s)
    A = r * r * e2s * (rho - p) / 4.0
    X = r * r * e2s * (rho + p) / 4.0
    return A, X * math.exp(2.0 * eta), X * math.exp(-2.0 * eta), p, Fp


@njit(cache=True)
def _fluid_solve(rN, sN, phN, AE, AW, AS, BuE, BuW, BuS, BvE, BvW, BvS, R1, R2, lT0, et0, aL, aH, kg2, vv, T_floor, th):
    """Newton for (ln T, eta_f) at N in the Preissmann-weighted box (th = 1/2 is the plain box; th > 1/2 damps its checkerboard mode):
       (A_N + A_W - A_E - A_S) + 2 th (Bu_N - Bu_W) + 2 (1 - th)(Bu_E - Bu_S) = R1,
       2 th (A_N - A_W) + 2 (1 - th)(A_E - A_S) + (Bv_N + Bv_W - Bv_E - Bv_S) = R2."""
    t1 = R1 - (AW - AE - AS) + 2.0 * th * BuW - 2.0 * (1.0 - th) * (BuE - BuS)
    t2 = R2 + 2.0 * th * AW - 2.0 * (1.0 - th) * (AE - AS) - (BvW - BvE - BvS)
    lT, et = lT0, et0
    lfl = math.log(T_floor)
    for it in range(30):
        A, Bu, Bv, p, Fp = _fluid_vars(rN, sN, math.exp(lT), et, phN, aL, aH, kg2, vv)
        f1 = A + 2.0 * th * Bu - t1
        f2 = 2.0 * th * A + Bv - t2
        d = 1e-7
        A1, Bu1, Bv1, p1, F1 = _fluid_vars(rN, sN, math.exp(lT + d), et, phN, aL, aH, kg2, vv)
        A2, Bu2, Bv2, p2, F2 = _fluid_vars(rN, sN, math.exp(lT), et + d, phN, aL, aH, kg2, vv)
        j11 = (A1 - A + 2.0 * th * (Bu1 - Bu)) / d
        j12 = (A2 - A + 2.0 * th * (Bu2 - Bu)) / d
        j21 = (2.0 * th * (A1 - A) + Bv1 - Bv) / d
        j22 = (2.0 * th * (A2 - A) + Bv2 - Bv) / d
        det = j11 * j22 - j12 * j21
        if det == 0.0 or not math.isfinite(det):
            break
        dl = (f1 * j22 - f2 * j12) / det
        de = (j11 * f2 - j21 * f1) / det
        if dl > 1.0:
            dl = 1.0
        if dl < -1.0:
            dl = -1.0
        if de > 1.0:
            de = 1.0
        if de < -1.0:
            de = -1.0
        lT -= dl
        et -= de
        if lT < lfl:
            lT = lfl
        if abs(dl) + abs(de) < 1e-12:
            break
    if lT <= lfl + 1e-12:              # vacuum (floor): no fluid velocity information; keep the rapidity bounded
        et = et0
    if et > 12.0:
        et = 12.0
    if et < -12.0:
        et = -12.0
    return lT, et


@njit(cache=True)
def _diamond(X, h, G, lam, vv, rho_f, lam_f, clip, aL, aH, kg2, Gam, T_floor, th, n_iter):
    """X[q, c], q = r, s, phi, eta_c, D, T_c, lnT, eta_f; c = E, W, S.  Returns the eight values at N."""
    NQ = 8
    out = np.empty(NQ)
    for q in range(NQ):
        out[q] = X[q, 0] + X[q, 1] - X[q, 2]
    c = np.empty(NQ)
    for it in range(n_iter):
        for q in range(NQ):
            c[q] = 0.25 * (out[q] + X[q, 0] + X[q, 1] + X[q, 2])
        ru = (out[0] + X[0, 1] - X[0, 0] - X[0, 2]) / (2.0 * h)
        rv = (out[0] + X[0, 0] - X[0, 1] - X[0, 2]) / (2.0 * h)
        su = (out[1] + X[1, 1] - X[1, 0] - X[1, 2]) / (2.0 * h)
        sv = (out[1] + X[1, 0] - X[1, 1] - X[1, 2]) / (2.0 * h)
        pu = (out[2] + X[2, 1] - X[2, 0] - X[2, 2]) / (2.0 * h)
        pv = (out[2] + X[2, 0] - X[2, 1] - X[2, 2]) / (2.0 * h)
        rc, sc, phc = c[0], c[1], c[2]
        a2 = math.exp(2.0 * sc)
        rho_eff = rho_f + lam_f * c[5]
        # plasma at the centre
        Tc_ = math.exp(c[6])
        Af, Buf, Bvf, pf, Fph = _fluid_vars(rc, sc, Tc_, c[7], phc, aL, aH, kg2, vv)
        uu = math.exp(-sc - c[7])
        uv = math.exp(-sc + c[7])
        fric = Gam * (uu * pu + uv * pv)
        # gravity
        Lam = 8.0 * math.pi * G * _Uc(phc, lam, vv, rho_eff, clip)
        T2uv = 0.25 * c[4] * math.exp(c[1] + c[3])
        As = ru * rv / (rc * rc) + a2 / (4.0 * rc * rc) - 4.0 * math.pi * G * pu * pv - G * T2uv / (rc * rc) \
            - 4.0 * math.pi * G * Af / (rc * rc) - 2.0 * math.pi * G * a2 * pf
        Fr = -ru * rv / rc + a2 * (Lam * rc * rc - 1.0) / (4.0 * rc) + G * T2uv / rc + 4.0 * math.pi * G * Af / rc
        out[0] = X[0, 0] + X[0, 1] - X[0, 2] + h * h * Fr
        out[1] = X[1, 0] + X[1, 1] - X[1, 2] + h * h * As
        # scalar: transport and friction explicit, potential + thermal force implicit
        Fp0 = -(ru * pv + rv * pu) / rc - 0.25 * a2 * fric
        base = X[2, 0] + X[2, 1] - X[2, 2] + h * h * Fp0
        S3 = X[2, 0] + X[2, 1] + X[2, 2]
        pN = out[2]
        for nk in range(6):
            pcn = 0.25 * (pN + S3)
            rh, pp, Fpn, Fppn = _eos(Tc_, pcn, aL, aH, kg2, vv)
            g = pN - base + 0.25 * h * h * a2 * (_dUc(pcn, lam, vv, rho_eff, clip) + Fpn)
            gp = 1.0 + h * h * a2 * (_d2Uc(pcn, lam, vv, rho_eff, clip) + Fppn) / 16.0
            if gp < 0.25:
                gp = 0.25
            pN = pN - g / gp
        out[2] = pN
        # condensate (geodesic flow, charge with tilt source, clock)
        lmb = math.exp(2.0 * c[3])
        sN = out[1]
        eE, eW, eS = X[3, 0], X[3, 1], X[3, 2]
        sE, sW, sS = X[1, 0], X[1, 1], X[1, 2]
        out[3] = -((lmb - 1.0) * sN + (eW - sW) - (eE - sE) - (eS - sS) + lmb * ((eE + sE) - (eW + sW) - (eS + sS))) / (1.0 + lmb)
        DE, DW, DS = X[4, 0], X[4, 1], X[4, 2]
        src = -4.0 * math.pi * a2 * rc * rc * lam_f * _chi(phc, vv, clip)
        out[4] = (DE - DW + DS - DE * math.exp(2.0 * eE) + DW * math.exp(2.0 * eW) + DS * math.exp(2.0 * eS) + 2.0 * h * src) \
            / (1.0 + math.exp(2.0 * out[3]))
        TE_, TW, TS = X[5, 0], X[5, 1], X[5, 2]
        out[5] = (2.0 * h * math.exp(c[1] + c[3]) - TW + TE_ + TS - lmb * (TE_ - TW - TS)) / (1.0 + lmb)
        # plasma at N (box scheme, Newton)
        Qu = -(Fph + fric) * pu
        Qv = -(Fph + fric) * pv
        S1 = 2.0 * su * Af - rc * ru * pf * a2 - 0.5 * a2 * rc * rc * Qu
        S2 = 2.0 * sv * Af - rc * rv * pf * a2 - 0.5 * a2 * rc * rc * Qv
        AE, BuE, BvE, p_, F_ = _fluid_vars(X[0, 0], X[1, 0], math.exp(X[6, 0]), X[7, 0], X[2, 0], aL, aH, kg2, vv)
        AW, BuW, BvW, p_, F_ = _fluid_vars(X[0, 1], X[1, 1], math.exp(X[6, 1]), X[7, 1], X[2, 1], aL, aH, kg2, vv)
        AS, BuS, BvS, p_, F_ = _fluid_vars(X[0, 2], X[1, 2], math.exp(X[6, 2]), X[7, 2], X[2, 2], aL, aH, kg2, vv)
        lT, et = _fluid_solve(out[0], out[1], out[2], AE, AW, AS, BuE, BuW, BuS, BvE, BvW, BvS, 2.0 * h * S1, 2.0 * h * S2,
                              out[6], out[7], aL, aH, kg2, vv, T_floor, th)
        out[6] = lT
        out[7] = et
    return out


@njit(cache=True)
def _evolve(A, B, h, G, lam, vv, rho_f, lam_f, clip, aL, aH, kg2, Gam, T_floor, th, n_iter, r_exc, r_inf, phi_max):
    """Corner evolution: A[q, j] on u = u0 (outgoing), B[q, i] on v = v0 (ingoing)."""
    NQ = A.shape[0]
    N = A.shape[1] - 1
    M = B.shape[1] - 1
    Z = np.full((NQ, M + 1, N + 1), np.nan)
    end = np.zeros(M + 1, dtype=np.int64)
    Z[:, 0, :] = A
    X = np.empty((NQ, 3))
    for i in range(1, M + 1):
        Z[:, i, 0] = B[:, i]
        for j in range(1, N + 1):
            if not (math.isfinite(Z[0, i - 1, j]) and math.isfinite(Z[0, i - 1, j - 1]) and math.isfinite(Z[0, i, j - 1])):
                break
            for q in range(NQ):
                X[q, 0] = Z[q, i - 1, j]
                X[q, 1] = Z[q, i, j - 1]
                X[q, 2] = Z[q, i - 1, j - 1]
            o = _diamond(X, h, G, lam, vv, rho_f, lam_f, clip, aL, aH, kg2, Gam, T_floor, th, n_iter)
            ok = True
            for q in range(NQ):
                if not math.isfinite(o[q]):
                    ok = False
            if not ok or abs(o[1]) > 500.0 or abs(o[2]) > phi_max or abs(o[3]) > 30.0:
                end[i] = 3
                break
            if o[0] < r_exc:
                end[i] = 1
                break
            if o[0] > r_inf:
                end[i] = 2
                break
            for q in range(NQ):
                Z[q, i, j] = o[q]
    return Z, end


class Plasma:
    """Thermal sector, parameterised by g_star (total), f_H (heavy fraction) and x_c = k g^2 (2v)^2 / T_c^2 (how strongly the heavy
    species are suppressed in the true vacuum at T_c).  T_c is the exact degeneracy temperature of this free energy:
    a_H T_c^4 (1 - e^{-x_c}) = rho_f.  Gamma in units of the false-vacuum mass; vacuum floor T_floor."""

    def __init__(self, rho_f, vphi, lam, g_star=10.0, f_H=0.95, x_c=5.0, Gamma_over_m=0.3, T_floor_over_Tc=0.02):
        a = np.pi**2 / 30 * g_star
        self.aL, self.aH = a * (1 - f_H), a * f_H
        self.Tc = (rho_f / (self.aH * (1 - np.exp(-x_c)))) ** 0.25
        self.kg2 = x_c * self.Tc**2 / (4 * vphi**2)
        m_f = np.sqrt(max(2 * lam * vphi**2 - 1.5 * rho_f / vphi**2, 1e-30))
        self.Gamma = Gamma_over_m * m_f
        self.T_floor = T_floor_over_Tc * self.Tc
        self.v, self.rho_f, self.lam, self.x_c, self.f_H = vphi, rho_f, lam, x_c, f_H

    def eos(self, T, ph):
        return _eos(float(T), float(ph), self.aL, self.aH, self.kg2, self.v)

    def args(self):
        return self.aL, self.aH, self.kg2, self.Gamma, self.T_floor

    def spinodal(self):
        """Lowest T/T_c at which the true-vacuum basin no longer exists (the field must leave +v)."""
        xs = np.linspace(-0.74, 1.3, 6000) * self.v
        dU = self.lam * xs * (xs**2 - self.v**2) - self.rho_f * (0.75 - 0.75 * (xs / self.v) ** 2) / self.v
        for r in np.linspace(1.0, 4.0, 1201):
            T = r * self.Tc
            f = dU + 2 * self.aH * self.kg2 * T * T * np.exp(-self.kg2 * (xs + self.v) ** 2 / T**2) * (xs + self.v)
            if not np.any((f[:-1] < 0) & (f[1:] >= 0)):
                return float(r)
        return None


def evolve_corner(A, B, h, pot, plasma, lam_f=0.0, clip=1, G=1.0, n_iter=3, r_exc=1e-3, r_inf=1e4, phi_max=None, theta=0.6):
    aL, aH, kg2, Gam, Tfl = plasma.args()
    Z, end = _evolve(np.ascontiguousarray(A, float), np.ascontiguousarray(B, float), h, G, pot.lam, pot.v, pot.rho_f, lam_f, int(clip),
                     aL, aH, kg2, Gam, Tfl, theta, n_iter, r_exc, r_inf, phi_max if phi_max is not None else 10 * pot.v)
    return dict(r=Z[0], s=Z[1], phi=Z[2], eta=Z[3], D=Z[4], T=Z[5], lnTf=Z[6], etaf=Z[7], end=end, h=h, G=G, Q=0.0, c=0.0)


def corner_data(r0, ru0, rv0, h, Nv, Nu, plasma, phi_v=None, phi_u=None, T_v=None, T_u=None, n_v=None, n_u=None, extra_Nv=None, G=1.0):
    """Constraint-satisfying corner data (gauge s = 0 on both rays), fluid and condensate at rest in the corner frame (eta = 0):
    4 pi r^2 T_vv = pi r^2 (rho + p) on u = 0 and likewise T_uu on v = 0.  Returns (A, B, info) for ``evolve_corner``."""
    from scipy.integrate import cumulative_trapezoid
    from . import throat_evolution as TE
    zero = lambda x: np.zeros_like(np.asarray(x, float))
    phi_v, phi_u = phi_v or zero, phi_u or zero
    T_v = T_v or (lambda x: plasma.T_floor + zero(x))
    T_u = T_u or (lambda x: plasma.T_floor + zero(x))
    n_v, n_u = n_v or zero, n_u or zero
    v = np.arange(Nv + 1) * h
    u = np.arange(Nu + 1) * h
    eos = np.vectorize(lambda T, ph: sum(plasma.eos(T, ph)[:2]))     # rho + p
    wv = eos(T_v(v), phi_v(v))
    wu = eos(T_u(u), phi_u(u))
    ra, rb = r0 + rv0 * v, r0 + ru0 * u
    for _ in range(3):
        Tv = np.pi * ra**2 * (wv + n_v(v))
        Tu = np.pi * rb**2 * (wu + n_u(u))
        Nvf = (lambda x, Tv=Tv: np.interp(x, v, Tv)) if extra_Nv is None else (lambda x, Tv=Tv: np.interp(x, v, Tv) + extra_Nv(x))
        a, b, info = TE.characteristic_data(0.0, 0.0, r0, ru0, rv0, h, Nv, Nu, phi_v=phi_v, phi_u=phi_u, N_u=lambda x, Tu=Tu: np.interp(x, u, Tu),
                                            N_v=Nvf, G=G)
        ra, rb = a[0], b[0]
    # final data from the profiles used in the last constraint integration (exactly consistent up to the r-iteration)
    Da = 4 * np.pi * ra**2 * n_v(v)
    Db = 4 * np.pi * rb**2 * n_u(u)
    A = np.vstack([a[0], a[1], a[2], 0 * v, Da, cumulative_trapezoid(0.5 + 0 * v, v, initial=0.0), np.log(T_v(v)), 0 * v])
    B = np.vstack([b[0], b[1], b[2], 0 * u, Db, cumulative_trapezoid(0.5 + 0 * u, u, initial=0.0), np.log(T_u(u)), 0 * u])
    return A, B, info


def constraint_residual(res, plasma, extra_Nv=None, r_lo=0.0, r_hi=np.inf, s_max=np.inf):
    """Median relative and max residuals of the two null constraints with every source (scalar, plasma, condensate, null shell)."""
    R, S, P, h, G = res["r"], res["s"], res["phi"], res["h"], res["G"]
    T, ef, ec, D = np.exp(res["lnTf"]), res["etaf"], res["eta"], res["D"]
    w = np.vectorize(lambda t, ph: sum(plasma.eos(t, ph)[:2]) if np.isfinite(t) and np.isfinite(ph) else np.nan)(T, P)
    base_f = np.pi * R**2 * w * np.exp(2 * S)
    T2uu = base_f * np.exp(2 * ef) + 0.25 * D * np.exp(S) * np.exp(3 * ec)
    T2vv = base_f * np.exp(-2 * ef) + 0.25 * D * np.exp(S) * np.exp(-ec)
    v = np.arange(R.shape[1]) * h
    if extra_Nv is not None:
        T2vv = T2vv + extra_Nv(v)[None, :]
    out = {}
    for ax, key, T2 in ((1, "vv", T2vv), (0, "uu", T2uu)):
        sl = lambda a, o: np.take(a, np.arange(1 + o, a.shape[ax] - 1 + o), axis=ax)
        R0 = sl(R, 0)
        rxx = (sl(R, 1) - 2 * R0 + sl(R, -1)) / h**2
        rx = (sl(R, 1) - sl(R, -1)) / (2 * h)
        sx = (sl(S, 1) - sl(S, -1)) / (2 * h)
        px = (sl(P, 1) - sl(P, -1)) / (2 * h)
        t = sl(T2, 0)
        E = rxx - 2 * sx * rx + 4 * np.pi * G * R0 * px**2 + G / R0 * t
        sc = np.abs(rxx) + np.abs(2 * sx * rx) + np.abs(4 * np.pi * G * R0 * px**2) + np.abs(G / R0 * t) + 1e-300
        m = (R0 > r_lo) & (R0 < r_hi) & (np.abs(sl(S, 0)) < s_max)
        out["E" + key] = float(np.nanmax(np.abs(np.where(m, E, np.nan))))
        out["E" + key + "_rel"] = float(np.nanmedian(np.abs(np.where(m, E, np.nan)) / sc))
    return out


# ------------------------------------------------------------------ Stage A: regular centre (scalar carried as psi = r phi)
@njit(cache=True)
def _diamond_c(X, Pc, h, G, lam, vv, rho_f, lam_f, clip, aL, aH, kg2, Gam, T_floor, th, n_iter, near_axis):
    """Centre-stage diamond.  X[q, c], q = r, s, psi, eta_c, D, T_c, lnT, eta_f (psi = r phi); Pc[c] = phi at the corners (used only in
    the plasma fluxes, which vanish on the axis).  Near the axis the combination (r_u r_v + e^{2s}/4)/r is replaced by its regular limit
    (2 pi G/3) e^{2s} rho_tot r."""
    NQ = 8
    out = np.empty(NQ)
    for q in range(NQ):
        out[q] = X[q, 0] + X[q, 1] - X[q, 2]
    c = np.empty(NQ)
    for it in range(n_iter):
        for q in range(NQ):
            c[q] = 0.25 * (out[q] + X[q, 0] + X[q, 1] + X[q, 2])
        ru = (out[0] + X[0, 1] - X[0, 0] - X[0, 2]) / (2.0 * h)
        rv = (out[0] + X[0, 0] - X[0, 1] - X[0, 2]) / (2.0 * h)
        su = (out[1] + X[1, 1] - X[1, 0] - X[1, 2]) / (2.0 * h)
        sv = (out[1] + X[1, 0] - X[1, 1] - X[1, 2]) / (2.0 * h)
        qu = (out[2] + X[2, 1] - X[2, 0] - X[2, 2]) / (2.0 * h)
        qv = (out[2] + X[2, 0] - X[2, 1] - X[2, 2]) / (2.0 * h)
        rc, sc = c[0], c[1]
        phc = c[2] / rc
        pu = (qu - phc * ru) / rc
        pv = (qv - phc * rv) / rc
        a2 = math.exp(2.0 * sc)
        rho_eff = rho_f + lam_f * c[5]
        Tc_ = math.exp(c[6])
        Af, Buf, Bvf, pf, Fph = _fluid_vars(rc, sc, Tc_, c[7], phc, aL, aH, kg2, vv)
        rhof, pf2, Fph2, Fpp2 = _eos(Tc_, phc, aL, aH, kg2, vv)
        uu = math.exp(-sc - c[7])
        uv = math.exp(-sc + c[7])
        fric = Gam * (uu * pu + uv * pv)
        Uc = _Uc(phc, lam, vv, rho_eff, clip)
        Lam = 8.0 * math.pi * G * Uc
        T2uv = 0.25 * c[4] * math.exp(c[1] + c[3])
        if near_axis:
            ncond = c[4] * math.exp(c[3] - sc) / (4.0 * math.pi * rc * rc)
            rho_tot = Uc + (pu * pu + pv * pv) / a2 + rhof + ncond
            reg_r = -(2.0 * math.pi * G / 3.0) * a2 * rho_tot * rc
            reg_s = (2.0 * math.pi * G / 3.0) * a2 * rho_tot
        else:
            reg_r = -(ru * rv + 0.25 * a2) / rc
            reg_s = (ru * rv + 0.25 * a2) / (rc * rc)
        Fr = reg_r + a2 * Lam * rc / 4.0 + G * T2uv / rc + 4.0 * math.pi * G * Af / rc
        As = reg_s - 4.0 * math.pi * G * pu * pv - G * T2uv / (rc * rc) - 4.0 * math.pi * G * Af / (rc * rc) - 2.0 * math.pi * G * a2 * pf
        out[0] = X[0, 0] + X[0, 1] - X[0, 2] + h * h * Fr
        out[1] = X[1, 0] + X[1, 1] - X[1, 2] + h * h * As
        # psi_uv = phi r_uv - (e^{2s}/4) r (U' + F_phi + friction); potential and thermal force implicit in psi_N
        rcn = 0.25 * (out[0] + X[0, 0] + X[0, 1] + X[0, 2])
        base = X[2, 0] + X[2, 1] - X[2, 2] - 0.25 * h * h * a2 * rcn * fric
        S3 = X[2, 0] + X[2, 1] + X[2, 2]
        qN = out[2]
        for nk in range(6):
            pcn = 0.25 * (qN + S3) / rcn
            rh, pp, Fpn, Fppn = _eos(Tc_, pcn, aL, aH, kg2, vv)
            g = qN - base - h * h * pcn * Fr + 0.25 * h * h * a2 * rcn * (_dUc(pcn, lam, vv, rho_eff, clip) + Fpn)
            gp = 1.0 - h * h * Fr / (4.0 * rcn) + h * h * a2 * (_d2Uc(pcn, lam, vv, rho_eff, clip) + Fppn) / 16.0
            if gp < 0.25:
                gp = 0.25
            qN = qN - g / gp
        out[2] = qN
        phN = qN / out[0] if out[0] > 0.0 else phc
        # condensate and clock
        lmb = math.exp(2.0 * c[3])
        sN = out[1]
        eE, eW, eS = X[3, 0], X[3, 1], X[3, 2]
        sE, sW, sS = X[1, 0], X[1, 1], X[1, 2]
        out[3] = -((lmb - 1.0) * sN + (eW - sW) - (eE - sE) - (eS - sS) + lmb * ((eE + sE) - (eW + sW) - (eS + sS))) / (1.0 + lmb)
        DE, DW, DS = X[4, 0], X[4, 1], X[4, 2]
        src = -4.0 * math.pi * a2 * rc * rc * lam_f * _chi(phc, vv, clip)
        out[4] = (DE - DW + DS - DE * math.exp(2.0 * eE) + DW * math.exp(2.0 * eW) + DS * math.exp(2.0 * eS) + 2.0 * h * src) \
            / (1.0 + math.exp(2.0 * out[3]))
        TE_, TW, TS = X[5, 0], X[5, 1], X[5, 2]
        out[5] = (2.0 * h * math.exp(c[1] + c[3]) - TW + TE_ + TS - lmb * (TE_ - TW - TS)) / (1.0 + lmb)
        # plasma
        Qu = -(Fph + fric) * pu
        Qv = -(Fph + fric) * pv
        S1 = 2.0 * su * Af - rc * ru * pf * a2 - 0.5 * a2 * rc * rc * Qu
        S2 = 2.0 * sv * Af - rc * rv * pf * a2 - 0.5 * a2 * rc * rc * Qv
        AE, BuE, BvE, p_, F_ = _fluid_vars(X[0, 0], X[1, 0], math.exp(X[6, 0]), X[7, 0], Pc[0], aL, aH, kg2, vv)
        AW, BuW, BvW, p_, F_ = _fluid_vars(X[0, 1], X[1, 1], math.exp(X[6, 1]), X[7, 1], Pc[1], aL, aH, kg2, vv)
        AS, BuS, BvS, p_, F_ = _fluid_vars(X[0, 2], X[1, 2], math.exp(X[6, 2]), X[7, 2], Pc[2], aL, aH, kg2, vv)
        lT, et = _fluid_solve(out[0], out[1], phN, AE, AW, AS, BuE, BuW, BuS, BvE, BvW, BvS, 2.0 * h * S1, 2.0 * h * S2,
                              out[6], out[7], aL, aH, kg2, vv, T_floor, th)
        out[6] = lT
        out[7] = et
    return out


@njit(cache=True)
def _evolve_c(Z0, h, n_rays, G, lam, vv, rho_f, lam_f, clip, aL, aH, kg2, Gam, T_floor, th, n_iter, k_axis, r_exc, r_inf, phi_max):
    """Regular-centre evolution from the initial outgoing cone Z0[q, j] (u = 0); the axis of ray i is at j = i (u = v there)."""
    NQ = Z0.shape[0]
    N = Z0.shape[1] - 1
    Z = np.full((NQ, n_rays + 1, N + 1), np.nan)
    Z[:, 0, :] = Z0
    X = np.empty((NQ, 3))
    Pc = np.empty(3)
    Fs_prev = np.zeros(N + 1)
    for i in range(1, n_rays + 1):
        if i >= N or not (math.isfinite(Z[0, i - 1, i]) and math.isfinite(Z[0, i - 1, i + 1])):
            break
        # axis point (i, i)
        Z[0, i, i] = 0.0
        Z[2, i, i] = 0.0
        Z[3, i, i] = 0.0
        Z[4, i, i] = 0.0
        Z[7, i, i] = 0.0
        Z[6, i, i] = Z[6, i - 1, i]
        s_ax = 2.0 * Z[1, i - 1, i] - Z[1, i - 1, i - 1] if i >= 2 else Z[1, 0, 1]
        Z[1, i, i] = s_ax
        Z[5, i, i] = Z[5, i - 1, i - 1] + h * math.exp(0.5 * (s_ax + Z[1, i - 1, i - 1]))   # clock first (refined below)
        for corr in range(2):
            for q in range(NQ):
                X[q, 0] = Z[q, i - 1, i + 1]
                X[q, 1] = Z[q, i, i]
                X[q, 2] = Z[q, i - 1, i]
            Pc[0] = Z[2, i - 1, i + 1] / Z[0, i - 1, i + 1]
            Pc[1] = 0.0
            Pc[2] = Z[2, i - 1, i] / Z[0, i - 1, i] if Z[0, i - 1, i] > 0.0 else Pc[0]
            o = _diamond_c(X, Pc, h, G, lam, vv, rho_f, lam_f, clip, aL, aH, kg2, Gam, T_floor, th, n_iter, True)
            rv_axis = (o[0] + Z[0, i - 1, i]) / (2.0 * h)
            if rv_axis > 0.0:
                Z[1, i, i] = math.log(2.0 * rv_axis)
        Z[5, i, i] = Z[5, i - 1, i - 1] + h * math.exp(0.5 * (Z[1, i, i] + Z[1, i - 1, i - 1]))
        for j in range(i + 1, N + 1):
            if not (math.isfinite(Z[0, i - 1, j]) and math.isfinite(Z[0, i - 1, j - 1]) and math.isfinite(Z[0, i, j - 1])):
                break
            for q in range(NQ):
                X[q, 0] = Z[q, i - 1, j]
                X[q, 1] = Z[q, i, j - 1]
                X[q, 2] = Z[q, i - 1, j - 1]
            Pc[0] = X[2, 0] / X[0, 0] if X[0, 0] > 0.0 else 0.0
            Pc[1] = X[2, 1] / X[0, 1] if X[0, 1] > 0.0 else 0.0
            Pc[2] = X[2, 2] / X[0, 2] if X[0, 2] > 0.0 else 0.0
            o = _diamond_c(X, Pc, h, G, lam, vv, rho_f, lam_f, clip, aL, aH, kg2, Gam, T_floor, th, n_iter, j - i <= k_axis)
            if j == i + 1 and math.isfinite(o[6]):
                # T is even in r at a regular centre: the first off-axis value on this ray fixes the axis value to O(h^2); redo the cell
                Z[6, i, i] = o[6]
                X[6, 1] = o[6]
                o = _diamond_c(X, Pc, h, G, lam, vv, rho_f, lam_f, clip, aL, aH, kg2, Gam, T_floor, th, n_iter, True)
            ok = True
            for q in range(NQ):
                if not math.isfinite(o[q]):
                    ok = False
            if ok and (abs(o[1]) > 500.0 or abs(o[3]) > 30.0 or (o[0] > 0.0 and abs(o[2] / o[0]) > phi_max)):
                ok = False
            if ok and o[0] < r_exc and o[0] < Z[0, i, j - 1]:
                ok = False
            if ok and o[0] > r_inf:
                ok = False
            if not ok:
                break
            for q in range(NQ):
                Z[q, i, j] = o[q]
    return Z


def centre_cone_data(h, v_max, plasma, phi_r, T_r, n_r=None, N_v=None, G=1.0):
    """Initial outgoing cone from a regular centre, gauge alpha = 1 (s = 0, r(0) = 0, r_v(0) = 1/2); profiles in areal radius:
    phi(r), plasma temperature T(r) (at rest, eta_f = 0), condensate density n(r) (at rest, eta_c = 0); ingoing null shell N_v(v)
    (= 4 pi r^2 T_vv).  The constraint is r_vv = -4 pi G r T_vv,total."""
    from scipy.integrate import solve_ivp
    n_r = n_r or (lambda r: 0.0 * r)
    N_v = N_v or (lambda v: 0.0 * v)
    N = int(round(v_max / h))
    v = np.arange(N + 1) * h
    dr = 1e-6
    dphi = lambda r: (phi_r(r + dr) - phi_r(r - dr)) / (2 * dr)
    w = lambda r: sum(plasma.eos(float(T_r(r)), float(phi_r(r)))[:2])

    def rhs(x, y):
        r, rv = y
        r = max(r, 0.0)
        Tvv = dphi(r) ** 2 * rv * rv + 0.25 * (w(r) + float(n_r(r))) + (N_v(x) / (4 * np.pi * r * r) if r > 1e-9 else 0.0)
        return [rv, -4 * np.pi * G * r * Tvv]
    sol = solve_ivp(rhs, (0.0, N * h), [0.0, 0.5], t_eval=v, method="DOP853", rtol=1e-11, atol=1e-13, max_step=h / 2)
    r = sol.y[0]
    if len(r) < N + 1:
        raise ValueError("cone integration failed")
    phi = phi_r(r)
    Z0 = np.vstack([r, 0 * v, r * phi, 0 * v, 4 * np.pi * r**2 * n_r(r), v / 2, np.log(np.vectorize(T_r)(r)), 0 * v])
    return Z0


def evolve_centre(Z0, h, pot, plasma, n_rays=None, lam_f=0.0, clip=1, G=1.0, n_iter=3, k_axis=2, r_inf=1e4, phi_max=None, theta=0.6):
    aL, aH, kg2, Gam, Tfl = plasma.args()
    N = Z0.shape[1] - 1
    n_rays = N - 1 if n_rays is None else min(N - 1, n_rays)
    Z = _evolve_c(np.ascontiguousarray(Z0, float), h, n_rays, G, pot.lam, pot.v, pot.rho_f, lam_f, int(clip), aL, aH, kg2, Gam, Tfl,
                  theta, n_iter, k_axis, 2 * h, r_inf, phi_max if phi_max is not None else 10 * pot.v)
    with np.errstate(invalid="ignore", divide="ignore"):
        phi = np.where(Z[0] > 0, Z[2] / Z[0], np.nan)
    return dict(r=Z[0], s=Z[1], psi=Z[2], phi=phi, eta=Z[3], D=Z[4], T=Z[5], lnTf=Z[6], etaf=Z[7], h=h, G=G, Q=0.0, c=0.0)


def handoff(resA, i1, j1, Lu, Lv=None):
    """Exact handoff from the centre stage to a corner stage at the grid point (i1, j1) (u1 = i1 h, v1 = j1 h): the outgoing ray u = u1
    for v >= v1 and the ingoing ray v = v1 for u1 <= u <= u1 + Lu, all eight fields, phi = psi / r.  Same grid and step."""
    h = resA["h"]
    nu = int(round(Lu / h))
    if i1 + nu >= j1:
        raise ValueError("the ingoing data ray would reach the axis; reduce Lu")
    jmax = resA["r"].shape[1] - 1 if Lv is None else j1 + int(round(Lv / h))
    fields = ("r", "s", "phi", "eta", "D", "T", "lnTf", "etaf")
    A = np.vstack([resA[f][i1, j1:jmax + 1] for f in fields])
    B = np.vstack([resA[f][i1:i1 + nu + 1, j1] for f in fields])
    if not (np.all(np.isfinite(A)) and np.all(np.isfinite(B))):
        raise ValueError("handoff rays are not fully computed in the centre stage")
    return A, B


def overlap_check(resA, resB, i1, j1, fields=("r", "s", "phi", "lnTf")):
    """Largest differences between the two stages where both are computed (the corner stage at (i - i1, j - j1))."""
    out = {}
    MB, NB = resB["r"].shape
    for f in fields:
        a = resA[f][i1:i1 + MB, j1:j1 + NB]
        b = resB[f][: a.shape[0], : a.shape[1]]
        d = np.abs(a - b)
        scale = np.nanmax(np.abs(a)) if f != "phi" else 1.0
        out[f] = float(np.nanmax(d) / (scale if scale > 0 else 1.0))
        out[f + "_points"] = int(np.sum(np.isfinite(d)))
    return out


# ------------------------------------------------------------------ item 1: can heating make a persistent, normal false-vacuum ball?
def restores_homogeneous(E, plasma, pot, t_max=150.0):
    """Homogeneous test of thermal restoration (flat, no gravity): the field starts at rest in the true vacuum +v and shares the total
    energy density E with the plasma (T from energy conservation, friction Gamma).  Returns (restored, time): restored iff the field
    ends in the false vacuum (phi < -0.9 v); time of first arrival there."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    v = pot.v

    def Tof(p, pd):
        return brentq(lambda T: plasma.eos(T, p)[0] - (E - 0.5 * pd * pd - pot.U(p)), 1e-8, 20.0)

    def dU(p):
        return pot.lam * p * (p * p - v * v) - pot.rho_f * (0.75 - 0.75 * (p / v) ** 2) / v
    try:
        s = solve_ivp(lambda t, y: [y[1], -plasma.Gamma * y[1] - dU(y[0]) - plasma.eos(Tof(y[0], y[1]), y[0])[2]], (0, t_max), [v, 0.0],
                      max_step=0.02, rtol=1e-8)
    except ValueError:
        return False, None
    k = np.where(s.y[0] < -0.9 * v)[0]
    return bool(s.y[0][-1] < -0.9 * v), (float(s.t[k[0]]) if len(k) else None)


def restoration_threshold(plasma, pot, lo=1.0, hi=6.0, n=11):
    """Smallest total energy density (in units of rho_f) that restores the false vacuum homogeneously (bisection; None if hi fails)."""
    if not restores_homogeneous(hi * pot.rho_f, plasma, pot)[0]:
        return None
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if restores_homogeneous(mid * pot.rho_f, plasma, pot)[0]:
            hi = mid
        else:
            lo = mid
    return float(hi)


def largest_normal_ball(plasma, T_hot, R_grid, h=0.05, w=0.5, G=1.0):
    """Largest radius R (from R_grid) of a hot ball (temperature T_hot inside, the floor outside, field in the true vacuum, all at rest)
    whose outgoing cone from a regular centre is normal (r_v > 0 everywhere on the cone).  None if no radius is normal."""
    best = None
    for R in R_grid:
        ball = lambda r, R=R: plasma.T_floor + (T_hot - plasma.T_floor) * 0.5 * (1 - np.tanh((np.asarray(r, float) - R) / w))
        Z0 = centre_cone_data(h, 3 * R + 20, plasma, phi_r=lambda r: plasma.v + 0 * np.asarray(r, float), T_r=ball, G=G)
        if np.min(np.gradient(Z0[0], h)) > 0:
            best = float(R)
    return best


def centre_history(res, times, vphi):
    """Along the slices t = (u + v)/2 of a centre-stage run: the largest areal radius still in the false vacuum (phi < -v/2) and whether
    any trapped point (r_u < 0, r_v < 0) exists on the slice."""
    h, R, P = res["h"], res["r"], res["phi"]
    ru = np.gradient(R, h, axis=0)
    rv = np.gradient(R, h, axis=1)
    out = []
    for t in times:
        k = int(round(2 * t / h))
        ii = np.arange(R.shape[0])
        jj = k - ii
        m = (jj >= ii) & (jj < R.shape[1])
        ii, jj = ii[m], jj[m]
        r_, p_ = R[ii, jj], P[ii, jj] / vphi
        ok = np.isfinite(r_)
        fv = ok & (p_ < -0.5)
        out.append(dict(t=float(t), r_false_vacuum=float(r_[fv].max()) if fv.any() else 0.0,
                        trapped=bool(np.any((ru[ii, jj] < 0) & (rv[ii, jj] < 0) & ok)), computed=bool(ok.any())))
    return out
