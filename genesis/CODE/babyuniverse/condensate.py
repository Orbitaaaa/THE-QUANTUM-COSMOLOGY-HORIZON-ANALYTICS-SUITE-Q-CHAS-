"""Round 7 (problem 1): the ghost condensate as a dynamical source in the double-null solver of ``throat_evolution``.

Infrared limit of L = M^4 P(X) - (alpha/2) M^2 (box phi)^2 around the condensate point X0 = 1/2 (phi normalised to proper time).  With
Y = X - X0 the stress tensor is a cosmological constant M^4 P(X0) (absorbed in the vacuum energy) plus a pressureless fluid
    T_ab = n u_a u_b,   n = M^4 P'(X) sqrt(2X) ~ M^4 P'' Y  (either sign; n < 0 is the NEC-violating phase),
    u_a = -d_a phi / sqrt(2X)  (geodesic: the flow lines are the condensate's clock),   div(n u) = alpha M^2 box(theta)  [shift charge];
pressure enters at O(Y^2) and alpha only through the gradient (k^2, k^4) terms that set the stability window (``ghost_condensate``).
Here the alpha source is dropped (Theorem 7-1 shows it only renormalises G for the l = 0 response) and n is conserved along the flow.

Double-null form (metric -e^{2s} du dv + r^2 dOmega^2): u_u = -(e^s/2) e^{eta}, u_v = -(e^s/2) e^{-eta}; the flow has dv/du = e^{2 eta}.
With D = 4 pi r^2 n e^{s - eta}:
    geodesic:   (eta - s)_u + e^{2 eta} (eta + s)_v = 0,
    continuity: D_u + (D e^{2 eta})_v = 0,
    4 pi r^2 T_uu = (D e^{s}/4) e^{3 eta},  4 pi r^2 T_vv = (D e^{s}/4) e^{-eta},  4 pi r^2 T_uv = (D e^{s}/4) e^{eta},  T_thth = 0,
so r_uv gains G T2_uv / r and s_uv gains -G T2_uv / r^2.  Both transport equations use the box (diamond) scheme, which is conservative
for D and unconditionally stable; the gravity update is the three-pass diamond of ``throat_evolution``.
"""

from __future__ import annotations

import math

import numpy as np

from . import throat_evolution as TE
from .throat_evolution import njit
from .throat_mol import _line


@njit(cache=True)
def _chi(p, vv, clip=0):
    """False-vacuum indicator (1 at phi = -v, 0 at +v): the tilt shape of the Round 5 potential.  clip = 1 holds it at 1 for phi < -v and
    at 0 for phi > +v (C^1), so that a destabilised false vacuum can only roll to the true vacuum."""
    x = p / vv
    if clip == 1:
        if x <= -1.0:
            return 1.0
        if x >= 1.0:
            return 0.0
    return 0.5 - 0.75 * x + 0.25 * x ** 3


@njit(cache=True)
def _dchi(p, vv, clip=0):
    x = p / vv
    if clip == 1 and (x <= -1.0 or x >= 1.0):
        return 0.0
    return (-0.75 + 0.75 * x * x) / vv


@njit(cache=True)
def _d2Uc(p, lam, vv, rho_f, clip):
    x = p / vv
    tilt = 0.0 if (clip == 1 and (x <= -1.0 or x >= 1.0)) else 1.5 * x / (vv * vv)
    return lam * (3.0 * p * p - vv * vv) + rho_f * tilt


@njit(cache=True)
def _Uc(p, lam, vv, rho_f, clip):
    return lam / 4.0 * (p * p - vv * vv) ** 2 + rho_f * _chi(p, vv, clip)


@njit(cache=True)
def _dUc(p, lam, vv, rho_f, clip):
    return lam * p * (p * p - vv * vv) + rho_f * _dchi(p, vv, clip)


@njit(cache=True)
def _rhs(r, ru, rv, s, su, sv, p, pu, pv, T2uv, G, Q2, k, lam, vv, rho_f, clip):
    a2 = math.exp(2.0 * s)
    Lam = 8.0 * math.pi * G * _Uc(p, lam, vv, rho_f, clip)
    W = Lam * r * r + Q2 / (r * r) - 1.0
    As = ru * rv / (r * r) + a2 / (4.0 * r * r) - 4.0 * math.pi * G * pu * pv - a2 * Q2 / (2.0 * r ** 4) - G * T2uv / (r * r)
    Fs = As / (1.0 - k / (r * r))
    Fr = -ru * rv / r + a2 * W / (4.0 * r) - (k / r) * Fs + G * T2uv / r
    Fp = -(ru * pv + rv * pu) / r - 0.25 * a2 * _dUc(p, lam, vv, rho_f, clip)
    return Fr, Fs, Fp


@njit(cache=True)
def _diamond(X, h, G, Q2, k, lam, vv, rho_f, lam_f, lam_s, n_iter, clip):
    """X[q, c]: q = r, s, phi, eta, D, T (condensate clock); c = E, W, S.  Returns the six values at N.  lam_f: tilt of the condensate
    potential V = lam_f chi(phi) T in the false vacuum: the false-vacuum energy is rho_f + lam_f T and div(n u) = -lam_f chi(phi)."""
    out = np.empty(6)
    for q in range(6):
        out[q] = X[q, 0] + X[q, 1] - X[q, 2]
    for it in range(n_iter):
        c = np.empty(6)
        for q in range(6):
            c[q] = 0.25 * (out[q] + X[q, 0] + X[q, 1] + X[q, 2])
        ru = (out[0] + X[0, 1] - X[0, 0] - X[0, 2]) / (2.0 * h)
        rv = (out[0] + X[0, 0] - X[0, 1] - X[0, 2]) / (2.0 * h)
        su = (out[1] + X[1, 1] - X[1, 0] - X[1, 2]) / (2.0 * h)
        sv = (out[1] + X[1, 0] - X[1, 1] - X[1, 2]) / (2.0 * h)
        pu = (out[2] + X[2, 1] - X[2, 0] - X[2, 2]) / (2.0 * h)
        pv = (out[2] + X[2, 0] - X[2, 1] - X[2, 2]) / (2.0 * h)
        T2uv = 0.25 * c[4] * math.exp(c[1] + c[3])
        rho_eff = rho_f + lam_f * c[5]
        Fr, Fs, Fp = _rhs(c[0], ru, rv, c[1], su, sv, c[2], pu, pv, T2uv, G, Q2, k, lam, vv, rho_eff, clip)
        out[0] = X[0, 0] + X[0, 1] - X[0, 2] + h * h * Fr
        out[1] = X[1, 0] + X[1, 1] - X[1, 2] + h * h * Fs
        # scalar: the potential force is implicit in phi_N (Newton), since the explicit form is unstable once m h e^s / 2 > 1
        a2c = math.exp(2.0 * c[1])
        Fp0 = Fp + 0.25 * a2c * _dUc(c[2], lam, vv, rho_eff, clip)
        base = X[2, 0] + X[2, 1] - X[2, 2] + h * h * Fp0
        S3 = X[2, 0] + X[2, 1] + X[2, 2]
        pN = out[2]
        for nk in range(6):
            pcn = 0.25 * (pN + S3)
            g = pN - base + 0.25 * h * h * a2c * _dUc(pcn, lam, vv, rho_eff, clip)
            gp = 1.0 + h * h * a2c * _d2Uc(pcn, lam, vv, rho_eff, clip) / 16.0
            if gp < 0.25:
                gp = 0.25
            pN = pN - g / gp
        out[2] = pN
        lmb = math.exp(2.0 * c[3])
        sN = out[1]
        eE, eW, eS = X[3, 0], X[3, 1], X[3, 2]
        sE, sW, sS = X[1, 0], X[1, 1], X[1, 2]
        out[3] = -((lmb - 1.0) * sN + (eW - sW) - (eE - sE) - (eS - sS) + lmb * ((eE + sE) - (eW + sW) - (eS + sS))) / (1.0 + lmb)
        DE, DW, DS = X[4, 0], X[4, 1], X[4, 2]
        FE, FW, FS = DE * math.exp(2.0 * eE), DW * math.exp(2.0 * eW), DS * math.exp(2.0 * eS)
        src = -4.0 * math.pi * math.exp(2.0 * c[1]) * c[0] * c[0] * lam_s * _chi(c[2], vv, clip)
        out[4] = (DE - DW + DS - FE + FW + FS + 2.0 * h * src) / (1.0 + math.exp(2.0 * out[3]))
        TE_, TW, TS = X[5, 0], X[5, 1], X[5, 2]
        out[5] = (2.0 * h * math.exp(c[1] + c[3]) - TW + TE_ + TS - lmb * (TE_ - TW - TS)) / (1.0 + lmb)
    return out


@njit(cache=True)
def _evolve(A, B, h, G, Q2, k, lam, vv, rho_f, lam_f, lam_s, n_iter, r_exc, r_inf, phi_max, clip):
    """A[q, j] on u = 0, B[q, i] on v = 0 (q = r, s, phi, eta, D, T).  Returns Z[q, i, j] (NaN beyond breakdown) and end codes."""
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
            o = _diamond(X, h, G, Q2, k, lam, vv, rho_f, lam_f, lam_s, n_iter, clip)
            ok = True
            for q in range(NQ):
                if not math.isfinite(o[q]):
                    ok = False
            if not ok or abs(o[1]) > 500.0 or abs(o[2]) > phi_max or abs(o[3]) > 30.0:
                end[i] = 3
                break
            if o[0] < r_exc or o[0] * o[0] <= 1.05 * k:
                end[i] = 1
                break
            if o[0] > r_inf:
                end[i] = 2
                break
            for q in range(NQ):
                Z[q, i, j] = o[q]
    return Z, end


def evolve(A, B, Q, c, pot, h, G=1.0, n_iter=3, r_exc=1e-3, r_inf=1e6, phi_max=None, lam_f=0.0, lam_s=None, clip=0):
    """lam_f: condensate tilt (drift of the false-vacuum energy); lam_s: its charge source (= lam_f; differ only in diagnostics)."""
    k = G * c / (12 * np.pi)
    Z, end = _evolve(np.ascontiguousarray(A, float), np.ascontiguousarray(B, float), h, G, Q * Q, k, pot.lam, pot.v, pot.rho_f, lam_f,
                     lam_f if lam_s is None else lam_s, n_iter, r_exc, r_inf, phi_max if phi_max is not None else 1e300, int(clip))
    return dict(r=Z[0], s=Z[1], phi=Z[2], eta=Z[3], D=Z[4], T=Z[5], end=end, h=h, Q=Q, c=c, k=k, G=G, lam_f=lam_f)


def stress(res):
    """4 pi r^2 T_ab of the condensate (T2_uu, T2_vv, T2_uv) and its proper density n."""
    D, s, e, R = res["D"], res["s"], res["eta"], res["r"]
    base = 0.25 * D * np.exp(s)
    return dict(T2uu=base * np.exp(3 * e), T2vv=base * np.exp(-e), T2uv=base * np.exp(e), n=D * np.exp(e - s) / (4 * np.pi * R**2))


def constraint_residual(res, r_lo=0.0, r_hi=np.inf, s_max=np.inf):
    """Residuals of r_vv - 2 s_v r_v + 4 pi G r phi_v^2 + (G/r) T2_vv (and the uu one) over the interior points with r_lo < r < r_hi and
    |s| < s_max (to keep away from the singularity edge): max absolute and median relative."""
    R, S, P, h, G = res["r"], res["s"], res["phi"], res["h"], res["G"]
    T = stress(res)
    out = {}
    for ax, key in ((1, "vv"), (0, "uu")):
        sl = lambda a, o: np.take(a, np.arange(1 + o, a.shape[ax] - 1 + o), axis=ax)
        Rpp, Rm, R0 = sl(R, 1), sl(R, -1), sl(R, 0)
        rxx = (Rpp - 2 * R0 + Rm) / h**2
        rx = (Rpp - Rm) / (2 * h)
        sx = (sl(S, 1) - sl(S, -1)) / (2 * h)
        px = (sl(P, 1) - sl(P, -1)) / (2 * h)
        t = sl(T["T2" + key], 0)
        E = rxx - 2 * sx * rx + 4 * np.pi * G * R0 * px**2 + G / R0 * t
        sc = np.abs(rxx) + np.abs(2 * sx * rx) + np.abs(4 * np.pi * G * R0 * px**2) + np.abs(G / R0 * t) + 1e-300
        mask = (R0 > r_lo) & (R0 < r_hi) & (np.abs(sl(S, 0)) < s_max)
        E = np.where(mask, E, np.nan)
        out["E" + key] = float(np.nanmax(np.abs(E)))
        out["E" + key + "_rel"] = float(np.nanmedian(np.abs(E) / sc))
    return out


# ------------------------------------------------------------------ characteristic data with a condensate
def _static_rapidity(r_x, r_y, sign):
    """eta of the static (r = const) frame, 0.5 ln(-r_u / r_v), clipped where the frame does not exist."""
    num, den = (-r_x, r_y) if sign > 0 else (-r_y, r_x)
    ok = (num > 0) & (den > 0)
    e = np.where(ok, 0.5 * np.log(np.where(ok, num, 1.0) / np.where(ok, den, 1.0)), np.nan)
    e = np.clip(e, -4.0, 4.0)
    idx = np.where(np.isfinite(e))[0]           # hold the last static value beyond the edge (no condensate charge there)
    return np.interp(np.arange(len(e)), idx, e[idx]) if len(idx) else np.zeros_like(e)


def condensate_data(Q, c, r0, ru0, rv0, h, Nv, Nu, n_v, n_u, phi_v=None, phi_u=None, G=1.0, pot=None, frame="static", n_pass=4,
                    extra_Nv=None):
    """Constraint-satisfying data (gauge s = 0 on both rays) with a condensate of proper density n_v(v) on u = 0 and n_u(u) on v = 0, at
    rest in the static frame (frame = "static": eta = 0.5 ln(-r_u/r_v), r_u on u = 0 and r_v on v = 0 from the transport equations) or
    in the corner frame (frame = "corner": eta = 0).  (r, eta) are iterated; the last pass is exactly consistent: 4 pi r^2 T_vv of the
    data equals the flux used in the constraint.  Returns (A, B, info) for ``evolve``."""
    v = np.arange(Nv + 1) * h
    u = np.arange(Nu + 1) * h
    phi_v = phi_v or (lambda x: np.zeros_like(np.asarray(x, float)))
    phi_u = phi_u or (lambda x: np.zeros_like(np.asarray(x, float)))
    pot = pot or TE.Potential()
    k = G * c / (12 * np.pi)
    dphi = lambda f, x: (f(x + 1e-6 * h) - f(x - 1e-6 * h)) / (2e-6 * h)
    ea, eb = np.zeros(Nv + 1), np.zeros(Nu + 1)
    ra, rb = r0 + rv0 * v, r0 + ru0 * u
    for it in range(n_pass if frame == "static" else 2):
        Tv = np.pi * ra**2 * n_v(v) * np.exp(-2 * ea)            # 4 pi r^2 T_vv = pi r^2 n e^{-2 eta}  (s = 0)
        Tu = np.pi * rb**2 * n_u(u) * np.exp(2 * eb)
        Nvf = (lambda x, Tv=Tv: np.interp(x, v, Tv)) if extra_Nv is None else (lambda x, Tv=Tv: np.interp(x, v, Tv) + extra_Nv(x))
        Nuf = lambda x, Tu=Tu: np.interp(x, u, Tu)
        a, b, info = TE.characteristic_data(Q, c, r0, ru0, rv0, h, Nv, Nu, phi_v=phi_v, phi_u=phi_u, N_u=Nuf, N_v=Nvf, G=G)
        ea_used, eb_used, Tv_used, Tu_used = ea, eb, Tv, Tu
        ra, rb = a[0], b[0]
        if frame == "static":
            RU = _line(a[0], a[1], a[2], len(a[0]), ru0, 0.0, float(dphi(phi_u, 0.0)), h, G, Q * Q, k, pot.lam, pot.v, pot.rho_f)[0]
            RV = _line(b[0], b[1], b[2], len(b[0]), rv0, 0.0, float(dphi(phi_v, 0.0)), h, G, Q * Q, k, pot.lam, pot.v, pot.rho_f)[0]
            ea = _static_rapidity(RU, info["rv_a"], +1)
            eb = _static_rapidity(info["ru_b"], RV, +1)
    from scipy.integrate import cumulative_trapezoid
    Ta = cumulative_trapezoid(0.5 * np.exp(-ea_used), v, initial=0.0)
    Tb = cumulative_trapezoid(0.5 * np.exp(eb_used), u, initial=0.0)
    A = np.vstack([a[0], a[1], a[2], ea_used, 4 * Tv_used * np.exp(ea_used), Ta])
    B = np.vstack([b[0], b[1], b[2], eb_used, 4 * Tu_used * np.exp(-3 * eb_used), Tb])
    info.update(T2vv=Tv_used, T2uu=Tu_used, v=v, u=u)
    return A, B, info


def nariai_condensate(n0, l_H=30.0, delta=1e-4, vphi=0.0015, eps=4.0, wall_frac=0.95, h=0.04, Lu=200.0, Lv=200.0, u_end_frac=0.6,
                      w_end=3.0, v_end=None, frame="static", profile="fill", lump_width=4.0, r_inf=1e4, lam_f=0.0, lam_s=None,
                      dphi0=0.0, clip=0):
    """The near-Nariai light-wall configuration of Round 6d (``TE.nariai_false_vacuum`` with delta = 1e-4, eps = 4, vphi = 1.5e-3) with
    NO prescribed flux: instead the ghost condensate, at rest in the static frame, carries the proper NEC-phase density n0 (< 0).
    profile = "fill": n0 throughout the static false-vacuum region on the data (v < v_wall on u = 0, u < u_end_frac * u_hole on v = 0,
    smooth edges of width w_end); profile = "lump": a Gaussian of width lump_width (in u and v) at the corner.  Returns the diagnostics of
    ``TE.nariai_false_vacuum`` plus the condensate's own fluxes and weighted deficits.  lam_f > 0: the condensate's potential is tilted
    in the false vacuum (V = lam_f chi(phi) T), so the condensate makes its own NEC-phase charge there (n0 = 0 allowed).
    Note: unlike the Round 6d driver, the corner sits at r_N exactly and the parent mass is M_N (1 - delta).""" 
    from scipy.optimize import brentq
    Lam = 3 / l_H**2
    rho_f = Lam / (8 * np.pi)
    lam = 4 * rho_f / (eps * vphi**4)
    w = np.sqrt(2) / (np.sqrt(lam) * vphi)
    pot = TE.Potential(lam, vphi, rho_f)
    rN = 1 / np.sqrt(Lam)
    MN = 1 / (3 * np.sqrt(Lam))
    M = MN * (1 - delta)
    f = lambda r: 1 - 2 * M / r - Lam * r**2 / 3
    f0 = f(rN)
    rc = brentq(f, rN, 3 * rN)
    rgrid = np.linspace(0.3 * rN, rN, 4001)
    sgn = np.where(np.diff(np.sign(f(rgrid))) != 0)[0]
    rb = brentq(f, rgrid[sgn[-1]], rgrid[sgn[-1] + 1])
    rv0 = np.sqrt(f0) / 2
    v2 = wall_frac * (rc - rN) / rv0
    u_hole = (rN - rb) / rv0
    # dphi0: displacement off the false vacuum (for eps > 16/3 the false vacuum is a hilltop and the daughter exits inflation classically)
    phi_v = lambda vv: vphi * np.tanh((np.asarray(vv, float) - v2) / w) + dphi0 * 0.5 * (1 - np.tanh((np.asarray(vv, float) - v2) / w))
    phi_u = lambda uu: -vphi + dphi0 + 0 * np.asarray(uu, float)
    ve = v2 if v_end is None else v_end
    ue = u_end_frac * u_hole
    if profile == "fill":
        n_v = lambda x: n0 * 0.5 * (1 - np.tanh((np.asarray(x, float) - ve) / w_end))
        n_u = lambda x: n0 * 0.5 * (1 - np.tanh((np.asarray(x, float) - ue) / w_end))
    else:
        n_v = lambda x: n0 * np.exp(-(np.asarray(x, float) / lump_width) ** 2)
        n_u = lambda x: n0 * np.exp(-(np.asarray(x, float) / lump_width) ** 2)
    Nv, Nu = int(Lv / h), int(Lu / h)
    A, B, info = condensate_data(0.0, 0.0, rN, -rv0, rv0, h, Nv, Nu, n_v, n_u, phi_v=phi_v, phi_u=phi_u, pot=pot, frame=frame)
    res = evolve(A, B, 0.0, 0.0, pot, h, r_inf=r_inf, phi_max=10 * vphi, lam_f=lam_f, lam_s=lam_s, clip=clip)
    cl = TE.classify(res)
    g = TE.inflated(res, cl, rN)
    P, R, Cc = res["phi"], res["r"], cl["cls"]
    # true vacuum = the field sits at the true minimum (|phi/v - 1| < 0.5).  phi > 0 alone also flags the grid-scale checkerboard of the
    # unresolved, extremely blue-shifted wall (phi ~ 3v alternating with -v at s ~ 5, F42); it is reported separately as phi_positive_flags.
    tv = (np.abs(P / vphi - 1) < 0.5) & np.isfinite(R)
    tv_any = (P > 0) & np.isfinite(R)
    tau = 3.2 * rN
    u, v = info["u"], info["v"]
    W_u = float(-np.trapezoid(np.minimum(info["T2uu"], 0) * np.exp(-u / tau), u))
    W_v = float(-np.trapezoid(np.minimum(info["T2vv"], 0) * np.exp(-v / tau), v))
    T = stress(res)
    n_corner = T["n"][:, 0]
    return dict(res=res, pot=pot, cls=Cc, grew=g["grew"], n_big=g["n_big"], r_max=g["r_max_antifalse"], r_N=float(rN), M_N=float(MN), f0=float(f0),
                r_c=float(rc), r_b=float(rb), u_hole=float(u_hole), v_wall=float(v2), W_u=W_u, W_v=W_v, n0=float(n0),
                neg_energy_u=float(-np.trapezoid(np.minimum(info["T2uu"], 0), u)), true_vacuum_untrapped_large_r=int(np.sum(tv & (Cc >= 0) & (R > 5 * rN))),
                true_vacuum_trapped=int(np.sum(tv & (Cc == -1))),
                phi_positive_flags_untrapped_large_r=int(np.sum(tv_any & (Cc >= 0) & (R > 5 * rN))), constraint=constraint_residual(res) if h >= 0.04 and Lu <= 60 else None,
                n_ingoing_ray=n_corner, info=info)


def route_condensate(lam_f=0.0, n0=0.0, E=7.258540773391722, l_H=30.0, r0=15.0, eps=4.0, vphi=0.0015, h=0.04, Lu=30.0, Lv=600.0,
                     v_shell=4.5, w_shell=0.25, wall_frac=0.6, delta_target=4e-3, r_inf=1e4, clip=0, dphi0=0.0):
    """The complete route of Round 6d (``preparation.build``: corner inside a false-vacuum de Sitter ball, an imploding shell of ordinary
    radiation of local energy E, a light wall, normal true vacuum; no black hole and no bridge on the data), with the push replaced by the
    ghost condensate: n0 = initial NEC-phase density in the false vacuum (corner frame), lam_f = tilt that makes the charge dynamically.
    clip = 1 with dphi0 > 0 (the false vacuum displaced by dphi0 inside the wall) lets the tilted false vacuum exit (Round 6f, item 4)."""
    from scipy.optimize import brentq
    Lam = 3 / l_H**2
    rho_f = Lam / (8 * np.pi)
    lam = 4 * rho_f / (eps * vphi**4)
    w = np.sqrt(2) / (np.sqrt(lam) * vphi)
    pot = TE.Potential(lam, vphi, rho_f)
    rN, MN = 1 / np.sqrt(Lam), 1 / (3 * np.sqrt(Lam))
    rv0 = np.sqrt(1 - Lam * r0**2 / 3) / 2
    M = MN * (1 - delta_target)
    rc = brentq(lambda r: 1 - 2 * M / r - Lam * r**2 / 3, rN, 3 * rN)
    shell = lambda v: E / (np.sqrt(np.pi) * w_shell) * np.exp(-((np.asarray(v, float) - v_shell) / w_shell) ** 2)
    Nv, Nu = int(Lv / h), int(Lu / h)
    a0, _, _ = TE.characteristic_data(0.0, 0.0, r0, -rv0, rv0, h, Nv, Nu, phi_v=lambda v: -vphi + 0 * np.asarray(v, float),
                                      phi_u=lambda u: -vphi + 0 * np.asarray(u, float), N_v=shell)
    j = np.where(a0[0] >= rN + wall_frac * (rc - rN))[0]
    v_wall = j[0] * h
    phi_v = lambda v: vphi * np.tanh((np.asarray(v, float) - v_wall) / w) + dphi0 * 0.5 * (1 - np.tanh((np.asarray(v, float) - v_wall) / w))
    phi_u = lambda u: -vphi + dphi0 + 0 * np.asarray(u, float)
    n_v = lambda x: n0 * 0.5 * (1 - np.tanh((np.asarray(x, float) - v_wall) / w))
    n_u = lambda x: n0 + 0 * np.asarray(x, float)
    A, B, info = condensate_data(0.0, 0.0, r0, -rv0, rv0, h, Nv, Nu, n_v, n_u, phi_v=phi_v, phi_u=phi_u, pot=pot, frame="corner",
                                 extra_Nv=shell)
    res = evolve(A, B, 0.0, 0.0, pot, h, r_inf=r_inf, phi_max=10 * vphi, lam_f=lam_f, clip=clip)
    cl = TE.classify(res)
    g = TE.inflated(res, cl, rN)
    R = res["r"]
    rows = []
    for i in (50, 80):
        js = np.where(np.isfinite(R[i]))[0]
        rows.append(float(R[i, js.max()]))
    return dict(res=res, pot=pot, cls=cl["cls"], grew=g["grew"], n_big=g["n_big"], r_max=g["r_max_antifalse"], r_N=float(rN), M_N=float(MN),
                v_wall=float(v_wall), normal_on_data=bool(np.all(info["ru_b"] < 0) and np.all(info["rv_a"] > 0)),
                r_at_end_of_rows_u2_u3p2=rows, lam_f=lam_f, n0=n0)


def hiding_report(o, vphi=0.0015):
    """Is the daughter behind the parent's event horizon?  u_EH: last outgoing row that reaches the outer edge v = L_v with r increasing in
    true vacuum (``gsl_solution.parent_event_horizon``); a daughter is hidden iff none of its rows reaches the outer edge.  Also counts
    true-vacuum points (|phi/v - 1| < 0.5) that are locally untrapped at large r after u_EH, and how their rows end (1 = singularity)."""
    from . import gsl_solution as GS
    R, P, Cc, end = o["res"]["r"], o["res"]["phi"], o["cls"], o["res"]["end"]
    h, rN = o["res"]["h"], o["r_N"]
    eh = GS.parent_event_horizon(dict(res=o["res"], cls=Cc, r_N=rN))
    i_eh = int(round((eh.get("u_EH") or 0.0) / h))
    dau = (Cc == 1) & (P < 0) & np.isfinite(R) & (R > 5 * rN)
    rows_d = np.unique(np.argwhere(dau)[:, 0])
    tv = (np.abs(P / vphi - 1) < 0.5) & np.isfinite(R) & (Cc >= 0) & (R > 5 * rN)
    tv[: i_eh + 1] = False
    rows_tv = np.unique(np.argwhere(tv)[:, 0])
    return dict(u_EH=eh.get("u_EH"), parent_deficit=eh.get("deficit"), first_daughter_u=float(rows_d.min() * h) if len(rows_d) else None,
                daughter_rows=int(len(rows_d)), daughter_rows_reaching_outer_edge=int(sum(np.isfinite(R[i, -1]) for i in rows_d)),
                hidden=bool(len(rows_d) > 0 and all(not np.isfinite(R[i, -1]) for i in rows_d) and (eh.get("u_EH") is not None)
                            and rows_d.min() > i_eh),
                untrapped_true_vacuum_large_r_inside=int(tv.sum()), their_rows_reaching_outer_edge=int(sum(np.isfinite(R[i, -1]) for i in rows_tv)))
