"""Semiclassical tunnelling between thin-wall branches (TUNNELLING.md).

Canonical (reduced phase space, clock = exterior Killing time, H = M):

    P = (R/G) [atanh(Rdot/A_out) - atanh(Rdot/A_in)]                (symbolic.py)

Under the barrier Rdot = i R', R'^2 = V(R) > 0 and P = i P_E with

    P_E(R) = (R/G) [Theta_out(R) - Theta_in(R)],  Theta = angle of (A, sqrt V).

The angle is defined only up to the branch of arctan.  Two prescriptions:

  * "analytic": Theta = atan2(sqrt V, A) in [0, pi], the analytic continuation
    of arctan(sqrt V / A) through A = 0 (A is a rational function of R that
    changes sign at the eps-flip radius while V stays positive there);
  * "principal": Theta = arctan(sqrt V / A) in (-pi/2, pi/2), which jumps by
    pi wherever A changes sign.

B = 2 int_{R1}^{R2} P_E dR  (full Euclidean bounce, there and back).

Results implemented and verified here (RESULTS.md):

  R1  M = 0 closed form: B_an = S_E[O(4) instanton] + (Sum_max A - Sum_min A)/4G
      over the extremal 2-spheres of the time-symmetric equatorial slice.
  R2  First law:  dBt = Vol_in/(8 pi G) dLambda_in + A_W dsigma - dt_W dM,
      Bt = B_an - (endpoint area terms), dt_W = exterior Euclidean Killing time
      elapsed along the wall.
  R3  Smarr identity (Euler homogeneity + R2):
      Bt = -(Lambda_in/8 pi G) Vol_in - (sigma/2) A_W - (M/2) dt_W
           + (Lambda_out/8 pi G) W_out       (W_out: exterior volume excised, eta > 0).
  R4  Embedding obstruction: if eps_out (eps_in) changes sign in the barrier,
      the exterior (interior) Euclidean region is not embedded (pseudomanifold).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.integrate import quad, solve_ivp

from .model import ThinWallModel


# ----------------------------------------------------------- canonical momentum
def theta(A, sqrtV, prescription="analytic"):
    if prescription == "analytic":
        return np.arctan2(sqrtV, A)
    if prescription == "principal":
        return np.arctan(sqrtV / A)
    raise ValueError(prescription)


def euclidean_momentum(model: ThinWallModel, R, prescription="analytic"):
    R = np.asarray(R, dtype=float)
    y = np.sqrt(np.maximum(model.V(R), 0.0))
    return R / model.G * (theta(model.A_out(R), y, prescription) - theta(model.A_in(R), y, prescription))


def B_arccos(model: ThinWallModel) -> float:
    """Independent formula (Phase II): on the analytic branch
    0 < Theta_out - Theta_in < pi in the barrier (cot Theta = A/sqrt V and
    A_out < A_in), hence Theta_out - Theta_in = arccos C with
    C = (f_in + f_out - kappa^2 R^2)/(2 sqrt(f_in f_out)), and
    B = (2/G) int R arccos C dR.  Uses no angle or branch logic."""
    Ra, Rb = barrier(model)

    def C(R):
        fi, fo = float(model.f_in(R)), float(model.f_out(R))
        return (fi + fo - model.kappa**2 * R**2) / (2 * np.sqrt(fi * fo))

    lo = Ra if Ra > 0 else 1e-14 * Rb
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        val = quad(lambda R: R * np.arccos(np.clip(C(R), -1.0, 1.0)), lo, Rb, limit=800,
                   epsabs=1e-13, epsrel=1e-12)[0]
    return 2.0 * val / model.G


def theta_D_critical(gamma: float) -> float:
    """|theta_D|/pi in the limit M -> M_cr (regime B, eta = 0), closed form:
    F(t) = 3 sqrt3 (t+1)^2 / (4 (t+2) sqrt(t^2+2)),  t = z_*^3 in (1, 2].
    F > 1 for t > 1 since 27(t+1)^4 - 16(t+2)^2(t^2+2) = (t-1)(11t^3+55t^2+121t+101)."""
    t = 0.5 * (-gamma + np.sqrt(gamma**2 + 8))
    return 3 * np.sqrt(3) * (t + 1) ** 2 / (4 * (t + 2) * np.sqrt(t**2 + 2))


def barrier(model: ThinWallModel):
    """(R_a, R_b): classically forbidden interval between the two branches."""
    tps = model.turning_points()
    if model.M == 0:
        return 0.0, tps[0]
    if len(tps) < 2:
        raise ValueError("no barrier: M >= M_critical")
    return tps[0], tps[-1]


def B_canonical(model: ThinWallModel, prescription="analytic") -> float:
    Ra, Rb = barrier(model)
    pts = [x for x in (model.flip_radius_out(), model.flip_radius_in()) if x and Ra < x < Rb]

    def f(R):
        return float(euclidean_momentum(model, R, prescription))

    lo = Ra if Ra > 0 else 1e-14 * Rb
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        val, err = quad(f, lo, Rb, points=pts or None, limit=800, epsabs=1e-13, epsrel=1e-12)
    return 2.0 * val


def endpoint_index(model: ThinWallModel, R) -> int:
    """n_out - n_in at a turning point, n = 1 if eps = -1 else 0.

    P_E(R_turn) = pi R n / G with n in {-1, 0, 1}.  n = +1: the wall is a local
    maximum of the area on that time-symmetric slice; n = -1: local minimum."""
    no = 1 if model.A_out(R) < 0 else 0
    ni = 1 if model.A_in(R) < 0 else 0
    return no - ni


def endpoint_area_term(model: ThinWallModel) -> float:
    """(1/4G)[n(R_b) A(R_b) - n(R_a) A(R_a)] + horizon terms that are
    M-independent: the interior de Sitter horizon area when eps_in(R_b) = -1
    (the de Sitter equator is then a maximal sphere of the final slice)."""
    Ra, Rb = barrier(model)
    G = model.G
    t = endpoint_index(model, Rb) * np.pi * Rb**2 / G
    if Ra > 0:
        t -= endpoint_index(model, Ra) * np.pi * Ra**2 / G
    if model.A_in(Rb) < 0 and model.chi > 0:
        t += np.pi / (G * model.chi**2)
    if Ra > 0 and model.A_in(Ra) < 0 and model.chi > 0:
        t -= np.pi / (G * model.chi**2)
    if model.M == 0 and model.eta > 0 and model.A_out(Rb) > 0:
        t += np.pi / (G * model.eta**2)
    return t


# ------------------------------------------------------- Euclidean trajectory
@dataclass
class EuclideanWall:
    s: np.ndarray            # parameter in [0, pi]: R = Ra + (Rb-Ra)(1-cos s)/2
    R: np.ndarray
    tau: np.ndarray          # Euclidean proper time along the wall
    t_out: np.ndarray        # exterior Euclidean Killing time
    t_in: np.ndarray         # interior Euclidean Killing time
    theta_out: Optional[np.ndarray]   # exterior angle around the Euclidean horizon
    theta_in: Optional[np.ndarray]    # interior angle around the de Sitter horizon


def _sqrtV_from_distances(model, R, d_a, d_b, Ra, Rb):
    """sqrt(V) in the barrier from the distances d_a = R - Ra, d_b = Rb - R
    (computed without cancellation by the caller).  Within a relative distance
    1e-6 of a turning point V is replaced by its linear expansion, so round-off
    in V near R_turn cannot produce spurious 0/0."""
    span = Rb - Ra
    if Ra > 0 and d_a < 1e-6 * span:
        return np.sqrt(max(float(model.dV(Ra)) * d_a, 0.0))
    if d_b < 1e-6 * span:
        return np.sqrt(max(-float(model.dV(Rb)) * d_b, 0.0))
    return np.sqrt(max(float(model.V(R)), 0.0))


def _quad_s(model, fun, Ra, Rb, s0=0.0, s1=np.pi):
    """int_{R(s0)}^{R(s1)} fun(R) dR / sqrt(V) with R = Ra + (Rb-Ra) sin^2(s/2):
    dR/sqrt(V) is smooth in s and finite at the turning points, where it tends
    to sqrt((Rb-Ra)/|V'(R_turn)|)."""
    span = Rb - Ra

    def g(s):
        d_a = span * np.sin(s / 2) ** 2
        d_b = span * np.cos(s / 2) ** 2
        R = Ra + d_a
        if Ra > 0 and d_a < 1e-6 * span:
            return fun(Ra + d_a) * np.sqrt(span / float(model.dV(Ra))) * np.cos(s / 2)
        if d_b < 1e-6 * span:
            return fun(Rb - d_b) * np.sqrt(span / -float(model.dV(Rb))) * np.sin(s / 2)
        dR = span * np.sin(s) / 2
        return fun(R) * dR / max(_sqrtV_from_distances(model, R, d_a, d_b, Ra, Rb), 1e-300)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return quad(g, s0, s1, limit=800, epsabs=1e-13, epsrel=1e-12)[0]


def euclidean_wall(model: ThinWallModel, n: int = 1500) -> EuclideanWall:
    """Wall trajectory through the barrier (half bounce), R from Ra to Rb.

    dt_E/dtau_E = A/f on each side (Euclidean continuation of f tdot = A);
    the substitution R(s) removes the integrable 1/sqrt(V) endpoint
    singularities."""
    Ra, Rb = barrier(model)
    s = np.linspace(0, np.pi, n)
    R = Ra + (Rb - Ra) * (1 - np.cos(s)) / 2

    def cumulative(fun):
        out = np.zeros_like(s)
        for i in range(1, n):
            out[i] = out[i - 1] + _quad_s(model, fun, Ra, Rb, s[i - 1], s[i])
        return out

    tau = cumulative(lambda r: 1.0)
    t_out = cumulative(lambda r: model.A_out(r) / model.f_out(r))
    t_in = cumulative(lambda r: model.A_in(r) / model.f_in(r)) if model.chi > 0 or True else None
    th_out = None
    hz = model.horizons_out()
    if model.M > 0 and model.eta == 0:
        th_out = t_out / (4 * model.G * model.M)
    elif model.M == 0 and model.eta > 0:
        th_out = model.eta * t_out
    th_in = model.chi * t_in if model.chi > 0 else None
    del hz
    return EuclideanWall(s, R, tau, t_out, t_in, th_out, th_in)


def euclidean_pieces(model: ThinWallModel) -> dict:
    """Full-bounce Euclidean quantities entering the Smarr identity.

    Vol_in uses Stokes with the 1-form (r^3 - r_h^3)/3 dt (regular at the
    de Sitter Euclidean horizon r_h = 1/chi) plus the axis contribution; the
    final slice passes through the horizon iff eps_in(R_b) = -1, which adds
    half a period pi/chi to the axis extent."""
    Ra, Rb = barrier(model)
    G = model.G
    D_out = _quad_s(model, lambda r: model.A_out(r) / model.f_out(r), Ra, Rb)
    area_half = 4 * np.pi * _quad_s(model, lambda r: r**2, Ra, Rb)
    out = dict(dt_W=2 * D_out, A_W=2 * area_half)
    if model.chi > 0:
        rh = 1.0 / model.chi
        D_in = _quad_s(model, lambda r: model.A_in(r) / model.f_in(r), Ra, Rb)
        W = _quad_s(model, lambda r: (r**3 - rh**3) / 3 * model.A_in(r) / model.f_in(r), Ra, Rb)
        ax = D_in
        if model.A_in(Rb) < 0:
            ax += np.pi / model.chi
        if Ra > 0 and model.A_in(Ra) < 0:
            ax -= np.pi / model.chi
        vol_half = 4 * np.pi * (W + rh**3 / 3 * ax)
        out.update(Vol_in=2 * vol_half, dt_W_in=2 * D_in)
    else:
        out.update(Vol_in=None, dt_W_in=None)
    out["S_in"] = -(model.Lambda_in / (8 * np.pi * G)) * out["Vol_in"] if out["Vol_in"] is not None else 0.0
    out["S_wall"] = -(model.sigma / 2) * out["A_W"]
    out["S_ext"] = -(model.M / 2) * out["dt_W"]
    # exterior vacuum energy (eta > 0): 4-volume of the region r < R excised from the
    # exterior static chart along the wall, W_out = 4 pi int (R^3/3) dt_E^out (full bounce);
    # it enters with the opposite sign to the interior volume (Phase II extension).
    if model.eta > 0:
        Wh = 4 * np.pi * _quad_s(model, lambda r: r**3 / 3 * model.A_out(r) / model.f_out(r), Ra, Rb)
        out["W_out"] = 2 * Wh
        out["S_out_bulk"] = (model.Lambda_out / (8 * np.pi * G)) * out["W_out"]
    else:
        out["W_out"] = 0.0
        out["S_out_bulk"] = 0.0
    del G
    return out


def smarr_check(model: ThinWallModel) -> dict:
    """Compare B_an - endpoint terms with the covariant Smarr expression."""
    B = B_canonical(model, "analytic")
    Bt = B - endpoint_area_term(model)
    pcs = euclidean_pieces(model)
    cov = pcs["S_in"] + pcs["S_wall"] + pcs["S_ext"] + pcs["S_out_bulk"]
    return dict(B=B, B_tilde=Bt, covariant=cov, residual=Bt - cov, pieces=pcs,
                endpoint=endpoint_area_term(model))


def first_law_check(model: ThinWallModel, h: float = 1e-5) -> dict:
    """Finite-difference derivatives of Bt versus the Euclidean pieces."""
    def Bt(m):
        return B_canonical(m) - endpoint_area_term(m)
    c, k, M, e, G = model.chi, model.kappa, model.M, model.eta, model.G
    mk = lambda c=c, k=k, M=M: ThinWallModel(c, k, M, e, G)  # noqa: E731
    dM = (Bt(mk(M=M * (1 + h))) - Bt(mk(M=M * (1 - h)))) / (2 * h * M)
    dc = (Bt(mk(c=c * (1 + h))) - Bt(mk(c=c * (1 - h)))) / (2 * h * c)
    dk = (Bt(mk(k=k * (1 + h))) - Bt(mk(k=k * (1 - h)))) / (2 * h * k)
    pcs = euclidean_pieces(model)
    out = dict(dBt_dM=dM, minus_dt_W=-pcs["dt_W"],
               dBt_dLambda=dc / (6 * c), Vol_over_8piG=pcs["Vol_in"] / (8 * np.pi * G),
               dBt_dsigma=dk * 4 * np.pi * G, A_W=pcs["A_W"])
    if e > 0:
        mke = lambda ee: ThinWallModel(c, k, M, ee, G)  # noqa: E731
        de = (Bt(mke(e * (1 + h))) - Bt(mke(e * (1 - h)))) / (2 * h * e)
        out.update(dBt_dLambda_out=de / (6 * e), minus_Wout_over_8piG=-pcs["W_out"] / (8 * np.pi * G))
    return out


# ------------------------------------------------------------- M = 0 closed form
def _stable_cap_cosines(model: ThinWallModel):
    """p = -A_out(R0), q = A_in(R0) and 1 +- p, 1 +- q without cancellation.

    p = (chi^2+kappa^2-eta^2)/L^2, q = (kappa^2+eta^2-chi^2)/L^2 and the
    identities L^4 - (chi^2+kappa^2-eta^2)^2 = 4 kappa^2 eta^2,
    L^4 - (kappa^2+eta^2-chi^2)^2 = 4 kappa^2 chi^2 give 1 -+ p, 1 -+ q as
    ratios of positive quantities when they are small."""
    c2, e2, k2 = model.chi**2, model.eta**2, model.kappa**2
    L2 = model.L2
    cp = c2 + k2 - e2
    cq = k2 + e2 - c2
    p, q = cp / L2, cq / L2

    def one_minus(cx, other):   # 1 - x with x = cx/L2, L^4 - cx^2 = 4 kappa^2 other
        return (4 * k2 * other / (L2 + cx)) / L2 if cx > 0 else (L2 - cx) / L2

    def one_plus(cx, other):
        return (4 * k2 * other / (L2 - cx)) / L2 if cx < 0 else (L2 + cx) / L2

    return dict(p=p, q=q, one_plus_p=one_plus(cp, e2), one_minus_p=one_minus(cp, e2),
                one_plus_q=one_plus(cq, c2), one_minus_q=one_minus(cq, c2))


def instanton_M0(model: ThinWallModel) -> dict:
    """Compact O(4) instanton (M = 0) in closed form.

    p = -A_out(R0) and q = A_in(R0) are the signed cosines of the polar angles
    of the two 4-sphere caps (eta = 0: the exterior cap is a flat 4-ball), and
    p + q = kappa R0 > 0, so at most one of them is negative.

      S_E      = -(pi R0^2/2G) [1/(1+p) + 1/(1+q)]
      B_an     = (pi R0^2/2G) (p + q + 2 p_+ q_+) / ((1+|p|)(1+|q|))     (= S_E + Morse)
      B_down   = (pi R0^2/2G) [1/(1-q) - 1/(1+p)]   (CdL decay of the interior vacuum)
      B_up     = (pi R0^2/2G) [1/(1-p) - 1/(1+q)]   (CdL decay of the exterior vacuum,
                                                     = upward transition from eta, eta > 0)
    all evaluated in cancellation-free form.
    """
    if model.M != 0:
        raise ValueError("M must be 0")
    G = model.G
    R0 = 2 * model.kappa / model.L2
    cc = _stable_cap_cosines(model)
    p, q = cc["p"], cc["q"]
    pre = np.pi * R0**2 / (2 * G)
    S_E = -pre * (1 / cc["one_plus_p"] + 1 / cc["one_plus_q"])
    pp, qp = max(p, 0.0), max(q, 0.0)
    absP = cc["one_plus_p"] if p >= 0 else cc["one_minus_p"]
    absQ = cc["one_plus_q"] if q >= 0 else cc["one_minus_q"]
    B_closed = pre * (p + q + 2 * pp * qp) / (absP * absQ)
    B_down = pre * (1 / cc["one_minus_q"] - 1 / cc["one_plus_p"])
    B_up = pre * (1 / cc["one_minus_p"] - 1 / cc["one_plus_q"]) if model.eta > 0 else np.inf
    A_w = 4 * np.pi * R0**2
    maxima, minima = [], []
    if q < 0:
        maxima.append(4 * np.pi / model.chi**2)
    if p < 0:
        maxima.append(4 * np.pi / model.eta**2)
    if q >= 0 and p >= 0:
        maxima.append(A_w)
    morse = (sum(maxima) - sum(minima)) / (4 * G)
    return dict(R0=R0, p=p, q=q, S_E=S_E, B_closed=B_closed, morse=morse,
                B_CdL_down=B_down, B_up=B_up, maxima=maxima, minima=minima)


# ------------------------------------------------- embedding / multiplicity
def boundary_loop(model: ThinWallModel, side: str = "out", n_wall: int = 800):
    """Closed boundary loop of the half-bounce Euclidean region on one side,
    in a plane chart of the (r, t_E) section centred on the Euclidean horizon.

    Exterior (Schwarzschild): rho = sqrt(r/r_h - 1), angle = t_E / 4GM.
    Interior (de Sitter):     rho = sqrt(1 - chi r) (rho = 1 is the axis r = 0),
                              angle = chi t_E.
    In both charts the region lies to the RIGHT of the wall's direction of
    motion (exterior: n is T rotated clockwise; interior: the chart reverses
    the radial direction and the region lies opposite to n).  Slices are
    radial; a slice whose side sign is -1 passes through the horizon point,
    and orientation fixes the lift: final outer ray at theta_D + pi, initial
    outer ray at -pi.  The loop is closed at infinity (exterior) or along the
    axis (interior) from the final to the initial angle in the lifted angle.
    """
    ew = euclidean_wall(model, n=n_wall)
    if side == "out":
        if not (model.M > 0 and model.eta == 0):
            raise NotImplementedError("exterior map implemented for Schwarzschild")
        rh = 2 * model.G * model.M
        rho = lambda r: np.sqrt(np.maximum(r / rh - 1.0, 0.0))  # noqa: E731
        th = ew.theta_out
        A = model.A_out
        far = float(rho(rh * 400.0))
    else:
        rho = lambda r: np.sqrt(np.maximum(1.0 - model.chi * r, 0.0))  # noqa: E731
        th = ew.theta_in
        A = model.A_in
        far = 1.0
    e1, e2 = int(np.sign(A(ew.R[0]))), int(np.sign(A(ew.R[-1])))
    pts = []

    def seg(r0, r1, ang, k=150):
        for x in np.linspace(r0, r1, k):
            pts.append((x, ang))

    r1, r2 = float(rho(ew.R[0])), float(rho(ew.R[-1]))
    # initial slice, traversed from the far boundary to the wall
    if e1 > 0:
        phi_i = 0.0
        seg(far, r1, 0.0)
    else:
        phi_i = -np.pi
        seg(far, 0.0, -np.pi)
        seg(0.0, r1, 0.0)
    i_wall0 = len(pts)
    for R, a in zip(ew.R, th):
        pts.append((float(rho(R)), float(a)))
    i_wall1 = len(pts)
    thD = float(th[-1])
    if e2 > 0:
        phi_f = thD
        seg(r2, far, thD)
    else:
        phi_f = thD + np.pi
        seg(r2, 0.0, thD)
        seg(0.0, far, thD + np.pi)
    for a in np.linspace(phi_f, phi_i, 600):
        pts.append((far, a))
    L = np.array(pts)
    return dict(x=L[:, 0] * np.cos(L[:, 1]), y=L[:, 0] * np.sin(L[:, 1]), far=far,
                phi_i=phi_i, phi_f=phi_f, eps=(e1, e2), theta=th,
                rho_wall_max=float(L[i_wall0:i_wall1, 0].max()))


def multiplicity_map(model: ThinWallModel, side: str = "out", ngrid: int = 161, extent: float = None):
    """Covering multiplicity of the half-bounce region (degree theory):
    m(x) = -winding(loop, x) since the loop has the region on its right."""
    lp = boundary_loop(model, side)
    xs, ys = lp["x"], lp["y"]
    if extent is None:
        if side == "out":
            # the map must contain the whole wall: beyond the wall's largest radius
            # the multiplicity can change only across the two radial rays
            extent = min(lp["far"], 1.25 * lp["rho_wall_max"] + 0.2)
        else:
            extent = 1.0
    ext = extent
    g = np.linspace(-ext, ext, ngrid)
    GX, GY = np.meshgrid(g, g)
    wind = np.zeros_like(GX)
    for i in range(len(xs) - 1):
        a1 = np.arctan2(ys[i] - GY, xs[i] - GX)
        a2 = np.arctan2(ys[i + 1] - GY, xs[i + 1] - GX)
        wind += np.angle(np.exp(1j * (a2 - a1)))
    m = -np.rint(wind / (2 * np.pi))
    if side == "in":
        m = np.where(GX**2 + GY**2 <= 1.0, m, 0.0)
    return GX, GY, m, lp


def embedding_report(model: ThinWallModel, ngrid: int = 121) -> dict:
    """Embeddability of the half-bounce Euclidean regions (R4).

    Analytic criterion (theorem in NO_GO_THEOREMS.md): a sign change of eps
    along the barrier forbids an embedding.  Independently we compute the
    covering multiplicity by degree theory; the region is a manifold region
    iff the multiplicity takes values in {0, 1}."""
    ew = euclidean_wall(model, n=600)
    rep = dict(eps_out=(int(np.sign(model.A_out(ew.R[0]))), int(np.sign(model.A_out(ew.R[-1])))),
               eps_in=(int(np.sign(model.A_in(ew.R[0]))), int(np.sign(model.A_in(ew.R[-1])))))
    sides = []
    if model.M > 0 and model.eta == 0:
        sides.append(("exterior", "out", ew.theta_out))
    if model.chi > 0:
        sides.append(("interior", "in", ew.theta_in))
    for name, sd, th in sides:
        Av = model.A_out(ew.R) if sd == "out" else model.A_in(ew.R)
        flip = bool(np.any(np.sign(Av[1:]) != np.sign(Av[:-1])))
        _, _, m, lp = multiplicity_map(model, sd, ngrid=ngrid)
        mmin, mmax = float(m.min()), float(m.max())
        rep.update({f"{name}_flip": flip, f"theta_{sd}_D": float(th[-1]),
                    f"theta_{sd}_max": float(th.max()), f"theta_{sd}_min": float(th.min()),
                    f"{name}_multiplicity_range": (mmin, mmax),
                    f"{name}_embedded": (mmin >= 0 and mmax <= 1)})
    # backwards-compatible keys
    if "theta_out_D" in rep:
        rep["theta_D"], rep["theta_max"] = rep["theta_out_D"], rep["theta_out_max"]
    return rep


def winding_multiplicity(model: ThinWallModel, ngrid: int = 241, rho_max: float = 3.0):
    """Exterior multiplicity map (kept for the figure scripts)."""
    GX, GY, m, lp = multiplicity_map(model, "out", ngrid=ngrid, extent=rho_max)
    return GX, GY, m, (lp["x"], lp["y"])
