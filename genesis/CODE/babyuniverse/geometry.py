"""Explicit global trajectories of the wall and Penrose compactification.

For a 2D metric written in double-null form ds^2 = -Omega^2 dp dq (p increasing
to the future-right, q to the future-left) and a wall whose exterior lies on
the side of increasing p - q, the unit normal points to the right and

    pdot = (Rdot + A) / (2 dr/dp),     qdot = (Rdot - A) / (2 dr/dq),

with A = eps sqrt(f + Rdot^2) the on-shell signed root (model.A_in / A_out).
With s_+ = Rdot + A, s_- = A - Rdot (s_+ s_- = f on shell) every chart has a
second, algebraically equivalent form; we switch to whichever has the larger
denominator so that horizon crossings (where s_+ or s_- -> 0 together with one
null coordinate) are integrated without loss of precision:

Kruskal (Schwarzschild), UV = (1 - r/2GM) e^{r/2GM}, k = (4 G^2 M^2/r) e^{-r/2GM}:
    Udot = -U / (4GM s_+)  =  s_- / (2 k V)
    Vdot =  V / (4GM s_-)  = -s_+ / (2 k U)
global de Sitter, r = sin psi / (H cos T), p = T + psi, q = T - psi:
    pdot = H cos p / s_-  =  H s_+ cos^2 T / cos q
    qdot = H cos q / s_+  =  H s_- cos^2 T / cos p
Minkowski, u = t - r, v = t + r:   vdot = s_+,  udot = s_-.

Interior charts are drawn with the interior on the LEFT of the wall, exterior
charts with the exterior on the RIGHT (normal pointing right in both).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp

from .model import ThinWallModel


@dataclass
class Trajectory:
    side: str
    chart: str
    tau: np.ndarray
    R: np.ndarray
    Rdot: np.ndarray
    p: np.ndarray        # V (Kruskal), p (dS), v (Minkowski)
    q: np.ndarray        # U (Kruskal), q (dS), u (Minkowski)
    X: np.ndarray        # compactified horizontal coordinate
    Y: np.ndarray        # compactified time coordinate
    eps: np.ndarray


def _chart_for(model: ThinWallModel, side: str) -> str:
    if side == "in":
        return "dS" if model.chi > 0 else "Mink"
    if model.M > 0 and model.eta == 0:
        return "Kruskal"
    if model.M == 0 and model.eta > 0:
        return "dS"
    if model.M == 0 and model.eta == 0:
        return "Mink"
    raise NotImplementedError("explicit chart for Schwarzschild-de Sitter exterior: use itinerary")


def _A(model, side, R):
    return model.A_in(R) if side == "in" else model.A_out(R)


def _f(model, side, R):
    return model.f_in(R) if side == "in" else model.f_out(R)


def _H(model, side):
    return model.chi if side == "in" else model.eta


def _null_rates(model, side, chart, R, Rd, p, q):
    A = float(_A(model, side, R))
    sp_, sm_ = Rd + A, A - Rd
    if chart == "Kruskal":
        GM = model.G * model.M
        V, U = p, q
        if abs(sm_) >= abs(sp_):
            Vd = V / (4 * GM * sm_)
        else:
            k = 4 * GM**2 / R * np.exp(-R / (2 * GM))
            Vd = -sp_ / (2 * k * U)
        if abs(sp_) >= abs(sm_):
            Ud = -U / (4 * GM * sp_)
        else:
            k = 4 * GM**2 / R * np.exp(-R / (2 * GM))
            Ud = sm_ / (2 * k * V)
        return Vd, Ud
    if chart == "dS":
        H = _H(model, side)
        T = 0.5 * (p + q)
        if abs(sm_) >= abs(sp_):
            pd = H * np.cos(p) / sm_
        else:
            pd = H * sp_ * np.cos(T) ** 2 / np.cos(q)
        if abs(sp_) >= abs(sm_):
            qd = H * np.cos(q) / sp_
        else:
            qd = H * sm_ * np.cos(T) ** 2 / np.cos(p)
        return pd, qd
    if chart == "Mink":
        return sp_, sm_
    raise ValueError(chart)


def _initial_null(model, side, chart, R0, eps0):
    """Null coordinates of a point at radius R0 in a static region, placed on
    the time-symmetric slice (Killing time 0) by the boost isometry."""
    if chart == "Kruskal":
        GM = model.G * model.M
        X = np.sqrt((R0 / (2 * GM) - 1) * np.exp(R0 / (2 * GM)))
        return (X, -X) if eps0 > 0 else (-X, X)   # (V, U): region I or III
    if chart == "dS":
        H = _H(model, side)
        s = np.arcsin(min(1.0, H * R0))
        psi = s if eps0 > 0 else np.pi - s
        return (psi, -psi)                          # T = 0
    if chart == "Mink":
        return (R0, -R0) if eps0 > 0 else (-R0, R0)  # eps=-1: point on the 'left' copy
    raise ValueError(chart)


def _compactify(chart, p, q):
    if chart == "Kruskal" or chart == "Mink":
        P, Q = np.arctan(p), np.arctan(q)
        return 0.5 * (P - Q), 0.5 * (P + Q)          # (X, Y)
    if chart == "dS":
        return 0.5 * (p - q), 0.5 * (p + q)          # (psi, T)
    raise ValueError(chart)


def integrate_wall(model: ThinWallModel, side: str, R0: float, Rdot0: float,
                   tau_span=(-40.0, 40.0), R_max: Optional[float] = None,
                   R_min: Optional[float] = None, n_out: int = 2000) -> Trajectory:
    """Integrate the wall forward and backward from (R0, Rdot0) on one side.

    R0 must lie in a static region of that side (f(R0) > 0); the starting point
    is placed at Killing time zero.  Rdot obeys Rddot = -V'(R)/2.
    """
    chart = _chart_for(model, side)
    if _f(model, side, R0) <= 0:
        raise ValueError("reference point must be in a static region")
    eps0 = int(np.sign(_A(model, side, R0)))
    p0, q0 = _initial_null(model, side, chart, R0, eps0)
    scale = max([R0] + [h for h in model.horizons_out() if h] + ([1 / model.chi] if model.chi else []))
    R_max = R_max or 60.0 * scale
    R_min = R_min or 1e-4 * scale

    def rhs(t, y):
        R, Rd, p, q = y
        pd, qd = _null_rates(model, side, chart, R, Rd, p, q)
        return [Rd, -0.5 * float(model.dV(R)), pd, qd]

    def ev_big(t, y):
        return y[0] - R_max
    ev_big.terminal = True

    def ev_small(t, y):
        return y[0] - R_min
    ev_small.terminal = True

    events = [ev_big, ev_small]
    if chart == "Kruskal":
        def ev_null(t, y):
            return 1e8 - max(abs(y[2]), abs(y[3]))   # arctan within 1e-8 of pi/2: at null infinity
        ev_null.terminal = True
        events.append(ev_null)
    if chart == "dS":
        def ev_scri(t, y):
            return np.cos(0.5 * (y[2] + y[3])) - 1e-9
        ev_scri.terminal = True
        events.append(ev_scri)

    pieces = []
    for t_end in (tau_span[1], tau_span[0]):
        sol = solve_ivp(rhs, (0.0, t_end), [R0, Rdot0, p0, q0], method="DOP853",
                        rtol=1e-11, atol=1e-13, events=events, dense_output=True,
                        max_step=abs(t_end) / 400)
        pieces.append(sol)
    fw, bw = pieces
    # uniform resampling of the dense output (independent finite-difference checks
    # of the normalisation and of the junction condition rely on it)
    t_fw = np.linspace(0.0, fw.t[-1], n_out)
    t_bw = np.linspace(bw.t[-1], 0.0, n_out)
    Y = np.concatenate([bw.sol(t_bw)[:, :-1], fw.sol(t_fw)], axis=1)
    tau = np.concatenate([t_bw[:-1], t_fw])
    R, Rd, p, q = Y
    X, T = _compactify(chart, p, q)
    eps = np.sign(_A(model, side, R))
    return Trajectory(side, chart, tau, R, Rd, p, q, X, T, eps)


def reference_point(model: ThinWallModel, branch, side: str):
    """(R0, Rdot0) on the branch lying in a static region of `side`.
    Prefers a turning point; otherwise samples the branch."""
    f = (lambda R: _f(model, side, R))
    for R in branch.turning_points:
        if f(R) > 0:
            return R, 0.0
    Rs = np.geomspace(1e-3, 1e3, 4000) * max(model.horizons_out() + [1.0 / model.chi if model.chi else 1.0])
    ok = [R for R in Rs if f(R) > 0 and model.V(R) < 0]
    if not ok:
        raise ValueError("branch never enters a static region of this side")
    R0 = ok[len(ok) // 2]
    d = branch.legs[0][2] if branch.legs else 1
    return R0, d * float(np.sqrt(-model.V(R0)))


def branch_trajectories(model: ThinWallModel, branch, tau_span=(-60.0, 60.0)):
    out = {}
    for side in ("in", "out"):
        R0, Rd0 = reference_point(model, branch, side)
        out[side] = integrate_wall(model, side, R0, Rd0, tau_span=tau_span)
    return out


# ---------------------------------------------------------------- invariants
def kruskal_region(U, V):
    """Conventional Kruskal region of a point (U, V)."""
    if U < 0 < V:
        return "I"
    if U > 0 and V > 0:
        return "II"
    if U > 0 > V:
        return "III"
    if U < 0 and V < 0:
        return "IV"
    return "horizon"


def ds_region(p, q):
    T, psi = 0.5 * (p + q), 0.5 * (p - q)
    if psi < np.pi / 2 - abs(T):
        return "N"
    if psi > np.pi / 2 + abs(T):
        return "S"
    return "F" if T > 0 else "P"


def regions_along(traj: Trajectory):
    labels = []
    for p, q in zip(traj.p, traj.q):
        if traj.chart == "Kruskal":
            labels.append(kruskal_region(q, p))
        elif traj.chart == "dS":
            labels.append(ds_region(p, q))
        else:
            labels.append("Mink")
    comp = [labels[0]]
    for lab in labels[1:]:
        if lab != comp[-1] and lab != "horizon":
            comp.append(lab)
    return [c for c in comp if c != "horizon"]
