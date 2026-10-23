"""Round 7: radiation from an explicit quantum field model (GENESIS OBSERVABLES).

Model.  A massless scalar in the s-wave, 1+1 approximation (no angular barrier, no greybody factors), in the vacuum defined on past null
infinity.  Background: exterior Schwarzschild, interior de Sitter static patch, thin wall on the regime-A bounded branch (expanding leg from the
past horizon to R1 at tau = t = 0, then the collapse leg to the future horizon).  An outgoing ray reaching I+ at retarded time u is traced back:
it left the wall at event tau_o with interior outgoing coordinate U = t_in - r*_in(R); inside, it was the ingoing ray V = U reflected at the centre
(r*_in(0) = 0); it entered the wall at tau_i with V_w(tau_i) = U_w(tau_o); outside, it came from I- at v = p(u).  Along the wall
    du/dtau = (A_out - Rdot)/f_out,  dv/dtau = (A_out + Rdot)/f_out,  dU/dtau = (A_in - Rdot)/f_in,  dV/dtau = (A_in + Rdot)/f_in,
so p'(u) = [(dv/dtau)/(dV/dtau)]_{tau_i} [(dU/dtau)/(du/dtau)]_{tau_o} exactly.  The outgoing energy flux (Davies-Fulling) is
    <T_uu> = (1/24 pi) [ (3/2)(p''/p')^2 - p'''/p' ] = (1/24 pi) [ y'^2/2 - y'' ],   y = ln p',
which tends to the 1+1 Hawking flux kappa^2/(48 pi), kappa = 1/4GM.

Daughter branch.  After the transition the wall is in region III and the region 2GM < r < R1 of the time-symmetric slice is new exterior.  The
semiclassical theory does not fix the field's state there.  Assumption D1 (stated, not derived): the state regular on the Euclidean section
(Hartle-Hawking-like near the horizon), whose outgoing flux is kappa^2/(48 pi) from u_0 on.  Any state regular at the horizon gives the same
late-time flux; the branches differ by the collapse transient.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

from .model import ThinWallModel


def _wall_history(model: ThinWallModel, n: int = 3000, eps: float = 1e-10):
    """tau, R, t_out, t_in along the collapse leg (tau >= 0) from R1 to 2GM(1 + eps); the expanding leg is its time reverse."""
    R1 = model.turning_points()[0]
    rs = 2 * model.G * model.M
    s = np.linspace(0, 1, n)
    R = R1 - (R1 - rs * (1 + eps)) * s**2
    tau, tout, tin = np.zeros(n), np.zeros(n), np.zeros(n)
    sq = lambda r: np.sqrt(max(-float(model.V(r)), 1e-300))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i in range(1, n):
            a, b = R[i], R[i - 1]
            tau[i] = tau[i - 1] + quad(lambda r: 1 / sq(r), a, b, limit=200)[0]
            tout[i] = tout[i - 1] + quad(lambda r: float(model.A_out(r) / model.f_out(r)) / sq(r), a, b, limit=200)[0]
            tin[i] = tin[i - 1] + quad(lambda r: float(model.A_in(r) / model.f_in(r)) / sq(r), a, b, limit=200)[0]
    return tau, R, tout, tin


def _rstar_out(model, r):
    rs = 2 * model.G * model.M
    return r + rs * np.log(r / rs - 1)


def _rstar_in(model, r):
    c = model.chi
    return np.log((1 + c * r) / (1 - c * r)) / (2 * c)


def ray_map(model: ThinWallModel, n: int = 3000) -> dict:
    """p(u) and y(u) = ln p'(u) for outgoing rays leaving the wall on the collapse leg (tau_o >= 0)."""
    tau, R, tout, tin = _wall_history(model, n)
    # full worldline: expanding leg (time reverse) + collapse leg
    T = np.concatenate([-tau[:0:-1], tau])
    RR = np.concatenate([R[:0:-1], R])
    TO = np.concatenate([-tout[:0:-1], tout])
    TI = np.concatenate([-tin[:0:-1], tin])
    sgn = np.concatenate([np.ones(n - 1), -np.ones(n)])           # Rdot sign: + expanding, - collapsing
    Rdot = sgn * np.sqrt(np.maximum(-np.array([float(model.V(r)) for r in RR]), 0.0))
    Ao = np.array([float(model.A_out(r)) for r in RR]); fo = np.array([float(model.f_out(r)) for r in RR])
    Ai = np.array([float(model.A_in(r)) for r in RR]); fi = np.array([float(model.f_in(r)) for r in RR])
    u_w = TO - _rstar_out(model, RR)
    v_w = TO + _rstar_out(model, RR)
    U_w = TI - _rstar_in(model, RR)
    V_w = TI + _rstar_in(model, RR)
    dv_dt, du_dt = (Ao + Rdot) / fo, (Ao - Rdot) / fo
    dV_dt, dU_dt = (Ai + Rdot) / fi, (Ai - Rdot) / fi
    # exits on the collapse leg; entry events by inverting V_w (monotone in tau)
    k0 = n - 1
    exits = np.arange(k0, len(T) - 1)
    order = np.argsort(V_w)
    tau_of_V = CubicSpline(V_w[order], T[order])
    rows = []
    for k in exits:
        Ue = U_w[k]
        if Ue < V_w.min() or Ue > V_w.max():
            continue
        ti = float(tau_of_V(Ue))
        j = np.searchsorted(T, ti)
        j = min(max(j, 1), len(T) - 1)
        w = (ti - T[j - 1]) / (T[j] - T[j - 1])
        lerp = lambda a: (1 - w) * a[j - 1] + w * a[j]
        p = lerp(v_w)
        dp = (lerp(dv_dt) / lerp(dV_dt)) * (dU_dt[k] / du_dt[k])
        rows.append((u_w[k], p, dp, ti, T[k]))
    rows = np.array(rows)
    return dict(u=rows[:, 0], p=rows[:, 1], dp=rows[:, 2], tau_in=rows[:, 3], tau_out=rows[:, 4], u0=float(u_w[k0]),
                kappa=1 / (4 * model.G * model.M))


def flux(model: ThinWallModel, n: int = 3000) -> dict:
    """Davies-Fulling flux <T_uu>(u) on the collapse branch, its late-time Hawking value, and the transient energy relative to the
    Hawking flux, int (F - F_H) du from u_0 (the turning-point signal) to the end of the computed range."""
    rm = ray_map(model, n)
    u, y = rm["u"], np.log(rm["dp"])
    keep = np.concatenate([[True], np.diff(u) > 1e-12])
    u, y = u[keep], y[keep]
    cs = CubicSpline(u, y)
    ug = np.linspace(u[0], u[-1], 4000)
    y1, y2 = cs(ug, 1), cs(ug, 2)
    F = (0.5 * y1**2 - y2) / (24 * np.pi)
    FH = rm["kappa"] ** 2 / (48 * np.pi)
    late = slice(int(0.9 * len(ug)), len(ug))
    i90 = int(np.argmax(F >= 0.9 * FH)) if np.any(F >= 0.9 * FH) else None
    transient = float(np.trapezoid(F - FH, ug))
    return dict(u=ug, F=F, F_hawking=float(FH), late_ratio=float(np.median(F[late]) / FH), u0=rm["u0"], kappa=rm["kappa"],
                switch_on_u=None if i90 is None else float(ug[i90] - rm["u0"]), transient_energy=transient,
                transient_in_hawking_units=float(transient / (FH / rm["kappa"])), min_F_over_FH=float(F.min() / FH))


def flux_lab(model: ThinWallModel, n: int = 4000) -> dict:
    """Laboratory state: time-symmetric data at R1 (t = 0) with the field in the conformal vacuum of the interior null coordinates on
    that slice (no horizon lies on the slice, so this state is regular).  Rays reaching I+ after the turning-point signal u_0 left the
    wall on the collapse leg and, traced back, meet the t = 0 slice inside the bubble as interior outgoing rays; the relevant map is
    u -> U = U_w(tau_o), with p'(u) = (dU/dtau)/(du/dtau) at the exit.  Returns the Davies-Fulling flux on the collapse branch for
    u > u_0, its late-time Hawking value and the energy emitted in excess of Hawking."""
    tau, R, tout, tin = _wall_history(model, n)
    Rdot = -np.sqrt(np.maximum(-np.array([float(model.V(r)) for r in R]), 0.0))
    Ao = np.array([float(model.A_out(r)) for r in R]); fo = np.array([float(model.f_out(r)) for r in R])
    Ai = np.array([float(model.A_in(r)) for r in R]); fi = np.array([float(model.f_in(r)) for r in R])
    u = tout - _rstar_out(model, R)
    y = np.log(((Ai - Rdot) / fi) / ((Ao - Rdot) / fo))
    keep = np.concatenate([[True], np.diff(u) > 1e-12])
    u, y = u[keep], y[keep]
    cs = CubicSpline(u, y)
    kappa = 1 / (4 * model.G * model.M)
    ug = np.linspace(u[0], u[-1], 6000)
    F = (0.5 * cs(ug, 1) ** 2 - cs(ug, 2)) / (24 * np.pi)
    FH = kappa**2 / (48 * np.pi)
    late = slice(int(0.9 * len(ug)), len(ug))
    i90 = int(np.argmax(F >= 0.9 * FH)) if np.any(F >= 0.9 * FH) else None
    excess = float(np.trapezoid(F - FH, ug))
    return dict(u=ug - u[0], F=F, F_hawking=float(FH), kappa=float(kappa), late_ratio=float(np.median(F[late]) / FH),
                switch_on_in_4GM=None if i90 is None else float((ug[i90] - u[0]) * kappa),
                excess_energy=excess, excess_in_thermal_units=float(excess / (FH / kappa)),
                F_at_u0_over_FH=float(F[0] / FH), min_F_over_FH=float(F.min() / FH))
