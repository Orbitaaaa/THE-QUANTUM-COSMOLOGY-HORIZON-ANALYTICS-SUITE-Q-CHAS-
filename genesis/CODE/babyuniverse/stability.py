"""Aspherical stability of a collapsing false-vacuum bubble wall (Phase IV), probe approximation.

A wall of tension sigma bounding a false-vacuum interior (energy density rho_f) in a flat
ambient spacetime is a timelike hypersurface of constant mean curvature.  Its spherical
motion R(tau) (proper time tau) conserves

    E = 4 pi sigma R^2 gamma + (4 pi/3) rho_f R^3,   gamma = sqrt(1 + Rdot^2),

so a bubble released from rest at R_0 collapses.  Normal deformations xi Y_lm obey the Jacobi
equation of the worldsheet (the linearisation of the constant-mean-curvature condition),

    box_h xi + (K_ab K^ab + Ric(n,n)) xi = 0,
    xi'' + 2 (Rdot/R) xi' + [ l(l+1)/R^2 - K_ab K^ab - Ric(n,n) ] xi = 0,

with K^theta_theta = gamma/R, K^tau_tau = Rddot/gamma and Ric = 0 in flat space.  The l = 1
translation mode xi = gamma and the l = 0 time-translation mode xi = Rdot are exact solutions;
they are used as checks.  The shape distortion is delta_l = xi_l / R.

Scope: the probe approximation neglects the wall's own gravity and the gravity of the
interior, so it is quantitative only while G M / R << 1.  It is used to estimate how much
asphericity grows before the self-gravitating bubble reaches its gravitational radius.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp


def cmc_trajectory(sigma: float, rho_f: float, R0: float, R_end: float, n: int = 4000):
    """Spherical collapse from rest at R0 down to R_end in flat space.
    Returns tau, R, Rdot, Rddot, gamma."""
    E = 4 * np.pi * sigma * R0**2 + 4 * np.pi / 3 * rho_f * R0**3

    def gam(R):
        return (E - 4 * np.pi / 3 * rho_f * R**3) / (4 * np.pi * sigma * R**2)

    # substitution R = R0 - s^2 removes the turning-point singularity of dtau = dR / Rdot
    sg = np.linspace(0.0, np.sqrt(R0 - R_end), n)
    Rs = R0 - sg**2
    g = gam(Rs)
    Rdot = -np.sqrt(np.clip(g**2 - 1, 0, None))
    # Rddot = (1/2) d(Rdot^2)/dR = gamma dgamma/dR
    num = E - 4 * np.pi / 3 * rho_f * Rs**3
    dg = (-4 * np.pi * rho_f * Rs**2) / (4 * np.pi * sigma * Rs**2) - 2 * num / (4 * np.pi * sigma * Rs**3)
    Rddot = g * dg
    with np.errstate(divide="ignore", invalid="ignore"):
        integrand = 2 * sg / np.abs(Rdot)
    integrand[0] = 2 / np.sqrt(2 * abs(Rddot[0]))      # limit s -> 0: Rdot^2 ~ 2 |Rddot_0| s^2
    tau = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(sg))])
    return dict(tau=tau, R=Rs, Rdot=Rdot, Rddot=Rddot, gamma=g, E=E)


def jacobi_evolve(traj: dict, l: int, ric_nn: float = 0.0, xi0=(1.0, 0.0)):
    """Integrate the Jacobi equation for mode l along the trajectory (interpolated in tau)."""
    tau, R, Rd, Rdd, g = traj["tau"], traj["R"], traj["Rdot"], traj["Rddot"], traj["gamma"]

    def coef(t):
        Rt = np.interp(t, tau, R); Rdt = np.interp(t, tau, Rd); Rddt = np.interp(t, tau, Rdd)
        gt = np.sqrt(1 + Rdt**2)
        K2 = (Rddt / gt) ** 2 + 2 * gt**2 / Rt**2
        return Rt, Rdt, l * (l + 1) / Rt**2 - K2 - ric_nn

    def rhs(t, y):
        Rt, Rdt, c = coef(t)
        return [y[1], -2 * Rdt / Rt * y[1] - c * y[0]]

    sol = solve_ivp(rhs, (tau[0], tau[-1]), list(xi0), t_eval=tau, method="DOP853", rtol=1e-10, atol=1e-12)
    return sol.t, sol.y[0], sol.y[1]


def zero_mode_check(traj: dict) -> dict:
    """The l = 1 translation mode xi = gamma and the l = 0 mode xi = Rdot must solve the
    Jacobi equation; returns the relative deviation of the numerical solutions started
    with the zero-mode initial data."""
    tau, g, Rd, Rdd = traj["tau"], traj["gamma"], traj["Rdot"], traj["Rddot"]
    g_dot = Rd * Rdd / g
    t1, x1, _ = jacobi_evolve(traj, 1, xi0=(g[0], g_dot[0]))
    i_end = len(t1) * 9 // 10
    dev1 = np.max(np.abs(x1[:i_end] - g[:i_end]) / g[:i_end])
    # l = 0: xi = Rdot, xi' = Rddot
    t0, x0, _ = jacobi_evolve(traj, 0, xi0=(Rd[0], Rdd[0]))
    dev0 = np.max(np.abs(x0[:i_end] - Rd[:i_end])) / np.max(np.abs(Rd[:i_end]))
    return dict(l1_translation=float(dev1), l0_time_translation=float(dev0))


def amplification(traj: dict, ls=(2, 3, 4, 6, 10), R_stop: float = None) -> dict:
    """Growth of the shape distortion delta_l = xi_l / R between release (R0, at rest) and
    R_stop: the largest singular value of the map (delta, delta') -> delta(R_stop), using the
    two fundamental solutions, with delta' measured in units of 1/R0."""
    tau, R = traj["tau"], traj["R"]
    k = len(R) - 1 if R_stop is None else int(np.searchsorted(-R, -R_stop))
    out = {}
    for l in ls:
        _, xa, _ = jacobi_evolve(traj, l, xi0=(R[0], 0.0))     # delta = 1, delta' = 0
        _, xb, _ = jacobi_evolve(traj, l, xi0=(0.0, 1.0))      # delta = 0, delta' = 1/R0
        out[int(l)] = float(np.hypot(xa[k] / R[k], xb[k] / R[k]))
    return out
