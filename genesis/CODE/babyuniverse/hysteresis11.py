"""Round 11 (G11-08): does cosmological hysteresis survive beyond homogeneity?

Background: the closed (S^{D-1}) dust + scalar interior of ``openprob10.scalar_from_bounce`` through repeated model-U bounces.
Perturbations: scalar-field inhomogeneities expanded in hyperspherical harmonics, -Laplacian Y_l = l (l + D - 2) Y_l / a^2, evolved as test
fields on the bouncing background,
    d^2 dphi_l/dt^2 + (D-1) H d dphi_l/dt + [l (l + D - 2)/a^2 + V''(phi)] dphi_l = 0,
with vacuum-sized initial amplitudes.  Metric perturbations are not included (first-order test-field analysis); a mode is dangerous if
its gradient + kinetic energy grows to a sizeable fraction of the homogeneous field energy before the interior escapes into inflation.
Energy accounting: the Friedmann constraint h(H^2 + 1/a^2) = rho is monitored along the solution, and the work done by the field per
cycle, W = -oint p_phi dV, is recorded (the hysteresis loop integral of Sahni & Toporensky).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from . import openprob10 as O
from . import unified6l as U


def run(D=10, VF_over_Lc=0.05, eps=1e-3, phi_b=1.0, sgn=1.0, ells=(1, 2, 4, 8, 16, 40), amp0=1e-6, t_max=3000.0, kind="hayward"):
    C_d0 = 3.0 * U.extremal_mu(D, 1.0, kind)
    Lc = U.vacuum_threshold(C_d0, D, 1.0, kind)["Lam_c"]
    VF = VF_over_Lc * Lc
    A = 16.0 * VF
    a_b = U.closed_dust(C_d0, D, 1.0, kind)["a_b"]
    rho_b = C_d0 / a_b ** (D - 1)
    K0 = eps * rho_b
    C_d = (1.0 - eps) * C_d0 - O._V(phi_b, VF, A) * a_b ** (D - 1)
    nl = len(ells)
    lam = np.array([l * (l + D - 2.0) for l in ells])

    def d2V(phi):
        s = np.clip(phi, 0.0, 1.0)
        dd = A * (12 * phi**2 - 12 * phi + 2)
        if 0.0 < phi < 1.0:
            dd += VF * (6.0 - 12.0 * s)
        return dd

    def rho(a, phi, dphi):
        return C_d / a ** (D - 1) + 0.5 * dphi**2 + O._V(phi, VF, A)

    def rhs(t, y):
        a, H, phi, dphi = y[:4]
        r_ = rho(a, phi, dphi)
        psi = float(U.h_inv(r_, 1.0, kind))
        x = max(1.0 - psi, 1e-14)
        h1 = 1.0 / x**2 if kind == "hayward" else 1.0 / x
        out = [a * H, 1.0 / a**2 - (D - 1) * (C_d / a ** (D - 1) + dphi**2) / (2.0 * h1), dphi, -(D - 1) * H * dphi - O._dV(phi, VF, A)]
        dp = y[4:4 + nl]; dq = y[4 + nl:]
        m2 = lam / a**2 + d2V(phi)
        return np.concatenate([out, dq, -(D - 1) * H * dq - m2 * dp])

    y0 = np.concatenate([[a_b, 1e-9, phi_b, sgn * np.sqrt(2 * K0)], np.full(nl, amp0), np.zeros(nl)])
    ev_big = lambda t, y: y[0] - 1e3 * C_d0 ** (1.0 / (D - 3))
    ev_big.terminal = True
    sol = solve_ivp(rhs, (0.0, t_max), y0, method="DOP853", rtol=1e-9, atol=1e-14, events=[ev_big], max_step=0.5, dense_output=False)
    a, H, phi, dphi = sol.y[:4]
    dp = sol.y[4:4 + nl]; dq = sol.y[4 + nl:]
    # constraint residual
    res = np.array([float(U.h_inv(rho(a[k], phi[k], dphi[k]), 1.0, kind)) - (H[k] ** 2 + 1.0 / a[k] ** 2) for k in range(0, len(a), max(1, len(a) // 400))])
    E_hom = 0.5 * dphi**2 + np.array([O._V(p, VF, A) for p in phi])
    E_l = 0.5 * dq**2 + 0.5 * (lam[:, None] / a**2) * dp**2
    turn = np.where((H[:-1] > 0) & (H[1:] <= 0))[0]
    amp = np.max(np.abs(dp), axis=1) / amp0
    frac = np.max(E_l / np.maximum(E_hom, 1e-300)[None, :], axis=1)
    escaped = bool(len(sol.t_events[0]) > 0)
    # work per cycle: W = -oint p dV with V ~ a^{D-1}, p = dphi^2/2 - V
    p = 0.5 * dphi**2 - np.array([O._V(x, VF, A) for x in phi])
    Vol = a ** (D - 1)
    W = []
    bounds = [0] + list(turn)
    for i0, i1 in zip(bounds[:-1], bounds[1:]):
        W.append(float(-np.trapezoid(p[i0:i1 + 1], Vol[i0:i1 + 1])))
    return dict(D=D, VF_over_Lc=VF_over_Lc, eps=eps, escaped=escaped, n_cycles=int(len(turn)), cycle_maxima=[float(a[i]) for i in turn][:20],
                ells=list(ells), max_amplification=amp.tolist(), max_energy_fraction=frac.tolist(),
                constraint_max_abs=float(np.max(np.abs(res))), work_per_cycle=W[:20], t_end=float(sol.t[-1]))
