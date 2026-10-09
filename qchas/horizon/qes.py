"""Quantum extremal surface of an evaporating 4D Schwarzschild black hole, [BH] Sec. 9.

Matter is an s-wave CFT of central charge c in the in-vacuum (vacuum in the Kruskal
coordinate U and in advanced time v). The island is I = [0, a], the radiation
R = [b, inf) collected at radius r_b and retarded time u_b. With metric
ds^2 = -e^{2 rho} dU dv, e^{2 rho} = (r_s/r) e^{-r/r_s} e^{kappa v},

    S_gen(a) = pi r_a^2 / G + (c/6) ln[(U_a - U_b)(v_b - v_a) e^{rho_a + rho_b} / eps^2].

We extremize over the island endpoint (r_a, v_a) by Newton/hybrid root finding on the
analytic gradient. Planck units (G = hbar = 1), cutoff eps = l_P, r_s = 2M.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

import numpy as np
from scipy.optimize import root
from scipy.special import lambertw

from ..constants import QES_DEPTH_COEFF, QES_LAG_OFFSET_RS


@dataclass
class QESResult:
    mass: float
    c: float
    rb_over_rs: float
    s_bh: float
    depth: float  # (r_s - r_QES)/r_s
    lag: float  # u_b - v_QES
    t_scr: float  # (1/2 pi T) ln(S_BH / c)
    s_gen_minus_s_bh: float
    converged: bool

    def as_dict(self) -> dict:
        return asdict(self)


def _r_of_uv(rs: float, U: float, v: float) -> float:
    """Areal radius from Kruskal U and advanced time v: (r/rs - 1) e^{r/rs} = -kappa U e^{kappa v}."""
    kappa = 1.0 / (2.0 * rs)
    X = -kappa * U * math.exp(kappa * v)
    return rs * (1.0 + float(np.real(lambertw(X / math.e))))


def solve_qes(mass: float, c: float = 1.0, rb_over_rs: float = 20.0, u_b: float = 0.0) -> QESResult:
    rs = 2.0 * mass
    kappa = 1.0 / (2.0 * rs)
    s_bh = math.pi * rs**2  # A/4G = 4 pi M^2
    rb = rb_over_rs * rs
    rstar_b = rb + rs * math.log(rb / rs - 1.0)
    v_b = u_b + 2.0 * rstar_b
    U_b = -math.exp(-kappa * u_b) / kappa
    rho_b = 0.5 * (math.log(rs / rb) - rb / rs + kappa * v_b)

    ln_abs_Ub = -kappa * u_b - math.log(kappa)

    def ln_U(delta: float, v: float) -> float:
        # r_a = rs (1 - delta) inside the horizon: (x - 1) e^x = -kappa U e^{kappa v}, U_a > 0
        return math.log(delta) + 1.0 - delta - kappa * v - math.log(kappa)

    def s_gen(delta: float, v: float) -> float:
        ra = rs * (1.0 - delta)
        rho_a = 0.5 * (math.log(rs / ra) - ra / rs + kappa * v)
        ln_dU = float(np.logaddexp(ln_U(delta, v), ln_abs_Ub))  # ln(U_a - U_b), U_b < 0
        return math.pi * ra**2 + c / 6.0 * (ln_dU + math.log(v_b - v) + rho_a + rho_b)

    def grad(p):
        log_delta, v = p
        delta = math.exp(min(log_delta, 0.0))
        v = min(v, v_b - 1e-9 * rs)
        # w = U_a / (U_a - U_b), computed stably in log space
        w = 1.0 / (1.0 + math.exp(min(700.0, ln_abs_Ub - ln_U(delta, v))))
        d_area = -2.0 * math.pi * rs**2 * (1.0 - delta)
        d_rho = 0.5 * (1.0 / (1.0 - delta) + 1.0)
        # delta * dS/d delta  and  r_s * dS/dv, both scaled by 6/c
        g_delta = delta * d_area * 6.0 / c + (1.0 - delta) * w + delta * d_rho
        g_v = rs * (-kappa * w - 1.0 / (v_b - v) + 0.5 * kappa)
        return [g_delta, g_v]

    # initial guess from the fit laws of [BH] Sec. 9
    t_scr = (1.0 / kappa) * math.log(s_bh / c)
    guess = [math.log(QES_DEPTH_COEFF * c / s_bh), u_b - (t_scr + QES_LAG_OFFSET_RS * rs)]
    best = None
    for dv in (0.0, -5 * rs, 5 * rs, -20 * rs, 20 * rs):
        sol = root(grad, [guess[0], guess[1] + dv], method="hybr", options={"xtol": 1e-13})
        if sol.success and (best is None or np.linalg.norm(sol.fun) < np.linalg.norm(best.fun)):
            best = sol
    assert best is not None
    delta, v_a = math.exp(best.x[0]), float(best.x[1])
    return QESResult(
        mass=mass, c=c, rb_over_rs=rb_over_rs, s_bh=s_bh, depth=delta, lag=u_b - v_a, t_scr=t_scr,
        s_gen_minus_s_bh=s_gen(delta, v_a) - s_bh, converged=bool(best.success),
    )


def depth_fit(c: float, s_bh: float) -> float:
    """[BH] Sec. 9 law: (r_s - r_QES)/r_s = 0.0399 c / S_BH."""
    return QES_DEPTH_COEFF * c / s_bh


def lag_fit(mass: float, c: float) -> float:
    """[BH] Sec. 9 law: u_b - v_QES = (1/2 pi T) ln(S_BH/c) + 4.24 r_s."""
    rs = 2.0 * mass
    return 2.0 * rs * math.log(4.0 * math.pi * mass**2 / c) + QES_LAG_OFFSET_RS * rs


#: the table published in [BH] Sec. 9, used for validation and shown in the dashboard
PAPER_TABLE = [
    # M, c, rb/rs, depth, lag, t_scr, Sgen - SBH
    (10, 1, 20, 3.11e-05, 370.16, 285.45, 0.98), (10, 1, 60, 3.23e-05, 371.49, 285.45, 1.12),
    (10, 10, 20, 3.10e-04, 277.86, 193.34, 11.56), (10, 10, 60, 3.22e-04, 279.34, 193.34, 13.03),
    (30, 1, 20, 3.47e-06, 1374.59, 1120.01, 1.17), (30, 1, 60, 3.59e-06, 1378.24, 1120.01, 1.30),
    (30, 10, 20, 3.45e-05, 1097.81, 843.70, 13.51), (30, 10, 60, 3.58e-05, 1101.83, 843.70, 14.92),
    (100, 1, 20, 3.14e-07, 5546.60, 4696.55, 1.38), (100, 1, 60, 3.23e-07, 5557.64, 4696.55, 1.51),
    (100, 10, 20, 3.12e-06, 4624.19, 3775.51, 15.64), (100, 10, 60, 3.23e-06, 4636.30, 3775.51, 16.98),
]
