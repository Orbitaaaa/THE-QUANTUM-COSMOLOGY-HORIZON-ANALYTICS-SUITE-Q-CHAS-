"""Horizon capacity from the gravitational crossed product, [BH] Theorem G, [QG] Theorems C, U.

The exterior field algebra A of a black hole is type III_1: no trace, no entropy.
Dressing it to a clock at infinity gives the crossed product A x_sigma R, a type II_inf
factor whose trace weight is e^{beta X} = e^{S_BH(M) - S_BH(M*)}. Restricting the ADM
energy to a window [E1, E2] cuts it to a type II_1 factor with finite trace, and every
state on it has entropy at most ln tr(Pi). Everything here is in Planck units.
"""

from __future__ import annotations

import math

import numpy as np

from ..constants import BARBERO_IMMIRZI_SU2, C_LIGHT, L_PLANCK, MPC, PLANCK2018

ALGEBRA_LADDER = [
    {"type": "III_1", "label": "exterior QFT algebra A",
     "trace": "none", "entropy": "undefined", "source": "Buchholz-D'Antoni-Fredenhagen; [BH] Sec. 12"},
    {"type": "II_inf", "label": "dressed algebra A x_sigma R (clock at infinity)",
     "trace": "semifinite, weight e^{beta X}", "entropy": "generalized entropy A/4G + S_out",
     "source": "Takesaki; [BH] Theorem G(ii)-(iii)"},
    {"type": "II_1", "label": "energy window Pi (A x R) Pi",
     "trace": "finite, tr Pi = int e^{beta X} dX", "entropy": "S <= ln tr Pi ~ A/4G",
     "source": "[BH] Theorem G(iv)-(v)"},
]


def schwarzschild_entropy(m_planck):
    """S_BH = 4 pi M^2 (Planck units)."""
    return 4.0 * math.pi * np.asarray(m_planck, dtype=float) ** 2


def capacity_window(m_star: float, width: float) -> dict:
    """Maximal entropy of the type II_1 algebra for the ADM-energy window [M*, M* + width].

    Exact: S_max = ln int_0^width e^{beta X} dX with beta = 8 pi M* (first law at M*), i.e.
    ln[(e^{beta width} - 1)/beta]. [BH] Theorem G(v) evaluates the same integral as
    A(E2)/4G - S_BH(M*) + ln[(1 - e^{-beta (E2 - E1)})/beta].
    """
    beta = 8.0 * math.pi * m_star
    x = beta * width
    # log((e^x - 1)/beta) without overflow
    s_max = (x + math.log1p(-math.exp(-x)) if x > 30 else math.log(math.expm1(x))) - math.log(beta)
    s_formula = (schwarzschild_entropy(m_star + width) - schwarzschild_entropy(m_star)
                 + math.log(-math.expm1(-x) / beta))
    return {
        "m_star": m_star, "width": width, "beta": beta,
        "ln_trace": s_max,
        "theorem_G_v": float(s_formula),
        "delta_S_BH": float(schwarzschild_entropy(m_star + width) - schwarzschild_entropy(m_star)),
    }


def capacity_curve(m_star: float, widths) -> dict:
    rows = [capacity_window(m_star, float(w)) for w in widths]
    return {k: [r[k] for r in rows] for k in rows[0]}


def de_sitter_entropy(h0_kms_mpc: float = PLANCK2018["H0"], omega_lambda: float = 0.6847) -> dict:
    """[QG] Theorem C: S = 3 pi / (G Lambda) for the measured cosmological constant."""
    h0 = h0_kms_mpc * 1e3 / MPC
    lam = 3.0 * omega_lambda * h0**2 / C_LIGHT**2  # m^-2
    s = 3.0 * math.pi / (lam * L_PLANCK**2)
    return {"Lambda_m2": lam, "entropy_nats": s, "log10_entropy": math.log10(s),
            "horizon_radius_Gpc": math.sqrt(3.0 / lam) / (1e3 * MPC), "omega_lambda": omega_lambda}


def minimal_area(gamma: float = BARBERO_IMMIRZI_SU2) -> float:
    """[QG] Theorem U: area gap a0 = 4 sqrt(3) pi gamma l_P^2 (5.97 l_P^2 for SU(2) counting)."""
    return 4.0 * math.sqrt(3.0) * math.pi * gamma


def bits_per_area(area_m2: float) -> float:
    """Capacity of a surface: one quarter nat per Planck area, in bits."""
    return area_m2 / (4.0 * L_PLANCK**2) / math.log(2.0)

