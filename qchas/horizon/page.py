"""The Page transfer law and its companions, [BH] Theorems C, D, E and F.

Time ``t`` is measured in units of the evaporation lifetime and entropies in units of
the initial Bekenstein-Hawking entropy S_BH(M0), unless an absolute S0 is passed.
The black hole loses mass as dM/dt = -alpha/M^2, so m(t) = (1 - t)^(1/3).
"""

from __future__ import annotations

import math

import numpy as np

from ..constants import BETA_RADIATION, N_QUANTA_COEFF


def mass_fraction(t):
    t = np.clip(np.asarray(t, dtype=float), 0.0, 1.0)
    return (1.0 - t) ** (1.0 / 3.0)


def bh_entropy(t):
    """S_BH(t) / S_BH(M0) = m^2: the capacity A/4G."""
    return mass_fraction(t) ** 2


def hawking_entropy(t, beta: float = BETA_RADIATION):
    """Fine-grained radiation entropy of field theory on a fixed background (Theorems B, C).

    Each emitted mode is thermal, so S_H(t) = beta [S_BH(M0) - S_BH(M(t))]; it ends at
    beta S_BH(M0) = 1.855 S_BH(M0) for a massless scalar.
    """
    return beta * (1.0 - bh_entropy(t))


def page_curve(t, beta: float = BETA_RADIATION):
    """Theorem E(i): S(R) <= min{S_H, A/4G}, with equality up to e^{-|S_H - A/4G|}."""
    return np.minimum(hawking_entropy(t, beta), bh_entropy(t))


def page_curve_lower(t, s0: float, beta: float = BETA_RADIATION):
    """Theorem E(ii) lower bound for 2-design scrambling, in absolute nats (S0 = S_BH(M0))."""
    sh, sb = s0 * hawking_entropy(t, beta), s0 * bh_entropy(t)
    return np.maximum(0.0, np.minimum(sh, sb) - np.exp(-np.abs(sh - sb)))


def page_time_fraction(beta: float = BETA_RADIATION) -> float:
    """Page time t_P / t_evap where beta (1 - m^2) = m^2. 0.476 for a massless scalar."""
    m2 = beta / (1.0 + beta)
    return 1.0 - m2**1.5


def page_mass_fraction(beta: float = BETA_RADIATION) -> float:
    """Remaining mass fraction at the Page time (0.806 for a massless scalar)."""
    return math.sqrt(beta / (1.0 + beta))


def island_min_entropy(t, beta: float = BETA_RADIATION):
    """Theorem D: the island's minimal entropy [beta (1 - m^2) - m^2] S_BH(M0), when positive."""
    m2 = bh_entropy(t)
    return np.maximum(0.0, beta * (1.0 - m2) - m2)


def page_transfer_slope(t):
    """Theorem E(iii): after the Page time dS(R)/du = (1/4G) dA/du, i.e. d(m^2)/dt."""
    t = np.asarray(t, dtype=float)
    m = mass_fraction(t)
    with np.errstate(divide="ignore"):
        return -2.0 / 3.0 / np.maximum(m, 1e-300)


def interior_volume(t):
    """Theorem F: Christodoulou-Rovelli volume in units of M0^2 t_evap.

    dV/dv = 3 sqrt(3) pi M(v)^2 gives V = 3 sqrt(3) pi (3/5)(1 - (1 - t)^(5/3)), reaching
    9 sqrt(3) pi / 5 at the endpoint while the capacity m^2 goes to zero.
    """
    t = np.clip(np.asarray(t, dtype=float), 0.0, 1.0)
    return 3.0 * math.sqrt(3.0) * math.pi * 0.6 * (1.0 - (1.0 - t) ** (5.0 / 3.0))


def total_quanta(m0_planck: float) -> float:
    """Theorem C(ii): N ~ 4.470 M0^2 / m_P^2 quanta emitted in total (massless scalar)."""
    return N_QUANTA_COEFF * m0_planck**2


def decoding_infidelity_bound(k: int, s) -> np.ndarray:
    """[QG] Theorem R: E[1 - F_dec] <= 2^(1 + k - s) for a k-qubit diary after s emissions."""
    return np.minimum(1.0, 2.0 ** (1.0 + k - np.asarray(s, dtype=float)))


def hayden_preskill_emissions(k: int, eps: float) -> float:
    """[BH] Theorem E(iv): about k + log2(1/eps) further emissions decode a k-qudit diary."""
    return k + math.log2(1.0 / eps)


def curves(n: int = 401, beta: float = BETA_RADIATION) -> dict:
    """All [BH] Figure-1 curves on a common time grid, for export."""
    t = np.linspace(0.0, 1.0, n)
    return {
        "t": t.tolist(),
        "hawking": hawking_entropy(t, beta).tolist(),
        "bh_entropy": bh_entropy(t).tolist(),
        "page": page_curve(t, beta).tolist(),
        "island_min": island_min_entropy(t, beta).tolist(),
        "interior_volume": (interior_volume(t) / interior_volume(1.0)).tolist(),
        "page_time": page_time_fraction(beta),
        "page_mass_fraction": page_mass_fraction(beta),
        "beta": beta,
    }
