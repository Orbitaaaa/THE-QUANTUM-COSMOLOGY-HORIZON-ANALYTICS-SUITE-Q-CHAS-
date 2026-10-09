"""Kerr horizons: thermodynamics, remnants of binary mergers, ringdown and echo delays.

Geometric units (G = c = 1) internally; helpers convert solar masses to seconds and to
Planck units. Used to put real LIGO-Virgo remnants into the horizon-capacity tools.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..constants import C_LIGHT, G_NEWTON, L_PLANCK, M_PLANCK, M_SUN, YEAR

MSUN_S = G_NEWTON * M_SUN / C_LIGHT**3  # 4.9255e-6 s
MSUN_M = G_NEWTON * M_SUN / C_LIGHT**2  # 1476.6 m
MSUN_PLANCK = M_SUN / M_PLANCK  # 9.14e37

#: total Hawking power of a massless scalar, P = alpha / M^2 (Planck units), [BH] Sec. 17
ALPHA_SCALAR = 7.440e-5


@dataclass(frozen=True)
class KerrHorizon:
    """Kerr black hole of mass ``m_sun`` (solar masses) and dimensionless spin ``chi``."""

    m_sun: float
    chi: float = 0.0

    @property
    def m_planck(self) -> float:
        return self.m_sun * MSUN_PLANCK

    @property
    def root(self) -> float:
        return math.sqrt(max(0.0, 1.0 - self.chi**2))

    @property
    def area_planck(self) -> float:
        """Horizon area A = 8 pi M^2 (1 + sqrt(1 - chi^2)) in Planck areas l_P^2."""
        return 8.0 * math.pi * self.m_planck**2 * (1.0 + self.root)

    @property
    def entropy(self) -> float:
        """Bekenstein-Hawking entropy A / 4G in nats: one quarter nat per Planck area."""
        return self.area_planck / 4.0

    @property
    def log10_entropy(self) -> float:
        return math.log10(self.entropy)

    @property
    def hawking_temperature_K(self) -> float:
        from ..constants import HBAR, K_BOLTZMANN

        m_kg = self.m_sun * M_SUN
        kappa = C_LIGHT**4 * self.root / (2 * G_NEWTON * m_kg * (1 + self.root))
        return HBAR * kappa / (2 * math.pi * C_LIGHT * K_BOLTZMANN)

    @property
    def omega_h(self) -> float:
        """Horizon angular velocity in units of 1/M."""
        return self.chi / (2.0 * (1.0 + self.root))

    @property
    def r_plus_km(self) -> float:
        return self.m_sun * MSUN_M * (1 + self.root) / 1e3

    def evaporation_time_years(self, alpha: float = ALPHA_SCALAR) -> float:
        """t_evap = M^3 / (3 alpha) for dM/dt = -alpha / M^2 (Schwarzschild, massless scalar)."""
        from ..constants import T_PLANCK

        return self.m_planck**3 / (3.0 * alpha) * T_PLANCK / YEAR

    def page_time_years(self, alpha: float = ALPHA_SCALAR) -> float:
        from .page import page_time_fraction

        return page_time_fraction() * self.evaporation_time_years(alpha)


def remnant_aligned(m1: float, m2: float, a1: float = 0.0, a2: float = 0.0) -> tuple[float, float]:
    """Final mass and spin of an aligned-spin binary black hole merger.

    Spin: Rezzolla et al. 2008 (ApJ 674, L29) aligned-spin fit; radiated energy:
    Barausse, Morozova & Rezzolla 2012 (ApJ 758, 63). Accurate to ~1 % in mass and
    ~0.02 in spin for the moderate spins of the events used here.
    """
    if m2 > m1:
        m1, m2, a1, a2 = m2, m1, a2, a1
    q = m2 / m1
    nu = q / (1 + q) ** 2
    at = (a1 + a2 * q**2) / (1 + q**2)
    s4, s5, t0, t2, t3 = -0.1229, 0.4537, -2.8904, -3.5171, 2.5763
    chi_f = at + s4 * at**2 * nu + s5 * at * nu**2 + t0 * at * nu + 2 * math.sqrt(3) * nu + t2 * nu**2 + t3 * nu**3
    chi_f = max(-0.998, min(0.998, chi_f))
    a_tot = (a1 + a2 * q**2) / (1 + q) ** 2
    e_isco = _isco_specific_energy(a_tot)
    p0, p1 = 0.04827, 0.01707
    mf = 1 - (1 - e_isco) * nu - 4 * nu**2 * (4 * p0 + 16 * p1 * a_tot * (a_tot + 1) + e_isco - 1)
    return (m1 + m2) * mf, abs(chi_f)


def _isco_specific_energy(chi: float) -> float:
    z1 = 1 + (1 - chi**2) ** (1 / 3) * ((1 + chi) ** (1 / 3) + (1 - chi) ** (1 / 3))
    z2 = math.sqrt(3 * chi**2 + z1**2)
    r_isco = 3 + z2 - math.copysign(math.sqrt((3 - z1) * (3 + z1 + 2 * z2)), chi)
    return math.sqrt(1 - 2 / (3 * r_isco))


#: Berti, Cardoso & Will 2006 fits (PRD 73, 064030), Table VIII: (f1, f2, f3, q1, q2, q3)
_BCW = {
    (2, 2, 0): (1.5251, -1.1568, 0.1292, 0.7000, 1.4187, -0.4990),
    (2, 2, 1): (1.3673, -1.0260, 0.1628, 0.1000, 0.5436, -0.4731),
    (3, 3, 0): (1.8956, -1.3043, 0.1818, 0.9000, 2.3430, -0.4810),
}


def qnm(m_sun: float, chi: float, mode: tuple[int, int, int] = (2, 2, 0)) -> tuple[float, float]:
    """Kerr quasi-normal frequency (Hz) and damping time (s), Berti-Cardoso-Will fit (~1 %)."""
    f1, f2, f3, q1, q2, q3 = _BCW[mode]
    momega = f1 + f2 * (1 - chi) ** f3
    quality = q1 + q2 * (1 - chi) ** q3
    m_s = m_sun * MSUN_S
    omega = momega / m_s
    return omega / (2 * math.pi), 2 * quality / omega


def echo_delay(m_sun: float, chi: float) -> float:
    """Echo delay for a reflecting surface one Planck length outside the horizon, [QG] Sec. 19.

    dt = 4 M (1 + 1/sqrt(1 - chi^2)) ln(M / l_P), in seconds.
    """
    m_s = m_sun * MSUN_S
    return 4 * m_s * (1 + 1 / math.sqrt(1 - chi**2)) * math.log(m_sun * MSUN_M / L_PLANCK)


def echo_reflectivity_bound(entropy_nats: float) -> float:
    """log10 of the echo power bound |R|^2 <= e^{-S}, [QG] Sec. 18 prediction 2."""
    return -entropy_nats / math.log(10)
