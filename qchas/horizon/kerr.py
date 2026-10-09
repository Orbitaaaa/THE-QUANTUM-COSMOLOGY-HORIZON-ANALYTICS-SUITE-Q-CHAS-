"""Kerr horizons: thermodynamics, remnants of binary mergers, ringdown and echo delays, [QG] Secs. 17-19.

Geometric units (G = c = 1) internally; helpers convert solar masses to seconds and to
Planck units. Used to put real LIGO-Virgo remnants into the horizon-capacity tools. The
remnant, QNM and echo-delay formulas reproduce the nine rows of the [QG] Sec. 19.1 table to
about 1 %; the Teukolsky spin-down table of [QG] Sec. 17 is quoted with consistency checks.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

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


# --- Hawking spin-down with exact Teukolsky greybody factors, quoted from [QG] Sec. 17 ------

#: [QG] Sec. 17 table: field, a_*, f = -M^2 dM/dt, g = -(M/a_*) dJ/dt, d ln a_*/d ln M,
#: GSL ratio dS_rad/(-dS_BH) (None where the horizon area grows), sign-law violations
SPIN_DOWN_TABLE = [
    ("scalar", 0.001, 7.38e-5, None, None, 1.78, 0), ("scalar", 0.200, 7.34e-5, 9.45e-5, -0.713, 1.84, 0),
    ("scalar", 0.400, 7.11e-5, 1.13e-4, -0.414, 1.85, 0), ("scalar", 0.555, 6.97e-5, 1.39e-4, -0.001, 1.89, 0),
    ("scalar", 0.600, 6.97e-5, 1.50e-4, 0.151, 1.91, 0), ("scalar", 0.800, 7.91e-5, 2.29e-4, 0.899, 2.25, 0),
    ("scalar", 0.900, 1.04e-4, 3.19e-4, 1.075, 3.98, 0), ("scalar", 0.990, 2.24e-4, 6.06e-4, 0.699, None, 0),
    ("photon", 0.001, 3.36e-5, None, None, 1.50, 0), ("photon", 0.200, 4.27e-5, 5.21e-4, 10.204, 1.56, 0),
    ("photon", 0.400, 7.57e-5, 6.60e-4, 6.718, 1.73, 0), ("photon", 0.600, 1.59e-4, 9.67e-4, 4.067, 2.09, 0),
    ("photon", 0.800, 4.01e-4, 1.71e-3, 2.249, 3.50, 0), ("photon", 0.900, 7.52e-4, 2.63e-3, 1.500, 27.25, 0),
    ("photon", 0.990, 2.13e-3, 5.78e-3, 0.712, None, 0),
    ("graviton", 0.001, 3.84e-6, None, None, 1.35, 0), ("graviton", 0.300, 1.49e-5, 2.23e-4, 12.938, 1.49, 0),
    ("graviton", 0.500, 6.82e-5, 6.24e-4, 7.142, 1.77, 0), ("graviton", 0.700, 3.91e-4, 2.32e-3, 3.939, 2.79, 0),
    ("graviton", 0.900, 4.10e-3, 1.54e-2, 1.760, None, 0), ("graviton", 0.990, 3.56e-2, 9.50e-2, 0.672, None, 0),
]
SPIN_DOWN_FACTS = {
    "fixed_point_scalar": 0.555,
    "page1976_f0": {"scalar": [7.38e-5, 7.44e-5], "graviton": [3.84e-6, 3.8e-6]},
    "superradiant_gain_a0999": {"scalar_l1m1": 0.0035, "photon_l1m1": 0.0439, "graviton_l2m2": 1.02},
    "born_at_0.99": {"a_0.5_mass_fraction": 0.77, "a_0.1_mass_fraction": 0.66, "fields": "photons + gravitons"},
    "source": "[QG] Sec. 17, Theorem K'",
}


def spin_down_rate(field: str, a_star: float) -> tuple[float, float]:
    """(f, g) of the [QG] Sec. 17 table, log-interpolated in a_* within the tabulated range."""
    rows = [r for r in SPIN_DOWN_TABLE if r[0] == field and r[3] is not None]
    a = [r[1] for r in rows]
    if not a[0] <= a_star <= a[-1]:
        raise ValueError(f"a_* = {a_star} outside the tabulated range {a[0]}-{a[-1]} for {field}")
    f = math.exp(float(np.interp(a_star, a, [math.log(r[2]) for r in rows])))
    g = math.exp(float(np.interp(a_star, a, [math.log(r[3]) for r in rows])))
    return f, g


def spin_down_mass_fraction(a_start: float = 0.99, a_end: float = 0.5,
                            fields=("photon", "graviton"), n: int = 2000) -> float:
    """M(a_end)/M(a_start) from d ln a_*/d ln M = g/f - 2 summed over fields (table-interpolated)."""
    grid = np.linspace(a_start, a_end, n)
    rate = []
    for a in grid:
        f = g = 0.0
        for fld in fields:
            fi, gi = spin_down_rate(fld, float(a))
            f, g = f + fi, g + gi
        rate.append(1.0 / (a * (g / f - 2.0)))  # d ln M / d a_*
    rate = np.asarray(rate)
    return float(math.exp(np.sum(0.5 * (rate[1:] + rate[:-1]) * np.diff(grid))))


def to_json() -> dict:
    """The "kerr_spin_down" part of horizon.json: quoted table plus consistency checks."""
    rows = []
    for fld, a, f, g, dlna, gsl, viol in SPIN_DOWN_TABLE:
        rows.append({"field": fld, "a_star": a, "f": f, "g": g, "dlna_dlnM": dlna, "gsl_ratio": gsl,
                     "area_grows": gsl is None, "sign_violations": viol,
                     "dlna_dlnM_from_f_g": None if g is None else g / f - 2.0})
    return {
        "fixed_point_scalar": SPIN_DOWN_FACTS["fixed_point_scalar"],
        "rows": rows,
        "facts": SPIN_DOWN_FACTS,
        "spin_down_check": {"a_start": 0.99, "a_end": 0.5, "mass_fraction": spin_down_mass_fraction(0.99, 0.5),
                            "paper": 0.77, "computed": True,
                            "note": "integrates the tabulated photon + graviton rates, interpolated in a_*"},
        "source": "[QG] Sec. 17",
    }
