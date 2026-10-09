"""Planck 2018 best-fit LambdaCDM from CAMB, shared by the transfer, covariance and stats modules."""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from ..constants import PLANCK2018


@lru_cache(maxsize=4)
def camb_results(lmax: int = 200, accuracy_boost: float = 2.0, transfer: bool = False):
    """CAMB results for the [TOP] Sec. 14 parameters (boosted accuracy as in the paper)."""
    import camb

    pars = camb.set_params(
        H0=PLANCK2018["H0"], ombh2=PLANCK2018["ombh2"], omch2=PLANCK2018["omch2"],
        tau=PLANCK2018["tau"], As=PLANCK2018["As"], ns=PLANCK2018["ns"], mnu=PLANCK2018["mnu"],
        lmax=lmax, WantTransfer=transfer,
    )
    pars.set_accuracy(AccuracyBoost=accuracy_boost, lSampleBoost=50.0, lAccuracyBoost=accuracy_boost)
    return pars, camb.get_results(pars)


def lcdm_cl(lmax: int = 64, lensed: bool = False) -> np.ndarray:
    """TT C_l in uK^2 for l = 0..lmax (raw C_l, not D_l)."""
    _, res = camb_results(max(lmax, 200))
    spec = res.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True, lmax=lmax)
    key = "total" if lensed else "unlensed_scalar"
    return np.asarray(spec[key][: lmax + 1, 0])


def lcdm_dl(lmax: int = 64, lensed: bool = False) -> np.ndarray:
    ell = np.arange(lmax + 1)
    return ell * (ell + 1) * lcdm_cl(lmax, lensed) / (2 * np.pi)


def chi_star_mpc() -> float:
    """Comoving distance to the last-scattering surface z_* (paper: 13 870 Mpc)."""
    _, res = camb_results()
    zstar = res.get_derived_params()["zstar"]
    return float(res.comoving_radial_distance(zstar))


def eta0_mpc() -> float:
    """Conformal time today (particle horizon), Mpc."""
    _, res = camb_results()
    return float(res.conformal_time(0.0))


def primordial_power(k_mpc) -> np.ndarray:
    """Dimensionless curvature power P_R(k) = A_s (k / 0.05 Mpc^-1)^(n_s - 1)."""
    k = np.asarray(k_mpc, dtype=float)
    return PLANCK2018["As"] * (k / 0.05) ** (PLANCK2018["ns"] - 1.0)
