"""Round 5: from field theory to the bubble's material properties (GENESIS FIELD THEORY).

Theory F1 (the false vacuum).  One real scalar with the tilted quartic of CLASSICAL_MODEL Sec. 7 (thickwall.QuarticPotential),
    U(phi) = (lam/4)(phi^2 - v^2)^2 + rho_f (1 - s(phi/v)),   s(x) = 1/2 + 3x/4 - x^3/4,   rho_f = eps lam v^4 / 4,
minimally coupled to gravity.  Everything the thin-wall analysis used as input is derived: the false-vacuum energy rho_f, the wall
tension sigma (= int sqrt(2U) dphi), the wall width w = sqrt2/(sqrt(lam) v), the scalar mass m = sqrt(U''(v)), and the CdL exponent
(flat thin-wall estimate and the full gravitational O(4) bounce of thickwall.cdl_bounce).  -v remains a minimum for eps < 16/3.

Proposition F1 (what the frontier demands of a field theory).  In the thin-wall regime sigma ~ (2 sqrt2/3) sqrt(lam) v^3, so
    B_CdL = 27 pi^2 sigma^4 / (2 rho^3) ~ 6745 / (eps^3 lam),       lam = 0.790 (4 rho/eps)^3 / sigma^4  for given (sigma, rho).
The semiclassical frontier of Corollary S7 used B_CdL = 10 at rho^(1/4) = 2.7e17 GeV, sigma^(1/3) = 2.2e17 GeV: that needs
lam ~ 690/eps^3 (2.5e4 at eps = 0.3), where the one-loop correction ~ lam/(4 pi^2) exceeds the tree potential.  With a perturbative
coupling the frontier moves to B_lab >= 1.671 sqrt(6745/(eps^3 lam) * M R) (Planck units), e.g. e^-1100 for lam = 4 pi, eps = 1/2.

Theory F2 (the NEC-violating sector of Route G).  Galilean genesis (Creminelli-Nicolis-Trincherini 2010),
    L = f^2 e^{2 pi} (d pi)^2 + (f^3/Lambda^3)(d pi)^2 Box pi + (f^3/2 Lambda^3)(d pi)^4,
with H_0^2 ~ Lambda^3/f and h ~ f^2/(9 M_P^2 H_0^2) (reduced Planck mass; O(1) factors not fixed here).
"""

from __future__ import annotations

import warnings

import numpy as np

from .thickwall import QuarticPotential, cdl_bounce

M_PL_GEV = 1.22e19                    # Planck mass (G = 1/M_pl^2)
M_RED_GEV = M_PL_GEV / np.sqrt(8 * np.pi)
HBARC_GEV_M = 1.97327e-16


def theory(lam: float, v: float, eps: float) -> QuarticPotential:
    """Planck units: v in m_P, rho_f = eps lam v^4 / 4."""
    return QuarticPotential(lam=lam, v=v, rho_f=eps * lam * v**4 / 4)


def derived_parameters(lam: float, v: float, eps: float) -> dict:
    """Material properties of the bubble derived from (lam, v, eps), all in Planck units, plus the validity diagnostics."""
    pot = theory(lam, v, eps)
    sigma = pot.sigma0                      # thin-wall tension in the degenerate limit (the bounce agrees to 3-17%, gravitational_bounce)
    rho = pot.rho_f
    h = 1e-6 * v
    m = np.sqrt((pot.U(v + h) - 2 * pot.U(v) + pot.U(v - h)) / h**2)
    m_false = np.sqrt(max((pot.U(-v + h) - 2 * pot.U(-v) + pot.U(-v - h)) / h**2, 0.0))
    R0 = 3 * sigma / rho
    chi = np.sqrt(8 * np.pi * rho / 3)
    kappa = 4 * np.pi * sigma
    return dict(sigma=float(sigma), sigma_tilted=float(pot.sigma_thin()), rho=float(rho), width=float(pot.width), m_true=float(m), m_false=float(m_false),
                R0_flat=float(R0), width_over_R0=float(pot.width / R0), kappa_over_chi=float(kappa / chi),
                B_flat_thin=float(27 * np.pi**2 * sigma**4 / (2 * rho**3)),
                one_loop_relative=float(lam / (4 * np.pi**2)), planck_operator_relative=float(8 * np.pi * v**2),
                E_vac_GeV=float(rho**0.25 * M_PL_GEV), E_wall_GeV=float(sigma ** (1 / 3) * M_PL_GEV))


def gravitational_bounce(lam: float, v: float, eps: float) -> float:
    """Full O(4) CdL exponent (false vacuum de Sitter -> true vacuum Minkowski), thickwall.cdl_bounce, G = 1."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(cdl_bounce(theory(lam, v, eps), G=1.0)["B"])


def couplings_for(E_vac_GeV: float, E_wall_GeV: float, eps: float) -> dict:
    """The quartic theory with tilt ratio eps that reproduces a prescribed (rho, sigma) in the thin-wall regime."""
    rho = (E_vac_GeV / M_PL_GEV) ** 4
    sigma = (E_wall_GeV / M_PL_GEV) ** 3
    lam = (64 / 81) * (4 * rho / eps) ** 3 / sigma**4
    v = (4 * rho / (eps * lam)) ** 0.25
    return dict(lam=float(lam), v_planck=float(v), v_GeV=float(v * M_PL_GEV), one_loop_relative=float(lam / (4 * np.pi**2)),
                perturbative=bool(lam / (4 * np.pi**2) < 1))


def perturbative_frontier(lam_max: float, eps: float, M_planck: float = 1.0, R_lP: float = 100.0) -> dict:
    """Frontier of Corollary S7 when sigma and rho come from the quartic theory with lam <= lam_max at tilt ratio eps (thin wall):
    B_CdL >= 6745/(eps^3 lam_max), B_lab >= 1.671 sqrt(B_CdL M R)."""
    from .bounds import FRONTIER_CONST
    B_cdl = 6745 / (eps**3 * lam_max)
    B = FRONTIER_CONST * np.sqrt(B_cdl * M_planck * R_lP)
    return dict(B_CdL_min=float(B_cdl), B_lab_min=float(B), log10_P_max=float(-B / np.log(10)))


def genesis_scales(f_GeV: float, Lambda_GeV: float) -> dict:
    """Theory F2 scales (O(1) factors not fixed): H_0 = sqrt(2 Lambda^3 / 3 f), h = f^2/(9 M_P^2 H_0^2) in GeV^-2 and m^2;
    EFT consistency needs H_0 << Lambda (i.e. Lambda << f)."""
    H0 = np.sqrt(2 * Lambda_GeV**3 / (3 * f_GeV))
    h = f_GeV**2 / (9 * M_RED_GEV**2 * H0**2)
    return dict(H0_GeV=float(H0), h_GeV2=float(h), h_m2=float(h * HBARC_GEV_M**2), H0_over_Lambda=float(H0 / Lambda_GeV),
                min_seed_size_m=float(np.sqrt(h) * HBARC_GEV_M * 10))


def laboratory_exponent_from_theory(lam: float, v: float, eps: float, M_planck: float = 1.0) -> dict:
    """End-to-end chain for Theory F1: couplings -> (sigma, rho) -> E = (R_CdL/R_M)^2 -> B_lab = B_CdL Phi(E) (Theorem S3, light walls)
    with B_CdL from the gravitational bounce, R_M = (3M / 4 pi rho)^(1/3)."""
    from .bounds import Phi
    d = derived_parameters(lam, v, eps)
    B_cdl = gravitational_bounce(lam, v, eps)
    R_M = (3 * M_planck / (4 * np.pi * d["rho"])) ** (1 / 3)
    E = (d["R0_flat"] / R_M) ** 2
    B_lab = B_cdl * Phi(E)
    return dict(B_CdL=float(B_cdl), R_M_lP=float(R_M), E=float(E), Phi=float(Phi(E)), B_lab=float(B_lab),
                log10_P=float(-B_lab / np.log(10)), kappa_over_chi=d["kappa_over_chi"], width_over_R0=d["width_over_R0"])
