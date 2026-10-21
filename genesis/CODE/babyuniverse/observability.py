"""Phase VI: what could be observed.

* Laboratory: lower bound on the exponent, B_lab >= 0.56 B(M = 0) (T5, Proposition-N) with B(0) the
  closed form T1, over a grid of vacuum and wall energy scales.  Also the WKB validity parameter.
* Characteristic masses: the regime A/C boundary M_cr(rho, sigma) and the mass whose Schwarzschild
  radius equals the false-vacuum Hubble radius, M_H = 1/(2 G chi).  Both scale as 1/(G chi) for
  sigma^(1/3) ~ rho^(1/4) << m_Pl, i.e. as m_Pl^3 / E_vac^2.
* Hawking lifetimes and the primordial-black-hole mass windows these masses fall into.

Energy scales in GeV; Planck units internally (G = hbar = c = 1, m_Pl = 1.22e19 GeV).
"""

from __future__ import annotations

import numpy as np

from .model import ThinWallModel
from . import tunnelling as T

M_PL_GEV = 1.22e19
M_PL_G = 2.176e-5        # grams
T_PL_S = 5.391e-44


def planck_params(E_vac_GeV: float, E_wall_GeV: float) -> dict:
    rho = (E_vac_GeV / M_PL_GEV) ** 4
    sigma = (E_wall_GeV / M_PL_GEV) ** 3
    chi = np.sqrt(8 * np.pi * rho / 3)
    kappa = 4 * np.pi * sigma
    return dict(rho=rho, sigma=sigma, chi=chi, kappa=kappa, k=kappa / chi)


def B_zero(E_vac_GeV: float, E_wall_GeV: float) -> float:
    """B(M -> 0) = b(kappa/chi) / chi^2 (T1 closed form via the canonical integral)."""
    p = planck_params(E_vac_GeV, E_wall_GeV)
    return float(T.B_canonical(ThinWallModel(1.0, p["k"], 0.0)) / p["chi"] ** 2)


def B_zero_flat(E_vac_GeV: float, E_wall_GeV: float) -> float:
    """Flat-space thin-wall CdL exponent 27 pi^2 sigma^4 / (2 rho^3) (small kappa/chi limit)."""
    p = planck_params(E_vac_GeV, E_wall_GeV)
    return float(27 * np.pi**2 * p["sigma"] ** 4 / (2 * p["rho"] ** 3))


def masses(E_vac_GeV: float, E_wall_GeV: float) -> dict:
    """Characteristic masses in grams: M_cr (regime boundary) and M_H = 1/(2 G chi)."""
    p = planck_params(E_vac_GeV, E_wall_GeV)
    L = np.sqrt(1 + p["k"] ** 2)                      # in units of chi
    with np.errstate(all="ignore"):
        Mcr_dimless = ThinWallModel(1.0, p["k"], 1e-9).M_critical if p["k"] > 1e-3 else np.nan
    if not np.isfinite(Mcr_dimless):
        # kappa/chi -> 0: g_min ~ 4 kappa^2/L^2, so M_cr -> 1/(2 G L) (PHASE_V Sec. 2; checked numerically)
        Mcr_dimless = 1 / (2 * L)
    Mcr = Mcr_dimless / p["chi"]
    MH = 1 / (2 * p["chi"])
    return dict(M_cr_g=float(Mcr * M_PL_G), M_H_g=float(MH * M_PL_G))


def hawking_lifetime_s(M_g: float) -> float:
    m = M_g / M_PL_G
    return float(5120 * np.pi * m**3 * T_PL_S)


# observational primordial-black-hole windows (orders of magnitude; see PHASE_VI for sources/caveats)
PBH_WINDOWS = [
    ("evaporated before BBN (no relic)", 0.0, 1e9),
    ("BBN / CMB constrained (evaporating)", 1e9, 1e13),
    ("extragalactic gamma-ray constrained", 1e13, 1e17),
    ("asteroid-mass window (largely unconstrained)", 1e17, 1e22),
    ("microlensing constrained", 1e22, 1e33),
    ("stellar mass (GW observatories)", 1e33, 1e35),
]


def window_of(M_g: float) -> str:
    for name, lo, hi in PBH_WINDOWS:
        if lo <= M_g < hi:
            return name
    return "supermassive / beyond"


def laboratory_case(E_vac_GeV: float, E_wall_GeV: float, M_planck: float) -> dict:
    """One laboratory scenario: false vacuum rho = E_vac^4, wall sigma = E_wall^3, bubble mass M (Planck
    masses).  Returns the exponents (B_an, sector-blind B_sq), the CdL exponent of the false vacuum,
    the survival exponent ln(Gamma-volume) needed for the bubble not to decay before it collapses, the
    WKB parameter and validity flags."""
    from . import contours as Cc
    from .quantum_system import wkb_parameter
    p = planck_params(E_vac_GeV, E_wall_GeV)
    chi, k = p["chi"], p["k"]
    Mcr = ThinWallModel(1.0, k, 1e-9).M_critical if k > 1e-3 else 1 / (2 * np.sqrt(1 + k**2))
    mu = M_planck * chi                          # G M chi
    out = dict(E_vac=E_vac_GeV, E_wall=E_wall_GeV, M_planck=M_planck, kappa_over_chi=k, M_over_Mcr=mu / Mcr)
    if not (0 < mu < 0.98 * Mcr):
        out.update(valid=False, reason="no barrier (M >= M_cr)")
        return out
    m = ThinWallModel(1.0, k, mu)
    try:
        R1, R2 = m.turning_points()[0], m.turning_points()[-1]
        B_an = T.B_canonical(m) / chi**2
        B_sq = Cc.B_sector_blind(m) / chi**2
        B0 = T.B_canonical(ThinWallModel(1.0, k, 0.0)) / chi**2
    except Exception as e:                        # pragma: no cover
        out.update(valid=False, reason=str(e))
        return out
    R1_pl = R1 / chi
    # bubble survival: Gamma ~ rho e^{-B_CdL} per 4-volume, 4-volume ~ (4 pi/3) R1^3 x (R1/c)
    ln_GV = np.log(p["rho"] * 4 * np.pi / 3 * R1_pl**4)
    w = wkb_parameter(m, chi)
    B_min = min(B_an, B_sq)
    out.update(valid=True, B_an=float(B_an), B_sq=float(B_sq), B_best=float(B_min), B_CdL=float(B0),
               survival_margin=float(B0 - ln_GV), metastable=bool(B0 > ln_GV + 3),
               wkb_eps=float(w["eps_physical"]), wkb_ok=bool(w["eps_physical"] < 0.1),
               R1_lP=float(R1_pl), R2_lP=float(R2 / chi), two_GM_lP=float(2 * M_planck),
               log10_P_best=float(-B_min / np.log(10)))
    return out


def best_laboratory_case(E_vac_grid, ratio_grid, M_planck: float = 10.0) -> dict:
    """Minimise the exponent over vacuum scale and wall/vacuum ratio, subject to metastability of the
    false vacuum over the bubble's life, WKB validity, and a semiclassical black hole (M >= M_planck)."""
    rows = []
    for ev in E_vac_grid:
        for r in ratio_grid:
            rows.append(laboratory_case(ev, r * ev, M_planck))
    ok = [x for x in rows if x.get("valid") and x["metastable"] and x["wkb_ok"]]
    best = min(ok, key=lambda x: x["B_best"]) if ok else None
    return dict(M_planck=M_planck, n_scanned=len(rows), n_admissible=len(ok), best=best, rows=rows)
