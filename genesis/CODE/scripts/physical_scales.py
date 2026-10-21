"""Physical scales: GUT-like false vacuum, scale-invariant reduction.

B(chi, kappa, M, G) = b(kappa/chi, G M chi) / (G chi^2) (dimensional analysis),
so we compute b with chi = G = 1 and restore units.  Planck units:
m_Pl = 1.22e19 GeV = 2.18e-8 kg, l_Pl = 1.62e-35 m."""
import warnings

import numpy as np

from _common import dump
from babyuniverse.model import ThinWallModel
from babyuniverse import tunnelling as T

warnings.simplefilter("ignore")
mPl_GeV, mPl_kg, lPl_m = 1.22e19, 2.176e-8, 1.616e-35
rows = []
for E_vac, E_wall in ((1e16, 1e16), (1e16, 3e16), (1e16, 1e17), (1e15, 1e16), (1e18, 1e18)):
    rho = (E_vac / mPl_GeV) ** 4           # Planck units, G = 1
    sigma = (E_wall / mPl_GeV) ** 3
    chi = np.sqrt(8 * np.pi * rho / 3)
    kappa = 4 * np.pi * sigma
    k = kappa / chi
    base = ThinWallModel(1.0, k, 1e-9)
    Mcr_dimless = base.M_critical                     # in units 1/(G chi)
    Mcr_kg = Mcr_dimless / chi * mPl_kg
    B0 = T.B_canonical(ThinWallModel(1.0, k, 0.0)) / chi**2
    d0 = T.instanton_M0(ThinWallModel(1.0, k, 0.0))
    Bs = []
    for frac in (1e-3, 0.1, 0.5, 0.9):
        mm = ThinWallModel(1.0, k, frac * Mcr_dimless)
        R1, R2 = mm.turning_points()
        Bs.append(dict(M_over_Mcr=frac, M_kg=frac * Mcr_kg, B=T.B_canonical(mm) / chi**2,
                       R1_m=R1 / chi * lPl_m, R2_m=R2 / chi * lPl_m, two_GM_m=2 * mm.M / chi * lPl_m))
    rows.append(dict(vacuum_scale_GeV=E_vac, wall_scale_GeV=E_wall, chi_inv_m=lPl_m / chi,
                     kappa_over_chi=k, gamma=base.gamma, M_critical_kg=Mcr_kg,
                     B_M0=B0, B_CdL_down=d0["B_CdL_down"] / chi**2,
                     flat_space_CdL=27 * np.pi**2 * sigma**4 / (2 * rho**3),
                     dS_entropy=np.pi / chi**2, samples=Bs))
    print(f"vac {E_vac:.0e} GeV, wall {E_wall:.0e} GeV: kappa/chi={k:.2e}, M_cr={Mcr_kg:.3g} kg, "
          f"B(M->0)={B0:.4g} (CdL {d0['B_CdL_down']/chi**2:.4g}, flat {27*np.pi**2*sigma**4/(2*rho**3):.4g}),"
          f" S_dS={np.pi/chi**2:.3g}")
dump("physical_scales.json", rows)
