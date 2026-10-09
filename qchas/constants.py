"""Physical constants and the fixed numbers Q-CHAS takes from its three source papers.

Paper keys used throughout the code base:

* ``[TOP]``  "Selecting the shape of flat space" (cosmic topology research note, 27 Sep 2026)
* ``[QG]``   "Quantum gravity from the crossed product" (research note v5, 29 Sep 2026)
* ``[BH]``   "What field theory forbids, and what the horizon must do" (28 Sep 2026)

Every number below that is *quoted* rather than computed carries the section it comes
from, so the dashboard can show provenance next to it.
"""

from __future__ import annotations

import math

# --- SI / astronomical units -------------------------------------------------------
C_LIGHT = 299_792_458.0  # m s^-1
G_NEWTON = 6.674_30e-11  # m^3 kg^-1 s^-2 (CODATA 2018)
HBAR = 1.054_571_817e-34  # J s
K_BOLTZMANN = 1.380_649e-23  # J K^-1
M_SUN = 1.988_47e30  # kg
MPC = 3.085_677_581_491_367e22  # m
GPC = 1e3 * MPC
YEAR = 365.25 * 86400.0  # s

L_PLANCK = math.sqrt(HBAR * G_NEWTON / C_LIGHT**3)  # m
M_PLANCK = math.sqrt(HBAR * C_LIGHT / G_NEWTON)  # kg
T_PLANCK = L_PLANCK / C_LIGHT  # s

# --- Cosmology: Planck 2018 best fit as used in [TOP] Sec. 14 -------------------------
PLANCK2018 = dict(
    H0=67.32,
    ombh2=0.022383,
    omch2=0.12011,
    tau=0.0543,
    As=2.101e-9,
    ns=0.96605,
    mnu=0.06,
)
T_CMB_K = 2.7255
T_CMB_UK = T_CMB_K * 1e6

#: comoving radius of the last-scattering sphere, [TOP] Sec. 7 (CAMB, Planck 2018)
CHI_STAR_MPC = 13_870.0
#: last-scattering diameter D = 2 chi_*, [TOP] Sec. 7
D_LSS_GPC = 2 * CHI_STAR_MPC / 1e3  # 27.74
#: particle horizon eta_0, [TOP] Prop. 10.1
ETA0_GPC = 14.15

# --- Standard Model content entering the Casimir energy, [TOP] Sec. 6 ---------------
N_BOSONIC_DOF = 30  # 24 gauge + 4 Higgs + 2 graviton
N_FERMIONIC_DOF = 90
JARLSKOG_J = 3.12e-5  # PDG 2024, [TOP] P1

# --- Topology results quoted from [TOP] ----------------------------------------------
#: circle-search lower bound on the shortest loop (COMPACT 2023), in units of D
CIRCLE_SEARCH_BOUND = 0.985
#: 95 % exclusion of [TOP] Sec. 9 (Planck 2018 likelihood), in units of D
LIKELIHOOD_BOUNDS = {"FCC": 1.04, "BCC": 1.00, "cubic": 0.985}
#: certified Epstein zeta values at unit covolume, [TOP] Prop. E.2
Z2_CERTIFIED = {"SC": 16.532315959761670, "BCC": 15.971880245588027, "FCC": 15.962131481337693}
#: certified Hessian eigenvalues in the 5-dimensional shape space, [TOP] Prop. E.2, as printed.
#: The printed BCC value 3.3953 violates the exact identity tr H = 2 Z / 3; our analytic and
#: finite-difference Hessians agree on 3.39504 (see qchas.topology.lattice). Keep the quote as printed.
HESSIAN_CERTIFIED = {
    "FCC": [0.5661, 0.5661, 3.1697, 3.1697, 3.1697],
    "BCC": [0.2314, 0.2314, 3.3953, 3.3953, 3.3953],
    "SC": [-1.8348, -1.8348, -1.8348, 8.2630, 8.2630],
}
#: best-fit (not significant) FCC lattice-axis directions at 1.05 D, Galactic (l, b), [TOP] Sec. 9
FCC_PREREGISTERED_AXES_GAL = [
    (280.5, 62.7), (219.4, -14.0), (185.5, 36.8), (277.5, 2.8), (95.0, 57.2), (337.0, 16.9),
]

# --- Horizon results quoted from [BH] ------------------------------------------------
#: Hawking radiation entropy / S_BH(M0) for a massless scalar, [BH] Theorem C
BETA_RADIATION = 1.855
#: total number of emitted quanta N ~ 4.470 M0^2 / m_P^2, [BH] Theorem C
N_QUANTA_COEFF = 4.470
#: QES depth law and scrambling-lag offset, [BH] Sec. 9
QES_DEPTH_COEFF = 0.0399
QES_LAG_OFFSET_RS = 4.24
#: universal atmosphere profile f(g) in cutoff lengths (continuum limit), [BH] Sec. 15
ATMOSPHERE_CONTINUUM = {"area_coeff": 0.2954, 1: 0.1483, 2: 0.0305, 3: 0.0133}

# --- Crossed-product cosmology quoted from [QG] -------------------------------------
#: bounce scale <a>_min in units of ell (ell^2 = 8G/pi, ell ~ 1.60 l_P), [QG] Theorem Q
ELL_OVER_LP = math.sqrt(8.0 / math.pi)  # 1.596
BOUNCE_AMIN_ISOTROPIC = {2: 0.164, 4: 0.162, 8: 0.159}
BOUNCE_AMIN_MIXMASTER = 0.256
#: minimal non-zero area from the capacity A/4G with SU(2) counting, [QG] Theorem U
MIN_AREA_LP2 = 5.97
BARBERO_IMMIRZI_SU2 = 0.2741
