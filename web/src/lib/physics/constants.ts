/**
 * Constants used by the live in-browser physics: a mirror of ``qchas/constants.py``.
 *
 * Numbers quoted from the papers carry their section in a `Quoted` record so panels can show
 * provenance next to them. Keep this file in sync with the Python module.
 */

export interface Quoted<T = number> {
  value: T;
  source: string;
}

const q = <T>(value: T, source: string): Quoted<T> => ({ value, source });

// --- SI / astronomical units --------------------------------------------------------
export const C_LIGHT = 299_792_458.0;
export const G_NEWTON = 6.6743e-11;
export const HBAR = 1.054_571_817e-34;
export const K_BOLTZMANN = 1.380_649e-23;
export const M_SUN = 1.988_47e30;
export const MPC = 3.085_677_581_491_367e22;
export const GPC = 1e3 * MPC;
export const YEAR = 365.25 * 86400.0;
export const L_PLANCK = Math.sqrt((HBAR * G_NEWTON) / C_LIGHT ** 3);
export const M_PLANCK = Math.sqrt((HBAR * C_LIGHT) / G_NEWTON);

// --- Cosmology, [TOP] Sec. 14 -------------------------------------------------------
export const PLANCK2018 = { H0: 67.32, ombh2: 0.022383, omch2: 0.12011, tau: 0.0543, As: 2.101e-9, ns: 0.96605 };
export const CHI_STAR_MPC = q(13_870.0, '[TOP] Sec. 7');
export const D_LSS_GPC = q((2 * 13_870.0) / 1e3, '[TOP] Sec. 7');
export const ETA0_GPC = q(14.15, '[TOP] Prop. 10.1');

// --- Topology, [TOP] -------------------------------------------------------------------
export const CIRCLE_SEARCH_BOUND = q(0.985, '[TOP] Sec. 9 (COMPACT 2023 circle search)');
export const LIKELIHOOD_BOUNDS = q({ FCC: 1.04, BCC: 1.0, cubic: 0.985 }, '[TOP] Sec. 9');
export const Z2_CERTIFIED = q(
  { SC: 16.53231595976167, BCC: 15.971880245588027, FCC: 15.962131481337693 },
  '[TOP] Prop. E.2',
);
export const HESSIAN_CERTIFIED = q(
  {
    FCC: [0.5661, 0.5661, 3.1697, 3.1697, 3.1697],
    BCC: [0.2314, 0.2314, 3.3953, 3.3953, 3.3953],
    SC: [-1.8348, -1.8348, -1.8348, 8.263, 8.263],
  },
  '[TOP] Prop. E.2',
);
/** Best-fit (not significant) FCC lattice axes at 1.05 D, Galactic (l, b) in degrees. */
export const FCC_PREREGISTERED_AXES_GAL = q<[number, number][]>(
  [
    [280.5, 62.7],
    [219.4, -14.0],
    [185.5, 36.8],
    [277.5, 2.8],
    [95.0, 57.2],
    [337.0, 16.9],
  ],
  '[TOP] Sec. 9 (best fit, not significant)',
);
export const FIG4_PAPER = q({ d_over_D: 0.8, radius_deg: 36.9, pairs: { SC: 3, BCC: 7, FCC: 6 } }, '[TOP] Fig. 4');

// --- Horizon, [BH] -------------------------------------------------------------------
export const BETA_RADIATION = q(1.855, '[BH] Theorem C');
export const N_QUANTA_COEFF = q(4.47, '[BH] Theorem C');
export const PAGE_TIME_PAPER = q(0.476, '[BH] Fig. 1');
export const PAGE_MASS_FRACTION_PAPER = q(0.806, '[BH] Fig. 1');
export const QES_DEPTH_COEFF = q(0.0399, '[BH] Sec. 9');
export const QES_LAG_OFFSET_RS = q(4.24, '[BH] Sec. 9');

// --- Crossed product, [QG] ----------------------------------------------------------
export const ELL_OVER_LP = q(Math.sqrt(8 / Math.PI), '[QG] Theorem Q');
export const BOUNCE_AMIN_MIXMASTER = q(0.256, '[QG] Theorem N');
export const MIN_AREA_LP2 = q(5.97, '[QG] Theorem U');
export const BARBERO_IMMIRZI_SU2 = q(0.2741, '[QG] Theorem U');
