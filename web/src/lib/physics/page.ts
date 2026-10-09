/**
 * The Page transfer law and its companions, [BH] Theorems C, D, E and F: a port of
 * ``qchas/horizon/page.py`` (same formulas, same defaults).
 *
 * Time ``t`` is in units of the evaporation lifetime and entropies in units of the initial
 * Bekenstein-Hawking entropy S_BH(M0) unless an absolute S0 is passed. dM/dt = -alpha/M^2,
 * so the remaining mass fraction is m(t) = (1 - t)^(1/3).
 */

import { BETA_RADIATION, N_QUANTA_COEFF } from './constants';

const BETA = BETA_RADIATION.value;
const clip01 = (t: number): number => Math.min(1, Math.max(0, t));

/** Remaining mass fraction m(t) = (1 - t)^(1/3). */
export const massFraction = (t: number): number => Math.cbrt(1 - clip01(t));

/** S_BH(t) / S_BH(M0) = m^2: the capacity A/4G. */
export const bhEntropy = (t: number): number => massFraction(t) ** 2;

/** Field theory on a fixed background (Theorems B, C): S_H = beta (1 - m^2), ends at beta. */
export const hawkingEntropy = (t: number, beta = BETA): number => beta * (1 - bhEntropy(t));

/** Theorem E(i): S(R) <= min{S_H, A/4G}. */
export const pageCurve = (t: number, beta = BETA): number => Math.min(hawkingEntropy(t, beta), bhEntropy(t));

/** Theorem E(ii) lower bound for 2-design scrambling, in absolute nats (s0 = S_BH(M0)). */
export function pageCurveLower(t: number, s0: number, beta = BETA): number {
  const sh = s0 * hawkingEntropy(t, beta);
  const sb = s0 * bhEntropy(t);
  return Math.max(0, Math.min(sh, sb) - Math.exp(-Math.abs(sh - sb)));
}

/** Page time t_P / t_evap where beta (1 - m^2) = m^2 (0.476 for a massless scalar). */
export function pageTimeFraction(beta = BETA): number {
  const m2 = beta / (1 + beta);
  return 1 - m2 ** 1.5;
}

/** Remaining mass fraction at the Page time (0.806 for a massless scalar). */
export const pageMassFraction = (beta = BETA): number => Math.sqrt(beta / (1 + beta));

/** Theorem D: the island's minimal entropy [beta (1 - m^2) - m^2] S_BH(M0), when positive. */
export function islandMinEntropy(t: number, beta = BETA): number {
  const m2 = bhEntropy(t);
  return Math.max(0, beta * (1 - m2) - m2);
}

/** Theorem E(iii): after the Page time dS(R)/dt = d(m^2)/dt = -(2/3)/m. */
export const pageTransferSlope = (t: number): number => -2 / 3 / Math.max(massFraction(t), 1e-300);

/**
 * Theorem F: Christodoulou-Rovelli interior volume in units of M0^2 t_evap,
 * V = 3 sqrt(3) pi (3/5)(1 - (1 - t)^(5/3)), reaching 9 sqrt(3) pi / 5 at the endpoint.
 */
export const interiorVolume = (t: number): number =>
  3 * Math.sqrt(3) * Math.PI * 0.6 * (1 - (1 - clip01(t)) ** (5 / 3));

/** Theorem C(ii): N ~ 4.470 M0^2 / m_P^2 quanta in total (massless scalar). */
export const totalQuanta = (m0Planck: number): number => N_QUANTA_COEFF.value * m0Planck ** 2;

/** [QG] Theorem R: E[1 - F_dec] <= 2^(1 + k - s) for a k-qubit diary after s emissions. */
export const decodingInfidelityBound = (k: number, s: number): number => Math.min(1, 2 ** (1 + k - s));

/** [BH] Theorem E(iv): about k + log2(1/eps) further emissions decode a k-qudit diary. */
export const haydenPreskillEmissions = (k: number, eps: number): number => k + Math.log2(1 / eps);

export interface PageCurves {
  t: number[];
  hawking: number[];
  bh_entropy: number[];
  page: number[];
  island_min: number[];
  interior_volume: number[];
  page_time: number;
  page_mass_fraction: number;
  beta: number;
}

/** All [BH] Figure-1 curves on a common grid (same keys as ``horizon.json.page``). */
export function curves(n = 401, beta = BETA): PageCurves {
  const t = Array.from({ length: n }, (_, i) => i / (n - 1));
  const vEnd = interiorVolume(1);
  return {
    t,
    hawking: t.map((x) => hawkingEntropy(x, beta)),
    bh_entropy: t.map(bhEntropy),
    page: t.map((x) => pageCurve(x, beta)),
    island_min: t.map((x) => islandMinEntropy(x, beta)),
    interior_volume: t.map((x) => interiorVolume(x) / vEnd),
    page_time: pageTimeFraction(beta),
    page_mass_fraction: pageMassFraction(beta),
    beta,
  };
}
