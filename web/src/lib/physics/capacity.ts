/**
 * Horizon capacity from the gravitational crossed product, [BH] Theorem G, [QG] Theorems C, U:
 * a port of ``qchas/horizon/capacity.py``. Planck units throughout.
 *
 * The exterior algebra is type III_1 (no trace). Dressing it to a clock at infinity gives the
 * type II_inf crossed product with trace weight e^{beta X}; an ADM-energy window cuts it to a
 * type II_1 factor whose states have entropy at most ln tr(Pi).
 */

import { BARBERO_IMMIRZI_SU2, C_LIGHT, L_PLANCK, MPC, PLANCK2018 } from './constants';

export interface LadderRung {
  type: string;
  label: string;
  trace: string;
  entropy: string;
  source: string;
}

export const ALGEBRA_LADDER: LadderRung[] = [
  {
    type: 'III_1',
    label: 'exterior QFT algebra A',
    trace: 'none',
    entropy: 'undefined',
    source: "Buchholz-D'Antoni-Fredenhagen; [BH] Sec. 12",
  },
  {
    type: 'II_inf',
    label: 'dressed algebra A x_sigma R (clock at infinity)',
    trace: 'semifinite, weight e^{beta X}',
    entropy: 'generalized entropy A/4G + S_out',
    source: 'Takesaki; [BH] Theorem G(ii)-(iii)',
  },
  {
    type: 'II_1',
    label: 'energy window Pi (A x R) Pi',
    trace: 'finite, tr Pi = int e^{beta X} dX',
    entropy: 'S <= ln tr Pi ~ A/4G',
    source: '[BH] Theorem G(iv)-(v)',
  },
];

/** S_BH = 4 pi M^2 (Planck units). */
export const schwarzschildEntropy = (m: number): number => 4 * Math.PI * m * m;

export interface CapacityWindow {
  m_star: number;
  width: number;
  beta: number;
  /** exact ln int_0^width e^{beta X} dX */
  ln_trace: number;
  /** the [BH] Theorem G(v) expression A(E2)/4G - S_BH(M*) + ln[(1 - e^{-beta width})/beta] */
  theorem_G_v: number;
  delta_S_BH: number;
}

/**
 * Maximal entropy of the type II_1 algebra for the ADM window [M*, M* + width]:
 * ln[(e^{beta width} - 1)/beta] with beta = 8 pi M* (first law at M*).
 */
export function capacityWindow(mStar: number, width: number): CapacityWindow {
  const beta = 8 * Math.PI * mStar;
  const x = beta * width;
  const lnTrace = (x > 30 ? x + Math.log1p(-Math.exp(-x)) : Math.log(Math.expm1(x))) - Math.log(beta);
  const dS = schwarzschildEntropy(mStar + width) - schwarzschildEntropy(mStar);
  return {
    m_star: mStar,
    width,
    beta,
    ln_trace: lnTrace,
    theorem_G_v: dS + Math.log(-Math.expm1(-x) / beta),
    delta_S_BH: dS,
  };
}

export function capacityCurve(mStar: number, widths: number[]): Record<keyof CapacityWindow, number[]> {
  const rows = widths.map((w) => capacityWindow(mStar, w));
  const keys = ['m_star', 'width', 'beta', 'ln_trace', 'theorem_G_v', 'delta_S_BH'] as const;
  return Object.fromEntries(keys.map((k) => [k, rows.map((r) => r[k])])) as Record<keyof CapacityWindow, number[]>;
}

export interface DeSitterEntropy {
  Lambda_m2: number;
  entropy_nats: number;
  log10_entropy: number;
  horizon_radius_Gpc: number;
  omega_lambda: number;
}

/** [QG] Theorem C: S = 3 pi / (G Lambda) for the measured cosmological constant. */
export function deSitterEntropy(h0KmsMpc = PLANCK2018.H0, omegaLambda = 0.6847): DeSitterEntropy {
  const h0 = (h0KmsMpc * 1e3) / MPC;
  const lam = (3 * omegaLambda * h0 * h0) / C_LIGHT ** 2;
  const s = (3 * Math.PI) / (lam * L_PLANCK ** 2);
  return {
    Lambda_m2: lam,
    entropy_nats: s,
    log10_entropy: Math.log10(s),
    horizon_radius_Gpc: Math.sqrt(3 / lam) / (1e3 * MPC),
    omega_lambda: omegaLambda,
  };
}

/** [QG] Theorem U: area gap a0 = 4 sqrt(3) pi gamma l_P^2 (5.97 l_P^2 for SU(2) counting). */
export const minimalArea = (gamma = BARBERO_IMMIRZI_SU2.value): number => 4 * Math.sqrt(3) * Math.PI * gamma;

/** Capacity of a surface of area A (m^2): one quarter nat per Planck area, in bits. */
export const bitsPerArea = (areaM2: number): number => areaM2 / (4 * L_PLANCK ** 2) / Math.LN2;
