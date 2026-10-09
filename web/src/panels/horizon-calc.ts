/**
 * Pure helpers for the Module 2 (horizon) panels: Kerr horizon thermodynamics in SI units,
 * the massless-scalar lifetime of [BH] Sec. 17, the [BH] Theorem F endpoint volume, the
 * decoding matrix of [QG] Theorem R and the LIGO strain alignment / ringdown projection.
 *
 * No DOM access here, so everything is unit-tested in node (tests/horizon-calc.test.ts).
 * The formulas mirror qchas/horizon/kerr.py and qchas/horizon/page.py.
 */

import { capacityWindow, minimalArea } from '../lib/physics/capacity';
import { C_LIGHT, K_BOLTZMANN, L_PLANCK, M_PLANCK, M_SUN, YEAR } from '../lib/physics/constants';
import { pageTimeFraction } from '../lib/physics/page';

/** Total Hawking power of a massless scalar, P = alpha / M^2 in Planck units ([BH] Sec. 17). */
export const ALPHA_SCALAR = 7.44e-5;
export const T_PLANCK = L_PLANCK / C_LIGHT;
/** Planck temperature m_P c^2 / k_B in kelvin. */
export const T_PLANCK_K = (M_PLANCK * C_LIGHT ** 2) / K_BOLTZMANN;
export const MSUN_PLANCK = M_SUN / M_PLANCK;

export interface BlackHoleNumbers {
  massSun: number;
  spin: number;
  massPlanck: number;
  /** outer horizon radius r_+ in km */
  rPlusKm: number;
  areaM2: number;
  /** A / l_P^2 */
  areaLp2: number;
  /** S = A / 4G in nats (one quarter nat per Planck area) */
  entropyNats: number;
  log10EntropyNats: number;
  log10EntropyBits: number;
  temperatureK: number;
  /** A / a0 with the area gap a0 of [QG] Theorem U */
  areaQuanta: number;
  log10AreaQuanta: number;
  /** log10 of the echo-power bound |R|^2 <= e^{-S} ([QG] Sec. 18) */
  log10EchoBound: number;
  /** Schwarzschild, massless scalar, isolated hole: t_evap = M^3 / (3 alpha) */
  log10TEvapYr: number;
  log10PageTimeYr: number;
}

/**
 * Kerr black hole of `massSun` solar masses and dimensionless spin `spin` (0 <= a* < 1).
 * Area A = 8 pi M^2 (1 + sqrt(1 - a*^2)), T_H = sqrt(1 - a*^2) / (4 pi r_+) (Planck units).
 * The lifetime ignores spin, as qchas.horizon.kerr does ([BH] Sec. 17 massless scalar).
 */
export function blackHole(massSun: number, spin = 0, a0 = minimalArea(), beta?: number): BlackHoleNumbers {
  const a = Math.min(Math.max(spin, 0), 0.999999);
  const root = Math.sqrt(1 - a * a);
  const m = massSun * MSUN_PLANCK;
  const rPlus = m * (1 + root);
  const areaLp2 = 8 * Math.PI * m * rPlus;
  const entropy = areaLp2 / 4;
  const log10TEvap = 3 * Math.log10(m) - Math.log10(3 * ALPHA_SCALAR) + Math.log10(T_PLANCK / YEAR);
  return {
    massSun,
    spin: a,
    massPlanck: m,
    rPlusKm: (rPlus * L_PLANCK) / 1e3,
    areaM2: areaLp2 * L_PLANCK ** 2,
    areaLp2,
    entropyNats: entropy,
    log10EntropyNats: Math.log10(entropy),
    log10EntropyBits: Math.log10(entropy / Math.LN2),
    temperatureK: (root / (4 * Math.PI * rPlus)) * T_PLANCK_K,
    areaQuanta: areaLp2 / a0,
    log10AreaQuanta: Math.log10(areaLp2 / a0),
    log10EchoBound: -entropy / Math.LN10,
    log10TEvapYr: log10TEvap,
    log10PageTimeYr: log10TEvap + Math.log10(pageTimeFraction(beta)),
  };
}

/**
 * [BH] Theorem F: Christodoulou-Rovelli interior volume at the endpoint,
 * V_end = 3 sqrt(3) pi M0^5 / (5 alpha) l_P^3; returned as log10 of cubic metres.
 */
export function log10InteriorVolumeEndM3(massPlanck: number, alpha = ALPHA_SCALAR): number {
  return 5 * Math.log10(massPlanck) + Math.log10((3 * Math.sqrt(3) * Math.PI) / (5 * alpha)) + 3 * Math.log10(L_PLANCK);
}

/** Log-spaced grid of n points from lo to hi (both > 0). */
export function logspace(lo: number, hi: number, n: number): number[] {
  const a = Math.log10(lo);
  const b = Math.log10(hi);
  return Array.from({ length: n }, (_, i) => 10 ** (a + ((b - a) * i) / Math.max(1, n - 1)));
}

export interface WindowCurve {
  /** beta * DeltaE */
  x: number[];
  /** exact ln tr Pi relative to S_BH(M*) */
  lnTrace: number[];
  /** [BH] Theorem G(v) first-order expression */
  theoremGv: number[];
  /** A(E2)/4G - S_BH(M*) */
  deltaS: number[];
}

/** Energy-window capacity of [BH] Theorem G(iv)-(v) as a function of x = beta DeltaE, beta = 8 pi M*. */
export function windowCurve(mStarPlanck: number, xs: number[]): WindowCurve {
  const beta = 8 * Math.PI * mStarPlanck;
  const lnTrace: number[] = [];
  const theoremGv: number[] = [];
  const deltaS: number[] = [];
  for (const x of xs) {
    const w = x / beta;
    // S_BH(M* + w) - S_BH(M*) = 4 pi w (2 M* + w), written without the cancellation that
    // (M* + w)^2 - M*^2 suffers for astrophysical M* (~1e39 in Planck units)
    const dS = 4 * Math.PI * w * (2 * mStarPlanck + w);
    lnTrace.push(capacityWindow(mStarPlanck, w).ln_trace);
    deltaS.push(dS);
    theoremGv.push(dS + Math.log(-Math.expm1(-x) / beta));
  }
  return { x: xs, lnTrace, theoremGv, deltaS };
}

// --- [QG] Theorem R decoding ---------------------------------------------------------------------

/** Values whose magnitude is below `floor` are numerical zeros (all qubits emitted): null. */
export function log10OrNull(v: number | null | undefined, floor = 1e-12): number | null {
  if (v === null || v === undefined || !Number.isFinite(v) || v < floor) return null;
  return Math.log10(v);
}

export const positiveOrNull = (v: number | null | undefined, floor = 1e-12): number | null =>
  v === null || v === undefined || !Number.isFinite(v) || v < floor ? null : v;

/** [QG] Theorem R: the number of emissions after which 2^{1+k-s} <= eps. */
export const emissionsForBound = (k: number, eps: number): number => Math.ceil(1 + k + Math.log2(1 / eps));

// --- LIGO strain --------------------------------------------------------------------------------------

/**
 * Align L1 on H1: with lag > 0 meaning L1 later than H1 and sign = -1 meaning L1 inverted,
 * H1(t) ~ sign * L1(t + lag), so L1 is drawn at t - lag with its sign flipped.
 */
export function alignL1(t: number[], l1: number[], lagMs: number, sign: number): { t: number[]; y: number[] } {
  const lag = lagMs / 1e3;
  const s = sign < 0 ? -1 : 1;
  return { t: t.map((x) => x - lag), y: l1.map((v) => s * v) };
}

export interface RingdownProjection {
  t: number[];
  y: number[];
  amplitude: number;
  phase: number;
  /** fraction of the window's variance the template explains */
  r2: number;
}

/**
 * Least-squares amplitude and phase of a damped sinusoid with fixed frequency `fHz` and
 * damping time `tauS`, started at `t0` and fitted on [t0, t0 + windowS] of (t, y):
 * y ~ e^{-(t - t0)/tau} [A cos(2 pi f (t - t0)) + B sin(2 pi f (t - t0))]. Linear in A and B.
 */
export function projectDampedSinusoid(
  t: number[],
  y: number[],
  fHz: number,
  tauS: number,
  t0: number,
  windowS: number,
): RingdownProjection | null {
  let scc = 0;
  let sss = 0;
  let scs = 0;
  let syc = 0;
  let sys = 0;
  let syy = 0;
  const idx: number[] = [];
  for (let i = 0; i < t.length; i++) {
    const dt = t[i]! - t0;
    if (dt < 0 || dt > windowS) continue;
    idx.push(i);
    const e = Math.exp(-dt / tauS);
    const c = e * Math.cos(2 * Math.PI * fHz * dt);
    const s = e * Math.sin(2 * Math.PI * fHz * dt);
    const v = y[i]!;
    scc += c * c;
    sss += s * s;
    scs += c * s;
    syc += v * c;
    sys += v * s;
    syy += v * v;
  }
  const det = scc * sss - scs * scs;
  if (idx.length < 3 || !(Math.abs(det) > 0)) return null;
  const A = (syc * sss - sys * scs) / det;
  const B = (sys * scc - syc * scs) / det;
  const tt: number[] = [];
  const yy: number[] = [];
  let res = 0;
  for (const i of idx) {
    const dt = t[i]! - t0;
    const e = Math.exp(-dt / tauS);
    const m = e * (A * Math.cos(2 * Math.PI * fHz * dt) + B * Math.sin(2 * Math.PI * fHz * dt));
    tt.push(t[i]!);
    yy.push(m);
    res += (y[i]! - m) ** 2;
  }
  return { t: tt, y: yy, amplitude: Math.hypot(A, B), phase: Math.atan2(-B, A), r2: syy > 0 ? 1 - res / syy : 0 };
}

// --- formatting of numbers with huge exponents -----------------------------------------------------------

const SUP: Record<string, string> = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };
const sup = (n: number): string => String(n).replace(/[-0-9]/g, (c) => SUP[c]!);

/**
 * 10^x as text: "1.82×10⁷²" while the mantissa is meaningful (|x| < 1e6), otherwise
 * "10^(−1.83×10⁸⁰)" (the exponent itself in scientific notation).
 */
export function fmtPow10(log10: number | null | undefined, sig = 3): string {
  if (log10 === null || log10 === undefined || !Number.isFinite(log10)) return '—';
  if (Math.abs(log10) < 1e6) {
    let e = Math.floor(log10);
    let m = 10 ** (log10 - e);
    if (Number(m.toPrecision(sig)) >= 10) {
      m /= 10;
      e += 1;
    }
    return `${m.toPrecision(sig)}×10${sup(e)}`;
  }
  const e = Math.floor(Math.log10(Math.abs(log10)));
  const m = log10 / 10 ** e;
  return `10^(${log10 < 0 ? '−' : ''}${Math.abs(m).toPrecision(sig)}×10${sup(e)})`;
}
