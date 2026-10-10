/**
 * Pure helpers of the Resonance Shell Inspector and Sky Statistics panels (no DOM):
 *
 * * dual-lattice shells of the cubic, BCC and FCC tori ([TOP] Sec. 11.6) with the same
 *   normalization as qchas.cmb.resonance: lattice scaled so its shortest vector is 1, dual
 *   vectors m = B^{-T} n (no 2 pi), |m|^2 = the bundle's `m2`, k chi* = pi |m| / (d/D);
 * * direction families of a shell in the cube frame (<100>, <110>, <111>, <221>, ...);
 * * histograms and one-sided Monte Carlo p-values of null distributions;
 * * the angle between two headless axes (quadrupole / octopole alignment, [TOP] Sec. 9.1).
 */

import { cubicBasis } from '../lib/physics/lattice';
import { type Mat3, type Vec3, inv3, norm, scaleMat3, transpose3 } from '../lib/physics/vec';

const DEG = 180 / Math.PI;

/** Lattice basis (rows) of 'cubic' | 'SC' | 'BCC' | 'FCC' scaled so the shortest vector has length 1. */
export function unitShortestBasis(shape: string): Mat3 {
  const B = cubicBasis(shape);
  const shortest = Math.min(...B.map((r) => norm(r)));
  return scaleMat3(B, 1 / shortest);
}

/** Dual basis (rows) g_i with b_i . g_j = delta_ij, i.e. B^{-T}, without the 2 pi. */
export function dualBasisUnit(shape: string): Mat3 {
  return transpose3(inv3(unitShortestBasis(shape)));
}

/** Every dual-lattice vector m with |m|^2 = m2 (relative tolerance 1e-6), sorted for stable output. */
export function shellVectors(shape: string, m2: number): Vec3[] {
  const B = unitShortestBasis(shape);
  const G = dualBasisUnit(shape);
  // n_i = m . b_i, so |n_i| <= |m| |b_i|
  const N = Math.ceil(Math.sqrt(m2) * Math.max(...B.map((r) => norm(r)))) + 1;
  const out: Vec3[] = [];
  const tol = 1e-6 * Math.max(1, m2);
  for (let i = -N; i <= N; i++)
    for (let j = -N; j <= N; j++)
      for (let k = -N; k <= N; k++) {
        if (!i && !j && !k) continue;
        const [g0, g1, g2] = G;
        const v: Vec3 = [i * g0[0] + j * g1[0] + k * g2[0], i * g0[1] + j * g1[1] + k * g2[1], i * g0[2] + j * g1[2] + k * g2[2]];
        const q = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        if (Math.abs(q - m2) <= tol) out.push(v);
      }
  return out.sort((a, b) => b[2] - a[2] || b[1] - a[1] || b[0] - a[0]);
}

/** k chi* of the shell |m|^2 = m2 for a torus of size d/D (chi* = D/2). */
export const shellKChi = (m2: number, dOverD: number): number => (Math.PI * Math.sqrt(m2)) / dOverD;

const gcd = (a: number, b: number): number => (b ? gcd(b, a % b) : a);

/** Direction family of a vector in the cube frame, e.g. [3, 0, 0] -> '100', [2, -2, 1] -> '221'. */
export function directionFamily(v: readonly number[]): string {
  const a = v.map((x) => Math.abs(x)).sort((x, y) => y - x);
  const nz = a.filter((x) => x > 1e-9);
  if (!nz.length) return '000';
  const min = Math.min(...nz);
  // scale to small integers: try multipliers until every component is integral
  for (let mult = 1; mult <= 12; mult++) {
    const s = a.map((x) => (x / min) * mult);
    if (s.every((x) => Math.abs(x - Math.round(x)) < 1e-6)) {
      const ints = s.map((x) => Math.round(x));
      const g = ints.reduce((p, x) => gcd(p, x), 0) || 1;
      return ints.map((x) => x / g).join('');
    }
  }
  return a.map((x) => x.toFixed(2)).join(',');
}

/** Counts per direction family, most populous first: [['100', 6], ['221', 24]]. */
export function familyCounts(vectors: readonly (readonly number[])[]): [string, number][] {
  const m = new Map<string, number>();
  for (const v of vectors) {
    const f = directionFamily(v);
    m.set(f, (m.get(f) ?? 0) + 1);
  }
  return [...m.entries()].sort((a, b) => a[0].length - b[0].length || a[0].localeCompare(b[0]) || b[1] - a[1]);
}

/** Azimuth (deg, in (-180, 180]) and elevation (deg) of a vector in the lattice frame. */
export function azEl(v: readonly number[]): [number, number] {
  const [x, y, z] = v as [number, number, number];
  const r = Math.hypot(x, y, z) || 1;
  const az = Math.atan2(y, x) * DEG;
  return [Math.abs(az + 180) < 1e-9 ? 180 : az, Math.asin(Math.max(-1, Math.min(1, z / r))) * DEG];
}

export interface Histogram {
  /** bin centres */
  x: number[];
  counts: number[];
  edges: number[];
  width: number;
}

/** Equal-width histogram over [lo, hi] (defaults: data range padded to whole bins). */
export function histogram(values: readonly number[], bins: number, lo?: number, hi?: number): Histogram {
  const finite = values.filter((v) => Number.isFinite(v));
  const a = lo ?? Math.min(...finite);
  let b = hi ?? Math.max(...finite);
  if (!(b > a)) b = a + 1;
  const width = (b - a) / bins;
  const counts = new Array<number>(bins).fill(0);
  for (const v of finite) {
    if (v < a || v > b) continue;
    counts[Math.min(bins - 1, Math.floor((v - a) / width))]! += 1;
  }
  const edges = Array.from({ length: bins + 1 }, (_, i) => a + i * width);
  return { x: counts.map((_, i) => a + (i + 0.5) * width), counts, edges, width };
}

/** Monte Carlo p-value: fraction of null values >= observed. */
export function pValueGE(nulls: readonly number[], observed: number): number {
  if (!nulls.length) return NaN;
  return nulls.filter((v) => v >= observed).length / nulls.length;
}

/** Median of a sample. */
export function median(values: readonly number[]): number {
  const s = [...values].sort((a, b) => a - b);
  const n = s.length;
  if (!n) return NaN;
  return n % 2 ? s[(n - 1) / 2]! : 0.5 * (s[n / 2 - 1]! + s[n / 2]!);
}

/** Unit vector of Galactic (l, b) in degrees. */
export function lbVec(l: number, b: number): Vec3 {
  const lr = l / DEG;
  const br = b / DEG;
  return [Math.cos(br) * Math.cos(lr), Math.cos(br) * Math.sin(lr), Math.sin(br)];
}

/** Angle in degrees between two headless axes (n and -n are the same axis): in [0, 90]. */
export function axisAngleDeg(a: readonly [number, number], b: readonly [number, number]): number {
  const u = lbVec(a[0], a[1]);
  const v = lbVec(b[0], b[1]);
  const c = Math.abs(u[0] * v[0] + u[1] * v[1] + u[2] * v[2]);
  return Math.acos(Math.min(1, c)) * DEG;
}

/** The antipode of a Galactic direction. */
export const antipode = (l: number, b: number): [number, number] => [(l + 180) % 360, -b];

/**
 * Probability that two isotropic random axes are at least as aligned as `angleDeg`:
 * P(|cos| >= cos angle) = 1 - cos angle.
 */
export const isotropicAlignmentP = (angleDeg: number): number => 1 - Math.cos(angleDeg / DEG);
