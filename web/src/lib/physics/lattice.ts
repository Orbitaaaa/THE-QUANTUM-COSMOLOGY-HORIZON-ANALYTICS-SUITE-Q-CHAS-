/**
 * Lattice geometry and the Casimir functional Z(2), [TOP] Sec. 6, Prop. E.2, Fig. 2:
 * a port of the Ewald evaluation in ``qchas/topology/lattice.py``.
 *
 * Bases are 3x3 arrays whose **rows** are basis vectors, g = n B for an integer row n.
 * Z(2) = sum_{g != 0} |g|^-4 at unit covolume uses the Ewald split at t = 1 of Prop. E.2:
 *
 *     Z = pi^2 [ sum_{g != 0} phi1(|g|^2) + (2 + sum_{k != 0} phi2(|k|^2)) / V - 1/2 ],
 *     phi1(x) = e^{-pi x} (1/(pi x) + 1/(pi x)^2),
 *     phi2(y) = 2 e^{-pi y} - 2 pi sqrt(y) erfc(sqrt(pi y)),
 *
 * with k running over the dual lattice inv(B)^T (no 2 pi) and V the covolume.
 */

import { Z2_CERTIFIED } from './constants';
import { type Mat3, type Vec3, clone3, det3, dot, identity3, inv3, matmul3, norm, scaleMat3, transpose3, vecMat } from './vec';

/** Ewald cutoff on |g| and |k| at unit covolume (Prop. E.2 sums |g|, |k| <= 4.5 exactly). */
export const EWALD_RADIUS = 5.0;

/** Integer boxes larger than this are first LLL-reduced (as in the Python). */
const BOX_REDUCE = 20 ** 3;
/** Enumerations larger than this many integer points are refused (Python `_BOX_MAX`). */
export const MAX_BOX_POINTS = 4_000_000;

export type CubicShape = 'SC' | 'BCC' | 'FCC';
export const SHAPES: readonly CubicShape[] = ['SC', 'BCC', 'FCC'];

const ALIASES: Record<string, CubicShape> = { SC: 'SC', CUBIC: 'SC', 'SIMPLE CUBIC': 'SC', BCC: 'BCC', FCC: 'FCC' };

/** Map 'cubic' / 'SC' / 'bcc' / ... to one of 'SC', 'BCC', 'FCC'. */
export function canonicalShape(shape: string): CubicShape {
  const s = ALIASES[shape.trim().toUpperCase()];
  if (!s) throw new Error(`unknown lattice shape ${shape}`);
  return s;
}

/** Primitive basis (rows) of SC, BCC or FCC in the conventional cube frame (cube side 1). */
export function cubicBasis(shape: string): Mat3 {
  const s = canonicalShape(shape);
  if (s === 'SC') return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  if (s === 'BCC') return [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]];
  return [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]];
}

/**
 * Three unit vectors with equal pairwise angle alpha, symmetric about (1,1,1).
 * alpha = 60 deg is FCC, 90 deg SC, arccos(-1/3) = 109.47 deg BCC; valid for 0 < alpha < 120.
 */
export function rhombohedralBasis(alphaDeg: number): Mat3 {
  const c = Math.cos((alphaDeg * Math.PI) / 180);
  if (!(c > -0.5 + 1e-12 && c < 1 - 1e-12)) throw new Error('rhombohedral angle must lie in (0, 120) degrees');
  const x = Math.sqrt(1 - c);
  const y = (Math.sqrt(1 + 2 * c) - x) / 3;
  return [
    [x + y, y, y],
    [y, x + y, y],
    [y, y, x + y],
  ];
}

export const gram = (B: Mat3): Mat3 => matmul3(B, transpose3(B));

export const covolume = (B: Mat3): number => Math.abs(det3(B));

/** Dual-lattice basis (rows): 2 pi inv(B)^T, or inv(B)^T with twoPi = false. */
export const dualBasis = (B: Mat3, twoPi = true): Mat3 => scaleMat3(transpose3(inv3(B)), twoPi ? 2 * Math.PI : 1);

function gso(B: Mat3): { Bs: Vec3[]; mu: number[][] } {
  const Bs: Vec3[] = [];
  const mu = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) {
    let v: Vec3 = [...B[i]!];
    for (let j = 0; j < i; j++) {
      const bj = Bs[j]!;
      const m = dot(B[i]!, bj) / dot(bj, bj);
      mu[i]![j] = m;
      v = [v[0] - m * bj[0], v[1] - m * bj[1], v[2] - m * bj[2]];
    }
    Bs.push(v);
  }
  return { Bs, mu };
}

/** Throw unless B is a finite basis of rank 3. */
export function checkBasis(B: Mat3): void {
  const flat = [...B[0], ...B[1], ...B[2]];
  if (!flat.every(Number.isFinite)) throw new Error('a lattice basis must be a finite 3x3 array');
  const scale = Math.max(...flat.map(Math.abs));
  if (!(scale > 0) || !(Math.abs(det3(B)) > 1e-13 * scale ** 3)) throw new Error('lattice basis is singular (rank < 3)');
}

/**
 * LLL reduction of a row basis. Returns the reduced basis and the unimodular U with
 * `B_red = U B` (same lattice), as `lll_reduce` in the Python.
 */
export function lllReduceWithTransform(B0: Mat3, delta = 0.99): { B: Mat3; U: Mat3 } {
  const B = clone3(B0);
  const U = identity3();
  let k = 1;
  let guard = 0;
  while (k < 3) {
    if (guard++ > 100_000) throw new Error('LLL reduction did not terminate (basis too ill-conditioned)');
    for (let j = k - 1; j >= 0; j--) {
      const q = Math.round(gso(B).mu[k]![j]!);
      if (q)
        for (let a = 0; a < 3; a++) {
          B[k]![a] = B[k]![a]! - q * B[j]![a]!;
          U[k]![a] = U[k]![a]! - q * U[j]![a]!;
        }
    }
    const { Bs, mu } = gso(B);
    const m = mu[k]![k - 1]!;
    if (dot(Bs[k]!, Bs[k]!) >= (delta - m * m) * dot(Bs[k - 1]!, Bs[k - 1]!)) k++;
    else {
      [B[k - 1], B[k]] = [B[k]!, B[k - 1]!];
      [U[k - 1], U[k]] = [U[k]!, U[k - 1]!];
      k = Math.max(k - 1, 1);
    }
  }
  return { B, U };
}

/** LLL reduction of a row basis; returns the reduced basis (same lattice). */
export const lllReduce = (B0: Mat3, delta = 0.99): Mat3 => lllReduceWithTransform(B0, delta).B;

/** Rescale to unit covolume ('covolume') or unit shortest vector ('shortest'). */
export function normalize(B: Mat3, to: 'covolume' | 'shortest' = 'covolume'): Mat3 {
  if (to === 'covolume') return scaleMat3(B, 1 / Math.cbrt(covolume(B)));
  return scaleMat3(B, 1 / shortestVector(B).length);
}

export interface LatticeVector {
  /** integer coefficients, g = n B */
  n: Vec3;
  g: Vec3;
  length: number;
}

function boxFor(radius: number, invDiag: readonly number[]): [number, number, number] {
  return invDiag.map((d) => Math.floor(radius * Math.sqrt(Math.max(d, 0)) + 1e-9)) as [number, number, number];
}

const boxPoints = (box: readonly number[]): number => box.reduce((p, n) => p * (2 * n + 1), 1);

function boxOf(B: Mat3, radius: number): [number, number, number] {
  const hinv = inv3(gram(B));
  return boxFor(radius, [hinv[0][0], hinv[1][1], hinv[2][2]]);
}

function refuse(points: number): never {
  throw new RangeError(
    `lattice enumeration needs ${points.toExponential(1)} integer points (limit ${MAX_BOX_POINTS.toExponential(0)}): ` +
      'the radius is too large for this lattice (e.g. d/D too small)',
  );
}

/**
 * All lattice vectors with 0 < |g| <= radius (both signs), sorted by length. A skewed basis is
 * LLL-reduced first when that shrinks the enumeration box; boxes above `MAX_BOX_POINTS` throw
 * a RangeError instead of freezing the page.
 */
export function latticeVectors(B: Mat3, radius: number, includeZero = false): LatticeVector[] {
  checkBasis(B);
  if (!(radius >= 0 && Number.isFinite(radius))) throw new RangeError(`invalid enumeration radius ${radius}`);
  let E = B;
  let U: Mat3 | null = null;
  let box = boxOf(B, radius);
  if (boxPoints(box) > BOX_REDUCE) {
    const red = lllReduceWithTransform(B);
    const boxR = boxOf(red.B, radius);
    if (boxPoints(boxR) < boxPoints(box)) [E, U, box] = [red.B, red.U, boxR];
  }
  if (boxPoints(box) > MAX_BOX_POINTS) refuse(boxPoints(box));
  const [n0, n1, n2] = box;
  const r2max = radius * radius * (1 + 1e-12);
  const out: LatticeVector[] = [];
  for (let a = -n0; a <= n0; a++)
    for (let b = -n1; b <= n1; b++)
      for (let c = -n2; c <= n2; c++) {
        if (!includeZero && a === 0 && b === 0 && c === 0) continue;
        const g = vecMat([a, b, c], E);
        const r2 = dot(g, g);
        if (r2 <= r2max) out.push({ n: U ? vecMat([a, b, c], U) : [a, b, c], g, length: Math.sqrt(r2) });
      }
  out.sort((x, y) => x.length - y.length);
  return out;
}

/** A shortest nonzero lattice vector and its length (`n` refers to the given basis B). */
export function shortestVector(B: Mat3): LatticeVector {
  checkBasis(B);
  const { B: R, U } = lllReduceWithTransform(B);
  const cand = latticeVectors(R, norm(R[0]) * (1 + 1e-9));
  const best = cand[0];
  if (!best) throw new Error('degenerate lattice');
  return { ...best, n: vecMat(best.n, U) };
}

/** int_1^inf t e^{-pi x t} dt. */
export function phi1(x: number): number {
  const a = Math.PI * x;
  return Math.exp(-a) * (1 / a + 1 / (a * a));
}

/** int_1^inf u^{-3/2} e^{-pi y u} du = 2 e^{-pi y} - 2 pi sqrt(y) erfc(sqrt(pi y)). */
export function phi2(y: number): number {
  return 2 * Math.exp(-Math.PI * y) - 2 * Math.PI * Math.sqrt(y) * erfc(Math.sqrt(Math.PI * y));
}

/**
 * Complementary error function to ~1e-16 absolute: the positive-term Taylor series of erf
 * below x = 2.5, the continued fraction (evaluated backwards) above.
 */
export function erfc(x: number): number {
  if (Number.isNaN(x)) return NaN;
  if (x < 0) return 2 - erfc(-x);
  if (x < 2.5) {
    const x2 = x * x;
    let term = x;
    let sum = x;
    for (let n = 1; n < 200; n++) {
      term *= (2 * x2) / (2 * n + 1);
      sum += term;
      if (term < 1e-17 * sum) break;
    }
    return 1 - ((2 / Math.sqrt(Math.PI)) * Math.exp(-x2) * sum);
  }
  if (x > 27) return 0;
  let f = x;
  for (let n = 120; n >= 1; n--) f = x + n / 2 / f;
  return Math.exp(-x * x) / (Math.sqrt(Math.PI) * f);
}

/**
 * Z(h) = sum_{n != 0} (n^T h n)^-2 for a Gram matrix h of any covolume (Ewald, t = 1).
 * The split is tuned to covolume ~ 1 (`radius` refers to unit covolume).
 */
export function epsteinZ2Gram(h: Mat3, radius = EWALD_RADIUS): number {
  const hinv = inv3(h);
  const V = Math.sqrt(det3(h));
  const R2 = radius * radius;
  const half = (box: [number, number, number], metric: Mat3, f: (x: number) => number): number => {
    let s = 0;
    const [N0, N1, N2] = box;
    for (let a = 0; a <= N0; a++)
      for (let b = -N1; b <= N1; b++)
        for (let c = -N2; c <= N2; c++) {
          // one of each +-pair: first nonzero coefficient positive
          if (a === 0 && (b < 0 || (b === 0 && c <= 0))) continue;
          const n: Vec3 = [a, b, c];
          const x = dot(n, [dot(metric[0], n), dot(metric[1], n), dot(metric[2], n)]);
          if (x <= R2) s += f(x);
        }
    return 2 * s;
  };
  const nbox = boxFor(radius, [hinv[0][0], hinv[1][1], hinv[2][2]]);
  const kbox = boxFor(radius, [h[0][0], h[1][1], h[2][2]]);
  if (boxPoints(nbox) + boxPoints(kbox) > MAX_BOX_POINTS) refuse(boxPoints(nbox) + boxPoints(kbox));
  const s1 = half(nbox, h, phi1);
  const s2 = half(kbox, hinv, phi2);
  return Math.PI ** 2 * (s1 + (2 + s2) / V - 0.5);
}

/**
 * Z(2) = sum_{g != 0} |g|^-4 of the lattice with row basis B; by default at unit covolume
 * as in [TOP] Prop. E.2 (certified values in `Z2_CERTIFIED`).
 */
export function epsteinZ2(B: Mat3, unitCovolume = true, radius = EWALD_RADIUS): number {
  checkBasis(B);
  const V = covolume(B);
  const Bu = lllReduce(scaleMat3(B, 1 / Math.cbrt(V)));
  const z = epsteinZ2Gram(gram(Bu), radius);
  return unitCovolume ? z : z * V ** (-4 / 3);
}

/** Z(2) along the rhombohedral family ([TOP] Fig. 2); alpha = 120 deg is the planar limit. */
export function rhombohedralCurve(alphaDeg: number[]): { alpha_deg: number[]; Z: number[] } {
  return { alpha_deg: alphaDeg, Z: alphaDeg.map((a) => epsteinZ2(rhombohedralBasis(a))) };
}

export const certifiedZ2 = Z2_CERTIFIED;
