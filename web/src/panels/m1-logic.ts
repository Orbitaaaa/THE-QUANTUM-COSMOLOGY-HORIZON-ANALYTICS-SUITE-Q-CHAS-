/**
 * Pure logic shared by the Module 1 panels (manifolds, lattice, bounce, ghost-twin): the
 * principle filter of the [TOP] Main Theorem, rotations for the torus orientation, the
 * matched-circle snap and twin of Corollary F.1, the comoving distance-redshift relation used to
 * label source distances, and small array helpers. No DOM access; unit-tested in
 * tests/m1-logic.test.ts.
 */

import { PLANCK2018 } from '../lib/physics/constants';
import { type CirclePair, CHI_STAR } from '../lib/physics/circles';
import { type Mat3, type Vec3, dot, matmul3, norm, scale } from '../lib/physics/vec';

// --- principles P1-P4 and the Main Theorem ([TOP] Secs. 3-4) ----------------------------------

export type PrincipleKey = 'P1' | 'P2' | 'P3' | 'P4';
export const PRINCIPLE_KEYS: readonly PrincipleKey[] = ['P1', 'P2', 'P3', 'P4'];

/** The part of a topology.json manifold entry the selection logic reads. */
export interface ManifoldPrinciples {
  label: string;
  principles: Record<PrincipleKey, string>;
}

/** A principle holds on a manifold only when its status is exactly "allowed". */
export const holds = (status: string | undefined): boolean => status === 'allowed';

/** Does manifold `m` satisfy every principle in `active`? (The empty set is satisfied by all.) */
export function satisfies(m: ManifoldPrinciples, active: Iterable<PrincipleKey>): boolean {
  for (const p of active) if (!holds(m.principles[p])) return false;
  return true;
}

/** Labels of the manifolds compatible with every principle in `active` (Python `survivors`). */
export function survivors(ms: ManifoldPrinciples[], active: Iterable<PrincipleKey>): string[] {
  const act = [...active];
  for (const p of act) if (!PRINCIPLE_KEYS.includes(p)) throw new Error(`unknown principle ${String(p)}`);
  return ms.filter((m) => satisfies(m, act)).map((m) => m.label);
}

/** Every principle the manifold fails ('excluded' or 'no spin structure'). */
export const failedPrinciples = (m: ManifoldPrinciples): PrincipleKey[] => PRINCIPLE_KEYS.filter((p) => !holds(m.principles[p]));

/** The active principles that exclude `m`. */
export const excludedBy = (m: ManifoldPrinciples, active: Iterable<PrincipleKey>): PrincipleKey[] => {
  const act = new Set(active);
  return failedPrinciples(m).filter((p) => act.has(p));
};

/** All subsets of P1-P4 (16), smallest first, in key order within a size. */
export function principleSubsets(): PrincipleKey[][] {
  const out: PrincipleKey[][] = [];
  for (let mask = 0; mask < 16; mask++) out.push(PRINCIPLE_KEYS.filter((_, i) => mask & (1 << i)));
  return out.sort((a, b) => a.length - b.length || PRINCIPLE_KEYS.indexOf(a[0]!) - PRINCIPLE_KEYS.indexOf(b[0]!));
}

/** Inclusion-minimal sets of principles whose survivors are exactly [target] (Python `minimal_selecting_sets`). */
export function minimalSelectingSets(ms: ManifoldPrinciples[], target = 'E1'): PrincipleKey[][] {
  const found: PrincipleKey[][] = [];
  for (const s of principleSubsets()) {
    const surv = survivors(ms, s);
    if (surv.length === 1 && surv[0] === target && !found.some((f) => f.every((p) => s.includes(p)))) found.push(s);
  }
  return found;
}

/**
 * Manifolds that fail exactly `n` principles. With n = 1 this tests the sentence of [TOP] Secs. 4
 * and 5.6, "every manifold except E1 fails at least two principles": the paper's own table gives
 * E11, E16 (no screw) and E18, which fail only P2.
 */
export const failingExactly = (ms: ManifoldPrinciples[], n: number): string[] =>
  ms.filter((m) => failedPrinciples(m).length === n).map((m) => m.label);

/** Format a set of principles: ['P2', 'P3'] -> "{P2, P3}". */
export const fmtSet = (s: readonly string[]): string => `{${s.join(', ')}}`;

// --- rotations ---------------------------------------------------------------------------------

const DEG = Math.PI / 180;

const rotZ = (t: number): Mat3 => [
  [Math.cos(t), -Math.sin(t), 0],
  [Math.sin(t), Math.cos(t), 0],
  [0, 0, 1],
];
const rotY = (t: number): Mat3 => [
  [Math.cos(t), 0, Math.sin(t)],
  [0, 1, 0],
  [-Math.sin(t), 0, Math.cos(t)],
];

/** R = Rz(alpha) Ry(beta) Rz(gamma), angles in degrees (lattice frame -> Galactic frame). */
export const eulerZYZ = (alphaDeg: number, betaDeg: number, gammaDeg: number): Mat3 =>
  matmul3(matmul3(rotZ(alphaDeg * DEG), rotY(betaDeg * DEG)), rotZ(gammaDeg * DEG));

const wrap360 = (x: number): number => ((x % 360) + 360) % 360;

/** ZYZ Euler angles (degrees; alpha, gamma in [0, 360), beta in [0, 180]) of a rotation matrix. */
export function eulerFromMatrix(R: Mat3): [number, number, number] {
  const cb = Math.max(-1, Math.min(1, R[2][2]));
  const beta = Math.acos(cb);
  if (Math.sin(beta) > 1e-9) {
    const alpha = Math.atan2(R[1][2], R[0][2]);
    const gamma = Math.atan2(R[2][1], -R[2][0]);
    return [wrap360(alpha / DEG), beta / DEG, wrap360(gamma / DEG)];
  }
  // gimbal lock: only alpha +- gamma is defined; put it all in alpha
  if (cb > 0) return [wrap360(Math.atan2(R[1][0], R[0][0]) / DEG), 0, 0];
  return [wrap360(Math.atan2(-R[1][0], -R[0][0]) / DEG), 180, 0];
}

/** Haar-random rotation from three uniforms in [0, 1) (Shoemake's quaternion method). */
export function haarRotation(rand: () => number = Math.random): Mat3 {
  const u1 = rand();
  const u2 = rand();
  const u3 = rand();
  const a = Math.sqrt(1 - u1);
  const b = Math.sqrt(u1);
  const [x, y, z, w] = [a * Math.sin(2 * Math.PI * u2), a * Math.cos(2 * Math.PI * u2), b * Math.sin(2 * Math.PI * u3), b * Math.cos(2 * Math.PI * u3)];
  return [
    [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
    [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
    [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
  ];
}

/** Coerce a 3x3 array (e.g. topology.circles.fcc_preregistered.R) to a Mat3, or null. */
export function asMat3(m: unknown): Mat3 | null {
  if (!Array.isArray(m) || m.length !== 3) return null;
  const rows = m.map((r) => (Array.isArray(r) && r.length === 3 && r.every((x) => typeof x === 'number' && Number.isFinite(x)) ? [r[0], r[1], r[2]] : null));
  return rows.every(Boolean) ? (rows as Mat3) : null;
}

// --- matched circles: snap a picked direction and find its twin (Corollary F.1) ---------------

export interface CircleHit {
  /** index of the pair in the `pairs` array */
  index: number;
  /** +1: the circle about -g_hat (twin n + g/chi*); -1: the circle about +g_hat (twin n - g/chi*) */
  side: 1 | -1;
  /** the picked direction projected onto the circle (unit vector) */
  point: Vec3;
  /** its matched partner n' = n +- g/chi* (unit vector, same azimuth on the antipodal circle) */
  twin: Vec3;
  /** angular distance of the pick from the circle, degrees */
  offsetDeg: number;
}

function perpUnit(n: Vec3, c: Vec3): Vec3 {
  let e: Vec3 = [n[0] - dot(n, c) * c[0], n[1] - dot(n, c) * c[1], n[2] - dot(n, c) * c[2]];
  if (norm(e) < 1e-12) {
    // the pick is at the circle centre: any azimuth will do
    const t: Vec3 = Math.abs(c[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    e = [t[0] - dot(t, c) * c[0], t[1] - dot(t, c) * c[1], t[2] - dot(t, c) * c[2]];
  }
  return scale(e, 1 / norm(e));
}

/**
 * The circle nearest to the direction `n` (unit vector) among all matched-circle pairs, if it lies
 * within `tolDeg`. The pick is snapped onto that circle and its twin n' = n + g/chi* computed
 * exactly, so |n'| = 1 holds to rounding (Corollary F.1: zero twist, antipodal circles).
 */
export function snapToCircle(n: Vec3, pairs: CirclePair[], tolDeg = 1.5, chiStar = CHI_STAR): CircleHit | null {
  let best: CircleHit | null = null;
  pairs.forEach((p, index) => {
    const alpha = p.radius_deg * DEG;
    for (const side of [1, -1] as const) {
      const c = side === 1 ? p.centre_minus : p.centre_plus;
      const theta = Math.acos(Math.max(-1, Math.min(1, dot(n, c))));
      const offsetDeg = Math.abs(theta - alpha) / DEG;
      if (offsetDeg > tolDeg || (best && best.offsetDeg <= offsetDeg)) continue;
      const e = perpUnit(n, c);
      const point: Vec3 = [
        Math.cos(alpha) * c[0] + Math.sin(alpha) * e[0],
        Math.cos(alpha) * c[1] + Math.sin(alpha) * e[1],
        Math.cos(alpha) * c[2] + Math.sin(alpha) * e[2],
      ];
      const g = scale(p.g_sky, side / chiStar);
      best = { index, side, point, twin: [point[0] + g[0], point[1] + g[1], point[2] + g[2]], offsetDeg };
    }
  });
  return best;
}

/**
 * Smallest source distance (fraction of chi*) at which some image of the source can lie inside
 * the last-scattering sphere: every image x + g has |x + g| >= |g| - |x| >= d - r, with equality
 * for x antiparallel to the shortest g, so images inside chi* need r >= d - chi*.
 * In units of chi* this is 2 d/D - 1 (0.97 at the circle-search bound d = 0.985 D).
 */
export const imageThresholdChi = (dOverD: number): number => Math.max(0, 2 * dOverD - 1);

// --- comoving distance and redshift (flat LCDM, Planck 2018 parameters) ------------------------

const C_KM_S = 299_792.458;
const OMEGA_GAMMA_H2 = 2.4728e-5; // photons at T_CMB = 2.7255 K
const NEFF = 3.046;
const OMEGA_NU_MASSIVE_H2 = 0.06 / 93.14; // one 0.06 eV neutrino, counted as matter

export interface Cosmology {
  h: number;
  omegaM: number;
  omegaR: number;
  omegaL: number;
}

/** Planck 2018 best fit ([TOP] Sec. 14) as a flat LCDM + radiation background. */
export function planck2018Cosmology(): Cosmology {
  const h = PLANCK2018.H0 / 100;
  const omegaM = (PLANCK2018.ombh2 + PLANCK2018.omch2 + OMEGA_NU_MASSIVE_H2) / (h * h);
  const omegaR = (OMEGA_GAMMA_H2 * (1 + 0.2271 * NEFF)) / (h * h);
  return { h, omegaM, omegaR, omegaL: 1 - omegaM - omegaR };
}

/**
 * Comoving distance to redshift z in Mpc: (c/H0) int_0^z dz'/E(z'), integrated by Simpson's rule
 * in s = ln(1 + z). Agrees with chi* = 13 870 Mpc ([TOP] Sec. 7, CAMB) to 0.05 % at z* = 1089.9.
 */
export function comovingDistanceMpc(z: number, cosmo: Cosmology = planck2018Cosmology(), n = 400): number {
  if (!(z >= 0) || !Number.isFinite(z)) throw new RangeError('redshift must be >= 0');
  if (z === 0) return 0;
  const { h, omegaM, omegaR, omegaL } = cosmo;
  const smax = Math.log1p(z);
  const f = (s: number): number => {
    const a1 = Math.exp(s); // 1 + z
    return a1 / Math.sqrt(omegaM * a1 ** 3 + omegaR * a1 ** 4 + omegaL);
  };
  const m = n % 2 ? n + 1 : n;
  const hs = smax / m;
  let sum = f(0) + f(smax);
  for (let i = 1; i < m; i++) sum += (i % 2 ? 4 : 2) * f(i * hs);
  return (C_KM_S / (100 * h)) * (sum * hs) / 3;
}

/** Inverse of comovingDistanceMpc by bisection in ln(1 + z) (z <= 1e5). */
export function redshiftForDistance(chiMpc: number, cosmo: Cosmology = planck2018Cosmology()): number {
  if (!(chiMpc >= 0)) throw new RangeError('distance must be >= 0');
  if (chiMpc === 0) return 0;
  let lo = 0;
  let hi = Math.log1p(1e5);
  if (comovingDistanceMpc(Math.expm1(hi), cosmo) < chiMpc) return Infinity;
  for (let i = 0; i < 60; i++) {
    const mid = (lo + hi) / 2;
    if (comovingDistanceMpc(Math.expm1(mid), cosmo, 200) < chiMpc) lo = mid;
    else hi = mid;
  }
  return Math.expm1((lo + hi) / 2);
}

// --- small numeric helpers ------------------------------------------------------------------------

/** Linear interpolation of (xs, ys) at x (xs ascending); null outside the range. */
export function interp(xs: readonly number[], ys: readonly number[], x: number): number | null {
  if (!xs.length || x < xs[0]! || x > xs[xs.length - 1]!) return null;
  let lo = 0;
  let hi = xs.length - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (xs[mid]! <= x) lo = mid;
    else hi = mid;
  }
  const x0 = xs[lo]!;
  const x1 = xs[hi]!;
  return x1 === x0 ? ys[lo]! : ys[lo]! + ((ys[hi]! - ys[lo]!) * (x - x0)) / (x1 - x0);
}

/** Trace identity of [TOP] Prop. E.2 at a critical point: tr H = 2 Z / 3. Returns both sides. */
export function hessianTraceCheck(Z: number, eigenvalues: readonly number[]): { trace: number; twoZ3: number; residual: number } {
  const trace = eigenvalues.reduce((s, v) => s + v, 0);
  const twoZ3 = (2 * Z) / 3;
  return { trace, twoZ3, residual: trace - twoZ3 };
}

/**
 * Concatenate trajectories into one line series with gaps (null) between them, plotting
 * Z - zRef; points where Z - zRef <= floor are dropped (converged to zRef within the bundle's
 * rounding), so the series can sit on a log axis.
 */
export function flowSeries(
  trajs: readonly { tau: number[]; Z: number[] }[],
  zRef: number,
  floor = 1e-7,
): { x: number[]; y: (number | null)[] } {
  const x: number[] = [];
  const y: (number | null)[] = [];
  trajs.forEach((t, k) => {
    if (k) {
      x.push(x[x.length - 1] ?? 0);
      y.push(null);
    }
    for (let i = 0; i < t.tau.length; i++) {
      const dz = (t.Z[i] ?? NaN) - zRef;
      x.push(t.tau[i]!);
      y.push(Number.isFinite(dz) && dz > floor ? dz : null);
    }
  });
  return { x, y };
}

/** Classical closed radiation universe a(T) = a_max |cos T| ([QG] Fig. 1 left, dotted). */
export const classicalScaleFactor = (aMax: number, T: number): number => aMax * Math.abs(Math.cos(T));

/** SPHEREx band key ('1'..'6') whose wavelength range matches `bandUm`, or null. */
export function bandKey(bandUm: readonly number[], bands: Record<string, readonly number[]>): string | null {
  for (const [k, v] of Object.entries(bands)) if (Math.abs(v[0]! - bandUm[0]!) < 0.011 && Math.abs(v[1]! - bandUm[1]!) < 0.011) return k;
  return null;
}
