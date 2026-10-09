/**
 * Matched circles and ghost images on a 3-torus, [TOP] Sec. 7, Theorem F, Corollary F.1:
 * a port of ``qchas/topology/circles.py``.
 *
 * Lengths are in units of the last-scattering diameter D = 2 chi_*, so chi_* = 1/2. A torus is
 * a row basis (rescaled so its shortest vector has length d/D) and a rotation R (lattice frame
 * -> Galactic frame). Every lattice vector g with 0 < |g| < 2 chi_* gives a circle pair centred
 * on -/+ R g_hat with angular radius arccos(|g| / 2 chi_*), matched with zero twist by
 * n' = n + g / chi_* (Corollary F.1). Directions are Galactic (l, b) in degrees with unit vector
 * (cos b cos l, cos b sin l, sin b).
 */

import { D_LSS_GPC } from './constants';
import { latticeVectors, shortestVector } from './lattice';
import { type Mat3, type Vec3, add, cross, dot, identity3, matVec, norm, scale, scaleMat3, unit } from './vec';

export const CHI_STAR = 0.5;
const DEG = Math.PI / 180;

/** Galactic (l, b) in degrees. */
export type LB = readonly [number, number];

export function lbToVec(l: number, b: number): Vec3 {
  const lr = l * DEG;
  const br = b * DEG;
  return [Math.cos(br) * Math.cos(lr), Math.cos(br) * Math.sin(lr), Math.sin(br)];
}

export function vecToLb(v: readonly number[]): [number, number] {
  const n = norm(v);
  const l = ((Math.atan2(v[1]!, v[0]!) / DEG) % 360 + 360) % 360;
  const b = Math.asin(Math.max(-1, Math.min(1, v[2]! / n))) / DEG;
  return [l, b];
}

const asDirection = (d: LB | Vec3): Vec3 => (d.length === 2 ? lbToVec(d[0], d[1]) : unit([d[0], d[1], d[2]!]));

/** Rescale a row basis so its shortest lattice vector has length d/D (d/D must be positive). */
export function scaledBasis(basis: Mat3, dOverD: number): Mat3 {
  if (!(Number.isFinite(dOverD) && dOverD > 0)) throw new RangeError('d_over_D must be a positive number');
  return scaleMat3(basis, dOverD / shortestVector(basis).length);
}

export interface CirclePair {
  n: Vec3;
  g_lattice: Vec3;
  g_sky: Vec3;
  length_over_D: number;
  radius_deg: number;
  centre_minus: Vec3;
  centre_plus: Vec3;
  centre_minus_lb: [number, number];
  centre_plus_lb: [number, number];
}

/** One entry per +-g with 0 < |g| < 2 chi_*: centres -/+ R g_hat, radius arccos(|g| / 2 chi_*). */
export function circlePairs(basis: Mat3, dOverD: number, R: Mat3 = identity3(), chiStar = CHI_STAR): CirclePair[] {
  const B = scaledBasis(basis, dOverD);
  const out: CirclePair[] = [];
  for (const v of latticeVectors(B, 2 * chiStar * (1 - 1e-12))) {
    // one representative of each +-g: first nonzero coefficient positive (as in the Python)
    if ((v.n.find((k) => k !== 0) ?? 0) < 0) continue;
    if (v.length >= 2 * chiStar) continue;
    const gs = matVec(R, v.g);
    const u = scale(gs, 1 / v.length);
    const minus = scale(u, -1);
    out.push({
      n: v.n,
      g_lattice: v.g,
      g_sky: gs,
      length_over_D: v.length,
      radius_deg: Math.acos(v.length / (2 * chiStar)) / DEG,
      centre_minus: minus,
      centre_plus: u,
      centre_minus_lb: vecToLb(minus),
      centre_plus_lb: vecToLb(u),
    });
  }
  out.sort((a, b) => a.length_over_D - b.length_over_D);
  return out;
}

/** Partner direction n' = n + g/chi_* of a direction n on the circle centred at -g_hat. */
export function matchedPoint(n: LB | Vec3, g: Vec3, chiStar = CHI_STAR, atol = 1e-9): Vec3 {
  const m = add(asDirection(n), scale(g, 1 / chiStar));
  if (Math.abs(norm(m) - 1) > atol) throw new Error('n does not lie on the circle C_g (|n + g/chi_*| != 1)');
  return m;
}

/** Points on the circle centred at -g_hat and their matched partners (same azimuth phi). */
export function circlePoints(g: Vec3, nPoints = 120, chiStar = CHI_STAR): { circle: Vec3[]; twin: Vec3[] } {
  const length = norm(g);
  const u = scale(g, 1 / length);
  const e1 = unit(cross(u, Math.abs(u[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0]));
  const e2 = cross(u, e1);
  const ca = length / (2 * chiStar);
  const sa = Math.sqrt(Math.max(0, 1 - ca * ca));
  const circle: Vec3[] = [];
  const twin: Vec3[] = [];
  for (let i = 0; i < nPoints; i++) {
    const phi = (2 * Math.PI * i) / nPoints;
    const ring = add(scale(e1, sa * Math.cos(phi)), scale(e2, sa * Math.sin(phi)));
    circle.push(add(scale(u, -ca), ring));
    twin.push(add(scale(u, ca), ring));
  }
  return { circle, twin };
}

export interface GhostImage {
  n: Vec3;
  is_source: boolean;
  l: number;
  b: number;
  distance_over_D: number;
  distance_Gpc: number;
  inside_horizon: boolean;
}

/**
 * Every image x + R g of a source at x (direction, distance in D) out to `maxDistanceOverD`
 * (default twice the horizon), sorted by distance.
 */
export function ghostImages(
  sourceDirection: LB | Vec3,
  sourceDistanceOverD: number,
  basis: Mat3,
  dOverD: number,
  R: Mat3 = identity3(),
  horizonOverD = CHI_STAR,
  maxDistanceOverD?: number,
): GhostImage[] {
  if (!(sourceDistanceOverD >= 0 && Number.isFinite(sourceDistanceOverD)))
    throw new RangeError('source_distance_over_D must be >= 0');
  const x = scale(asDirection(sourceDirection), sourceDistanceOverD);
  const rmax = maxDistanceOverD ?? 2 * horizonOverD;
  const B = scaledBasis(basis, dOverD);
  const out: GhostImage[] = [];
  for (const v of latticeVectors(B, rmax + norm(x), true)) {
    const y = add(x, matVec(R, v.g));
    const r = norm(y);
    if (r > rmax) continue;
    const [l, b] = r > 0 ? vecToLb(y) : [NaN, NaN];
    out.push({
      n: v.n,
      is_source: v.n.every((k) => k === 0),
      l,
      b,
      distance_over_D: r,
      distance_Gpc: r * D_LSS_GPC.value,
      inside_horizon: r <= horizonOverD,
    });
  }
  out.sort((a, b) => a.distance_over_D - b.distance_over_D);
  return out;
}

/** The six <110> directions of FCC (one of each antipodal pair), unit vectors. */
export function fccAxesLatticeFrame(): Vec3[] {
  const out: Vec3[] = [];
  for (const [i, j] of [[0, 1], [0, 2], [1, 2]] as const)
    for (const s of [1, -1]) {
      const v: Vec3 = [0, 0, 0];
      v[i] = 1;
      v[j] = s;
      out.push(scale(v, Math.SQRT1_2));
    }
  return out;
}

/** Angular separation in degrees between two directions. */
export function angularSeparation(a: LB | Vec3, b: LB | Vec3): number {
  const c = dot(asDirection(a), asDirection(b));
  return Math.acos(Math.max(-1, Math.min(1, c))) / DEG;
}
