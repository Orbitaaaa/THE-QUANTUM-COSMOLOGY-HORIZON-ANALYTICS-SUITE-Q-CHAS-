/**
 * Small dense linear algebra for 3-vectors and 3x3 matrices (row-major, rows = basis vectors),
 * shared by the lattice, circle and sky modules.
 */

export type Vec3 = [number, number, number];
export type Mat3 = [Vec3, Vec3, Vec3];

export const dot = (a: readonly number[], b: readonly number[]): number =>
  a[0]! * b[0]! + a[1]! * b[1]! + a[2]! * b[2]!;

export const norm = (a: readonly number[]): number => Math.sqrt(dot(a, a));

export const add = (a: Vec3, b: Vec3): Vec3 => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];

export const sub = (a: Vec3, b: Vec3): Vec3 => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];

export const scale = (a: Vec3, s: number): Vec3 => [a[0] * s, a[1] * s, a[2] * s];

export const cross = (a: Vec3, b: Vec3): Vec3 => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];

export function unit(a: Vec3): Vec3 {
  const n = norm(a);
  if (!(n > 0)) throw new Error('zero vector has no direction');
  return scale(a, 1 / n);
}

export const identity3 = (): Mat3 => [
  [1, 0, 0],
  [0, 1, 0],
  [0, 0, 1],
];

export const clone3 = (m: readonly (readonly number[])[]): Mat3 =>
  [0, 1, 2].map((i) => [m[i]![0]!, m[i]![1]!, m[i]![2]!]) as Mat3;

export const transpose3 = (m: Mat3): Mat3 => [
  [m[0][0], m[1][0], m[2][0]],
  [m[0][1], m[1][1], m[2][1]],
  [m[0][2], m[1][2], m[2][2]],
];

export function matmul3(a: Mat3, b: Mat3): Mat3 {
  const out = identity3();
  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++) out[i]![j] = a[i]![0] * b[0]![j]! + a[i]![1] * b[1]![j]! + a[i]![2] * b[2]![j]!;
  return out;
}

/** Matrix times column vector, ``M v``. */
export const matVec = (m: Mat3, v: readonly number[]): Vec3 => [dot(m[0], v), dot(m[1], v), dot(m[2], v)];

/** Row vector times matrix, ``n M`` (a lattice vector ``g = n B`` for a row basis B). */
export const vecMat = (n: readonly number[], m: Mat3): Vec3 => [
  n[0]! * m[0][0] + n[1]! * m[1][0] + n[2]! * m[2][0],
  n[0]! * m[0][1] + n[1]! * m[1][1] + n[2]! * m[2][1],
  n[0]! * m[0][2] + n[1]! * m[1][2] + n[2]! * m[2][2],
];

export const det3 = (m: Mat3): number =>
  m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
  m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
  m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

export function inv3(m: Mat3): Mat3 {
  const d = det3(m);
  if (d === 0 || !Number.isFinite(d)) throw new Error('singular 3x3 matrix');
  const [[a, b, c], [e, f, g], [h, i, j]] = m;
  return [
    [(f * j - g * i) / d, (c * i - b * j) / d, (b * g - c * f) / d],
    [(g * h - e * j) / d, (a * j - c * h) / d, (c * e - a * g) / d],
    [(e * i - f * h) / d, (b * h - a * i) / d, (a * f - b * e) / d],
  ];
}

export const scaleMat3 = (m: Mat3, s: number): Mat3 => m.map((r) => scale(r, s)) as Mat3;
