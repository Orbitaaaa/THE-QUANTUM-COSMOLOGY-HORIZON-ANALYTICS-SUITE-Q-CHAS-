import { describe, expect, it } from 'vitest';
import { Z2_CERTIFIED } from '../src/lib/physics/constants';
import {
  covolume,
  cubicBasis,
  dualBasis,
  epsteinZ2,
  erfc,
  latticeVectors,
  lllReduce,
  normalize,
  phi2,
  rhombohedralBasis,
  shortestVector,
} from '../src/lib/physics/lattice';
import { type Mat3, matmul3 } from '../src/lib/physics/vec';

describe('erfc', () => {
  it('matches scipy.special.erfc', () => {
    const ref: [number, number][] = [
      [0, 1], [0.1, 0.8875370839817152], [0.5, 0.4795001221869535], [1, 0.15729920705028516],
      [1.5, 0.03389485352468927], [2, 0.004677734981047266], [2.4999, 0.0004071699003334512],
      [2.5, 0.00040695201744495886], [3, 2.2090496998585445e-5], [4, 1.541725790028002e-8],
      [5, 1.5374597944280347e-12], [7, 4.183825607779414e-23], [8.9, 2.5053574980338316e-36],
      [-0.7, 1.6778011938374182],
    ];
    for (const [x, v] of ref) {
      expect(Math.abs(erfc(x) - v)).toBeLessThan(5e-16 + 1e-12 * v);
    }
  });

  it('phi2 is the integral it claims to be', () => {
    // int_1^inf u^-3/2 e^{-pi y u} du by a crude but converged trapezoid in t = u^-1/2
    const y = 0.8;
    let s = 0;
    const n = 200_000;
    for (let i = 1; i <= n; i++) {
      const t = (i - 0.5) / n;
      s += 2 * Math.exp((-Math.PI * y) / (t * t));
    }
    expect(phi2(y)).toBeCloseTo(s / n, 9);
  });
});

describe('Epstein zeta Z(2) by Ewald summation, [TOP] Prop. E.2', () => {
  it('reproduces the certified SC, BCC and FCC values to 1e-9', () => {
    expect(Math.abs(epsteinZ2(cubicBasis('SC')) - 16.53231595976167)).toBeLessThan(1e-9);
    expect(Math.abs(epsteinZ2(cubicBasis('BCC')) - 15.971880245588027)).toBeLessThan(1e-9);
    expect(Math.abs(epsteinZ2(cubicBasis('FCC')) - 15.962131481337693)).toBeLessThan(1e-9);
    for (const s of ['SC', 'BCC', 'FCC'] as const) expect(Math.abs(epsteinZ2(cubicBasis(s)) - Z2_CERTIFIED.value[s])).toBeLessThan(1e-9);
  });

  it('agrees with the Python Ewald sum along the rhombohedral family and for a skew basis', () => {
    const py: [number, number][] = [
      [60, 15.962131481337684], [75, 16.26032369339946], [90, 16.532315959761664],
      [100, 16.345204107257103], [109.47122063449069, 15.971880245588027],
    ];
    for (const [a, z] of py) expect(Math.abs(epsteinZ2(rhombohedralBasis(a)) - z)).toBeLessThan(1e-10);
    const B: Mat3 = [[1.0, 0.2, 0.1], [0.3, 1.1, -0.2], [0.05, 0.4, 0.9]];
    expect(Math.abs(epsteinZ2(B) - 16.307429499784238)).toBeLessThan(1e-10);
  });

  it('is invariant under rescaling, rotation and change of basis; scales as V^-4/3', () => {
    const fcc = cubicBasis('FCC');
    const z = epsteinZ2(fcc);
    const c = Math.cos(0.7);
    const s = Math.sin(0.7);
    const R: Mat3 = [[c, -s, 0], [s, c, 0], [0, 0, 1]];
    const U: Mat3 = [[1, 1, 0], [0, 1, 0], [2, 3, 1]]; // unimodular
    expect(epsteinZ2(matmul3(fcc, R))).toBeCloseTo(z, 11);
    expect(epsteinZ2(matmul3(U, fcc))).toBeCloseTo(z, 11);
    expect(epsteinZ2(fcc.map((r) => r.map((x) => 3 * x)) as Mat3)).toBeCloseTo(z, 11);
    const V = covolume(fcc);
    expect(epsteinZ2(fcc, false)).toBeCloseTo(z * V ** (-4 / 3), 9);
  });

  it('rhombohedral basis hits FCC, SC and BCC at 60, 90 and 109.47 degrees', () => {
    expect(epsteinZ2(rhombohedralBasis(60))).toBeCloseTo(Z2_CERTIFIED.value.FCC, 9);
    expect(epsteinZ2(rhombohedralBasis(90))).toBeCloseTo(Z2_CERTIFIED.value.SC, 9);
    expect(epsteinZ2(rhombohedralBasis((Math.acos(-1 / 3) * 180) / Math.PI))).toBeCloseTo(Z2_CERTIFIED.value.BCC, 9);
    expect(() => rhombohedralBasis(120)).toThrow();
  });
});

describe('lattice geometry', () => {
  it('finds shells and shortest vectors', () => {
    expect(latticeVectors(cubicBasis('SC'), 1.0001)).toHaveLength(6);
    expect(latticeVectors(cubicBasis('FCC'), Math.SQRT1_2 * 1.0001)).toHaveLength(12);
    expect(latticeVectors(cubicBasis('BCC'), (Math.sqrt(3) / 2) * 1.0001)).toHaveLength(8);
    expect(shortestVector(cubicBasis('FCC')).length).toBeCloseTo(Math.SQRT1_2, 14);
    const skew: Mat3 = [[1, 0, 0], [7, 1, 0], [3, 5, 1]];
    expect(shortestVector(skew).length).toBeCloseTo(1, 14);
    expect(covolume(lllReduce(skew))).toBeCloseTo(1, 12);
    expect(shortestVector(normalize(cubicBasis('BCC'), 'shortest')).length).toBeCloseTo(1, 14);
    expect(covolume(normalize(cubicBasis('BCC')))).toBeCloseTo(1, 14);
    const D = dualBasis(cubicBasis('SC'));
    expect(D[0][0]).toBeCloseTo(2 * Math.PI, 14);
  });
});
