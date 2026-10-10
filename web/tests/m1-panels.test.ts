/**
 * Statements the Lattice, Ghost-Twin and Quantum Big Bang panels compute live and print as
 * facts: no circles once d >= D, images of sources inside the sphere lie beyond it, the
 * nearest-image bound is attained, the BCC trace-identity discrepancy, the unit-covolume
 * rhombohedral family, and the classical a_max |cos T| reference.
 */

import { describe, expect, it } from 'vitest';
import { CHI_STAR, circlePairs, ghostImages, lbToVec } from '../src/lib/physics/circles';
import { CHI_STAR_MPC, HESSIAN_CERTIFIED, Z2_CERTIFIED } from '../src/lib/physics/constants';
import { cubicBasis, epsteinZ2, rhombohedralBasis, shortestVector } from '../src/lib/physics/lattice';
import { type Mat3, type Vec3, det3, matVec, scale, unit } from '../src/lib/physics/vec';
import { classicalScaleFactor, comovingDistanceMpc, haarRotation, hessianTraceCheck, imageThresholdChi } from '../src/panels/m1-logic';

/** Small deterministic PRNG (mulberry32) so the random orientations are reproducible. */
function rng(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const SHAPES = ['SC', 'BCC', 'FCC'] as const;

describe('ghost-twin live facts ([TOP] Theorem F, Corollary F.1)', () => {
  it('has no matched circles once d >= D, for every lattice and orientation', () => {
    const rand = rng(7);
    for (const s of SHAPES)
      for (const d of [1.0, 1.001, 1.05, 1.3, 1.5])
        for (let k = 0; k < 4; k++) expect(circlePairs(cubicBasis(s), d, haarRotation(rand))).toHaveLength(0);
  });

  it('just below D there are circles, all smaller than 9.94 degrees at the circle-search bound', () => {
    for (const s of SHAPES) {
      const pairs = circlePairs(cubicBasis(s), 0.985);
      expect(pairs.length).toBeGreaterThan(0);
      for (const p of pairs) expect(p.radius_deg).toBeLessThan(9.94);
    }
  });

  it('puts every image of a source inside the sphere beyond last scattering when d >= D', () => {
    const rand = rng(11);
    expect(imageThresholdChi(1)).toBe(1);
    for (const s of SHAPES)
      for (const d of [1.0, 1.04, 1.2]) {
        const R = haarRotation(rand);
        for (let k = 0; k < 6; k++) {
          const dir: Vec3 = lbToVec(360 * rand(), Math.asin(2 * rand() - 1) * (180 / Math.PI));
          const r = CHI_STAR * (0.05 + 0.95 * rand());
          const imgs = ghostImages(dir, r, cubicBasis(s), d, R).filter((g) => !g.is_source);
          for (const g of imgs) {
            expect(g.inside_horizon).toBe(false);
            expect(g.distance_over_D).toBeGreaterThanOrEqual(d - r - 1e-12);
          }
        }
      }
  });

  it('attains the nearest-image bound d - r for a source opposite the shortest loop', () => {
    for (const s of SHAPES) {
      const d = 1.05;
      const sv = shortestVector(cubicBasis(s));
      const g = unit(sv.g);
      const r = 0.3 * CHI_STAR;
      const imgs = ghostImages(scale(g, -1), r, cubicBasis(s), d).filter((x) => !x.is_source);
      expect(imgs[0]!.distance_over_D).toBeCloseTo(d - r, 10);
      // in units of chi*: (2 d/D - r/chi*), what the key-fact box prints
      expect(imgs[0]!.distance_over_D / CHI_STAR).toBeCloseTo(2 * d - r / CHI_STAR, 10);
    }
  });

  it('a z = 1 galaxy sits at about a quarter of chi*', () => {
    const f = comovingDistanceMpc(1) / CHI_STAR_MPC.value;
    expect(f).toBeGreaterThan(0.23);
    expect(f).toBeLessThan(0.26);
  });
});

describe('lattice panel ([TOP] Prop. E.2)', () => {
  it('printed BCC triple eigenvalue violates tr H = 2Z/3 beyond rounding; SC and FCC do not', () => {
    const tol = 5 * 0.5e-4; // five eigenvalues printed to 4 decimals
    const res = (s: 'SC' | 'BCC' | 'FCC') => Math.abs(hessianTraceCheck(Z2_CERTIFIED.value[s], HESSIAN_CERTIFIED.value[s]).residual);
    expect(res('SC')).toBeLessThan(tol);
    expect(res('FCC')).toBeLessThan(tol);
    expect(res('BCC')).toBeGreaterThan(tol);
    // our analytic value (lattice.json critical.BCC.hessian) satisfies it to rounding
    const ours = [0.23139477861794489, 0.23139477861795604, 3.3950435354964883, 3.3950435354964896, 3.39504353549649];
    expect(Math.abs(hessianTraceCheck(Z2_CERTIFIED.value.BCC, ours).residual)).toBeLessThan(1e-12);
  });

  it('the unit-covolume rhombohedral family reproduces the certified Z(2) at 60, 90 and 109.47 degrees', () => {
    const unitCov = (a: number): Mat3 => {
      const B = rhombohedralBasis(a);
      const s = Math.abs(det3(B)) ** (-1 / 3);
      return B.map((r) => scale(r, s)) as Mat3;
    };
    const bcc = (Math.acos(-1 / 3) * 180) / Math.PI;
    expect(Math.abs(det3(unitCov(73)))).toBeCloseTo(1, 12);
    expect(epsteinZ2(unitCov(60))).toBeCloseTo(Z2_CERTIFIED.value.FCC, 9);
    expect(epsteinZ2(unitCov(90))).toBeCloseTo(Z2_CERTIFIED.value.SC, 9);
    expect(epsteinZ2(unitCov(bcc))).toBeCloseTo(Z2_CERTIFIED.value.BCC, 9);
    // a rotation does not change Z(2)
    const R = haarRotation(rng(3));
    const B = unitCov(75);
    expect(epsteinZ2(B.map((r) => matVec(R, r)) as Mat3)).toBeCloseTo(epsteinZ2(B), 10);
  });
});

describe('bounce panel ([QG] Fig. 1)', () => {
  it('the classical universe reaches a = 0 at T = pi/2', () => {
    expect(classicalScaleFactor(2, 0)).toBe(2);
    expect(classicalScaleFactor(2, Math.PI / 2)).toBeLessThan(1e-12);
    expect(classicalScaleFactor(2, Math.PI)).toBeCloseTo(2, 12);
  });
});
