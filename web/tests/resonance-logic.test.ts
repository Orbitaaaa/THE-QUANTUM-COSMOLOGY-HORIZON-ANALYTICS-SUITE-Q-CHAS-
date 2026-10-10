import { describe, expect, it } from 'vitest';
import {
  axisAngleDeg,
  azEl,
  directionFamily,
  familyCounts,
  histogram,
  isotropicAlignmentP,
  median,
  pValueGE,
  shellKChi,
  shellVectors,
} from '../src/panels/resonance-logic';

/** The built bundle, when present (vite resolves the glob at transform time; {} without a build). */
const BUNDLES = import.meta.glob<{ default: unknown }>('../public/data/resonance.json', { eager: true });
const RES = (Object.values(BUNDLES)[0]?.default ?? null) as {
  shells: Record<string, { m2: number; k_chi: number; n_vectors: number }[]>;
  search: { S_max: number; S_avg: number; p_max: number; p_avg: number; null: { S_max: number[]; S_avg: number[] } } | null;
} | null;

describe('dual-lattice shells ([TOP] Sec. 11.6)', () => {
  it('cubic shells: 6 cube axes, 30 vectors on |m|^2 = 9, 6 on |m|^2 = 16', () => {
    expect(shellVectors('cubic', 1)).toHaveLength(6);
    expect(shellVectors('cubic', 9)).toHaveLength(30);
    expect(shellVectors('cubic', 16)).toHaveLength(6);
    expect(familyCounts(shellVectors('cubic', 9))).toEqual([
      ['100', 6],
      ['221', 24],
    ]);
  });

  it('FCC dual is BCC-like (8 <111> first), BCC dual is FCC-like (12 <110> first)', () => {
    const f = shellVectors('FCC', 1.5);
    expect(f).toHaveLength(8);
    expect(familyCounts(f)).toEqual([['111', 8]]);
    const b = shellVectors('BCC', 1.5);
    expect(b).toHaveLength(12);
    expect(familyCounts(b)).toEqual([['110', 12]]);
  });

  it('k chi* = pi |m| / (d/D) matches the bundle normalization', () => {
    expect(shellKChi(1, 1.2)).toBeCloseTo(2.6179938779914944, 12);
    expect(shellKChi(1.5, 1.05)).toBeCloseTo(3.664428086176755, 12);
  });

  it('agrees with every shell of resonance.json when the bundle is built', () => {
    if (!RES) return;
    for (const [key, shells] of Object.entries(RES.shells)) {
      const [shape, d] = key.split('_') as [string, string];
      for (const s of shells) {
        expect(shellVectors(shape, s.m2), `${key} m2=${s.m2}`).toHaveLength(s.n_vectors);
        expect(shellKChi(s.m2, Number(d))).toBeCloseTo(s.k_chi, 9);
      }
    }
  });
});

describe('direction helpers', () => {
  it('direction families', () => {
    expect(directionFamily([3, 0, 0])).toBe('100');
    expect(directionFamily([2, -2, 1])).toBe('221');
    expect(directionFamily([0.5, 0.5, -0.5])).toBe('111');
    expect(directionFamily([0, -1, 1])).toBe('110');
  });

  it('azimuth / elevation', () => {
    expect(azEl([1, 0, 0])).toEqual([0, 0]);
    expect(azEl([0, 0, 2])[1]).toBeCloseTo(90, 12);
    const [az, el] = azEl([1, 1, 1]);
    expect(az).toBeCloseTo(45, 12);
    expect(el).toBeCloseTo(35.26438968, 6);
    expect(azEl([-1, 0, 0])[0]).toBe(180);
  });
});

describe('null-distribution statistics', () => {
  it('histogram counts every in-range value once', () => {
    const v = [0, 0.1, 0.5, 0.9, 1];
    const hst = histogram(v, 2, 0, 1);
    expect(hst.counts).toEqual([2, 3]);
    expect(hst.x).toEqual([0.25, 0.75]);
  });

  it('p-value and median', () => {
    expect(pValueGE([1, 2, 3, 4], 3)).toBe(0.5);
    expect(median([3, 1, 2])).toBe(2);
    expect(median([4, 1, 2, 3])).toBe(2.5);
  });

  it('reproduces the search p-values of resonance.json', () => {
    const s = RES?.search;
    if (!s) return;
    expect(pValueGE(s.null.S_max, s.S_max)).toBeCloseTo(s.p_max, 9);
    expect(pValueGE(s.null.S_avg, s.S_avg)).toBeCloseTo(s.p_avg, 9);
  });
});

describe('low-l axis alignment ([TOP] Sec. 9.1)', () => {
  it("the paper's Planck axes are 14 deg apart with isotropic probability 0.0297", () => {
    const a = axisAngleDeg([242.8, 49.0], [239.5, 62.9]);
    expect(a).toBeCloseTo(14.0, 0);
    expect(isotropicAlignmentP(a)).toBeCloseTo(0.0297, 3);
  });

  it('axes are headless', () => {
    expect(axisAngleDeg([10, 20], [190, -20])).toBeCloseTo(0, 6);
    expect(axisAngleDeg([0, 0], [90, 0])).toBeCloseTo(90, 9);
  });
});
