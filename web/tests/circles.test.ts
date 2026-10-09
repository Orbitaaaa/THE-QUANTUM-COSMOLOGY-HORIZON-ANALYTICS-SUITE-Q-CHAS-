import { describe, expect, it } from 'vitest';
import { FCC_PREREGISTERED_AXES_GAL } from '../src/lib/physics/constants';
import {
  angularSeparation,
  circlePairs,
  circlePoints,
  fccAxesLatticeFrame,
  ghostImages,
  lbToVec,
  matchedPoint,
  vecToLb,
} from '../src/lib/physics/circles';
import { cubicBasis } from '../src/lib/physics/lattice';
import { type Mat3, norm } from '../src/lib/physics/vec';

/** fcc_preregistered_rotation()['R'] from qchas/topology/circles.py */
const R_PY: Mat3 = [
  [0.6563928569914044, 0.5892020253511204, -0.4711575008559236],
  [-0.6462630800555001, 0.11694117784179434, -0.7541013143352504],
  [-0.38922030862718304, 0.7994784138718889, 0.45753886950186934],
];

describe('matched circles, [TOP] Sec. 7 and Corollary F.1', () => {
  it('reproduces the circle table for 0.985 D and Fig. 4 at 0.8 D', () => {
    const count = (s: string, d: number) => circlePairs(cubicBasis(s), d).length;
    expect([count('SC', 0.985), count('BCC', 0.985), count('FCC', 0.985)]).toEqual([3, 4, 6]);
    expect([count('SC', 0.8), count('BCC', 0.8), count('FCC', 0.8)]).toEqual([3, 7, 6]);
    expect(circlePairs(cubicBasis('FCC'), 0.985)[0]!.radius_deg).toBeCloseTo(9.936367, 5);
    expect(circlePairs(cubicBasis('SC'), 0.8)[0]!.radius_deg).toBeCloseTo(36.869898, 5);
    expect(circlePairs(cubicBasis('FCC'), 1.0)).toHaveLength(0);
    expect(count('BCC', 0.5)).toBe(25);
    expect(count('FCC', 0.5)).toBe(21);
  });

  it('places FCC circle centres like the Python for the [TOP] Sec. 9 orientation', () => {
    const py: [number[], [number, number]][] = [
      [[0, 0, 1], [336.976687, 16.863721]],
      [[0, 1, -1], [219.401722, -13.99208]],
      [[0, 1, 0], [277.535143, 2.768952]],
      [[1, -1, 0], [95.031223, 57.196861]],
      [[1, 0, -1], [185.463117, 36.780365]],
      [[1, 0, 0], [280.49599, 62.728523]],
    ];
    const pairs = circlePairs(cubicBasis('FCC'), 0.8, R_PY);
    expect(pairs.map((p) => p.n)).toEqual(py.map((x) => x[0]));
    pairs.forEach((p, i) => {
      expect(p.centre_plus_lb[0]).toBeCloseTo(py[i]![1][0], 5);
      expect(p.centre_plus_lb[1]).toBeCloseTo(py[i]![1][1], 5);
      expect(p.radius_deg).toBeCloseTo(36.869898, 5);
    });
    // the rotated <110> axes land on the quoted best-fit axes (rms 0.03 deg in the Python fit)
    for (const lb of FCC_PREREGISTERED_AXES_GAL.value) {
      const best = Math.min(...pairs.map((p) => Math.min(angularSeparation(lb, p.centre_plus), angularSeparation(lb, p.centre_minus))));
      expect(best).toBeLessThan(0.1);
    }
  });

  it('twin points n + g/chi_* are unit vectors on the partner circle', () => {
    for (const p of circlePairs(cubicBasis('BCC'), 0.7)) {
      const { circle, twin } = circlePoints(p.g_sky, 24);
      circle.forEach((n, i) => {
        const m = matchedPoint(n, p.g_sky);
        expect(norm(m)).toBeCloseTo(1, 12);
        expect(m[0]).toBeCloseTo(twin[i]![0], 12);
        expect(angularSeparation(n, p.centre_minus)).toBeCloseTo(p.radius_deg, 9);
        expect(angularSeparation(m, p.centre_plus)).toBeCloseTo(p.radius_deg, 9);
      });
    }
    expect(() => matchedPoint([0, 0], [0.3, 0, 0])).toThrow();
  });

  it('ghost images of a source match the Python', () => {
    const g = ghostImages([30, 20], 0.3, cubicBasis('SC'), 0.7);
    expect(g).toHaveLength(13);
    const py: [number[], number, number, number][] = [
      [[0, 0, 0], 30.0, 20.0, 0.3],
      [[-1, 0, 0], 162.818189, 12.135912, 0.488062],
      [[0, -1, 0], 293.591227, 9.547689, 0.618599],
      [[0, 0, -1], 30.0, -64.737565, 0.660569],
      [[-1, -1, 0], 230.805305, 8.095575, 0.728608],
      [[-1, 0, -1], 162.818189, -51.38467, 0.764563],
    ];
    py.forEach(([n, l, b, d], i) => {
      expect(g[i]!.n).toEqual(n);
      expect(g[i]!.l).toBeCloseTo(l, 5);
      expect(g[i]!.b).toBeCloseTo(b, 5);
      expect(g[i]!.distance_over_D).toBeCloseTo(d, 6);
    });
    expect(g[0]!.is_source).toBe(true);
    expect(g[1]!.inside_horizon).toBe(true);
    expect(g[1]!.distance_Gpc).toBeCloseTo(0.488062 * 27.74, 4);
  });

  it('coordinate helpers and the FCC axes', () => {
    const [l, b] = vecToLb(lbToVec(359.5, -12));
    expect(l).toBeCloseTo(359.5, 12);
    expect(b).toBeCloseTo(-12, 12);
    expect(fccAxesLatticeFrame()).toHaveLength(6);
  });
});
