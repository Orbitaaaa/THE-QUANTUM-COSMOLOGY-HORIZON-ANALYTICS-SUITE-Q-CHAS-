import { describe, expect, it } from 'vitest';
import { CHI_STAR_MPC } from '../src/lib/physics/constants';
import { circlePairs, ghostImages, lbToVec, matchedPoint } from '../src/lib/physics/circles';
import { cubicBasis } from '../src/lib/physics/lattice';
import { type Mat3, det3, dot, matmul3, norm, transpose3 } from '../src/lib/physics/vec';
import {
  type ManifoldPrinciples,
  type PrincipleKey,
  asMat3,
  bandKey,
  classicalScaleFactor,
  comovingDistanceMpc,
  eulerFromMatrix,
  eulerZYZ,
  excludedBy,
  failedPrinciples,
  failingExactly,
  flowSeries,
  haarRotation,
  hessianTraceCheck,
  imageThresholdChi,
  interp,
  minimalSelectingSets,
  redshiftForDistance,
  snapToCircle,
  survivors,
} from '../src/panels/m1-logic';

/** [TOP] Sec. 4 status table, grouped rows expanded (as qchas/topology/manifolds.py derives it). */
const A = 'allowed';
const X = 'excluded';
const N = 'no spin structure';
const TABLE: Record<string, [string, string, string, string]> = {
  E1: [A, A, A, A],
  E2: [A, A, X, X],
  E3: [A, A, X, X],
  E4: [A, A, X, X],
  E5: [A, A, X, X],
  E6: [A, A, X, X],
  E7: [X, A, X, N],
  E8: [X, A, X, N],
  E9: [X, A, X, N],
  E10: [X, A, X, N],
  E11: [A, X, A, A],
  E12: [A, X, X, X],
  E13: [X, X, X, N],
  E14: [X, X, X, N],
  E15: [X, X, X, N],
  E16: [A, X, A, A],
  E17: [X, X, X, N],
  E18: [A, X, A, A],
};
const MS: ManifoldPrinciples[] = Object.entries(TABLE).map(([label, s]) => ({
  label,
  principles: { P1: s[0], P2: s[1], P3: s[2], P4: s[3] },
}));
const E = (...ns: number[]) => ns.map((n) => `E${n}`);

/** fcc_preregistered_rotation()['R'] from qchas/topology/circles.py */
const R_PY: Mat3 = [
  [0.6563928569914044, 0.5892020253511204, -0.4711575008559236],
  [-0.6462630800555001, 0.11694117784179434, -0.7541013143352504],
  [-0.38922030862718304, 0.7994784138718889, 0.45753886950186934],
];

const maxAbsDiff = (a: Mat3, b: Mat3): number => Math.max(...a.flatMap((r, i) => r.map((v, j) => Math.abs(v - b[i]![j]!))));

describe('principle filter, [TOP] Main Theorem', () => {
  it('selects the 3-torus from {P2, P3}, {P2, P4} and all four principles', () => {
    expect(survivors(MS, ['P2', 'P3'])).toEqual(['E1']);
    expect(survivors(MS, ['P2', 'P4'])).toEqual(['E1']);
    expect(survivors(MS, ['P1', 'P2', 'P3', 'P4'])).toEqual(['E1']);
  });

  it('leaves the right families for single principles (Theorems A, C′, B-loc, D)', () => {
    expect(survivors(MS, [])).toHaveLength(18);
    expect(survivors(MS, ['P1'])).toEqual(E(1, 2, 3, 4, 5, 6, 11, 12, 16, 18));
    expect(survivors(MS, ['P2'])).toEqual(E(1, 2, 3, 4, 5, 6, 7, 8, 9, 10));
    expect(survivors(MS, ['P3'])).toEqual(E(1, 11, 16, 18));
    expect(survivors(MS, ['P4'])).toEqual(E(1, 11, 16, 18));
    expect(survivors(MS, ['P1', 'P3', 'P4'])).toEqual(E(1, 11, 16, 18));
    expect(() => survivors(MS, ['P5' as PrincipleKey])).toThrow();
  });

  it('shows the independence examples of [TOP] Sec. 3', () => {
    const by = (l: string) => MS.find((m) => m.label === l)!;
    // infinite space E18 satisfies P1, P3, P4 but not P2
    expect(failedPrinciples(by('E18'))).toEqual(['P2']);
    // the half-turn space E2 satisfies P1 and P2 but not P3 or P4
    expect(failedPrinciples(by('E2'))).toEqual(['P3', 'P4']);
    // the Klein space E7 satisfies P2 but not P1, P3 or P4
    expect(failedPrinciples(by('E7'))).toEqual(['P1', 'P3', 'P4']);
    // P2 alone and P3 alone each leave more than the torus: both are needed
    expect(survivors(MS, ['P2']).length).toBeGreaterThan(1);
    expect(survivors(MS, ['P3']).length).toBeGreaterThan(1);
    expect(excludedBy(by('E7'), ['P2', 'P3'])).toEqual(['P3']);
  });

  it('finds the inclusion-minimal selecting sets {P2, P3} and {P2, P4}', () => {
    expect(minimalSelectingSets(MS)).toEqual([
      ['P2', 'P3'],
      ['P2', 'P4'],
    ]);
  });

  it('flags the false sentence: E11, E16 and E18 fail only P2', () => {
    expect(failingExactly(MS, 1)).toEqual(['E11', 'E16', 'E18']);
    expect(failingExactly(MS, 0)).toEqual(['E1']);
  });
});

describe('rotations', () => {
  it('round-trips ZYZ Euler angles, including the pre-registered FCC orientation', () => {
    const [a, b, g] = eulerFromMatrix(R_PY);
    expect(maxAbsDiff(eulerZYZ(a, b, g), R_PY)).toBeLessThan(1e-12);
    for (const ang of [
      [10, 20, 30],
      [350, 179, 5],
      [123, 45, 270],
    ] as const) {
      const R = eulerZYZ(ang[0], ang[1], ang[2]);
      expect(maxAbsDiff(eulerZYZ(...eulerFromMatrix(R)), R)).toBeLessThan(1e-12);
    }
    // gimbal lock keeps the rotation
    for (const R of [eulerZYZ(40, 0, 25), eulerZYZ(40, 180, 25)]) expect(maxAbsDiff(eulerZYZ(...eulerFromMatrix(R)), R)).toBeLessThan(1e-9);
  });

  it('draws proper rotations (orthogonal, det +1)', () => {
    let seed = 7;
    const rand = () => ((seed = (seed * 16807) % 2147483647) - 1) / 2147483646;
    for (let k = 0; k < 20; k++) {
      const R = haarRotation(rand);
      expect(maxAbsDiff(matmul3(R, transpose3(R)), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])).toBeLessThan(1e-12);
      expect(det3(R)).toBeCloseTo(1, 12);
    }
    expect(asMat3(R_PY)).toEqual(R_PY);
    expect(asMat3([[1, 2], [3]])).toBeNull();
  });
});

describe('matched-circle snap and twin, Corollary F.1', () => {
  it('snaps a pick near a circle onto it and returns the exact twin n + g/chi*', () => {
    const pairs = circlePairs(cubicBasis('FCC'), 0.8, R_PY);
    for (const [i, p] of pairs.entries()) {
      for (const side of [1, -1] as const) {
        const c = side === 1 ? p.centre_minus : p.centre_plus;
        // a direction 0.6 deg off the circle, at an arbitrary azimuth
        const alpha = ((p.radius_deg + 0.6) * Math.PI) / 180;
        const e = norm([c[1], -c[0], 0]) > 1e-6 ? [c[1], -c[0], 0] : [1, 0, 0];
        const en = Math.hypot(e[0]!, e[1]!, e[2]!);
        const n: [number, number, number] = [
          Math.cos(alpha) * c[0] + (Math.sin(alpha) * e[0]!) / en,
          Math.cos(alpha) * c[1] + (Math.sin(alpha) * e[1]!) / en,
          Math.cos(alpha) * c[2] + (Math.sin(alpha) * e[2]!) / en,
        ];
        const hit = snapToCircle(n, pairs, 1.0);
        expect(hit).not.toBeNull();
        expect(hit!.index).toBe(i);
        expect(hit!.side).toBe(side);
        expect(hit!.offsetDeg).toBeCloseTo(0.6, 6);
        expect(norm(hit!.twin)).toBeCloseTo(1, 12);
        const twin = matchedPoint(hit!.point, side === 1 ? p.g_sky : [-p.g_sky[0], -p.g_sky[1], -p.g_sky[2]]);
        expect(Math.max(...twin.map((v, k) => Math.abs(v - hit!.twin[k]!)))).toBeLessThan(1e-12);
        // zero twist: the twin sits on the antipodal circle at the same angular radius
        const cTwin = side === 1 ? p.centre_plus : p.centre_minus;
        expect((Math.acos(dot(hit!.twin, cTwin)) * 180) / Math.PI).toBeCloseTo(p.radius_deg, 9);
      }
    }
  });

  it('returns null away from every circle and when there are no circles (d >= D)', () => {
    const pairs = circlePairs(cubicBasis('SC'), 0.985);
    // the centre of a 9.9 deg circle is 9.9 deg away from it
    expect(snapToCircle(pairs[0]!.centre_plus, pairs, 1.5)).toBeNull();
    expect(snapToCircle(lbToVec(10, 10), circlePairs(cubicBasis('FCC'), 1.0), 5)).toBeNull();
  });

  it('bounds where images can appear: no image inside chi* for sources closer than d - chi*', () => {
    expect(imageThresholdChi(0.985)).toBeCloseTo(0.97, 12);
    expect(imageThresholdChi(0.4)).toBe(0);
    for (const shape of ['SC', 'BCC', 'FCC']) {
      for (const d of [0.985, 1.04, 1.3]) {
        const rD = 0.5 * imageThresholdChi(d) * 0.999; // just inside the bound, in units of D
        for (const lb of [[0, 0], [123, -40], [280, 62], [33, 80]] as const) {
          const imgs = ghostImages(lb, rD, cubicBasis(shape), d, R_PY);
          expect(imgs.filter((g) => !g.is_source && g.inside_horizon)).toHaveLength(0);
        }
      }
    }
    // just beyond the bound, a source along -g_hat of a shortest vector has an image inside
    const pairs = circlePairs(cubicBasis('SC'), 0.985, R_PY);
    const imgs = ghostImages(pairs[0]!.centre_minus, 0.5 * imageThresholdChi(0.985) * 1.01, cubicBasis('SC'), 0.985, R_PY);
    expect(imgs.some((g) => !g.is_source && g.inside_horizon)).toBe(true);
  });
});

describe('comoving distance, flat LCDM with Planck 2018 parameters', () => {
  it('reaches chi* = 13 870 Mpc at z* within 0.1 % and puts z = 1 at about 0.24 chi*', () => {
    expect(Math.abs(comovingDistanceMpc(1089.9) / CHI_STAR_MPC.value - 1)).toBeLessThan(1e-3);
    expect(comovingDistanceMpc(1) / CHI_STAR_MPC.value).toBeCloseTo(0.245, 2);
    expect(comovingDistanceMpc(0)).toBe(0);
  });

  it('inverts to the redshift', () => {
    for (const z of [0.1, 1, 3, 300, 1089.9]) expect(redshiftForDistance(comovingDistanceMpc(z)) / z).toBeCloseTo(1, 4);
  });
});

describe('helpers', () => {
  it('interpolates linearly and refuses to extrapolate', () => {
    expect(interp([0, 1, 2], [0, 10, 40], 1.5)).toBe(25);
    expect(interp([0, 1, 2], [0, 10, 40], 2.5)).toBeNull();
  });

  it('checks the Hessian trace identity tr H = 2Z/3 (BCC typo in [TOP] Prop. E.2)', () => {
    const Z = 15.971880245588027;
    const ours = hessianTraceCheck(Z, [0.23139477861794489, 0.23139477861795604, 3.3950435354964883, 3.3950435354964896, 3.39504353549649]);
    expect(Math.abs(ours.residual)).toBeLessThan(1e-10);
    const printed = hessianTraceCheck(Z, [0.2314, 0.2314, 3.3953, 3.3953, 3.3953]);
    expect(printed.residual).toBeGreaterThan(7e-4);
  });

  it('builds gap-separated flow series for a log axis', () => {
    const s = flowSeries(
      [
        { tau: [0, 1, 2], Z: [17, 16, 15.96213148] },
        { tau: [0, 1], Z: [18, 15.98] },
      ],
      15.962131481337693,
    );
    expect(s.x).toEqual([0, 1, 2, 2, 0, 1]);
    expect(s.y[2]).toBeNull();
    expect(s.y[3]).toBeNull();
    expect(s.y[5]).toBeCloseTo(0.0178685, 6);
  });

  it('classical scale factor and SPHEREx band lookup', () => {
    expect(classicalScaleFactor(6, Math.PI / 2)).toBeCloseTo(0, 12);
    expect(classicalScaleFactor(6, Math.PI)).toBeCloseTo(6, 12);
    const bands = { '1': [0.75, 1.09, 39], '2': [1.1, 1.62, 41] };
    expect(bandKey([1.1, 1.62], bands)).toBe('2');
    expect(bandKey([3, 4], bands)).toBeNull();
  });
});
