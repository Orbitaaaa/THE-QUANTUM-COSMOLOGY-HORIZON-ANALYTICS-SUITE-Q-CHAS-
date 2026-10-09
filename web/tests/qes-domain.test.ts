import { describe, expect, it } from 'vitest';
import { depthFit, lagFit, solveQes } from '../src/lib/physics/qes';

// The QES panel exposes M in [10, 1000] (Planck units), c in [1, 30] and r_b/r_s in [20, 60].
describe('QES solver over the panel slider domain', () => {
  it('converges at the corners and centre, close to the [BH] Sec. 9 fit laws, in a few ms', () => {
    const t0 = performance.now();
    let n = 0;
    for (const m of [10, 100, 1000])
      for (const c of [1, 10, 30])
        for (const rb of [20, 40, 60]) {
          const r = solveQes(m, c, rb);
          n++;
          expect(r.converged, `M=${m} c=${c} rb=${rb}`).toBe(true);
          expect(r.depth).toBeGreaterThan(0);
          expect(r.depth).toBeLessThan(0.01);
          // the fit laws hold to a few per cent for S_BH >> c (and 4.24 r_s is quoted at r_b = 20 r_s)
          if (4 * Math.PI * m * m > 100 * c) {
            expect(Math.abs(r.depth / depthFit(c, r.s_bh) - 1)).toBeLessThan(0.08);
            expect(Math.abs(r.lag - lagFit(m, c)) / (2 * m)).toBeLessThan(1.5);
          }
        }
    expect((performance.now() - t0) / n).toBeLessThan(20);
  });
});
