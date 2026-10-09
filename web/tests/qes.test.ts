import { describe, expect, it } from 'vitest';
import { PAPER_TABLE, depthFit, lagFit, solveQes } from '../src/lib/physics/qes';

/** solve_qes() of qchas/horizon/qes.py (scipy hybr) for every PAPER_TABLE row. */
const PYTHON: [number, number, number, number, number, number, number][] = [
  [10, 1, 20, 3.109978091791561e-5, 370.1580370970976, 285.4477773182953, 0.9768338425562888],
  [10, 1, 60, 3.2255518016341006e-5, 371.4927261207777, 285.4477773182953, 1.1168384615837113],
  [10, 10, 20, 0.0003095927245213878, 277.86170576139676, 193.34437359853348, 11.563842576010757],
  [10, 10, 60, 0.0003224546884293157, 279.34263730402085, 193.34437359853348, 13.034423897423494],
  [30, 1, 20, 3.4699696862965523e-6, 1374.5946918291634, 1120.0102812352322, 1.1709338374730578],
  [30, 1, 60, 3.586656674048383e-6, 1378.2414305942777, 1120.0102812352322, 1.3048551253632468],
  [30, 10, 20, 3.4548039435341164e-5, 1097.807298538388, 843.7000700759467, 13.512721551589493],
  [30, 10, 60, 3.583828348578498e-5, 1101.8298804258886, 843.7000700759467, 14.915835279238308],
  [100, 1, 20, 3.13558105693376e-7, 5546.604412669712, 4696.54584757819, 1.3828713838593103],
  [100, 1, 60, 3.230680521771845e-7, 5557.636481832695, 4696.54584757819, 1.5107450817304198],
  [100, 10, 20, 3.1235600085284014e-6, 4624.193892312529, 3775.511810380571, 15.639877135778079],
  [100, 10, 60, 3.2281108691189603e-6, 4636.296826276461, 3775.511810380571, 16.976315337102278],
];

/** Number of decimals printed in [BH] Sec. 9 for each column. */
const rel = (a: number, b: number): number => Math.abs(a - b) / Math.abs(b);
const roundSig = (x: number, sig: number): number => Number(x.toPrecision(sig));

describe('QES solver, [BH] Sec. 9', () => {
  it('converges for every row of the paper table', () => {
    for (const r of PAPER_TABLE.rows) {
      const q = solveQes(r.mass, r.c, r.rb_over_rs);
      expect(q.converged).toBe(true);
      expect(q.residual).toBeLessThan(1e-9);
    }
  });

  it('agrees with the Python solver (scipy hybr) to 1e-9 relative', () => {
    for (const [M, c, rb, depth, lag, tscr, sgen] of PYTHON) {
      const q = solveQes(M, c, rb);
      expect(rel(q.depth, depth)).toBeLessThan(1e-9);
      expect(rel(q.lag, lag)).toBeLessThan(1e-9);
      expect(rel(q.t_scr, tscr)).toBeLessThan(1e-12);
      expect(rel(q.s_gen_minus_s_bh, sgen)).toBeLessThan(1e-8);
    }
  });

  it('reproduces all 12 rows of the printed table: lag and t_scr to 1e-3 relative, every column to the printed digit', () => {
    for (const r of PAPER_TABLE.rows) {
      const q = solveQes(r.mass, r.c, r.rb_over_rs);
      expect(rel(q.lag, r.lag)).toBeLessThan(1e-3);
      expect(rel(q.t_scr, r.t_scr)).toBeLessThan(1e-3);
      // depth printed to 3 significant figures, lag and t_scr to 2 decimals, S_gen - S_BH to 2 decimals
      expect(roundSig(q.depth, 3)).toBeCloseTo(r.depth, 15);
      expect(Math.abs(q.lag - r.lag)).toBeLessThanOrEqual(0.005 + 1e-9);
      expect(Math.abs(q.t_scr - r.t_scr)).toBeLessThanOrEqual(0.005 + 1e-9);
      expect(Math.abs(q.s_gen_minus_s_bh - r.s_gen_minus_s_bh)).toBeLessThanOrEqual(0.005 + 1e-9);
      // and to 1e-3 relative wherever three printed significant figures allow it
      expect(rel(q.depth, r.depth)).toBeLessThan(1.6e-3);
    }
  });

  it('obeys the fit laws of the paper within their stated accuracy', () => {
    for (const r of PAPER_TABLE.rows) {
      const q = solveQes(r.mass, r.c, r.rb_over_rs);
      expect(rel(q.depth, depthFit(r.c, q.s_bh))).toBeLessThan(0.25);
      expect(rel(q.lag, lagFit(r.mass, r.c))).toBeLessThan(0.02);
    }
  });
});
