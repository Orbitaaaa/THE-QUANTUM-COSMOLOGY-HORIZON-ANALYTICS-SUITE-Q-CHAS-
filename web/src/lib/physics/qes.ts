/**
 * Quantum extremal surface of an evaporating 4D Schwarzschild black hole, [BH] Sec. 9:
 * a port of ``qchas/horizon/qes.py`` (same generalized entropy, same log-space gradient).
 *
 * Matter is an s-wave CFT of central charge c in the in-vacuum. The island is I = [0, a] and the
 * radiation R = [b, inf) is collected at radius r_b and retarded time u_b. With
 * e^{2 rho} = (r_s/r) e^{-r/r_s} e^{kappa v},
 *
 *     S_gen(a) = pi r_a^2 / G + (c/6) ln[(U_a - U_b)(v_b - v_a) e^{rho_a + rho_b} / eps^2].
 *
 * The island endpoint is parametrized by (ln delta, v) with r_a = r_s (1 - delta), and the
 * gradient is driven to zero by a 2D Levenberg-Marquardt iteration with the analytic Jacobian.
 * Planck units (G = hbar = 1), cutoff eps = l_P, r_s = 2M.
 */

import { QES_DEPTH_COEFF, QES_LAG_OFFSET_RS } from './constants';

export interface QesResult {
  mass: number;
  c: number;
  rb_over_rs: number;
  s_bh: number;
  /** (r_s - r_QES)/r_s */
  depth: number;
  /** u_b - v_QES */
  lag: number;
  /** (1/2 pi T) ln(S_BH / c) */
  t_scr: number;
  s_gen_minus_s_bh: number;
  converged: boolean;
  iterations: number;
  residual: number;
}

const logAddExp = (a: number, b: number): number => {
  const m = Math.max(a, b);
  return m + Math.log(Math.exp(a - m) + Math.exp(b - m));
};

export function solveQes(mass: number, c = 1, rbOverRs = 20, uB = 0): QesResult {
  const rs = 2 * mass;
  const kappa = 1 / (2 * rs);
  const sBh = Math.PI * rs * rs;
  const rb = rbOverRs * rs;
  const rstarB = rb + rs * Math.log(rb / rs - 1);
  const vB = uB + 2 * rstarB;
  const rhoB = 0.5 * (Math.log(rs / rb) - rb / rs + kappa * vB);
  const lnAbsUb = -kappa * uB - Math.log(kappa);

  const lnU = (delta: number, v: number): number => Math.log(delta) + 1 - delta - kappa * v - Math.log(kappa);

  const sGen = (delta: number, v: number): number => {
    const ra = rs * (1 - delta);
    const rhoA = 0.5 * (Math.log(rs / ra) - ra / rs + kappa * v);
    const lnDU = logAddExp(lnU(delta, v), lnAbsUb);
    return Math.PI * ra * ra + (c / 6) * (lnDU + Math.log(vB - v) + rhoA + rhoB);
  };

  // Residual (delta dS/d delta, r_s dS/dv) * 6/c and its Jacobian in (ln delta, v/r_s).
  const system = (logDelta: number, q: number) => {
    const delta = Math.exp(Math.min(logDelta, 0));
    const v = Math.min(q * rs, vB - 1e-9 * rs);
    const w = 1 / (1 + Math.exp(Math.min(700, lnAbsUb - lnU(delta, v))));
    const ww = w * (1 - w);
    const a = (12 * Math.PI * rs * rs) / c;
    const gap = vB - v;
    const gDelta = -a * delta * (1 - delta) + (1 - delta) * w + 0.5 * delta * (1 / (1 - delta) + 1);
    const gV = rs * (-kappa * w - 1 / gap + 0.5 * kappa);
    const dGdL = delta * (-a * (1 - 2 * delta) - w + 0.5 / (1 - delta) ** 2 + 0.5) + (1 - delta) ** 2 * ww;
    const dGdQ = -(1 - delta) * ww * kappa * rs;
    const dVdL = -rs * kappa * ww * (1 - delta);
    const dVdQ = rs * rs * (kappa * kappa * ww - 1 / (gap * gap));
    return { g: [gDelta, gV] as const, J: [dGdL, dGdQ, dVdL, dVdQ] as const };
  };

  const tScr = (1 / kappa) * Math.log(sBh / c);
  const guess: [number, number] = [Math.log((QES_DEPTH_COEFF.value * c) / sBh), (uB - (tScr + QES_LAG_OFFSET_RS.value * rs)) / rs];

  let best = { p: guess, res: Infinity, it: 0 };
  for (const dq of [0, -5, 5, -20, 20]) {
    const run = levenbergMarquardt(system, [guess[0], guess[1] + dq]);
    if (run.res < best.res) best = run;
    if (best.res < 1e-12) break;
  }
  const delta = Math.exp(best.p[0]);
  const vA = best.p[1] * rs;
  return {
    mass,
    c,
    rb_over_rs: rbOverRs,
    s_bh: sBh,
    depth: delta,
    lag: uB - vA,
    t_scr: tScr,
    s_gen_minus_s_bh: sGen(delta, vA) - sBh,
    converged: best.res < 1e-9,
    iterations: best.it,
    residual: best.res,
  };
}

type System = (x: number, y: number) => {
  g: readonly [number, number];
  J: readonly [number, number, number, number];
};

function levenbergMarquardt(f: System, p0: [number, number], maxIter = 200) {
  let p: [number, number] = [...p0];
  let cur = f(p[0], p[1]);
  let res = Math.hypot(cur.g[0], cur.g[1]);
  let lambda = 1e-6;
  let it = 0;
  for (; it < maxIter && res > 1e-14; it++) {
    const [a, b, cc, d] = cur.J;
    const [g0, g1] = cur.g;
    // normal equations (J^T J + lambda diag(J^T J)) dp = -J^T g
    const m00 = a * a + cc * cc;
    const m01 = a * b + cc * d;
    const m11 = b * b + d * d;
    const r0 = -(a * g0 + cc * g1);
    const r1 = -(b * g0 + d * g1);
    let accepted = false;
    for (let k = 0; k < 40 && !accepted; k++) {
      const A00 = m00 * (1 + lambda);
      const A11 = m11 * (1 + lambda);
      const det = A00 * A11 - m01 * m01;
      if (!(Math.abs(det) > 0)) {
        lambda *= 10;
        continue;
      }
      const dp0 = (r0 * A11 - m01 * r1) / det;
      const dp1 = (A00 * r1 - m01 * r0) / det;
      const trial: [number, number] = [p[0] + dp0, p[1] + dp1];
      const next = f(trial[0], trial[1]);
      const nres = Math.hypot(next.g[0], next.g[1]);
      if (Number.isFinite(nres) && nres < res) {
        p = trial;
        cur = next;
        res = nres;
        lambda = Math.max(lambda / 10, 1e-12);
        accepted = true;
        if (Math.abs(dp0) + Math.abs(dp1) < 1e-15 * (1 + Math.abs(p[0]) + Math.abs(p[1]))) it = maxIter;
      } else {
        lambda *= 10;
      }
    }
    if (!accepted) break;
  }
  return { p, res, it };
}

/** [BH] Sec. 9 law: (r_s - r_QES)/r_s = 0.0399 c / S_BH. */
export const depthFit = (c: number, sBh: number): number => (QES_DEPTH_COEFF.value * c) / sBh;

/** [BH] Sec. 9 law: u_b - v_QES = (1/2 pi T) ln(S_BH/c) + 4.24 r_s. */
export function lagFit(mass: number, c: number): number {
  const rs = 2 * mass;
  return 2 * rs * Math.log((4 * Math.PI * mass * mass) / c) + QES_LAG_OFFSET_RS.value * rs;
}

export interface QesPaperRow {
  mass: number;
  c: number;
  rb_over_rs: number;
  depth: number;
  lag: number;
  t_scr: number;
  s_gen_minus_s_bh: number;
}

/** The table printed in [BH] Sec. 9 (same rows as ``qes.PAPER_TABLE``). */
export const PAPER_TABLE: { rows: QesPaperRow[]; source: string } = {
  source: '[BH] Sec. 9',
  rows: (
    [
      [10, 1, 20, 3.11e-5, 370.16, 285.45, 0.98],
      [10, 1, 60, 3.23e-5, 371.49, 285.45, 1.12],
      [10, 10, 20, 3.1e-4, 277.86, 193.34, 11.56],
      [10, 10, 60, 3.22e-4, 279.34, 193.34, 13.03],
      [30, 1, 20, 3.47e-6, 1374.59, 1120.01, 1.17],
      [30, 1, 60, 3.59e-6, 1378.24, 1120.01, 1.3],
      [30, 10, 20, 3.45e-5, 1097.81, 843.7, 13.51],
      [30, 10, 60, 3.58e-5, 1101.83, 843.7, 14.92],
      [100, 1, 20, 3.14e-7, 5546.6, 4696.55, 1.38],
      [100, 1, 60, 3.23e-7, 5557.64, 4696.55, 1.51],
      [100, 10, 20, 3.12e-6, 4624.19, 3775.51, 15.64],
      [100, 10, 60, 3.23e-6, 4636.3, 3775.51, 16.98],
    ] as const
  ).map(([mass, c, rb_over_rs, depth, lag, t_scr, s_gen_minus_s_bh]) => ({
    mass,
    c,
    rb_over_rs,
    depth,
    lag,
    t_scr,
    s_gen_minus_s_bh,
  })),
};
