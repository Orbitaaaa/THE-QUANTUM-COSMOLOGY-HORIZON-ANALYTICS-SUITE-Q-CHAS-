import { describe, expect, it } from 'vitest';
import { interiorVolume } from '../src/lib/physics/page';
import {
  ALPHA_SCALAR,
  T_PLANCK,
  alignL1,
  blackHole,
  emissionsForBound,
  fmtPow10,
  log10InteriorVolumeEndM3,
  log10OrNull,
  logspace,
  positiveOrNull,
  projectDampedSinusoid,
  windowCurve,
} from '../src/panels/horizon-calc';
import { L_PLANCK, YEAR } from '../src/lib/physics/constants';

// Numbers written by qchas (web/public/data/gw.json, GW150914 remnant, detector frame).
const GW150914 = {
  mf: 67.7675585301165,
  af: 0.664500584185361,
  S_BH_nats: 4.2085735808423613e80,
  T_H_K: 7.787928322261557e-10,
  r_plus_km: 174.85160831951555,
  t_evap_log10_yr: 72.25927696561605,
  page_time_log10_yr: 71.93713101626399,
  echo_bound_log10: -1.8277602828436463e80,
};

describe('Kerr horizon numbers (mirror of qchas.horizon.kerr.KerrHorizon)', () => {
  const bh = blackHole(GW150914.mf, GW150914.af);
  it('reproduces the pipeline remnant thermodynamics of GW150914', () => {
    expect(bh.entropyNats / GW150914.S_BH_nats - 1).toBeLessThan(1e-9);
    expect(Math.abs(bh.temperatureK / GW150914.T_H_K - 1)).toBeLessThan(1e-9);
    expect(Math.abs(bh.rPlusKm / GW150914.r_plus_km - 1)).toBeLessThan(1e-9);
    expect(bh.log10TEvapYr).toBeCloseTo(GW150914.t_evap_log10_yr, 9);
    expect(bh.log10PageTimeYr).toBeCloseTo(GW150914.page_time_log10_yr, 9);
    expect(Math.abs(bh.log10EchoBound / GW150914.echo_bound_log10 - 1)).toBeLessThan(1e-9);
  });

  it('one quarter nat per Planck area, bits = nats / ln 2, area quanta = A / a0', () => {
    expect(bh.entropyNats).toBeCloseTo(bh.areaLp2 / 4, -60);
    expect(bh.log10EntropyBits - bh.log10EntropyNats).toBeCloseTo(-Math.log10(Math.LN2), 12);
    expect(bh.areaM2 / (bh.areaLp2 * L_PLANCK ** 2)).toBeCloseTo(1, 12);
    expect(bh.areaQuanta * 5.965948868839189).toBeCloseTo(bh.areaLp2, -60);
  });

  it('Schwarzschild limit: S = 4 pi M^2, T = 1/(8 pi M), r_s = 2GM/c^2', () => {
    const s = blackHole(1, 0);
    expect(s.entropyNats / (4 * Math.PI * s.massPlanck ** 2)).toBeCloseTo(1, 12);
    expect(s.temperatureK).toBeCloseTo(6.17e-8, 9); // 6.17e-8 K for one solar mass
    expect(s.rPlusKm).toBeCloseTo(2.953, 3);
  });

  it('echo bound for the 62 Msun check of [QG] Sec. 18 (log10 |R|^2 = -1.75e80)', () => {
    const v = blackHole(62, 0).log10EchoBound;
    expect(Math.abs(v / -1.7511546455939697e80 - 1)).toBeLessThan(1e-9); // gw.json echo_bound_check
    expect((v / 1e80).toFixed(2)).toBe('-1.75'); // printed to three digits
  });

  it('lifetime t = M^3/(3 alpha) t_P, independent of spin', () => {
    const m = 1e20; // in solar masses: any mass, no overflow in log space
    const a = blackHole(m, 0);
    const b = blackHole(m, 0.9);
    expect(a.log10TEvapYr).toBe(b.log10TEvapYr);
    const direct = Math.log10((blackHole(1).massPlanck ** 3 / (3 * ALPHA_SCALAR)) * (T_PLANCK / YEAR));
    expect(blackHole(1).log10TEvapYr).toBeCloseTo(direct, 10);
  });
});

describe('[BH] Theorem F endpoint volume', () => {
  it('V_end = (9 sqrt3 pi / 5) M0^2 t_evap = 3 sqrt3 pi M0^5 / (5 alpha)', () => {
    const m = 1e10;
    const tEvapPlanck = m ** 3 / (3 * ALPHA_SCALAR);
    const v = interiorVolume(1) * m * m * tEvapPlanck * L_PLANCK ** 3;
    expect(log10InteriorVolumeEndM3(m)).toBeCloseTo(Math.log10(v), 10);
  });
});

describe('[BH] Theorem G window capacity', () => {
  it('matches the pipeline window examples at M* = 10 (horizon.json capacity.window_examples)', () => {
    const c = windowCurve(10, [0.01, 1, 100]);
    expect(c.lnTrace[0]).toBeCloseTo(-10.126922539848179, 8);
    expect(c.lnTrace[1]).toBeCloseTo(-4.985431665910364, 8);
    expect(c.lnTrace[2]).toBeCloseTo(94.47324347947672, 8);
    expect(c.theoremGv[2]).toBeCloseTo(96.46268026812534, 6);
    expect(c.deltaS[2]).toBeCloseTo(101.98943678864862, 6);
  });
  it('ln tr Pi -> beta DeltaE - ln beta for wide windows', () => {
    const m = 1e39;
    const c = windowCurve(m, [500]);
    expect(c.lnTrace[0]).toBeCloseTo(500 - Math.log(8 * Math.PI * m), 8);
    // the area term stays accurate for astrophysical M* (no (M + w)^2 - M^2 cancellation)
    expect(c.deltaS[0]).toBeCloseTo(500, 8);
    expect(c.theoremGv[0]).toBeCloseTo(500 - Math.log(8 * Math.PI * m), 8);
  });
  it('logspace spans its ends', () => {
    const g = logspace(0.01, 100, 5);
    expect(g[0]).toBeCloseTo(0.01, 14);
    expect(g[2]).toBeCloseTo(1, 14);
    expect(g[4]).toBeCloseTo(100, 10);
  });
});

describe('[QG] Theorem R helpers', () => {
  it('numerical zeros become null on log axes', () => {
    expect(log10OrNull(1e-3)).toBeCloseTo(-3, 12);
    expect(log10OrNull(5.55e-16)).toBeNull();
    expect(log10OrNull(-1.85e-16)).toBeNull();
    expect(positiveOrNull(0.5)).toBe(0.5);
    expect(positiveOrNull(null)).toBeNull();
  });
  it('emissions needed for 2^{1+k-s} <= eps', () => {
    expect(emissionsForBound(2, 0.0625)).toBe(7);
    expect(2 ** (1 + 2 - emissionsForBound(2, 0.01))).toBeLessThanOrEqual(0.01);
    expect(2 ** (1 + 2 - (emissionsForBound(2, 0.01) - 1))).toBeGreaterThan(0.01);
  });
});

describe('LIGO strain helpers', () => {
  it('alignL1 shifts by the lag and applies the sign', () => {
    // a pulse at 0.010 s in H1 reaches L1 7 ms earlier and inverted
    const t = Array.from({ length: 41 }, (_, i) => i * 0.001);
    const h1 = t.map((x) => Math.exp(-(((x - 0.02) / 0.002) ** 2)));
    const l1 = t.map((x) => -Math.exp(-(((x - 0.013) / 0.002) ** 2)));
    const a = alignL1(t, l1, -7, -1);
    const iPeak = a.y.indexOf(Math.max(...a.y));
    expect(a.t[iPeak]).toBeCloseTo(0.02, 9);
    expect(a.y[iPeak]).toBeCloseTo(h1[20]!, 12);
  });

  it('projectDampedSinusoid recovers amplitude and phase of a clean ringdown', () => {
    const fs = 2048;
    const t = Array.from({ length: 400 }, (_, i) => -0.05 + i / fs);
    const f = 250;
    const tau = 0.004;
    const t0 = 0.003;
    const y = t.map((x) => (x < t0 ? 0 : 1.7 * Math.exp(-(x - t0) / tau) * Math.cos(2 * Math.PI * f * (x - t0) + 0.6)));
    const p = projectDampedSinusoid(t, y, f, tau, t0, 0.024);
    expect(p).not.toBeNull();
    expect(p!.amplitude).toBeCloseTo(1.7, 9);
    expect(p!.phase).toBeCloseTo(0.6, 9);
    expect(p!.r2).toBeCloseTo(1, 9);
    expect(p!.t[0]).toBeGreaterThanOrEqual(t0);
  });

  it('returns null for an empty window', () => {
    expect(projectDampedSinusoid([0, 1], [0, 0], 100, 0.01, 5, 0.01)).toBeNull();
  });
});

describe('fmtPow10', () => {
  it('formats ordinary and astronomically large exponents', () => {
    expect(fmtPow10(72.25927696561605)).toBe('1.82×10⁷²');
    expect(fmtPow10(2)).toBe('1.00×10²');
    expect(fmtPow10(-1.8277602828436463e80)).toBe('10^(−1.83×10⁸⁰)');
    expect(fmtPow10(null)).toBe('—');
    expect(fmtPow10(0.99999999)).toBe('1.00×10¹');
  });
});
