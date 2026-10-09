import { describe, expect, it } from 'vitest';
import { classifyStatus, fmtBytes, fmtNum, fmtSci, fmtUtc } from '../src/lib/ui';
import { DEFAULT_ROUTE, PANEL_LOADERS, ROUTES, parseHash } from '../src/routes';
import { equirectTexture, greatArcPoints, makeColorizer, pickValue, polygonPoints, robustLimit, smallCirclePoints } from '../src/lib/sky';
import { lbToVec } from '../src/lib/healpix';

const ARCHITECTURE_ROUTES = [
  'overview', 'manifolds', 'lattice', 'ghost-twin', 'resonance', 'sky-stats', 'bounce',
  'capacity', 'page-curve', 'qes', 'hayden-preskill', 'gw-horizons', 'claims', 'data',
];

describe('routes', () => {
  it('registers exactly the routes of docs/ARCHITECTURE.md, each with a panel loader', () => {
    expect(ROUTES.map((r) => r.id)).toEqual(ARCHITECTURE_ROUTES);
    for (const id of ARCHITECTURE_ROUTES) expect(typeof PANEL_LOADERS[id as keyof typeof PANEL_LOADERS]).toBe('function');
  });

  it('parses hashes', () => {
    expect(parseHash('#/qes')).toEqual({ id: 'qes', rest: [] });
    expect(parseHash('#/data/warnings')).toEqual({ id: 'data', rest: ['warnings'] });
    expect(parseHash('#nope')).toEqual({ id: DEFAULT_ROUTE, rest: [] });
    expect(parseHash('')).toEqual({ id: DEFAULT_ROUTE, rest: [] });
  });
});

describe('status classification and formatting', () => {
  it('maps the papers status strings to provenance kinds', () => {
    expect(classifyStatus('Theorem, unconditional')).toBe('proved');
    expect(classifyStatus('Theorem, computer-assisted with certified error bounds')).toBe('proved');
    expect(classifyStatus('Numerical, linear theory, Planck 2018 parameters, converged to 0.1%')).toBe('computed');
    expect(classifyStatus('Exact lattice computation, continuum-extrapolated')).toBe('computed');
    expect(classifyStatus('Computation')).toBe('computed');
    expect(classifyStatus('Pre-registered data analysis')).toBe('measured');
    expect(classifyStatus('Analysis of public data')).toBe('measured');
    expect(classifyStatus('Conditional: needs a pre-inflationary Casimir era')).toBe('assumed');
    expect(classifyStatus('Not used; announced in an unrefereed preprint')).toBe('assumed');
    expect(classifyStatus('hypothesis')).toBe('assumed');
    expect(classifyStatus('Must be measured; currently d ≥ 0.985D')).toBe('open');
    expect(classifyStatus('Theorem; <a^-2> computed')).toBe('proved');
  });

  it('formats numbers, sizes and timestamps', () => {
    expect(fmtNum(15.962131481337693, 8)).toBe('15.962131');
    expect(fmtNum(1234567)).toBe('1.235×10⁶');
    expect(fmtNum(null)).toBe('—');
    expect(fmtSci(3.11e-5, 3)).toBe('3.11×10⁻⁵');
    expect(fmtBytes(1040592)).toBe('1016.2 KB');
    expect(fmtBytes(2202)).toBe('2.15 KB');
    expect(fmtUtc('2026-10-09T12:00:00+00:00')).toBe('2026-10-09 12:00 UTC');
  });
});

describe('sky geometry and colouring', () => {
  const sep = (a: number[], b: number[]): number => (Math.acos(Math.min(1, a[0]! * b[0]! + a[1]! * b[1]! + a[2]! * b[2]!)) * 180) / Math.PI;

  it('small circles keep a constant angular radius about their centre', () => {
    const c = lbToVec(123, -40);
    for (const p of smallCirclePoints([123, -40], 9.936, 64)) expect(sep(p, c)).toBeCloseTo(9.936, 9);
  });

  it('great arcs and polygons stay on the unit sphere and close', () => {
    const arc = greatArcPoints(lbToVec(0, 0), lbToVec(90, 0));
    expect(arc[0]![0]).toBeCloseTo(1, 12);
    expect(arc[arc.length - 1]![1]).toBeCloseTo(1, 12);
    for (const p of arc) expect(Math.hypot(...p)).toBeCloseTo(1, 12);
    const poly = polygonPoints([[10, 10], [12, 10], [12, 12], [10, 12]]);
    expect(sep(poly[0]!, poly[poly.length - 1]!)).toBeCloseTo(0, 5);
  });

  it('colours a map with a symmetric limit, masks pixels and picks values', () => {
    const map = Float32Array.from({ length: 12 }, (_, i) => (i - 6) * 10);
    const mask = Uint8Array.from({ length: 12 }, (_, i) => (i === 3 ? 0 : 1));
    expect(robustLimit(map, mask, 1)).toBe(60);
    const col = makeColorizer({ map, mask, limit: 60 }, 'dark');
    const out = new Uint8ClampedArray(3);
    expect(col.color(3, out, 0)).toBe(false);
    expect(col.color(6, out, 0)).toBe(true);
    expect(Array.from(out)).toEqual([col.lut[3 * 128]!, col.lut[3 * 128 + 1]!, col.lut[3 * 128 + 2]!]);
    const tex = equirectTexture({ map, mask, limit: 60 }, 'light', 16, 8);
    expect(tex).toHaveLength(16 * 8 * 4);
    expect(tex[3]).toBe(255);
    const p = pickValue({ map, mask }, 45, 41.8);
    expect(p.pix).toBe(0);
    expect(p.value).toBe(-60);
    expect(pickValue({ map, mask }, 315, 41.8)).toMatchObject({ pix: 3, masked: true });
  });
});
