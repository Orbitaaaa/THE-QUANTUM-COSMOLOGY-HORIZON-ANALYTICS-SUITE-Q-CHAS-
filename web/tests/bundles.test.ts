/**
 * The built data bundles (web/public/data, written by `qchas build`) carry every field the panels
 * read. Each bundle is checked only when it is present, so a fresh checkout without a build still
 * passes; once built, a renamed or missing field fails here instead of as a blank panel.
 */

import { describe, expect, it } from 'vitest';
import { ROUTES } from '../src/routes';

// Node's fs through process.getBuiltinModule (Node >= 22), so the browser-only tsconfig needs no
// @types/node; the tests run in vitest's node environment.
interface Fs {
  existsSync(p: URL): boolean;
  readFileSync(p: URL, enc: 'utf8'): string;
  statSync(p: URL): { size: number };
}
const fs = (globalThis as unknown as { process: { getBuiltinModule(id: string): unknown } }).process.getBuiltinModule('node:fs') as Fs;
const DATA = new URL('../public/data/', import.meta.url);
const at = (name: string): URL => new URL(name, DATA);
const has = (name: string): boolean => fs.existsSync(at(name));
// eslint-disable-next-line @typescript-eslint/no-explicit-any
const load = (name: string): any => JSON.parse(fs.readFileSync(at(name), 'utf8'));

const isNum = (x: unknown): x is number => typeof x === 'number' && Number.isFinite(x);
const numArray = (x: unknown, len?: number): void => {
  expect(Array.isArray(x)).toBe(true);
  const a = x as unknown[];
  if (len !== undefined) expect(a.length).toBe(len);
  expect(a.length).toBeGreaterThan(0);
  expect(a.every((v) => v === null || isNum(v))).toBe(true);
};
const keys = (obj: unknown, ...names: string[]): void => {
  expect(obj, `expected an object with ${names.join(', ')}`).toBeTypeOf('object');
  for (const k of names) expect(obj, `missing key ${k}`).toHaveProperty(k);
};

describe.skipIf(!has('manifest.json'))('manifest.json', () => {
  const m = has('manifest.json') ? load('manifest.json') : null;
  it('names the CMB map, lists the datasets with provenance and the bundles written', () => {
    keys(m, 'generated_utc', 'qchas_version', 'cmb_map_source', 'datasets', 'bundles', 'warnings');
    expect(typeof m.cmb_map_source).toBe('string');
    for (const d of Object.values(m.datasets) as Record<string, unknown>[]) keys(d, 'filename', 'bytes', 'sha256', 'retrieved_utc', 'mission');
  });
  it('lists only bundles that exist on disk', () => {
    for (const b of m.bundles as string[]) expect(has(b), `${b} is in the manifest but not on disk`).toBe(true);
  });
  it('covers every bundle a route declares, so the navigation dots are truthful', () => {
    const built = new Set<string>([...(m.bundles as string[]), 'manifest.json']);
    for (const r of ROUTES) for (const b of r.bundles) expect(built.has(b), `${r.id} declares ${b}`).toBe(true);
  });
});

describe.skipIf(!has('topology.json'))('topology.json', () => {
  const t = has('topology.json') ? load('topology.json') : null;
  it('has the 18 flat manifolds with their four principle verdicts', () => {
    expect(t.manifolds).toHaveLength(18);
    for (const m of t.manifolds) {
      keys(m, 'label', 'name', 'closed', 'orientable', 'holonomy', 'translations_only', 'principles');
      keys(m.principles, 'P1', 'P2', 'P3', 'P4');
    }
    expect(t.principles).toHaveLength(4);
  });
  it('carries the quoted Planck tables used by the selector and the ghost-twin tracker', () => {
    for (const r of t.evidence_torus) keys(r, 'shape', 'd_over_D', 'max_dlnL', 'lnB', 'null95', 'p');
    for (const r of t.exclusion_torus) keys(r, 'shape', 'd_over_D', 'lnB_planck', 'median_if_true', 'q05_if_true', 'p_le_planck', 'verdict');
    expect(t.other_flat.length).toBeGreaterThan(0);
    expect(t.curved.length).toBeGreaterThan(0);
    expect(t.weeks.length).toBeGreaterThan(0);
  });
  it('has claims with a status and a source, circle tables and the vacuum grids', () => {
    for (const c of t.claims) keys(c, 'statement', 'status', 'source');
    keys(t.circles, 'tables', 'fcc_preregistered');
    expect(t.circles.fcc_preregistered.R).toHaveLength(3);
    keys(t.vacuum, 'E1', 'E2');
  });
});

describe.skipIf(!has('lattice.json'))('lattice.json', () => {
  const L = has('lattice.json') ? load('lattice.json') : null;
  it('has the rhombohedral curve and the three certified critical points', () => {
    numArray(L.rhombohedral.alpha_deg);
    numArray(L.rhombohedral.Z, L.rhombohedral.alpha_deg.length);
    for (const s of ['SC', 'BCC', 'FCC']) {
      keys(L.critical[s], 'Z', 'Z_paper', 'hessian', 'hessian_paper');
      numArray(L.critical[s].hessian, 5);
      numArray(L.critical[s].hessian_paper, 5);
      expect(L.bases[s]).toHaveLength(3);
    }
  });
  it('has the shape flows with counts, the paper counts and trajectories', () => {
    keys(L.flows, 'n', 'counts', 'paper_counts', 'trajectories');
    const total = Object.values(L.flows.counts as Record<string, number>).reduce((a, b) => a + b, 0);
    expect(total).toBe(L.flows.n);
    for (const tr of L.flows.trajectories) {
      keys(tr, 'tau', 'Z', 'end');
      expect(tr.tau.length).toBe(tr.Z.length);
    }
    keys(L.relaxation, 'kappa_cE_rate_fcc', 'TL_today', 'source');
  });
});

describe.skipIf(!has('sky_map.json'))('sky_map.json and its binary maps', () => {
  const S = has('sky_map.json') ? load('sky_map.json') : null;
  it('describes RING Galactic maps in μK whose files have 12 nside² float32 pixels', () => {
    keys(S, 'nside', 'ordering', 'frame', 'unit', 'maps', 'source');
    expect(S.ordering).toBe('RING');
    const npix = 12 * S.nside * S.nside;
    for (const m of S.maps) {
      keys(m, 'name', 'file', 'description');
      expect(fs.statSync(at(m.file)).size).toBe(4 * npix);
    }
    if (S.mask_file) expect(fs.statSync(at(S.mask_file)).size).toBe(npix);
  });
  it('names the same map as the manifest', () => {
    if (has('manifest.json')) expect(S.source).toBe(load('manifest.json').cmb_map_source);
  });
});

describe.skipIf(!has('cmb.json'))('cmb.json', () => {
  const C = has('cmb.json') ? load('cmb.json') : null;
  it('has spectra on one ℓ grid', () => {
    numArray(C.ell);
    numArray(C.map_D, C.ell.length);
    numArray(C.lcdm_D, C.ell.length);
    keys(C.observed, 'wmap');
    keys(C.observed.wmap, 'ell', 'D', 'err');
    for (const v of Object.values(C.torus_D)) numArray(v, C.ell.length);
  });
  it('has the anomaly statistics, axes and the Sec. 8 table', () => {
    keys(C.stats, 'D2_obs', 'S12_obs', 'S12_lcdm_median', 'p_C2', 'p_S12', 'D2', 'S12', 'C_theta', 'S12_lcdm_hist', 'paper', 'axes');
    keys(C.stats.paper, 'S12_planck', 'D2_planck', 'p_S12_planck', 'p_C2_planck');
    keys(C.stats.axes, 'angle_deg', 'p_isotropic');
    for (const r of C.sec8.rows) keys(r, 'paper');
  });
  it('has the KL table, the ceiling and the orientation-scanned likelihood', () => {
    numArray(C.kl.d_over_D);
    for (const s of ['cubic', 'BCC', 'FCC']) numArray(C.kl[s], C.kl.d_over_D.length);
    keys(C.ceiling, 'd_over_D', 'causal_ball', 'cmb_l30', 'source');
    expect(C.likelihood.length).toBeGreaterThan(0);
    for (const r of C.likelihood) {
      keys(r, 'shape', 'd_over_D', 'max_dlnL', 'lnB', 'null_q95', 'p', 'p_max', 'best_axes_gal', 'exclusion');
      keys(r.exclusion, 'verdict');
    }
    keys(C.likelihood_meta, 'map', 'n_null', 'grid', 'global', 'fcc_family', 'fcc_105_axes_vs_paper');
  });
  it('labels its own analysis with the map the manifest names, never as Planck by default', () => {
    if (has('manifest.json')) {
      const src = load('manifest.json').cmb_map_source;
      expect(C.source).toBe(src);
      expect(C.likelihood_meta.map).toBe(src);
    }
  });
});

describe.skipIf(!has('resonance.json'))('resonance.json', () => {
  const R = has('resonance.json') ? load('resonance.json') : null;
  it('has the KL table with paper values and the gains', () => {
    for (const r of R.kl_table) keys(r, 'shape', 'd_over_D', 'no_line', 'kl', 'paper');
    expect(R.kl_table.some((r: { gain?: number }) => isNum(r.gain))).toBe(true);
  });
  it('has 80 shells, the Woodbury benchmark, single-field limits and the 3D exits', () => {
    const n = Object.values(R.shells as Record<string, unknown[]>).reduce((a, s) => a + s.length, 0);
    expect(n).toBeGreaterThan(0);
    keys(R.woodbury_benchmark, 'ranks', 'direct_ms', 'woodbury_ms', 'max_abs_diff');
    keys(R.single_field, 'rows', 'n_shells', 'n_f_below_0.05');
    for (const r of R.exits_3d) keys(r, 'spectrum', 'd1.5', 'd2', 'd3');
  });
  it('has the search on the build map next to the paper search on Planck', () => {
    keys(R.search_paper, 'S_max', 'p_max', 'S_avg', 'p_avg');
    if (R.search) {
      keys(R.search, 'map', 'models', 'S_max', 'S_avg', 'p_max', 'p_avg', 'null', 'best');
      numArray(R.search.null.S_max);
      numArray(R.search.null.S_avg);
      if (has('manifest.json')) expect(R.search.map).toBe(load('manifest.json').cmb_map_source);
    }
  });
});

describe.skipIf(!has('spherex.json'))('spherex.json', () => {
  const X = has('spherex.json') ? load('spherex.json') : null;
  it('has frames with Galactic footprints and thumbnails whose images exist', () => {
    expect(X.frames.length).toBe(X.n_frames);
    for (const f of X.frames) {
      keys(f, 'obsid', 'l', 'b', 'band_um', 'corners_gal');
      expect(f.corners_gal).toHaveLength(4);
    }
    for (const t of X.thumbnails) {
      keys(t, 'obsid', 'file', 'l', 'b', 'band_um', 'unit', 'p1', 'p995');
      expect(has(t.file), t.file).toBe(true);
    }
  });
});

describe.skipIf(!has('horizon.json'))('horizon.json', () => {
  const H = has('horizon.json') ? load('horizon.json') : null;
  it('has the Page curves on one time grid', () => {
    numArray(H.page.t);
    for (const k of ['hawking', 'bh_entropy', 'page', 'island_min', 'interior_volume']) numArray(H.page[k], H.page.t.length);
    keys(H.page, 'page_time', 'page_mass_fraction', 'beta');
  });
  it('has capacity, QES rows with paper values, atmosphere, evaporation and decoding', () => {
    keys(H.capacity, 'ladder', 'de_sitter', 'min_area_lp2', 'window_examples');
    for (const r of H.qes.rows) keys(r, 'mass', 'c', 'rb_over_rs', 'depth', 'lag', 't_scr', 'paper');
    keys(H.atmosphere, 'gaps', 'share', 'paper_share', 'area_coeff');
    keys(H.evaporation, 'no_diary', 'early_diary', 'late_diary', 'paper');
    keys(H.decoding, 'haar_old', 'eth', 'haar_young', 'paper');
    keys(H.kerr_spin_down, 'fixed_point_scalar', 'rows');
    for (const c of H.claims) keys(c, 'statement', 'status', 'source');
  });
});

describe.skipIf(!has('bounce.json'))('bounce.json', () => {
  const B = has('bounce.json') ? load('bounce.json') : null;
  it('has exact runs with the printed Theorem Λ values and the paper-grid rerun', () => {
    for (const r of B.runs) {
      keys(r, 'a_max', 'T', 'a_mean', 'a_inv2', 'a_min', 'peak_a_inv2', 'paper');
      expect(r.T.length).toBe(r.a_mean.length);
    }
    expect(B.runs.some((r: { paper_grid?: unknown }) => r.paper_grid)).toBe(true);
  });
  it('has snapshots, the spectrum check, gravitons, mixmaster, the Planck test and units', () => {
    keys(B.snapshots, 'a', 'T', 'psi2');
    keys(B.spectrum_check, 'j', 'E_numeric', 'E_exact');
    keys(B.graviton, 'n', 'h2_peak_by_a_max');
    keys(B.mixmaster, 'a_min', 'source');
    keys(B.planck_bounce_test, 'templates', 'kc_excluded_above', 'source');
    keys(B.units, 'ell_over_lP');
    keys(B.limits, 'a_min', 'fwhm_times_amax');
  });
});

describe.skipIf(!has('gw.json'))('gw.json', () => {
  const G = has('gw.json') ? load('gw.json') : null;
  it('has events with strain, ringdown fits and remnant thermodynamics', () => {
    expect(G.events.length).toBeGreaterThan(0);
    for (const e of G.events) {
      keys(e, 'name', 'm1', 'm2', 'mf', 'af', 'f220_Hz', 'tau220_ms', 'echo_delay_ms', 'S_BH_log10', 'T_H_K', 't_evap_log10_yr', 'strain', 'ringdown_fit');
      keys(e.strain, 't', 'H1', 'L1', 'fs');
      expect(e.strain.H1.length).toBe(e.strain.t.length);
    }
  });
  it('has the quoted echo search of the paper', () => {
    for (const e of G.paper_events) keys(e, 'event', 'snr', 'mf', 'af', 'f220', 'dt_ms', 'R90');
    keys(G.paper_echo_summary, 'combined_R90', 'p_combined');
    expect('gwtc' in G).toBe(true);
  });
});

describe.skipIf(!has('validation.json'))('validation.json', () => {
  const V = has('validation.json') ? load('validation.json') : null;
  it('has rows with ours, paper and verdict, and a consistent summary', () => {
    for (const r of V.rows) {
      keys(r, 'group', 'quantity', 'ours', 'paper', 'rel_diff', 'verdict', 'source', 'note');
      expect(['match', 'differs']).toContain(r.verdict);
      expect(isNum(r.ours) && isNum(r.paper)).toBe(true);
    }
    expect(V.summary.total).toBe(V.rows.length);
    expect(V.summary.match + V.summary.differs).toBe(V.summary.total);
    expect(V.rows.filter((r: { verdict: string }) => r.verdict === 'match').length).toBe(V.summary.match);
  });
  it('states every disagreement in a note', () => {
    for (const r of V.rows.filter((x: { verdict: string }) => x.verdict !== 'match')) expect(r.note.length).toBeGreaterThan(0);
  });
});
