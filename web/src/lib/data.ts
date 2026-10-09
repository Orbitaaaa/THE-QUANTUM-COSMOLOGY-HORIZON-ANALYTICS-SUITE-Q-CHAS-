/**
 * Typed loaders for the data bundles of docs/ARCHITECTURE.md ("Data contract").
 *
 * Bundles are produced by `qchas build` into web/public/data/ and fetched at runtime from
 * ./data/. A missing bundle is reported as such (`{ ok: false, reason: 'missing' }`) and the
 * panel renders an explicit "not built yet" state; nothing here ever substitutes placeholder
 * values for data. Results are cached per URL for the lifetime of the page.
 */

/** Galactic (l, b) in degrees. */
export type LBPair = [number, number];
export type Matrix3 = [number[], number[], number[]];

/** A number the paper prints, shipped next to ours. */
export type PaperValue = number | string | null;

// --- manifest.json ---------------------------------------------------------------------

export interface DatasetRecord {
  filename: string;
  url: string;
  bytes: number;
  sha256: string;
  retrieved_utc: string;
  mission: string;
  description: string;
  citation: string;
}

export interface Manifest {
  generated_utc: string;
  qchas_version: string;
  /** e.g. "WMAP9-ILC" or "Planck-PR3-Commander" */
  cmb_map_source: string | null;
  datasets: Record<string, DatasetRecord>;
  bundles: string[];
  warnings: string[];
}

// --- topology.json -------------------------------------------------------------------

export type PrincipleStatus = 'allowed' | 'excluded' | 'no spin structure';

export interface ManifoldEntry {
  label: string;
  name: string;
  closed: boolean;
  orientable: boolean;
  holonomy: string;
  holonomy_order?: number | null;
  translations_only: boolean;
  periodic_spin_possible: boolean;
  principles: Record<'P1' | 'P2' | 'P3' | 'P4', PrincipleStatus | string>;
  note?: string;
}

export interface PrincipleEntry {
  key: string;
  name: string;
  status: string;
  statement: string;
  theorem: string;
  source?: string;
}

export interface CurvedSpaceEntry {
  name: string;
  order: number;
  shortest_loop: string;
  omega_k_range: [number, number];
  lnB_vs_S3: number;
  note: string;
  source: string;
}

export interface WeeksEntry {
  omega_k: number;
  Rc_Gpc: number;
  chi_over_Rc: number;
  shortest_loop_over_D: number;
  volume_ratio: number;
  circles: string;
  source?: string;
}

export interface EvidenceTorusEntry {
  shape: string;
  d_over_D: number;
  max_dlnL: number;
  lnB: number;
  null95: number;
  p: number;
  source?: string;
}

export interface ExclusionTorusEntry {
  shape: string;
  d_over_D: number;
  lnB_planck: number;
  median_if_true: number;
  q05_if_true: number;
  p_le_planck: string;
  sims: number;
  verdict: string;
  source?: string;
}

export interface OtherFlatEntry {
  space: string;
  orientable: boolean;
  samples: number;
  lnB: string;
  largest: number;
  source?: string;
}

export interface Claim {
  statement: string;
  status: string;
  source: string;
}

export interface TopologyBundle {
  manifolds: ManifoldEntry[];
  principles: PrincipleEntry[];
  curved: CurvedSpaceEntry[];
  weeks: WeeksEntry[];
  evidence_torus: EvidenceTorusEntry[];
  exclusion_torus: ExclusionTorusEntry[];
  other_flat: OtherFlatEntry[];
  claims: Claim[];
}

// --- lattice.json --------------------------------------------------------------------

export type ShapeKey = 'SC' | 'BCC' | 'FCC';

export interface CriticalPoint {
  Z: number;
  Z_paper: number;
  hessian: number[];
  hessian_paper: number[];
  gradient_norm?: number;
  computed?: boolean;
  source?: string;
}

export interface FlowTrajectory {
  tau: number[];
  Z: number[];
  end: string;
}

export interface LatticeBundle {
  bases: Record<ShapeKey, Matrix3>;
  rhombohedral: { alpha_deg: number[]; Z: number[]; angles_deg?: Record<string, number> };
  critical: Record<ShapeKey, CriticalPoint>;
  gaps?: { computed: Record<string, number>; paper: Record<string, number>; source: string };
  flows: {
    n: number;
    counts: Record<string, number>;
    paper_counts: Record<string, number>;
    trajectories: FlowTrajectory[];
    source?: string;
    [extra: string]: unknown;
  };
  relaxation: { kappa_cE_rate_fcc: number; TL_today: number; source: string; [extra: string]: unknown };
}

// --- sky_map.json + binaries -----------------------------------------------------------

export interface SkyMapEntry {
  name: string;
  file: string;
  description: string;
  min: number;
  max: number;
}

export interface SkyMapMeta {
  nside: number;
  ordering: 'RING';
  frame: 'galactic';
  unit: string;
  maps: SkyMapEntry[];
  mask_file: string | null;
  source: string;
}

// --- cmb.json ------------------------------------------------------------------------

export interface SpectrumSeries {
  ell: number[];
  D: number[];
  err: number[];
}

export type ShapeSeries = Record<string, number[]>;

export interface LikelihoodRow {
  shape: string;
  d_over_D: number;
  max_dlnL: number;
  lnB: number;
  null_q95: number;
  p: number;
  n_null: number;
  best_axes_gal: LBPair[];
  paper?: { max_dlnL?: PaperValue; lnB?: PaperValue; p?: PaperValue };
}

export interface CmbBundle {
  source: string;
  ell: number[];
  observed: { wmap: SpectrumSeries | null; planck: SpectrumSeries | null };
  map_D: number[];
  lcdm_D: number[];
  torus_D: ShapeSeries;
  stats: {
    D2_obs: number;
    S12_obs: number;
    S12_lcdm_median: number;
    p_C2: Record<string, number | null>;
    p_S12: Record<string, number | null>;
    [extra: string]: unknown;
  };
  kl: {
    d_over_D: number[];
    cubic: number[];
    BCC: number[];
    FCC: number[];
    paper: Record<string, Record<string, number>>;
    pinsker_bound: ShapeSeries;
  };
  ceiling: {
    d_over_D: number[];
    causal_ball: number[];
    boundary_sphere: number[];
    lss_ball: number[];
    cmb_l30: number[];
    source: string;
  };
  likelihood: LikelihoodRow[];
}

// --- resonance.json ------------------------------------------------------------------

export interface ResonanceKlRow {
  shape: string;
  d_over_D: number;
  no_line: number;
  shell: number;
  n_vectors: number;
  f: number;
  kl: number;
  gain: number;
  paper?: { no_line?: PaperValue; kl?: PaperValue; gain?: PaperValue };
}

export interface ResonanceShell {
  m2: number;
  k_chi: number;
  n_vectors: number;
}

export interface ResonanceSearchModel {
  shape: string;
  d_over_D: number;
  shell_m2: number;
  f: number;
  lnB: number;
}

export interface ResonanceBundle {
  kl_table: ResonanceKlRow[];
  shells: Record<string, ResonanceShell[]>;
  woodbury_benchmark: {
    n: number;
    ranks: number[];
    direct_ms: number[];
    woodbury_ms: number[];
    max_abs_diff: number;
  };
  search: {
    map: string;
    models: ResonanceSearchModel[];
    S_max: number;
    S_avg: number;
    null: { S_max: number[]; S_avg: number[] };
    p_max: number;
    p_avg: number;
    paper: { S_max: number; p_max: number; S_avg: number; p_avg: number };
  } | null;
  exits_3d: { spectrum: string; [d: string]: number | string }[];
  source: string;
}

// --- spherex.json --------------------------------------------------------------------

export interface SpherexFrame {
  obsid: string;
  detector: number;
  date: string;
  ra: number;
  dec: number;
  l: number;
  b: number;
  corners_gal: LBPair[];
  band_um: [number, number];
  key: string;
}

export interface SpherexThumbnail {
  obsid: string;
  detector: number;
  file: string;
  l: number;
  b: number;
  band_um: [number, number];
  unit: string;
  p1: number;
  p995: number;
}

export interface SpherexBundle {
  release: string;
  bucket: string;
  n_frames: number;
  /** band -> [lambda_min_um, lambda_max_um, n_channels] */
  bands: Record<string, [number, number, number]>;
  frames: SpherexFrame[];
  thumbnails: SpherexThumbnail[];
}

// --- horizon.json --------------------------------------------------------------------

export interface QesRow {
  mass: number;
  c: number;
  rb_over_rs: number;
  depth: number;
  lag: number;
  t_scr: number;
  s_gen_minus_s_bh: number;
  paper: Partial<Record<'depth' | 'lag' | 't_scr' | 's_gen_minus_s_bh', number>>;
}

/** Output of qchas.horizon.evaporation.evaporate(). */
export interface EvaporationRun {
  n_qubits: number;
  diary_step: number | null;
  seed: number;
  steps: number[];
  s_rad: number[];
  s_bh: number[];
  hawking: number[];
  envelope: number[];
  page_average: number[];
  diary_mi: number[];
}

/** Output of qchas.horizon.evaporation.decoding_infidelity(). */
export interface DecodingRun {
  n_bh: number;
  k: number;
  old: boolean;
  dynamics: string;
  t: number | null;
  s: number[];
  infidelity: number[];
  bound: number[];
  [extra: string]: unknown;
}

export interface HorizonBundle {
  page: {
    t: number[];
    hawking: number[];
    bh_entropy: number[];
    page: number[];
    island_min: number[];
    interior_volume: number[];
    page_time: number;
    page_mass_fraction: number;
    beta: number;
  };
  capacity: {
    ladder: { type: string; label: string; trace: string; entropy: string; source: string }[];
    de_sitter: { entropy_nats: number; log10_entropy: number; horizon_radius_Gpc: number; [extra: string]: number };
    min_area_lp2: number;
    window_examples: { m_star: number; width: number; ln_trace: number; delta_S_BH: number }[];
  };
  qes: { rows: QesRow[]; laws: Record<string, unknown> };
  atmosphere: {
    gaps: number[];
    share: number[];
    paper_share: number[];
    area_coeff: number;
    continuum: Record<string, number>;
  };
  evaporation: {
    no_diary: EvaporationRun;
    early_diary: EvaporationRun;
    late_diary: EvaporationRun;
    paper: { max_bits: number; max_dev: number; early_recovered_step: number; late_recovered_step: number };
  };
  decoding: {
    haar_old: DecodingRun;
    eth: Record<string, DecodingRun>;
    haar_young: DecodingRun;
    paper: Record<string, unknown>;
  };
  kerr_spin_down: { fixed_point_scalar: number; rows: Record<string, number | string>[] };
  claims: Claim[];
}

// --- bounce.json ---------------------------------------------------------------------

export interface BounceRun {
  a_max: number;
  T: number[];
  a_mean: number[];
  a_inv2: number[];
  a_min: number;
  peak_a_inv2: number;
  fwhm_T: number;
  paper: { a_min?: PaperValue; peak_a_inv2?: PaperValue; fwhm_T?: PaperValue };
}

export interface BounceBundle {
  runs: BounceRun[];
  spectrum_check: { j: number[]; E_numeric: number[]; E_exact: number[] };
  snapshots: { a: number[]; T: number[]; psi2: number[][] };
  graviton: { n: number[]; h2_peak: number[] };
  mixmaster: { a_min: number; source: string };
  planck_bounce_test: { dchi2_cutoff_max: number; kc_excluded_above: number; source: string };
  units: { ell_over_lP: number; [extra: string]: number | string };
  /** [QG] status rows for the bounce results (optional) */
  claims?: Claim[];
}

// --- gw.json -------------------------------------------------------------------------

export interface GwEvent {
  name: string;
  source: string;
  m1: number;
  m2: number;
  a1: number | null;
  a2: number | null;
  mf: number;
  af: number;
  f220_Hz: number;
  tau220_ms: number;
  echo_delay_ms: number;
  S_BH_log10: number;
  echo_bound_log10: number | null;
  T_H_K: number;
  t_evap_log10_yr: number;
  page_time_log10_yr: number;
  strain: { t: number[]; H1: number[]; L1: number[]; fs: number; band: [number, number] } | null;
  ringdown_fit: { f_Hz: number; tau_ms: number; t0_ms_after_peak: number } | null;
}

export interface GwPaperEvent {
  event: string;
  snr: number;
  mf: number;
  af: number;
  f220: number;
  dt_ms: number;
  R90: number;
}

export interface GwBundle {
  events: GwEvent[];
  paper_events: GwPaperEvent[];
  gwtc: Record<string, unknown>[] | null;
  combined_R90: number;
}

// --- registry and loaders --------------------------------------------------------------

export interface BundleMap {
  'manifest.json': Manifest;
  'topology.json': TopologyBundle;
  'lattice.json': LatticeBundle;
  'sky_map.json': SkyMapMeta;
  'cmb.json': CmbBundle;
  'resonance.json': ResonanceBundle;
  'spherex.json': SpherexBundle;
  'horizon.json': HorizonBundle;
  'bounce.json': BounceBundle;
  'gw.json': GwBundle;
}

export type BundleName = keyof BundleMap;

export const BUNDLE_NAMES: BundleName[] = [
  'manifest.json',
  'topology.json',
  'lattice.json',
  'sky_map.json',
  'cmb.json',
  'resonance.json',
  'spherex.json',
  'horizon.json',
  'bounce.json',
  'gw.json',
];

/** Top-level keys each bundle must carry to be accepted. */
const REQUIRED_KEYS: Record<BundleName, string[]> = {
  'manifest.json': ['generated_utc', 'datasets'],
  'topology.json': ['manifolds', 'claims'],
  'lattice.json': ['bases', 'critical'],
  'sky_map.json': ['nside', 'maps'],
  'cmb.json': ['ell'],
  'resonance.json': ['kl_table'],
  'spherex.json': ['frames'],
  'horizon.json': ['page'],
  'bounce.json': ['runs'],
  'gw.json': ['events'],
};

export type LoadFailure = { ok: false; reason: 'missing' | 'invalid' | 'network'; name: string; message: string };
export type LoadResult<T> = { ok: true; data: T; name: string } | LoadFailure;

/** The command that produces the bundles, shown in every "not built yet" state. */
export const BUILD_COMMAND = 'qchas build';

const DATA_DIR = `${import.meta.env?.BASE_URL ?? './'}data/`;

let fetcher: typeof fetch = (...args) => fetch(...args);
const cache = new Map<string, Promise<unknown>>();

/** Replace the fetch implementation (tests) and drop the cache. */
export function setFetcher(f: typeof fetch): void {
  fetcher = f;
  cache.clear();
}

export function clearCache(): void {
  cache.clear();
}

export const dataUrl = (file: string): string => DATA_DIR + file;

const looksLikeHtml = (res: Response, text?: string): boolean =>
  (res.headers.get('content-type') ?? '').includes('text/html') || (text ?? '').trimStart().startsWith('<');

async function fetchJson(name: string): Promise<LoadResult<unknown>> {
  let res: Response;
  try {
    res = await fetcher(dataUrl(name), { cache: 'no-cache' });
  } catch (err) {
    return { ok: false, reason: 'network', name, message: String(err) };
  }
  if (res.status === 404) return { ok: false, reason: 'missing', name, message: `${name} not found` };
  if (!res.ok) return { ok: false, reason: 'network', name, message: `HTTP ${res.status} for ${name}` };
  const text = await res.text();
  // dev servers answer unknown paths with index.html (SPA fallback): that is a missing bundle
  if (looksLikeHtml(res, text)) return { ok: false, reason: 'missing', name, message: `${name} not found` };
  try {
    return { ok: true, name, data: JSON.parse(text) };
  } catch (err) {
    return { ok: false, reason: 'invalid', name, message: `${name} is not valid JSON: ${String(err)}` };
  }
}

/** Load and cache one JSON bundle. Never throws. */
export function loadBundle<K extends BundleName>(name: K): Promise<LoadResult<BundleMap[K]>> {
  let p = cache.get(name) as Promise<LoadResult<BundleMap[K]>> | undefined;
  if (!p) {
    p = fetchJson(name).then((r) => {
      if (!r.ok) return r;
      const d = r.data;
      const missing = REQUIRED_KEYS[name].filter((k) => !(d && typeof d === 'object' && k in d));
      if (missing.length)
        return { ok: false, reason: 'invalid', name, message: `${name} lacks required keys: ${missing.join(', ')}` } as LoadFailure;
      return { ok: true, name, data: d as BundleMap[K] };
    });
    cache.set(name, p);
  }
  return p;
}

/** Load several bundles at once; resolves to a record keyed by bundle name. */
export async function loadBundles<K extends BundleName>(names: K[]): Promise<{ [P in K]: LoadResult<BundleMap[P]> }> {
  const results = await Promise.all(names.map((n) => loadBundle(n)));
  return Object.fromEntries(names.map((n, i) => [n, results[i]])) as { [P in K]: LoadResult<BundleMap[P]> };
}

/** Load a little-endian binary file from ./data/ as Float32 or UInt8, checking its length. */
export function loadBinary(file: string, kind: 'f32', expectedLength?: number): Promise<LoadResult<Float32Array>>;
export function loadBinary(file: string, kind: 'u8', expectedLength?: number): Promise<LoadResult<Uint8Array>>;
export function loadBinary(
  file: string,
  kind: 'f32' | 'u8',
  expectedLength?: number,
): Promise<LoadResult<Float32Array | Uint8Array>> {
  const key = `bin:${kind}:${file}`;
  let p = cache.get(key) as Promise<LoadResult<Float32Array | Uint8Array>> | undefined;
  if (!p) {
    p = (async (): Promise<LoadResult<Float32Array | Uint8Array>> => {
      let res: Response;
      try {
        res = await fetcher(dataUrl(file));
      } catch (err) {
        return { ok: false, reason: 'network', name: file, message: String(err) };
      }
      if (res.status === 404 || (res.ok && looksLikeHtml(res)))
        return { ok: false, reason: 'missing', name: file, message: `${file} not found` };
      if (!res.ok) return { ok: false, reason: 'network', name: file, message: `HTTP ${res.status} for ${file}` };
      const buf = await res.arrayBuffer();
      const width = kind === 'f32' ? 4 : 1;
      if (buf.byteLength % width !== 0)
        return { ok: false, reason: 'invalid', name: file, message: `${file}: ${buf.byteLength} bytes is not a whole number of values` };
      const n = buf.byteLength / width;
      if (expectedLength !== undefined && n !== expectedLength)
        return { ok: false, reason: 'invalid', name: file, message: `${file}: ${n} values, expected ${expectedLength}` };
      if (kind === 'u8') return { ok: true, name: file, data: new Uint8Array(buf) };
      const view = new DataView(buf);
      const out = new Float32Array(n);
      for (let i = 0; i < n; i++) out[i] = view.getFloat32(4 * i, true);
      return { ok: true, name: file, data: out };
    })();
    cache.set(key, p);
  }
  return p;
}

export interface SkyMapData {
  meta: SkyMapMeta;
  entry: SkyMapEntry;
  map: Float32Array;
  mask: Uint8Array | null;
}

/** Load sky_map.json and one of its maps (default: the first) plus the mask when present. */
export async function loadSkyMap(name?: string): Promise<LoadResult<SkyMapData>> {
  const meta = await loadBundle('sky_map.json');
  if (!meta.ok) return meta;
  const entry = name ? meta.data.maps.find((m) => m.name === name) : meta.data.maps[0];
  if (!entry) return { ok: false, reason: 'invalid', name: 'sky_map.json', message: `no map named ${name ?? '(any)'}` };
  const n = 12 * meta.data.nside * meta.data.nside;
  const map = await loadBinary(entry.file, 'f32', n);
  if (!map.ok) return map;
  let mask: Uint8Array | null = null;
  if (meta.data.mask_file) {
    const m = await loadBinary(meta.data.mask_file, 'u8', n);
    if (m.ok) mask = m.data;
  }
  return { ok: true, name: entry.file, data: { meta: meta.data, entry, map: map.data, mask } };
}
