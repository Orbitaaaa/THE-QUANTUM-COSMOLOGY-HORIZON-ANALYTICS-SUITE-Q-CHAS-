/**
 * Typed loaders for the data bundles of docs/ARCHITECTURE.md ("Data contract").
 *
 * Bundles are produced by `qchas build` into web/public/data/ and fetched at runtime from
 * ./data/. A missing bundle is reported as such (`{ ok: false, reason: 'missing' }`) and the
 * panel renders an explicit "not built yet" state; nothing here ever substitutes placeholder
 * values for data. Results are cached per URL for the lifetime of the page, except transient
 * network failures, which are retried on the next request.
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

/** A bundle whose builder raised during `qchas build` (the other bundles are still written). */
export interface BuildError {
  bundle: string;
  error: string;
  trace?: string;
}

export interface Manifest {
  generated_utc: string;
  qchas_version: string;
  /** e.g. "WMAP9-ILC" or "Planck-PR3-Commander"; null when no sky map could be prepared */
  cmb_map_source: string | null;
  datasets: Record<string, DatasetRecord>;
  /** files actually written (JSON bundles and .bin maps) */
  bundles: string[];
  warnings: string[];
  /** builders that failed (export.py); absent in older manifests */
  errors?: BuildError[];
  fast_build?: boolean;
  timings_s?: Record<string, number>;
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

/** One cubic-family circle pattern in the lattice frame (qchas.topology.circles.circle_table). */
export interface CircleTable {
  shape: string;
  d_over_D: number;
  n_pairs: number;
  radius_deg: number | null;
  angles_between_centres_deg: number[];
  shells: { length_over_D: number; radius_deg: number; n_pairs: number; angles_between_centres_deg: number[] }[];
  computed?: boolean;
}

/** Theorem B / B-loc visual: renormalized local variance on a z-slice (qchas.topology.vacuum). */
export interface VacuumSlice {
  manifold: string;
  m_L: number;
  L: number;
  x: number[];
  y: number[];
  z: number;
  excess: number[][];
  excess_minus_min: number[][];
  min: number;
  max: number;
  smoothing_sigma: number | null;
  cutoff_mL: number;
  computed?: boolean;
  source?: string;
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
  /** matched-circle tables at 0.985 D and 0.8 D and the [TOP] Sec. 9 FCC orientation (not in the contract yet) */
  circles?: {
    chi_star_over_D: number;
    D_Gpc: number;
    tables: CircleTable[];
    paper_table: { shape: string; pairs: number; centres: string; radius: string | null; angles_deg: number[] }[];
    paper_fig4: { d_over_D: number; radius_deg: number; pairs: Record<string, number> };
    fcc_preregistered: {
      R: Matrix3;
      assignment: number[];
      signs: number[];
      residual_rms_deg: number;
      residual_max_deg: number;
      axes_gal_paper: LBPair[];
      source: string;
    };
    source: string;
  };
  /** Theorem B / B-loc vacuum-variance slices (not in the contract yet) */
  vacuum?: {
    E1?: VacuumSlice;
    E2?: VacuumSlice;
    /** Gaussian-smoothed versions (sigma = 0.35 lattice units, as [TOP] Fig. 1) */
    E1_smoothed?: VacuumSlice;
    E2_smoothed?: VacuumSlice;
    /** (x, y) in lattice units where the E2 half-turn screw axes cross the slice */
    E2_axes?: [number, number][];
    note?: string;
  };
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
  /** tau = ln(conformal time), Theorem H.2 */
  tau: number[];
  /** Z(2) at unit covolume along the flow (rounded to 1e-8 in the bundle) */
  Z: number[];
  /** 'FCC' | 'BCC' | 'other' */
  end: string;
  converged?: boolean;
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
    /** a sample of the flows (the counts cover all n) */
    trajectories: FlowTrajectory[];
    /** 'gradient' (Theorem H.2 limit) or 'kinetic' */
    mode?: string;
    kappa?: number | null;
    /** Haar starts truncated at |t| <= t_max */
    t_max?: number;
    seed?: number;
    fcc_fraction?: number;
    fcc_fraction_err?: number;
    paper_fcc_fraction?: number;
    unconverged?: number;
    computed?: boolean;
    source?: string;
    [extra: string]: unknown;
  };
  /** [TOP] Prop. H.3: quoted values, with our recomputation in `ours` */
  relaxation: {
    kappa_cE_rate_fcc: number;
    TL_today: number;
    kappa_cE_coeff?: number;
    TL4_per_dtau?: number;
    source: string;
    ours?: {
      kappa_cE_coeff?: number;
      fcc_mu_min?: number;
      TL4_per_dtau?: number;
      TL_today?: number;
      L_today_inv_eV?: number;
    };
    computed?: boolean;
    [extra: string]: unknown;
  };
}

// --- sky_map.json + binaries -----------------------------------------------------------

export interface SkyMapEntry {
  name: string;
  file: string;
  description: string;
  /** extremes of the map (often Galactic-plane residuals; not a good colour limit) */
  min: number;
  max: number;
  /** symmetric display range [-r, r], r = 99.5th percentile of |value| (export.py) */
  range?: [number, number];
}

export interface SkyMapMeta {
  nside: number;
  ordering: 'RING';
  frame: 'galactic';
  unit: string;
  maps: SkyMapEntry[];
  mask_file: string | null;
  mask?: { name: string; fsky: number; threshold: number | null } | null;
  source: string;
  computed?: boolean;
}

// --- cmb.json ------------------------------------------------------------------------

export interface SpectrumSeries {
  ell: number[];
  D: number[];
  /** symmetric error (or the mean of the asymmetric pair) */
  err: number[];
  err_minus?: number[];
  err_plus?: number[];
  label?: string | null;
  estimator?: string | null;
  citation?: string | null;
  file?: string | null;
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
  /** null when no sky map could be prepared */
  source: string | null;
  ell: number[];
  observed: { wmap: SpectrumSeries | null; planck: SpectrumSeries | null };
  /** D_l of our own a_lm extraction; null without a sky map */
  map_D: number[] | null;
  lcdm_D: number[];
  torus_D: ShapeSeries;
  stats: {
    /** 'planck' | 'wmap' | 'map': the input D2_obs and S12_obs refer to */
    primary?: string;
    D2_obs: number | null;
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
    [extra: string]: unknown;
  };
  ceiling: {
    d_over_D: number[];
    causal_ball: number[];
    boundary_sphere: number[];
    lss_ball: number[];
    /** the paper leaves the last entry blank: null */
    cmb_l30: (number | null)[];
    source: string;
  };
  /** null until the orientation-scanned likelihood (qchas.cmb.likelihood) has run */
  likelihood: LikelihoodRow[] | null;
  [extra: string]: unknown;
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
  computed?: boolean;
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
  /** blank cells of the [TOP] Sec. 11.8 table are null */
  exits_3d: { spectrum: string; [d: string]: number | string | null }[];
  source: string;
  [extra: string]: unknown;
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
  /** (RA, Dec) corners in degrees */
  corners_eq?: LBPair[];
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
  /** display stretch: 1st and 99.5th percentiles of the frame, in `unit` */
  p1: number;
  p995: number;
  key?: string;
  /** observation time (UTC, ISO) */
  date?: string;
}

export interface SpherexBundle {
  release: string;
  bucket: string;
  n_frames: number;
  /** band -> [lambda_min_um, lambda_max_um, n_channels] */
  bands: Record<string, [number, number, number]>;
  frames: SpherexFrame[];
  thumbnails: SpherexThumbnail[];
  note?: string;
  source?: string;
}

// --- horizon.json --------------------------------------------------------------------

export interface QesRow {
  mass: number;
  c: number;
  rb_over_rs: number;
  s_bh?: number;
  depth: number;
  lag: number;
  t_scr: number;
  s_gen_minus_s_bh: number;
  converged?: boolean;
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
  s_ref?: number[];
  max_bits?: number;
  max_step?: number;
  max_dev?: number;
  max_dev_step?: number;
  final_bits?: number;
  recovered_step?: number | null;
  recovery_threshold_bits?: number;
  computed?: boolean;
}

/** Output of qchas.horizon.evaporation.decoding_infidelity(). */
export interface DecodingRun {
  n_bh: number;
  k: number;
  old: boolean;
  dynamics: string;
  t: number | null;
  s: number[];
  /** mean over `samples` draws; values near 1e-16 (or slightly negative) are numerical zeros */
  infidelity: number[];
  /** sample standard deviation per s */
  std?: number[];
  bound: number[];
  samples?: number;
  seed?: number;
  computed?: boolean;
  [extra: string]: unknown;
}

/** The [QG] Theorem R table (Sec. 15) as printed. */
export interface DecodingPaperTable {
  n_bh: number;
  k: number;
  samples: number;
  s: number[];
  bound: number[];
  haar_old: number[];
  eth: Record<string, number[]>;
  haar_young: number[];
  source: string;
}

/** One row of the quoted [QG] Sec. 17 (Theorem K') Kerr spin-down table. */
export interface KerrSpinDownRow {
  field: string;
  a_star: number;
  f: number;
  g: number | null;
  dlna_dlnM: number | null;
  /** null where the paper prints "area grows" */
  gsl_ratio: number | null;
  area_grows: boolean;
  sign_violations: number;
  /** computed consistency check from the quoted f and g */
  dlna_dlnM_from_f_g: number | null;
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
    de_sitter: {
      entropy_nats: number;
      log10_entropy: number;
      horizon_radius_Gpc: number;
      Lambda_m2?: number;
      omega_lambda?: number;
      [extra: string]: number | undefined;
    };
    min_area_lp2: number;
    window_examples: { m_star: number; width: number; beta?: number; ln_trace: number; theorem_G_v?: number; delta_S_BH: number }[];
    window_mass_planck?: number;
  };
  qes: {
    rows: QesRow[];
    laws: { depth_coeff?: number; depth_err?: number; lag_offset_rs?: number; source?: string; [extra: string]: unknown };
  };
  atmosphere: {
    gaps: number[];
    share: number[];
    /** null where the paper prints no value for that gap */
    paper_share: (number | null)[];
    area_coeff: number;
    /** quoted continuum values ([BH] Sec. 15: area_coeff, f1, f2, f3) plus their "source" string */
    continuum: { area_coeff?: number; f1?: number; f2?: number; f3?: number; source?: string; [extra: string]: number | string | undefined };
    n_sites?: number;
    n_ball?: number;
    ell_max?: number;
    entropy?: number;
  };
  evaporation: {
    no_diary: EvaporationRun;
    early_diary: EvaporationRun;
    late_diary: EvaporationRun;
    paper: {
      max_bits: number;
      max_dev: number;
      early_recovered_step: number;
      late_recovered_step: number;
      n_qubits?: number;
      max_step?: number;
      final_bits_with_diary?: number;
      early_diary_step?: number;
      late_diary_step?: number;
      recovery_threshold_bits?: number;
      source?: string;
      [extra: string]: number | string | undefined;
    };
    /** false when the run size differs from the paper's (fast builds) */
    paper_comparable?: boolean;
  };
  decoding: {
    haar_old: DecodingRun;
    eth: Record<string, DecodingRun>;
    haar_young: DecodingRun;
    paper: DecodingPaperTable;
  };
  kerr_spin_down: {
    fixed_point_scalar: number;
    rows: KerrSpinDownRow[];
    /** further quoted [QG] Sec. 17 numbers (Page 1976 checks, superradiant gains, spin-down) */
    facts?: Record<string, unknown>;
    /** computed: spin-down from a = 0.99 to 0.5 by integrating the quoted rates */
    spin_down_check?: { a_start: number; a_end: number; mass_fraction: number; paper: number; computed?: boolean; note?: string };
    source?: string;
    [extra: string]: unknown;
  };
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
  /** null when <a^-2>(T) has no resolved half maximum */
  fwhm_T: number | null;
  fwhm_times_amax?: number | null;
  t_bounce?: number;
  period?: number;
  norm_max_dev?: number;
  /** 'exact' (odd coherent-state combination) */
  method?: string;
  computed?: boolean;
  closed_form?: { a_min?: number; peak_a_inv2?: number };
  /** printed values where the paper tabulates this a_max (empty object otherwise) */
  paper: {
    a_min?: PaperValue;
    peak_a_inv2?: PaperValue;
    fwhm_T?: PaperValue;
    fwhm_times_amax?: PaperValue;
    /** <a>_min printed in Theorem Q(iii) / Theorem N text */
    a_min_theorem_Q?: PaperValue;
    a_min_theorem_N?: PaperValue;
    source?: string;
    [extra: string]: PaperValue | undefined;
  };
  /** the same run on the paper's finite-difference grid ([QG] Sec. 21) */
  paper_grid?: {
    a_max: number;
    grid: { L: number; n: number };
    t_bounce: number;
    bounce_delay: number;
    peak_a_inv2: number;
    a_min: number;
    fwhm_T: number | null;
    fwhm_times_amax: number | null;
    levels_used?: number;
    computed?: boolean;
    method?: string;
  };
  [extra: string]: unknown;
}

/** |Psi(a, T)|^2 at three times around the bounce ([QG] Fig. 1 right). */
export interface BounceSnapshots {
  a_max?: number;
  a: number[];
  T: number[];
  psi2: number[][];
  grid?: { L: number; n: number };
  method?: string;
  computed?: boolean;
  source?: string;
}

export interface BounceSpectrumCheck {
  j: number[];
  E_numeric: number[];
  E_exact: number[];
  abs_error?: number[];
  grid?: { L: number; n: number; h?: number };
  source?: string;
}

export interface MixmasterRun {
  run: string;
  a_min: number;
  a_at_pi: number;
  max_a_inv2: number;
  max_p_small: number;
  beta2_at_pi: number;
  norm_drift: number;
}

export interface BounceBundle {
  runs: BounceRun[];
  spectrum_check: BounceSpectrumCheck;
  spectrum_check_evolution_grid?: BounceSpectrumCheck;
  snapshots: BounceSnapshots;
  snapshots_paper_grid?: BounceSnapshots;
  snapshots_a_max_2?: BounceSnapshots;
  graviton: {
    n: number[];
    h2_peak: number[];
    occupation?: number;
    a_max?: number;
    /** a_max -> peak <h_n^2> for each n */
    h2_peak_by_a_max?: Record<string, number[]>;
    /** printed peak <a^-2> of Theorem Q'(ii) by a_max */
    paper_a_inv2_peak?: Record<string, number>;
    formula?: string;
    source?: string;
  };
  /** large-a_max limits of the exact solution */
  limits?: { a_min: number; fwhm_times_amax: number; paper_fwhm_times_amax?: number; computed?: boolean };
  /** quoted [QG] Theorem N / Theorem H tables */
  mixmaster: {
    a_min: number;
    source: string;
    a_min_range?: [number, number];
    a_mean_at_pi?: number;
    a_max?: number;
    norm_drift?: number;
    p_below_1e-3_max?: number;
    a_inv2_bound?: number;
    runs?: MixmasterRun[];
    saturation?: { channels: number; peak_H_beta: number; bound: number; peak_a_inv2: number; a_min: number; max_p_small: number; E: number }[];
    hardy_threshold?: number;
  };
  /** quoted [QG] Sec. 19.3 */
  planck_bounce_test: {
    dchi2_cutoff_max: number;
    kc_excluded_above: number;
    source: string;
    data?: string;
    chi2_no_feature?: number;
    templates?: { name: string; form: string; extra_params: number; dchi2_max: number; p_chance: number }[];
    kc_best_Mpc?: number;
    kc_exclusion_dchi2?: number;
    n_null_skies?: number;
    verdict?: string;
  };
  units: { ell_over_lP: number; a_min_lP?: number; definition?: string; source?: string; [extra: string]: number | string | undefined };
  notes?: string[];
  /** [QG] status rows for the bounce results (optional) */
  claims?: Claim[];
  [extra: string]: unknown;
}

// --- gw.json -------------------------------------------------------------------------

export interface GwEvent {
  name: string;
  source: string;
  /** GPS time of the event; strain.t is in seconds relative to it */
  tevent?: number;
  m1: number;
  m2: number;
  a1: number | null;
  a2: number | null;
  /** 'detector' (redshifted) or 'source' frame for m1, m2, mf */
  mass_frame?: string;
  approximant?: string;
  mf: number;
  af: number;
  mf_source?: number | null;
  redshift?: number | null;
  /** frame in which S_BH, T_H and the lifetimes are evaluated */
  thermo_frame?: string;
  f220_Hz: number;
  tau220_ms: number;
  echo_delay_ms: number;
  S_BH_nats?: number;
  S_BH_log10: number;
  S_BH_bits_log10?: number;
  /** log10 of the echo-power bound |R|^2 <= e^{-S} ([QG] Sec. 18) */
  echo_bound_log10: number | null;
  T_H_K: number;
  r_plus_km?: number;
  t_evap_log10_yr: number;
  page_time_log10_yr: number;
  /** the emission model behind t_evap (massless scalar, Schwarzschild) */
  lifetime_model?: string;
  strain: {
    /** seconds relative to tevent */
    t: number[];
    H1: number[];
    L1: number[];
    fs: number;
    band: [number, number];
    fs_export?: number;
    units?: string;
    /** positive: L1 later than H1 */
    lag_ms?: number;
    /** -1: L1 inverted relative to H1 */
    lag_sign?: number;
    lag_corr?: number;
    lag_convention?: string;
    lag_window_s?: [number, number];
    [extra: string]: unknown;
  } | null;
  ringdown_fit: {
    f_Hz: number;
    tau_ms: number;
    t0_ms_after_peak: number;
    f_err_Hz?: number;
    tau_err_ms?: number;
    f_sys_Hz?: number;
    tau_sys_ms?: number;
    f_stat_Hz?: number | null;
    tau_stat_ms?: number | null;
    snr?: number | null;
    p_noise?: number | null;
    sigma?: number | null;
    significant?: boolean;
    t0_range_ms?: [number, number];
    window_ms?: number;
    per_start?: { t0_ms: number[]; f_Hz: number[]; tau_ms: number[] };
    n_noise?: number;
    t_peak_gps?: number;
    /** band of the ringdown fit (wider than strain.band) */
    band?: [number, number];
    kerr_f220_Hz?: number;
    kerr_tau220_ms?: number;
    f_pull?: number | null;
    tau_pull?: number | null;
    [extra: string]: unknown;
  } | null;
  [extra: string]: unknown;
}

/** One row of the [QG] Sec. 19.1 nine-event echo table, with our recomputed f220 and delay. */
export interface GwPaperEvent {
  event: string;
  snr: number;
  mf: number;
  af: number;
  f220: number;
  dt_ms: number;
  gamma?: number;
  echo_snr?: number;
  p?: number;
  /** 90 % upper limit on R; null where the paper prints "> 1" (see R90_label) */
  R90: number | null;
  R90_label?: string;
  f220_ours?: number;
  dt_ms_ours?: number;
}

export interface GwBundle {
  events: GwEvent[];
  paper_events: GwPaperEvent[];
  /** GWOSC catalogue rows; null when the build host could not reach GWOSC */
  gwtc: (Record<string, unknown> & { possible_neutron_star?: boolean })[] | null;
  combined_R90: number;
  paper_echo_summary?: {
    combined_R90?: number;
    sum_snr2?: number;
    p_combined?: number;
    p_fisher?: number;
    A90_gamma09?: number;
    dt_gw150914_ms?: { ours: number; abedi: number };
    masses?: string;
    source?: string;
  };
  /** quoted [QG] Sec. 19.2 ringdown table; pairs are [value, uncertainty] */
  paper_ringdown?: { rows: Record<string, string | number | [number, number]>[]; source?: string };
  echo_bound_check?: { m_sun: number; log10_R2: number; log10_R2_ours?: number; source?: string; computed?: boolean };
  source?: string;
  [extra: string]: unknown;
}

// --- validation.json (written by qchas.validation; not in the contract yet) ------------

export interface ValidationRow {
  group: string;
  quantity: string;
  ours: number;
  paper: number;
  rel_diff: number;
  verdict: 'match' | 'differs';
  source: string;
  note: string;
}

export interface ValidationBundle {
  rows: ValidationRow[];
  summary: { total: number; match: number; differs: number };
  skipped: string[];
  note: string;
  computed?: boolean;
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
  'validation.json': ValidationBundle;
}

export type BundleName = keyof BundleMap;

/** The bundles of the docs/ARCHITECTURE.md data contract (validation.json is an optional extra). */
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
  'validation.json': ['rows'],
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
      if (!r.ok) {
        if (r.reason === 'network') cache.delete(name);
        return r;
      }
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
    })().then((r) => {
      if (!r.ok && r.reason === 'network') cache.delete(key);
      return r;
    });
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
