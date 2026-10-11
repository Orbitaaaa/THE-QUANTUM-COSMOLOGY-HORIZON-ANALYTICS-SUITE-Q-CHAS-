/**
 * All copy of the show in one place: act titles, beat captions, honesty chips, statements.
 * Every scientific statement carries a tag (docs/SHOW.md): ESTABLISHED, PROVED, COMPUTED,
 * HYPOTHESIS, ILLUSTRATIVE. The binding facts of SHOW.md are never contradicted here.
 */

export type Tag = 'ESTABLISHED' | 'PROVED' | 'COMPUTED' | 'HYPOTHESIS' | 'ILLUSTRATIVE' | 'DATA';

export interface Chip {
  tag: Tag;
  text: string;
  /** longer explanation, shown on hover/focus */
  more?: string;
}

export interface Beat {
  /** act time in seconds at which the beat starts */
  at: number;
  title: string;
  body: string;
  chips?: Chip[];
  /** a large one-sentence statement instead of a caption */
  statement?: boolean;
}

export type ActId = 'discovery' | 'experiment' | 'breakthrough' | 'bounce' | 'evidence';

export interface ActCopy {
  id: ActId;
  number: string;
  nav: string;
  title: string;
  /** scripted length used by the guided tour (seconds) */
  duration: number;
  beats: Beat[];
  /** always-visible truth panel */
  truth: Chip[];
}

export const TAG_INFO: Record<Tag, string> = {
  ESTABLISHED: 'Textbook physics or a measurement',
  PROVED: 'A theorem proved in the papers',
  COMPUTED: 'Computed by Q-CHAS; reproduces the papers’ numbers',
  HYPOTHESIS: 'A principle or assumption the result rests on',
  ILLUSTRATIVE: 'A visual exaggeration of scale',
  DATA: 'What the real data say',
};

export const BRAND = {
  name: 'Q-CHAS',
  long: 'Quantum Cosmology & Horizon Analytics Suite',
  event: 'NASA Space Apps Challenge',
  tagline: 'The shape of space, the fate of information and the first instant — tested against NASA’s sky.',
};

export const DISCOVERY_STATEMENT =
  'Q-CHAS turns a theorem about the shape of space into an exact map of where its seams would lie on NASA’s sky — and shows what the data can and cannot prove.';

export const EXPERIMENT_TAKEAWAY =
  'A narrow primordial line would make a torus up to 65× more visible — if we knew its orientation. Neither Planck nor WMAP shows one.';

export const ACTS: ActCopy[] = [
  {
    id: 'discovery',
    number: 'I',
    nav: 'Discovery',
    title: 'Space may be a crystal',
    duration: 19,
    truth: [
      { tag: 'DATA', text: 'Real WMAP 9-year sky (2 ≤ ℓ ≤ 30)', more: 'NASA WMAP 9-year ILC map, our extraction of the large-scale multipoles 2 ≤ ℓ ≤ 30.' },
      { tag: 'PROVED', text: '18 flat shapes → the 3-torus', more: '[TOP] Main Theorem: only the 3-torus satisfies P1–P3 (or P2 + P4).' },
      { tag: 'DATA', text: 'No seams seen: d ≥ 0.985 D', more: 'Circle searches find no matched circles; the torus must be at least 0.985 × the diameter D of the visible sky (FCC ≥ 1.04 D in the Planck likelihood).' },
    ],
    beats: [
      {
        at: 0,
        title: 'Inside the oldest light',
        body: 'The real sky measured by NASA’s WMAP: light released 13.8 billion years ago, the farthest we can see in every direction.',
        chips: [{ tag: 'ESTABLISHED', text: 'CMB: the afterglow of the Big Bang' }],
      },
      {
        at: 4.6,
        title: 'Light may wrap around',
        body: 'If space is finite and repeating — a 3-torus — light leaving one side returns through the other. We would see the same sky again, from copies of ourselves.',
        chips: [{ tag: 'ILLUSTRATIVE', text: 'size 0.8 D — already excluded by data' }],
      },
      {
        at: 8,
        title: 'Seams where the copies meet',
        body: 'Two copies of the sky overlap in an exact circle: plane x·ĝ = |g|/2, radius √(χ*² − |g|²/4). Q-CHAS computes every one.',
        chips: [{ tag: 'COMPUTED', text: 'exact circles ([TOP] Thm F)' }],
      },
      {
        at: 11,
        title: 'One shape survives',
        body: 'Of the 18 possible flat shapes of space, four physical principles leave exactly one: the 3-torus. Here its building block is a rhombic dodecahedron (FCC).',
        chips: [
          { tag: 'PROVED', text: 'Main Theorem: 18 → 1' },
          { tag: 'ESTABLISHED', text: 'P1 chirality' },
          { tag: 'HYPOTHESIS', text: 'P2 finite thermodynamics' },
          { tag: 'HYPOTHESIS', text: 'P3 local Copernican' },
          { tag: 'HYPOTHESIS', text: 'P4 supersymmetry · optional' },
        ],
      },
      {
        at: 14.4,
        title: 'The seams on our sky',
        body: 'On the real sky the seams come in matched pairs with identical temperatures. As space grows past the visible sky they shrink — and vanish at d = D.',
        chips: [{ tag: 'DATA', text: 'no matched circles in WMAP or Planck' }],
      },
      { at: 17.2, title: 'Discovery', body: DISCOVERY_STATEMENT, statement: true },
    ],
  },
  {
    id: 'experiment',
    number: 'II',
    nav: 'Experiment',
    title: 'Hunt for the seams',
    duration: 13,
    truth: [
      { tag: 'COMPUTED', text: 'Information K reproduced to ~1%', more: 'Kullback–Leibler information of the torus covariance ([TOP] Theorem P) and the rank-N resonance gains of Sec. 11.6, recomputed by Q-CHAS.' },
      { tag: 'HYPOTHESIS', text: 'A narrow primordial line', more: 'The resonance needs a narrow line in the primordial spectrum on one dual-lattice shell; nothing in standard inflation requires it.' },
      { tag: 'DATA', text: 'No resonance: p = 0.78 Planck · 0.79 WMAP', more: 'Pre-registered 160-model search: Planck 2018 (paper) S_max = 1.01, p = 0.78; our WMAP 9-yr re-analysis S_max = 0.86, p = 0.79.' },
    ],
    beats: [
      {
        at: 0,
        title: 'Two ways to look',
        body: 'Left: the real sky. Right: the only waves a torus allows — a crystal in wave-space. The standard search assumes a smooth spectrum and an unknown orientation.',
      },
      {
        at: 3.4,
        title: 'Add a resonance',
        body: 'Suppose the early universe rang at one pitch. One shell of the crystal flares, its pattern prints on the sky, and a 1.3 D torus becomes 65× more visible.',
        chips: [{ tag: 'ILLUSTRATIVE', text: 'hypothetical pattern — not in the data' }],
      },
      {
        at: 7.2,
        title: 'Same torus, smaller size',
        body: 'At 1.05 D the smooth spectrum already carries most of the signal: the line adds only ×3.',
      },
      {
        at: 10,
        title: 'The catch',
        body: 'We never know the orientation. Searching every orientation eats most of the gain — and neither Planck (p = 0.78) nor WMAP (p = 0.79) shows a line.',
      },
    ],
  },
  {
    id: 'breakthrough',
    number: 'III',
    nav: 'Black hole',
    title: 'A black hole gives its information back',
    duration: 19,
    truth: [
      { tag: 'ESTABLISHED', text: 'Black holes radiate; field theory forces S to keep rising', more: 'Hawking radiation; on a fixed background the radiation entropy keeps rising to 1.855 S_BH ([BH] Theorem B).' },
      { tag: 'COMPUTED', text: 'Exactly unitary 20-qubit evaporation follows the Page curve', more: '10.28 bits at the Page time, back to 1 bit; diaries recovered at steps 13 and 17 — reproduced ([BH] Sec. 11).' },
      { tag: 'HYPOTHESIS', text: 'Finite capacity A/4G, no independent interior (H*)', more: 'The capacity is derived from the crossed product ([BH] Theorem G); H* is proved for eternal black holes ([BH] Theorem K).' },
    ],
    beats: [
      {
        at: 0,
        title: 'A black hole gives its information back',
        body: 'An evaporating black hole bends the starlight around it. Its horizon stores A/4G: a quarter nat per Planck cell.',
        chips: [{ tag: 'ILLUSTRATIVE', text: 'Planck cells drawn enormously large' }],
      },
      {
        at: 2.4,
        title: 'Hawking pairs',
        body: 'Each quantum that escapes leaves an entangled partner inside, tied by a gold thread. Field theory alone says the threads pile up for ever.',
      },
      {
        at: 3.8,
        title: 'Throw in a diary',
        body: 'A diary (one qubit) falls in after 5 of 21 emissions. Where does its information go?',
        chips: [{ tag: 'COMPUTED', text: 'exact 20-qubit unitary run' }],
      },
      {
        at: 8.3,
        title: 'The Page time: an island ignites',
        body: 'At 47.6% of the lifetime (80.6% of the mass left) the horizon is full. An island appears and old threads re-route into the radiation itself.',
        chips: [{ tag: 'ILLUSTRATIVE', text: 'island depth; true depth 0.0399 c/S_BH of r_s' }],
      },
      {
        at: 10.9,
        title: 'The diary comes back',
        body: 'After the Page time the diary is decoded from the radiation — step 13. Thrown in late (step 14), it returns just 3 emissions later.',
      },
      {
        at: 15.2,
        title: 'Nothing is lost',
        body: 'The hole vanishes. Every thread now lives inside the radiation, and its entropy returns to zero: the Page curve.',
      },
    ],
  },
  {
    id: 'bounce',
    number: 'III·B',
    nav: 'Big Bang',
    title: 'The Big Bang bounces',
    duration: 9,
    truth: [
      { tag: 'PROVED', text: 'Unitary through the singularity', more: '[QG] Theorem Q: with Ψ(0) = 0 the Wheeler–DeWitt evolution is unitary; the universe is exactly periodic.' },
      { tag: 'COMPUTED', text: 'Bounce at ⟨a⟩ = 0.163 ℓ ≈ 0.26 ℓ_P (exact)', more: 'Exact solution; the printed table is grid-biased — both are in Evidence ([QG] Theorem Λ).' },
      { tag: 'HYPOTHESIS', text: 'Closed, radiation-filled universe; radiation as clock' },
      { tag: 'DATA', text: 'Planck shows no bounce feature (consistent)', more: 'The framework predicts no imprint (β_k = 0); Planck 2018 finds none (Δχ² ≤ 1.2).' },
    ],
    beats: [
      {
        at: 0,
        title: 'The Big Bang bounces',
        body: 'The exact quantum state |Ψ(a, T)|² of a closed universe of size a. Classically (red ghost) it crashes to zero size.',
      },
      {
        at: 3.6,
        title: 'A bounce, not a singularity',
        body: 'The quantum universe never reaches zero: its size bounces at 0.163 ℓ ≈ 0.26 Planck lengths and re-expands. Gravitons ride through unchanged.',
      },
      {
        at: 6.6,
        title: 'Unitary all the way',
        body: 'No information is lost at the Big Bang either: the evolution is unitary and the universe breathes with period π.',
      },
    ],
  },
  {
    id: 'evidence',
    number: 'IV',
    nav: 'Evidence',
    title: 'Evidence',
    duration: 5,
    truth: [],
    beats: [{ at: 0, title: 'Evidence', body: 'Real data, reproduced numbers, honest limits — and every tool behind the show.' }],
  },
];

export const actCopy = (id: ActId): ActCopy => ACTS.find((a) => a.id === id)!;

/** The Evidence landing (Act IV). */
export const EVIDENCE = {
  lead: 'Everything in the show is computed from public NASA and ESA data or from exact formulas in the three papers. Here is the ledger.',
  data: [
    { name: 'WMAP 9-year', what: 'ILC CMB map, TT spectrum, KQ75 mask', where: 'NASA LAMBDA' },
    { name: 'Planck PR3', what: 'Commander map, TT spectrum (GitHub build)', where: 'NASA/IPAC IRSA · ESA PLA' },
    { name: 'SPHEREx QR3', what: '240 Level-2 frames (pointing atlas)', where: 'NASA/IPAC IRSA S3' },
    { name: 'LIGO-Virgo', what: 'GW150914, LVT151012, GW151226, GW170104 strain', where: 'GWOSC' },
  ],
  methods: [
    'Torus covariance of the CMB from CAMB transfer functions; KL information (Theorem P)',
    'Rank-N Woodbury resonance likelihood, scanned over SO(3) by an exact Fourier series',
    'Exact circles-in-the-sky and ghost images (Theorem F); Casimir Z(2) by Ewald sums',
    'Quantum extremal surface, Srednicki atmosphere, exactly unitary 20-qubit evaporation',
    'Exact Wheeler–DeWitt bounce (odd coherent states) and graviton oscillators',
  ],
  results: [
    { k: 'Torus on the real sky', v: 'no evidence', d: 'smallest single-model p = 0.03 among 21 models — the look-elsewhere scatter the paper also reports (global p = 0.15)' },
    { k: 'Resonance search', v: 'p = 0.79', d: 'WMAP 9-yr, 160 models, 100 null skies (paper on Planck: p = 0.78)' },
    { k: 'Resonance gain', v: '65× (paper 64×)', d: 'at known orientation, cubic torus at 1.3 D' },
    { k: 'Page time', v: '0.476 · 80.6%', d: 'of the lifetime · of the mass left ([BH] Fig. 1)' },
    { k: 'Unitary evaporation', v: '10.28 bits', d: 'peak at step 11; diaries recovered at steps 13 and 17' },
    { k: 'Big Bang bounce', v: '0.163 ℓ', d: '≈ 0.26 ℓ_P, exact; printed table is grid-biased' },
  ],
  limits: [
    'No detection is claimed: the data show no torus and no resonance; small tori (d < 0.985 D) are excluded.',
    'A torus as large as 3 D could only be decided by a multi-field primordial line mapped in 3D (e.g. 21-cm tomography), not by the CMB alone.',
    'The island and capacity results rest on hypothesis H* (proved for eternal black holes only).',
    'The build sandbox uses WMAP; the GitHub Pages build runs the same pipeline on Planck PR3.',
    'Paper discrepancies we found are reported, not hidden: BCC Hessian 3.39504 (not 3.3953); exact bounce peak 270.7 (printed 262); E11/E16/E18 fail only P2.',
  ],
  doors: [
    { route: 'overview', code: '0.0', label: 'Overview' },
    { route: 'manifolds', code: '1.1', label: 'Manifold Selector' },
    { route: 'lattice', code: '1.2', label: 'Lattice Shape Optimizer' },
    { route: 'ghost-twin', code: '1.3', label: 'Ghost-Twin Tracker' },
    { route: 'resonance', code: '1.4', label: 'Resonance Shell Inspector' },
    { route: 'sky-stats', code: '1.5', label: 'Sky Statistics' },
    { route: 'bounce', code: '1.6', label: 'Quantum Big Bang' },
    { route: 'capacity', code: '2.1', label: 'Crossed-Product Capacity' },
    { route: 'page-curve', code: '2.2', label: 'Page Curve' },
    { route: 'qes', code: '2.3', label: 'QES & Atmosphere' },
    { route: 'hayden-preskill', code: '2.4', label: 'Hayden–Preskill Recovery' },
    { route: 'gw-horizons', code: '2.5', label: 'GW Horizons (LIGO)' },
    { route: 'claims', code: 'L.1', label: 'Claims Ledger' },
    { route: 'data', code: 'L.2', label: 'Data & Provenance' },
  ],
};
