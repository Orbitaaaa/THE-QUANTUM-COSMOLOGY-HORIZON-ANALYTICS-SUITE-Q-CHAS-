# Q-CHAS architecture and data contract

Q-CHAS has two halves joined by a set of JSON/binary bundles:

```
 public archives ──► qchas (Python analytics core) ──► web/public/data/*.json|bin ──► web/ (Vite + TS dashboard)
 (WMAP, Planck,        fetch → compute → export            the data contract below        static site, GitHub Pages
  LIGO, GWOSC,
  SPHEREx)
```

The Python package does every heavy or validated computation once, at build time, from
real data. The dashboard loads the bundles, renders them, and runs only cheap closed-form
physics live in the browser (circle geometry, Page law, QES Newton solve, lattice sums).

## Source papers (keys used in code comments)

| key | paper | PDF |
|---|---|---|
| `[TOP]` | *Selecting the shape of flat space* (27 Sep 2026) | `docs/papers/` (not committed; see README) |
| `[QG]`  | *Quantum gravity from the crossed product* (v5, 29 Sep 2026) | |
| `[BH]`  | *What field theory forbids, and what the horizon must do* (28 Sep 2026) | |

Framing rule: the code and UI state what each paper **proves**, **computes**, **assumes**
or **measures**, using the papers' own status tables ([TOP] Sec. 13, [QG] Sec. 20,
[BH] Sec. 16). Never claim a detection or a capability the papers do not claim. Notably:
the Planck sky shows **no evidence** for any torus or resonance; resonance gives up to 64x
information **at known orientation**; deciding a torus as large as 3D needs a
**multi-field primordial line mapped in 3D**, not the CMB alone.

## Python package layout (`qchas/`)

| module | owner content |
|---|---|
| `constants.py` | physical constants; numbers quoted from the papers with section refs |
| `data/registry.py`, `data/fetch.py` | dataset registry with mirrors; downloader with provenance (`data/raw/provenance.json`) |
| `data/spherex.py` | SPHEREx QR Level-2 frame footprints via S3 range reads; thumbnails |
| `data/gw.py` | LIGO strain loading (LOSC HDF5), whitening, band-passing, ringdown fit; GWOSC catalog parsing |
| `topology/manifolds.py` | E1–E18 catalogue, principles P1–P4, selection logic; spherical/hyperbolic candidates; quoted evidence tables |
| `topology/lattice.py` | lattice bases (SC/BCC/FCC/rhombohedral/general), dual lattice, Epstein zeta Z(2) by Ewald, gradient & Hessian in 5-dim shape space |
| `topology/shapeflow.py` | Theorem H.2 shape flow, LLL reduction, Haar-random initial lattices, basin statistics |
| `topology/circles.py` | Theorem F matched circles, Corollary F.1, ghost images of a source at any distance |
| `topology/vacuum.py` | Theorem B / B-loc visual: renormalized local variance on E1 vs E2 by image sums |
| `cmb/harmonics.py` | real spherical harmonics, real-basis Wigner rotation matrices for l <= 64 |
| `cmb/transfer.py` | CAMB transfer functions Delta_l(k), l <= 30 (cached `.npz`) |
| `cmb/covariance.py` | Eq. (2) torus covariance, LCDM covariance, KL divergence, Pinsker bound |
| `cmb/resonance.py` | primordial line on a dual-lattice shell; rank-N Woodbury update; KL gain; Bayes factor |
| `cmb/maps.py` | load WMAP ILC / Planck Commander, downgrade to Nside 64, beam/pixel-window, real a_lm for 2 <= l <= 30 |
| `cmb/stats.py` | D_l, quadrupole, S_1/2, spectra loaders (WMAP/Planck text files) |
| `cmb/likelihood.py` | orientation-scanned log-likelihood ratio and Bayes factor on real a_lm, null calibration |
| `horizon/page.py` | [BH] Theorems C–F closed forms |
| `horizon/capacity.py` | [BH] Theorem G capacity, [QG] Theorem C de Sitter entropy, Theorem U area gap |
| `horizon/qes.py` | [BH] Sec. 9 quantum extremal surface (validated to every printed digit) |
| `horizon/atmosphere.py` | [BH] Secs. 10/15 Srednicki lattice atmosphere (validated) |
| `horizon/evaporation.py` | [BH] Sec. 11 exact unitary evaporation; [QG] Theorem R decoding |
| `horizon/bounce.py` | [QG] Theorems Q, Q', Λ: quantum Big Bang spectral evolution |
| `horizon/kerr.py` | Kerr thermodynamics, remnant fits, BCW QNMs, echo delay |
| `export.py` | builds every bundle in the contract below |
| `cli.py` | `qchas fetch [--planck]`, `qchas build [--fast]`, `qchas info` |

Tests live in `tests/` (pytest). Anything slower than ~10 s is marked `@pytest.mark.slow`.

## Conventions

* Sky coordinates: **Galactic** (l, b) in degrees; HEALPix **RING** ordering; maps in **μK**.
* Unit vector of (l, b): `(cos b cos l, cos b sin l, sin b)`.
* Cosmology: Planck 2018 best fit (`constants.PLANCK2018`), χ* = 13 870 Mpc, D = 2χ* = 27.74 Gpc.
* Torus sizes are given as `d/D`, d = length of the shortest lattice vector.
* Lattice bases are 3x3 arrays whose **rows** are basis vectors, scaled so the shortest vector has length 1 unless stated.
* Horizon quantities in Planck units unless suffixed (`_K`, `_s`, `_years`, `_m`).
* Every number quoted from a paper is shipped with a `"source"` string such as `"[BH] Sec. 9"`; every number we compute is shipped with `"computed": true` and, where the paper prints the same quantity, the paper value next to it (`"paper": ...`).

## Data contract (`web/public/data/`)

All JSON files are UTF-8, keys in `snake_case`, arrays of numbers as plain JSON arrays.
Missing optional inputs (for example Planck when its archive is unreachable) produce
`null` fields plus an entry in `manifest.json → warnings`, never fabricated values.

### `manifest.json`
```json
{
  "generated_utc": "2026-10-09T12:00:00+00:00",
  "qchas_version": "0.1.0",
  "cmb_map_source": "WMAP9-ILC" ,
  "datasets": { "<key>": { "filename", "url", "bytes", "sha256", "retrieved_utc", "mission", "description", "citation" } },
  "bundles": ["topology.json", "..."],
  "warnings": ["Planck PR3 not reachable from build host; CMB analyses use WMAP 9-year ILC"]
}
```

### `topology.json`
```json
{
  "manifolds": [ { "label": "E1", "name": "3-torus", "closed": true, "orientable": true,
                   "holonomy": "trivial", "translations_only": true, "periodic_spin_possible": true,
                   "principles": { "P1": "allowed|excluded", "P2": "...", "P3": "...", "P4": "allowed|excluded|no spin structure" } } ],
  "principles": [ { "key": "P1", "name": "Uniform chiral physics", "status": "experimentally grounded", "statement": "...", "theorem": "Theorem A" } ],
  "curved": [ { "name": "Poincaré dodecahedral S^3/I*", "order": 120, "shortest_loop": "π/5", "omega_k_range": [-0.0104, -0.0010], "lnB_vs_S3": -0.05, "note": "...", "source": "[TOP] Sec. 9.3" } ],
  "weeks": [ { "omega_k": 0.0005, "Rc_Gpc": 198.3, "chi_over_Rc": 0.070, "shortest_loop_over_D": 4.18, "volume_ratio": 659.0, "circles": "no circles for any observer" } ],
  "evidence_torus": [ { "shape": "FCC", "d_over_D": 1.05, "max_dlnL": 7.20, "lnB": -1.20, "null95": 1.12, "p": 0.44 } ],
  "exclusion_torus": [ { "shape": "FCC", "d_over_D": 0.985, "lnB_planck": -22.14, "median_if_true": 17.99, "q05_if_true": 3.07, "p_le_planck": "<=0.016", "sims": 60, "verdict": "excluded (95%)" } ],
  "other_flat": [ { "space": "E2 half-turn", "orientable": true, "samples": 48, "lnB": "+0.20 ± 0.03", "largest": 0.89 } ],
  "claims": [ { "statement": "...", "status": "Theorem, unconditional", "source": "[TOP] Sec. 13" } ]
}
```

### `lattice.json`
```json
{
  "bases": { "SC": [[...],[...],[...]], "BCC": [...], "FCC": [...] },
  "rhombohedral": { "alpha_deg": [...], "Z": [...] },
  "critical": { "SC": { "Z": 16.5323, "Z_paper": 16.532315959761670, "hessian": [...5], "hessian_paper": [...5] }, "BCC": {...}, "FCC": {...} },
  "flows": { "n": 120, "counts": { "FCC": 95, "BCC": 25, "other": 0 }, "paper_counts": { "FCC": 190, "BCC": 50 },
             "trajectories": [ { "tau": [...], "Z": [...], "end": "FCC" } ] },
  "relaxation": { "kappa_cE_rate_fcc": 0.566, "TL_today": 1e30, "source": "[TOP] Prop. H.3" }
}
```

### `sky_map.json` + `sky_map_<name>.bin`
`sky_map.json`:
```json
{ "nside": 64, "ordering": "RING", "frame": "galactic", "unit": "uK",
  "maps": [ { "name": "full", "file": "sky_map_full.bin", "description": "WMAP 9-yr ILC at Nside 64", "min": -400, "max": 400 },
            { "name": "lowl", "file": "sky_map_lowl.bin", "description": "2 <= l <= 30 reconstruction from the extracted a_lm" } ],
  "mask_file": "sky_mask.bin", "source": "WMAP9-ILC" }
```
`.bin` files: little-endian Float32, length 12·nside², RING order. Mask is UInt8 (1 = kept).

### `cmb.json`
```json
{
  "source": "WMAP9-ILC",
  "ell": [2, ..., 30],
  "observed": { "wmap": { "ell": [...], "D": [...], "err": [...] }, "planck": null },
  "map_D": [...],                           // D_l of our own a_lm extraction (same ell)
  "lcdm_D": [...],                          // CAMB, Planck 2018 best fit
  "torus_D": { "FCC_0.985": [...], "cubic_0.985": [...], "BCC_0.985": [...], "FCC_0.70": [...] },
  "stats": { "D2_obs": 150.6, "S12_obs": 6694, "S12_lcdm_median": 34576, "p_C2": {...}, "p_S12": {...} },
  "kl": { "d_over_D": [...], "cubic": [...], "BCC": [...], "FCC": [...],
          "paper": { "cubic": {"1.0": 19.09, "1.05": 5.05, ...}, "BCC": {...}, "FCC": {...} },
          "pinsker_bound": { "cubic": [...], "BCC": [...], "FCC": [...] } },
  "ceiling": { "d_over_D": [1.05, 1.1, 1.2, 1.3, 1.5, 2], "causal_ball": [...], "boundary_sphere": [...], "lss_ball": [...], "cmb_l30": [...], "source": "[TOP] Sec. 11.2" },
  "likelihood": [ { "shape": "FCC", "d_over_D": 1.05, "max_dlnL": ..., "lnB": ..., "null_q95": ..., "p": ..., "n_null": 50,
                    "best_axes_gal": [[l, b], ...], "paper": { "max_dlnL": 7.20, "lnB": -1.20, "p": 0.44 } } ]
}
```

### `resonance.json`
```json
{
  "kl_table": [ { "shape": "cubic", "d_over_D": 1.3, "no_line": 0.12, "shell": 16, "n_vectors": 6, "f": 0.5, "kl": 7.77, "gain": 64.7,
                  "paper": { "no_line": 0.12, "kl": 7.77, "gain": 64 } } ],
  "shells": { "cubic_1.2": [ { "m2": 9, "k_chi": 2.62, "n_vectors": 30 } ] },
  "woodbury_benchmark": { "n": 957, "ranks": [6, 8, 12], "direct_ms": [...], "woodbury_ms": [...], "max_abs_diff": 1e-12 },
  "search": { "map": "WMAP9-ILC", "models": [ { "shape", "d_over_D", "shell_m2", "f", "lnB" } ], "S_max": ..., "S_avg": ..., "null": { "S_max": [...], "S_avg": [...] }, "p_max": ..., "p_avg": ..., "paper": { "S_max": 1.01, "p_max": 0.78, "S_avg": -0.03, "p_avg": 0.53 } },
  "exits_3d": [ { "spectrum": "smooth", "d1.5": 0.063, "d2": 0.006, "d3": 0.0003 } ],
  "source": "[TOP] Secs. 11.6-11.8"
}
```

### `spherex.json` + `spherex/<obsid>.png`
```json
{ "release": "qr3", "bucket": "nasa-irsa-spherex", "n_frames": 240,
  "bands": { "1": [0.75, 1.09, 39], ... },
  "frames": [ { "obsid", "detector", "date", "ra", "dec", "l", "b", "corners_gal": [[l, b] x4], "band_um": [lo, hi], "key" } ],
  "thumbnails": [ { "obsid", "detector", "file": "spherex/<obsid>_D1.png", "l", "b", "band_um", "unit", "p1", "p995" } ] }
```

### `horizon.json`
```json
{
  "page": { "t": [...], "hawking": [...], "bh_entropy": [...], "page": [...], "island_min": [...], "interior_volume": [...], "page_time": 0.476, "page_mass_fraction": 0.806, "beta": 1.855 },
  "capacity": { "ladder": [...], "de_sitter": { "entropy_nats", "log10_entropy", "horizon_radius_Gpc" }, "min_area_lp2": 5.97,
                "window_examples": [ { "m_star", "width", "ln_trace", "delta_S_BH" } ] },
  "qes": { "rows": [ { "mass", "c", "rb_over_rs", "depth", "lag", "t_scr", "s_gen_minus_s_bh", "paper": {...} } ], "laws": {...} },
  "atmosphere": { "gaps": [...], "share": [...], "paper_share": [...], "area_coeff": 0.2952, "continuum": {...} },
  "evaporation": { "no_diary": {...evaporate() output...}, "early_diary": {...}, "late_diary": {...}, "paper": { "max_bits": 10.28, "max_dev": 0.72, "early_recovered_step": 13, "late_recovered_step": 17 } },
  "decoding": { "haar_old": {...}, "eth": { "1": {...}, "3": {...}, "10": {...}, "100": {...} }, "haar_young": {...}, "paper": {...} },
  "kerr_spin_down": { "fixed_point_scalar": 0.555, "rows": [...quoted [QG] Sec. 17 table...] },
  "claims": [ { "statement", "status", "source" } ]
}
```

### `bounce.json`
```json
{
  "runs": [ { "a_max": 2, "T": [...], "a_mean": [...], "a_inv2": [...], "a_min": 0.1644, "peak_a_inv2": 262, "fwhm_T": 0.252,
              "paper": { "a_min": 0.1644, "peak_a_inv2": 262, "fwhm_T": 0.252 } } ],
  "spectrum_check": { "j": [...], "E_numeric": [...], "E_exact": [...] },
  "snapshots": { "a": [...], "T": [...], "psi2": [[...], ...] },
  "graviton": { "n": [3, 4, 5, 10], "h2_peak": [...] },
  "mixmaster": { "a_min": 0.256, "source": "[QG] Theorem N" },
  "planck_bounce_test": { "dchi2_cutoff_max": 1.2, "kc_excluded_above": 4.9e-4, "source": "[QG] Sec. 19.3" },
  "units": { "ell_over_lP": 1.596 }
}
```

### `gw.json`
```json
{
  "events": [ { "name": "GW150914", "source": "LOSC tutorial (GWOSC)", "m1", "m2", "a1", "a2", "mf", "af",
                "f220_Hz", "tau220_ms", "echo_delay_ms", "S_BH_log10", "echo_bound_log10", "T_H_K", "t_evap_log10_yr",
                "page_time_log10_yr", "strain": { "t": [...], "H1": [...], "L1": [...] , "fs": 4096, "band": [35, 350] },
                "ringdown_fit": { "f_Hz", "tau_ms", "t0_ms_after_peak" } } ],
  "paper_events": [ { "event": "GW150914", "snr": 23.8, "mf": 68.7, "af": 0.68, "f220": 248, "dt_ms": 293, "R90": 0.67 } ],
  "gwtc": null,
  "combined_R90": 0.28
}
```

## Web app layout (`web/`)

Vite + TypeScript, no framework. Dependencies: `three` (globe, 3D lattice),
`plotly.js-basic-dist-min` (charts). Structure:

```
web/
  index.html
  src/main.ts              router (hash-based), nav, theme, data preload
  src/styles.css           design tokens (dark mission-control theme + light), layout grid
  src/lib/data.ts          typed loaders for every bundle above (graceful "not built yet" state)
  src/lib/charts.ts        Plotly wrappers with the shared palette and layout
  src/lib/healpix.ts       ang2pix/pix2ang (RING), Mollweide projection helpers
  src/lib/sky.ts           Three.js globe with HEALPix texture, overlays (circles, points, lines, footprints)
  src/lib/physics/*.ts     live closed-form physics (circles, page, qes, lattice sums, capacity)
  src/lib/ui.ts            panel chrome, stat tiles, sliders, tables, provenance badges
  src/panels/<id>.ts       one file per tool; each exports `mount(el, ctx)` and `title`
```

Routes: `#/overview`, Module 1 — `#/manifolds`, `#/lattice`, `#/ghost-twin`, `#/resonance`,
`#/sky-stats`, `#/bounce`; Module 2 — `#/capacity`, `#/page-curve`, `#/qes`,
`#/hayden-preskill`, `#/gw-horizons`; plus `#/claims` and `#/data`.
