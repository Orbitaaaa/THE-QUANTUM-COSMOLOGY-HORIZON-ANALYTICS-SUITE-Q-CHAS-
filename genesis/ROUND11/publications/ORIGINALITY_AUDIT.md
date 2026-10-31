# Originality audit (G11-25 … G11-27)

**Format per item.** The published starting point, then what is new here, then how to reproduce it.

**Claims are graded:**
* **NEW**: not found in the verified literature below.
* **EXTENSION**: a known method applied to a new setting.
* **REPRODUCTION**: known result, re-derived as a check.

The search was by web search and the papers' abstracts (arXiv itself is blocked from this environment). A reviewer should repeat it with
full-text access before submission. Any NEW item below could turn out to exist in a paper we have not read.

## Verified starting points

| Ref. | What it establishes |
|---|---|
| Bueno, Cano, Hennigar, PLB 861 (2025), [arXiv:2403.04827](https://arxiv.org/abs/2403.04827) (BCH) | Infinite quasi-topological towers in D ≥ 5 have regular black holes (Hayward-, Bardeen-, Dymnikova-like) as unique spherically symmetric vacua |
| Bueno, Cano, Hennigar, Murcia, [arXiv:2412.02740](https://arxiv.org/abs/2412.02740), [arXiv:2412.02742](https://arxiv.org/abs/2412.02742) | Thin-shell collapse, Birkhoff theorem, junction conditions, bounces and white-hole explosions into a new universe |
| Bueno, Cano, Hennigar, Murcia, Vicente-Cano, [arXiv:2505.09680](https://arxiv.org/abs/2505.09680) | Oppenheimer–Snyder collapse; Friedmann equation h(Φ) ∝ ρ |
| Di Filippo, Kolář, Kubizňák, [arXiv:2404.07058](https://arxiv.org/abs/2404.07058) | Inner-extremal regular black holes from quasi-topological gravity (tuned mass) |
| Di Filippo, Kubizňák, Srinivasan, [arXiv:2604.27980](https://arxiv.org/abs/2604.27980) | No null thin shells in these theories, so inner-horizon stability is left open |
| Frolov, Zelnikov, [arXiv:2601.01861](https://arxiv.org/abs/2601.01861) | Significant mass inflation only for crossings very close to the horizon |
| Carballo-Rubio et al., [arXiv:2205.13556](https://arxiv.org/abs/2205.13556) | Inner-extremal regular black holes without classical mass inflation |
| McMaken, [arXiv:2303.03562](https://arxiv.org/abs/2303.03562) | Semiclassical instability of inner-extremal regular black holes |
| Konoplya, Zhidenko, [arXiv:2403.07848](https://arxiv.org/abs/2403.07848) | Quasinormal modes and late-time behaviour of D-dimensional regular holes from an infinite tower |
| Kanti, March-Russell, [hep-ph/0203223](https://arxiv.org/abs/hep-ph/0203223), [hep-ph/0212199](https://arxiv.org/abs/hep-ph/0212199) | Brane greybody factors for Tangherlini holes (scalars; fermions and gauge fields) |
| Sahni, Toporensky, [arXiv:1203.0395](https://arxiv.org/abs/1203.0395) | Cosmological hysteresis |
| Poisson–Israel (1990), Ori (1991); Dafermos–Luk | Mass inflation; weak null singularity; C⁰ extendibility (GR) |

## Paper A (inner-horizon instability)

| # | Claim | Grade | Starting point → new element | Reproduce |
|---|---|---|---|---|
| A1 | Double-null field equations of the tower with a scalar: constraints h′(ψ)(r_uu − 2σ_u r_u) = −rφ_u², first-law fluxes μ_u, μ_v, and σ_uv from the reduced action | EXTENSION | BCH give the static and thin-shell sectors; we derive the continuous-flux characteristic system (sympy; exact GR limit) | `doublenull10.derive_sigma_uv`; `tests/test_round10.py` |
| A2 | Hybrid (r, σ, μ, φ) scheme with ψ recovered from μ, keeping ψ < 1/α | NEW (method) | Standard characteristic codes (Hamadé–Stewart, Burko) evolve r, σ in GR → we carry μ by its exact fluxes, needed because ψ → 1/α near the Cauchy horizon | `doublenull10.run_hyb`; static-hole test to 10⁻⁷ |
| A3 | Continuous-flux mass inflation in a quasi-topological regular black hole (D = 5): K grows by > 10¹⁰ with r finite | NEW (to our knowledge) | 2604.27980 leaves stability open; 2601.01861 uses crossing-flux models → nonlinear scalar evolution | `scripts/round11_physics.py` (classification); R209 |
| A4 | Theorem 11-1: area control, I_T finite ⟺ bounded focusing, with R̃ = (D−2)φ_λ²/h′ | EXTENSION | Raychaudhuri/Tipler argument in GR → modified focusing with the tower factor h′ | One-line proof in Paper A |
| A5 | Tipler-weak null singularity (I_T = 0.0327, converged), Królak undetermined | NEW for this theory | Ori-type weakness known in GR → first demonstration in a pure-gravity regular black hole | `strength11.ray_integrals` |
| A6 | D = 10: mass ×3 with curvature at the core scale (converged to v ≈ 14.5) | NEW, limited range | — | `strength11.refined_run` |
| A7 | Inner-extremal tower evolved nonlinearly: no mass inflation for v ≤ 36 | EXTENSION | 2404.07058 constructs the static solution → dynamical test with a pulse | `stable11.run_ie` |

## Paper B (collider)

| # | Claim | Grade | Starting point → new element | Reproduce |
|---|---|---|---|---|
| B1 | Brane greybody factors on the model-U metric; gauge-boson emission suppressed to 7 % near extremality | NEW (to our knowledge for this metric) | Kanti–March-Russell for Tangherlini; Konoplya–Zhidenko QNMs for tower holes → greybody spectra of the regular metric near extremality | `greybody11.emission_tables`; `scripts/round11_greybody.py` |
| B2 | Collider events with a cold extremal remnant from the tower; exact U-vs-GR separation | NEW (application) | BlackMax/CHARYBDIS-type generators for GR holes → regular-hole generator with T → 0 | `pipeline10.signal`; `round10a.py` |
| B3 | Discovery/exclusion reach with a full uncertainty budget | NEW (application) | Standard profile likelihood (Cowan et al.) | `round11b.py`; `analysis_config.json` |
| B4 | 13 TeV ATLAS validation and reinterpretation | REPRODUCTION plus reinterpretation | ATLAS 2102.10874 | `round10.py` |

## Paper C (Round 12 draft: finite daughters inside the Cauchy development)

| # | Claim | Grade | Starting point → new element | Reproduce |
|---|---|---|---|---|
| C1 | The bounced interior lies in the maximal globally hyperbolic development up to the ingoing cone from the surface's exit point E | NEW (to our knowledge) | BCH-M-V give Oppenheimer–Snyder collapse and the extended geometry → identification of the daughter region inside the Cauchy development | Theorem 12-1; `daughter12.cdd` |
| C2 | de Sitter timing law: the Cauchy horizon is the core's de Sitter horizon (η_cone − η_b → π/2) | NEW (to our knowledge) | — | `daughter12.cdd_fast`; `tests/test_round12.py` |
| C3 | Maximum, recollapse and full cycle inside the determined region (D ≥ 7, near threshold; Dymnikova D = 10) | NEW (to our knowledge) | — | `scripts/round12.py` (maps) |

## Paper C (earlier outline)

The daughter claims are conditional on G11-02/04, so Paper C is not ready (see `PAPER_C_OUTLINE.md`).

## What we do not claim

* That Model U predicts a daughter universe for generic data. It does not, classically, in D = 5.
* That the tower is a consistent UV completion. Well-posedness and causality are open.
* That any collider observable tests the interior.
