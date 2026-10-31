# Round 11 — Does Model U support a persistent daughter? First hardware experiment and external evaluation

**Primary objective.** Establish whether Model U genuinely supports a physically persistent daughter spacetime, while preparing the
first funded detector-hardware experiment and external scientific evaluation.

## Answer

**Not generically, by classical evolution alone.** The acceptance criterion was *a mathematically justified continuation with a clear
global causal structure*.

| Case | Status |
|---|---|
| Generic D = 5 | **Fails.** The Cauchy horizon becomes a Tipler-weak null singularity (finite affine distance, finite area, converged Tipler integral; Królak integral growing). The continuation is at best C⁰ and is not determined by the field equations (Theorem 11-1, Proposition-N 11-1) |
| Inner-extremal tower (κ₋ = 0) | **Passes conditionally:** no mass inflation over v ≤ 36 l, and rays cross r₋. The family is fine-tuned and is detuned by accretion; it is semiclassically unstable in the literature |
| Generic D = 10 | **Open.** Curvature-bounded mass growth through v ≈ 14.5 l (converged); later times unresolved |
| Exterior | Every collider observable is exterior, so a discovery would test the tower but not the daughter |

## Work packages

| ID | Item | Status | Where |
|---|---|---|---|
| G11-01 | Singularity strength | Done: Tipler-weak; Królak undetermined | `ROUND11/physics/G11_PHYSICS.md`, Paper A |
| G11-02 | Stable inner-horizon geometry | Done: inner-extremal stable over the evolved range; no generic stable core | same |
| G11-03 | D = 10 nonlinear evolution | Partly done: converged to v ≈ 14.5 l | same |
| G11-04 | Post-bounce geometry, causal future | Done (verdict above); schematic diagrams | `FIGURES/round11_causal.png` |
| G11-05 | Field theory: hyperbolicity, causality, EFT | Assessment; well-posedness open | G11_PHYSICS.md |
| G11-06 | Rotation, nonspherical perturbations | Theorem 11-2 (exterior test fields); gravitational modes and rotation open | same |
| G11-07 | Unify formation, bounce, evaporation | Chain documented; greybody-corrected evaporation; single code open | same |
| G11-08 | Nonlinear hysteresis (parallel) | Test-field level done; metric perturbations open | same |
| G11-09 | MC reproducibility | Done: manifest, configuration, hashes | `DATA/round11_manifest.json`, `ROUND11/collider/` |
| G11-10 | Higher-order uncertainties | Done: scales, PDFs, K, EW estimate, signal PDF | `ROUND11/collider/G11_COLLIDER.md` |
| G11-11 | Detector systematics | Done for tail, pile-up and scale; timing and trigger noted | same |
| G11-12 | Greybody factors | Done: computed on the model-U metric; generator patch | `babyuniverse/greybody11.py` |
| G11-13 | Spinning holes | Envelope (Myers–Perry proxy); rotating solution open | same |
| G11-14 | Remnant phenomenology | Done: velocity, ionisation, stopping, TOF | same |
| G11-15 | Alternative hypotheses | Done: separations robust to the emission model | same |
| G11-16 | LHC open-data run | **Blocked**: open-data hosts refused by the network policy | same |
| G11-17 … 24 | Timing demonstrator V1 engineering | Done: requirements, sensors, front-end, DAQ, mechanics, BOM (prices UNRESOLVED), test protocol, safety | `ROUND11/engineering/` |
| G11-25 / 26 | Papers A and B | Drafts compiled | `PAPER/paperA_inner_horizon.pdf`, `PAPER/paperB_collider.pdf` |
| G11-27 | Paper C | Outline with a gate (not ready) | `ROUND11/publications/PAPER_C_OUTLINE.md` |
| G11-28 … 30 | Reproducibility, review dossier, collaboration, originality audit | Done | `ROUND11/publications/` |
| G11-31 … 37 | Prospectus, Proposals A and B, budget, funding strategy, organisation, milestones | Done; deadlines checked 31 Oct 2026 | `ROUND11/funding/` |

## Key numbers

| Area | Result |
|---|---|
| D = 5 strength | I_T = 0.0327 (converged); r → 1.0521; λ* = 1.03; I_K ≈ 410–520 (growing) |
| Inner-extremal tower | μ +1 %; K l⁴ ≤ 106; crosses r₋ at v = 27.4 l |
| D = 10 | μ ×3.15; ψα → 0.990; K l⁴ ≤ 400 through v ≈ 14.5 l |
| Greybody | spin-1 emission at 7 % of blackbody near extremality |
| Reach | discovery 33.4 (+0.8/−1.1) TeV; exclusion 35.0 (+1.2/−1.0) TeV (FCC-hh 100 TeV, 30 ab⁻¹) |
| Remnants | β median 0.11; if charged, ~43 % stop; +41 ns at 1.5 m |
| Timing | 30–35 ps per plane achievable (calibrated model); 1 m beam test separates p/π at 2 GeV/c by 6.8σ (σ_TOF = 50 ps) |

## Data and figures

* `DATA/round11_physics.json`, `round11_greybody.json`, `round11_collider.json`, `round11_manifest.json`, `round11_timing.json`.
* `FIGURES/round11_physics.png`, `round11_collider.png`, `round11_timing.png`, `round11_causal.png`.

## Blocked and unresolved

* G11-16 open data: opendata.cern.ch, eospublic.cern.ch, hepdata.net, zenodo.org and arxiv.org are refused by the environment's
  network policy.
* All BOM prices and component quotations are unresolved.
* Funding eligibility needs a host institution and PI.
