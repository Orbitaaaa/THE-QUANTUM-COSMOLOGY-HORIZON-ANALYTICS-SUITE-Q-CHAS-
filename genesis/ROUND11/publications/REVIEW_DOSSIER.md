# Independent review dossier (G11-29)

For an external reviewer: what to check, in what order, and where each claim lives.

## 1. Headline claims and their status

| Claim | Status | Evidence | Weakest point |
|---|---|---|---|
| Generic D = 5 model-U interiors end at a Tipler-weak null singularity | Proposition-N, converged at two resolutions | Paper A §3; `round11_physics.json` | Królak integral not converged; spherical symmetry; single scalar |
| Inner-extremal tower: no mass inflation within v ≤ 36 | Proposition-N | Paper A §4 | Fine-tuned; finite time; semiclassical instability (McMaken) |
| D = 10: curvature-bounded mass growth to v ≈ 14.5 | Proposition-N (range-limited) | Paper A §4 | Later times not converged |
| A persistent daughter exists generically | **Not established** | G11_PHYSICS.md (G11-04) | — |
| FCC-hh reach 33.4 (+0.8/−1.1) TeV | Proposition-N | Paper B §3 | LO backgrounds × K; parametric detector; no full simulation |
| Gauge-boson emission suppressed near extremality | Proposition-N (validated solver) | Paper B §1 | Brane fields only; no bulk gravitons |

## 2. Suggested checks, in order of value

1. **Theorem 11-1.** Re-derive the focusing equation from the constraint (5 lines).
2. **Convergence.** Run `strength11.evolve` at a third resolution (e.g. 8000×32000 in D = 5).
3. **σ_uv.** Re-derive it with `doublenull10.derive_sigma_uv` and compare with an independent derivation.
4. **Greybody solver.** Compare with any published Tangherlini brane greybody result for n = 6.
5. **Collider budget.** Rerun `round11b.py --recombine` with different exclusion choices. Check that the smoothing is not hiding a
   physical shape change (`smoothed_ratios` in the JSON).
6. **Literature.** Search for prior nonlinear evolutions of quasi-topological regular black holes with continuous fluxes. If one
   exists, A3/A5 are downgraded.

## 3. Known errors found and fixed (full list in `FAILED_APPROACHES.md`)

* **F66:** the (r, σ) formulation drifts at the Cauchy horizon.
* **F67:** linear nuisances inflated Z from 4.0 to 11.5.
* **F68 / F69 (this round):** the inner-extremal polynomial ansatz diverged; the mass normalisation was wrong (13.5 → 16).
* **F70:** the uncalibrated LGAD model was too optimistic (12 ps).
* **F71:** the D = 10 late-time spikes are not resolution-stable and are not interpreted.

## 4. Reviewer access

Repository branch `claude/bold-noether-qbcahe`; ZIP of `genesis/` attached to each round. The papers are drafts and have not been
submitted.
