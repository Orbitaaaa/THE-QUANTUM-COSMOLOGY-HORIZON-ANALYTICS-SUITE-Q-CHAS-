# Round 11 physics (G11-01 … G11-08): does Model U support a persistent daughter spacetime?

**Short answer.** Not generically, by classical evolution alone. Acceptance criterion (from the request): *a mathematically justified
continuation with a clear global causal structure*.

| Case | Continuation | Global causal structure | Criterion |
|---|---|---|---|
| Unperturbed hole (Round 6m) | exact, regular | daughter IV ∪ I′ outside J⁻(I⁺) of the parent | met (known; idealised) |
| Generic tower, D = 5, after a scalar pulse | Cauchy horizon becomes a **weak null singularity**: Tipler-weak, r finite, affine distance finite; Królak integral still growing | extension at most C⁰, not fixed by the field equations | **not met** |
| Inner-extremal tower (κ₋ = 0, tuned mass), D = 5 | regular; rays cross r₋ at finite v; no mass inflation (converged, v ≤ 36 l) | daughter reached by the field equations within the evolved range | **conditionally met**: measure-zero family, detuned by accretion; semiclassical instability in the literature |
| Generic tower, D = 10 | converged only to v ≈ 14.5 l: mass ×3.2, ψα → 0.990, curvature bounded (K l⁴ ≤ 400) | undecided beyond v ≈ 14.5 l | **open** |

Data: `DATA/round11_physics.json`, `DATA/round11_greybody.json`.
Figures: `FIGURES/round11_physics.png`, `FIGURES/round11_causal.png` (schematic).
Code: `babyuniverse/strength11.py`, `stable11.py`, `hysteresis11.py`, `greybody11.py`.
Scripts: `scripts/round11_physics.py`, `round11_greybody.py`, `round11_causal.py`.

---

## G11-01 Singularity strength

**Set-up.** Hybrid double-null evolution (R208), with Hayward tower α = 1, D = 5, μ₀ = 10 (ψ₋α = 0.887).
* Scalar pulse of amplitude 0.05. A variant adds a p = 3 power-law tail (amplitude 0.01) as a Price-type proxy.
* Two resolutions: 2000×8000 and 4000×16000.
* Outgoing interior ray at 70 % of the interior u-range.
* Affine parameter dλ = e^{2σ}dv. The constraint h′(r_vv − 2σ_v r_v) = −rφ_v² becomes

  r_λλ = −r R̃/(D−2), with R̃ ≡ (D−2) φ_λ²/h′(ψ)

  (R̃ = R_kk in the GR limit h′ = 1).
* Królak integral I_K = ∫R̃ dλ; Tipler integral I_T = ∫dλ′∫R̃ dλ″.

**Theorem 11-1 (area control).** Suppose r ∈ [r₁, r₂] with r₁ > 0 on [0, λ*), λ* < ∞. Then

  r₁ I_T(λ)/(D−2) ≤ r(0) + r_λ(0)λ − r(λ) ≤ r₂ I_T(λ)/(D−2).

So I_T(λ*) is finite iff r has a finite limit with finite total focusing. Likewise r_λ(λ*) is finite iff I_K(λ*) is finite.
*Proof:* integrate the focusing equation twice; R̃ ≥ 0 because h′ > 0.

**Proposition-N 11-1 (deformationally weak null singularity).** Results are at v = 26, for two resolutions and two data sets.
* λ* = 1.028 → 1.036 (finite; 0.8 % resolution change).
* r → 1.05210 / 1.05211 (pulse) and 1.05223 / 1.05224 (tail).
* I_T → 0.03263 / 0.03269 (pulse) and 0.03223 / 0.03230 (tail).
* I_K = 413 / 409 (pulse) and 520 / 515 (tail). Still growing, but decelerating.

By Theorem 11-1, finite I_T and finite r mean the singularity is **weak in Tipler's sense**: extended bodies are not crushed to zero
volume. The growth of I_K means r_λ may diverge at λ*. The data cannot decide whether I_K converges (Królak strength undetermined).

**Metric extension (Assessment).** Finite λ*, finite r and continuous transverse geometry give the Ori-type picture. A C⁰ extension
across the Cauchy horizon is plausible, as in GR (Dafermos–Luk for Kerr). If I_K diverges, the extension is not C¹, and the
Christoffel symbols are not locally square-integrable. The field equations then do not determine what lies beyond. Any "daughter" is
fixed by a prescription, not by Model U's dynamics. In higher-curvature gravity this is worse, because the equations involve more
derivatives of the metric.

**Classification (Proposition-N 11-2).** Six runs:
* **D = 5, ψ₋α = 0.887 / 0.958 / 0.990:** curvature growth ×6×10¹¹ / ×2×10¹⁰ / ×2×10¹¹ → *singular growth* for all three masses.
* **D = 10, ψ₋α = 0.896 / 0.965 / 0.990:** growth ×20 / ×2.2 / ×140; inner apparent horizon crossed at v = 13.7 / 10.4 / 8.3 →
  *regular core* within the evolved range.

---

## G11-02 Stable inner-horizon geometry

**Construction.** An inner-extremal tower (κ₋ = 0) following Di Filippo–Kolář–Kubizňák
([arXiv:2404.07058](https://arxiv.org/abs/2404.07058)).
* f(x) built directly in x = r²: a triple root at x₋ = 1.5 and a simple root at x₊ = 9, so f changes sign at r₋ with f′(r₋) = 0.
* Mass μ = 3x₋ + x₊ + p₃ = 16, so that h′ → 1 in the GR limit.
* ψ(μ/r⁴) is monotone and h′ > 0 (checks in `stable11.ie_checks`; error F69).

**Proposition-N 11-3.** Same pulse, two resolutions (2000×8000 and 4000×16000), v ≤ 36 l:
* μ grows only 16 → 16.17 (+1.0 %): no mass inflation;
* K l⁴ ≤ 106.8 / 106.0;
* the interior ray crosses r₋ at v = 27.36 / 27.38 and continues into the core (r = 1.203 < r₋ = 1.225).

**Proposition 11-4 (Ori-model argument).** In the Ori / Poisson–Israel model the inner mass function obeys m ∝ (outflux) × (influx tail) × e^{κ₋v}.
* With κ₋ = 0 and an influx tail ∝ v^{−p}, p > 1, m stays bounded.
* The pointwise blueshift is the same for every angular mode ℓ, so the argument does not depend on spherical symmetry at linear order.

**Caveats (Assessment).**
1. **Measure zero.** κ₋ = 0 holds at one mass for a given tower. Accretion detunes it: the run's own +1 % mass increase splits the
   triple root. For a generic perturbation δ the inner surface gravity becomes small but non-zero, so the blueshift returns on a long
   timescale ~1/κ₋. We have not evolved that long.
2. **Semiclassical instability.** McMaken ([arXiv:2303.03562](https://arxiv.org/abs/2303.03562), PRD 107, 125023) finds that an
   inner-horizon instability remains for generic inner-extremal regular black holes formed in collapse. Our classical result does not
   address that.
3. **Alternative towers.** Dymnikova shares Hayward's κ₋ ≠ 0 behaviour (R209). No other generic stable core was found.

---

## G11-03 D = 10 nonlinear evolution

**Proposition-N 11-5.** Run parameters: μ₀ = 3μ_ext (32.5), amplitude 0.02, 4000×16000 and 8000×32000 (refined event-horizon location).
* **Converged range:** through v ≈ 14.5 l.
  * μ: 32.5 → 102.6 / 102.7 (×3.15).
  * ψα: 0.9899 / 0.9899.
  * I_K: 4.773 / 4.767; I_T: 0.1268 / 0.1264.
  * K l⁴ between 150 and 400. This is the de Sitter-core scale; 2D(D−1) = 180.
* **Beyond v ≈ 14.5:** the two resolutions disagree.
  * Curvature spikes to 3×10⁸ at v = 14.8 (fine grid) and to 1.7×10⁵ at v = 16 (coarse grid, then back to ~180).
  * Mass jumps to 350 (fine) and 250 (coarse).
  * This is not converged, and we do not interpret it.

**Conjecture 11-1 (curvature-saturated mass inflation, D = 10).** Mass inflation proceeds (×3 by v ≈ 14), but ψ saturates near 1/α
and the curvature stays at the core scale. Status: supported only in the converged range. The late-time behaviour, and whether the
spikes are physical, need a compactified-v code with adaptive mesh refinement (P-R11-a).

**Why D = 5 and D = 10 differ (Assessment).** In D = 10 the static hole has a much larger μ_ext and a thinner trapped shell
(ψ₊α ≈ 0.44 at 3μ_ext). The evolved range reaches the inner apparent horizon before the blueshift has amplified the flux by many
e-folds. This is a statement about the evolved range, not a proof of stability.

---

## G11-04 Post-bounce geometry and causal future

See `FIGURES/round11_causal.png` (schematic conformal diagrams).

* **Proposition 6m-3 (Round 6m, unchanged).** In the unperturbed geometry the daughter IV ∪ I′ lies outside J⁻(I⁺) of the parent and
  has its own null infinity. Whatever happens beyond r₋ cannot signal the parent.

  **Corollary.** Every collider observable (G11-09 … G11-16) is an exterior observable: production, Hawking emission, a stable
  remnant. None of them tests the daughter's existence directly. The internal interpretation can be tested only through the theory
  (the tower) that the exterior fixes. The two layers are kept separate, as required.
* **Generic towers (D = 5).** The causal future of the parent's black-hole region ends at the Cauchy horizon, which becomes a weak
  null singularity (G11-01). The classical theory gives no unique continuation. *Model U does not, as currently formulated, predict a
  daughter for generic data.*
* **Inner-extremal towers.** Within the evolved range the solution is regular across r₋. Rays enter the untrapped core at finite v,
  as in panel (c) of the figure. This is the only case where the field equations themselves carry the solution towards the daughter,
  subject to the G11-02 caveats.
* **D = 10.** Open: regular through the inner apparent horizon in the converged range; later evolution unresolved.

**Verdict on the primary objective.** The physically persistent daughter is established only for (i) the idealised unperturbed
solution and (ii) a tuned, measure-zero inner-extremal family within a finite evolution time. For generic data the answer in D = 5 is
negative (classically undetermined continuation). In D = 10 it is undecided.

---

## G11-05 Field-theory consistency

* **Spherical sector.** The double-null system is a well-posed characteristic problem while h′(ψ) > 0. For Hayward
  h′ = 1/(1−ψ)² > 0 for ψ < 1/α, and the hybrid solver keeps ψ < 1/α (R208). Theorem-level for the reduced system.
* **Off spherical symmetry.** Quasi-topological theories have fourth-order field equations for general metrics (they reduce to second
  order on spherical/static ansätze; Myers–Robinson, [arXiv:1004.2055](https://arxiv.org/abs/1004.2055)). The tower contains
  Gauss–Bonnet with coefficient l² (R210). Known results:
  * Gauss–Bonnet can lose hyperbolicity near small black holes (Reall–Tanahashi–Way,
    [arXiv:1406.3379](https://arxiv.org/abs/1406.3379)).
  * Lovelock theories are not strongly hyperbolic in harmonic gauge (Papallo–Reall,
    [arXiv:1705.04370](https://arxiv.org/abs/1705.04370)).
  * They admit a well-posed formulation at weak coupling in modified harmonic gauge (Kovács–Reall, arXiv:2003.08398).
  * For the full quasi-topological tower no well-posed formulation is known. **Status: open (P-R11-b).**
* **EFT control (Proposition-N, computed).** On the outer horizon of collider holes (D = 10):
  * Hayward ψ₊α = 0.75 (μ = 1.01 μ_ext) down to 0.29 (10 μ_ext);
  * Dymnikova 0.90 down to 0.38.

  The tower is O(1) at the horizon, so no truncation is controlled. The model must be treated as a complete (resummed) theory. CEMZ
  causality then requires new higher-spin states at m ≲ 1/l (R210).
* **Modes and strong coupling (Assessment, to verify against Bueno–Cano–Hennigar).** On maximally symmetric backgrounds,
  quasi-topological theories propagate only the massless graviton, with an effective coupling ∝ 1/h′(ψ̄). In the core (ψ → 1/α)
  h′ → ∞: the geometry becomes stiff (matter back-reaction suppressed), not strongly coupled. The same factor appears in our
  constraints.

---

## G11-06 Rotation and nonspherical perturbations

**Theorem 11-2 (exterior test fields).** For the Hayward tower,

  f′(r) = 2μ r (r^{D−1}(D−3)/2 − μ) / (r^{D−1}+μ)²,

so f′ > 0 for r > r_+. Hence each of the following potentials is a sum of non-negative terms on (r_+, ∞):
* brane scalar, V = f[l(l+1)/r² + f′/r];
* brane Maxwell field, V = f l(l+1)/r²;
* bulk scalar, V = f[l(l+D−3)/r² + (D−2)(D−4)f/4r² + (D−2)f′/2r].

The brane Dirac field has partner potentials W² ± W′, so its operator is A†A ≥ 0. None of these test fields has an exterior mode
instability.
* Numerical check (`greybody11.exterior_positivity`): no negative value for Hayward or Dymnikova at μ/μ_ext = 1.01 … 100.
* Literature: quasinormal modes and late-time tails of D-dimensional Hayward/Bardeen holes from an infinite tower
  (Konoplya–Zhidenko, [arXiv:2403.07848](https://arxiv.org/abs/2403.07848)).

**Not done.** Gravitational nonspherical perturbations of the tower, which need its fourth-order perturbation equations.

**Interior.** The inner-horizon blueshift does not depend on ℓ (Proposition 11-4). Nonspherical modes are expected to inflate like the
spherical one for κ₋ ≠ 0. This is an Assessment, consistent with GR.

**Rotation.** No rotating solution of the quasi-topological tower is known to us. The collider analysis uses the Myers–Perry proxy
(R212) and the "U-spin" signal variant as an envelope. Kerr-like inner horizons are expected to share the instability (Assessment).
**Status: open (P-R11-d).**

---

## G11-07 Unifying formation, bounce and evaporation

**Current chain (one parameter set, explicit handoffs):**
1. Formation and interior: double-null code (R208, Round 11).
2. Bounce and daughter: closed-FRW interior (Propositions 6m-2/3).
3. Evaporation: model-U Vaidya with brane emission (Proposition-N 6o-1).

**New this round (G11-12 feeds G11-07).** The evaporation luminosity now uses the greybody factors of the model-U brane metric rather
than the geometric blackbody, with Standard Model degrees of freedom:

| μ/μ_ext | L / L_geo | Ṅ / Ṅ_geo | spin-1 share of quanta | neutrino share |
|---|---|---|---|---|
| 10 | 1.20 | 1.26 | 0.225 (blackbody 0.259) | 0.045 |
| 1.5 | 1.32 | 1.26 | 0.122 | 0.051 |
| 1.05 | 0.92 | 0.77 | **0.023** | 0.056 |
| (Tangherlini n = 6) | 1.11 | 1.20 | 0.263 | 0.043 |

**Proposition-N 11-6.** Near extremality (T r_+ = 0.047 at 1.05 μ_ext) brane gauge-boson emission falls to 7 % of the blackbody rate.
Scalars stay at 142 % and fermions at 101 %. Model-U evaporation ends with a soft, gauge-boson-poor tail. Proposition-N 6o-1 (T → 0,
remnant) is unchanged, because the luminosity still vanishes at μ_ext.

**Still open (P-R11-c).** A single code that evolves collapse, inner region and evaporation together. This needs a semiclassical
source (an RSET model or negative-energy influx) in the double-null code.

---

## G11-08 Nonlinear hysteresis (parallel track; does not block the objective)

**Proposition-N 11-7** (`hysteresis11.run`). Set-up: D = 10 interior, homogeneous scalar background, test inhomogeneities on S⁹
(eigenvalue ℓ(ℓ+8)/a²) for ℓ = 1, 2, 4, 8, 16, 40, initial amplitude 10⁻⁶.

| V_F/Λ_c | Cycles before de Sitter escape | Maximum amplification | Maximum energy fraction | Friedmann constraint |
|---|---|---|---|---|
| 0.5 | 11 | 5.2 (ℓ=1) | 5.8×10⁻⁹ | 2×10⁻¹² |
| 2 | 4 | 17 (ℓ=1) | 3.5×10⁻⁹ | 4×10⁻¹⁴ |
| 0.05 | 8 | 9.5×10³ (ℓ=2) | 2.3×10⁻⁸ | 1×10⁻⁷ |

* Net work per cycle is positive: the pumping.
* In the small-V_F case the work per cycle grows with the cycle number, as the bounce volume grows.

**Interpretation.**
* Hysteresis survives at the test-field level if the initial inhomogeneity is ≲ 10⁻³/amplification. That is ≲ 10⁻⁷ in the most
  amplifying case.
* High-ℓ modes are not amplified (gradient energy dominates) but carry most of the energy fraction.

**Not done.** Metric perturbations; coupling to the horizon region. For generic towers that region is itself singular in D = 5
(G11-01), which removes the boundary the hysteresis picture assumes.
