# Round 13, priority 1: the causal-domain result

## 1. Definitions

**Setting.** A spherically symmetric solution (M, g) of Model U in D = n + 4 dimensions with a perfect-fluid source (velocity u, areal
radius R), and a spacelike hypersurface Σ carrying regular, untrapped initial data with an asymptotically flat parent region. Write:
* D⁺(Σ) for the future domain of dependence;
* I⁺ for the parent's future null infinity;
* θ = ∇·u for the expansion of the fluid congruence.

**Definition 13-0 (daughter region and its cycle).**
1. The *bounce* of a fluid worldline γ is the first point where θ changes sign from negative to positive along γ.
2. The *daughter region* is

   𝒟 = { p ∈ D⁺(Σ) \ J⁻(I⁺) : p lies on a fluid worldline to the future of that worldline's bounce }.

   *Determined* means 𝒟 ⊂ D⁺(Σ). *Hidden* means 𝒟 ∩ J⁻(I⁺) = ∅.
3. The *stages* of the cycle, along a worldline γ, are:
   * S1, bounce;
   * S2, expansion (θ > 0);
   * S3, maximum (θ changes from + to −);
   * S4, recollapse (θ < 0);
   * S5, second bounce.

   The daughter *completes stage k* if there is an open set of fluid worldlines (a neighbourhood of the centre) that reach S_k at
   points of 𝒟.
4. *Size* of a stage: the fraction of the ball's proper volume (equivalently, of its rest mass) whose worldlines reach S_k inside 𝒟.

For homogeneous dust (closed FRW patch, metric a(η)²(−dη² + dχ² + sin²χ dΩ²)), θ = (D−1)a′/a³. The stages are then simultaneous in
the conformal time η: η_b (S1), η_max (S3), η_b2 (S5).

## 2. Theorem 13-1 (domain of dependence of the bounced interior)

**Set-up.** The Oppenheimer–Snyder solution of Model U:
* a closed FRW dust patch χ ≤ χ_s;
* glued along the dust surface to the static exterior f(r) = 1 − r²ψ(μ/r^{D−1}), with μ > μ_ext;
* Σ the moment of time symmetry, with the surface at rest at R₀ > r₊ (untrapped data).

Let E be the event where the surface, after its bounce, reaches r = r₋.

**Theorem 13-1.** An interior event p = (η, χ) to the future of Σ lies in D⁺(Σ) if and only if

  η + χ < η_E + χ_s.

Every interior event to the future of the first bounce lies outside J⁻(I⁺).

**Proof.**

*(i) The quotient picture.* By spherical symmetry, causal relations reduce to the 2D quotient (t, r), with the centre χ = 0 as a
reflecting boundary. Every past causal curve from p projects into the region bounded by the two past radial null rays from p. Hence
p ∈ D⁺(Σ) if and only if both rays, continued through the centre and the surface where necessary, meet Σ.

*(ii) Interior rays.* In conformal coordinates the radial null rays are η ± χ = const.
* The past ray *toward the centre*, η − χ = η_p − χ_p, reflects at χ = 0 into η + χ = η_p − χ_p.
* The past ray *toward the surface*, η + χ = η_p + χ_p, reaches the surface at η_S = η_p + χ_p − χ_s.
* The reflected ray meets the surface earlier, at η_p − χ_p − χ_s < η_S.

So everything turns on the surface point S at conformal time η_S.

*(iii) The surface has no turning point in the core (Lemma 12-0), and its exit point is at v = ∞.*
* Lemma 12-0: f decreases on (0, r₋]. Hayward: analytic; Dymnikova: numerical.
* In ingoing Eddington–Finkelstein coordinates (v, r), which cover the parent region, the black-hole region and the core
  r < r₋: v̇ = (ṙ + √(ṙ² + f))/f.
* After the bounce, ṙ → E > 0 as r → r₋ (where E² = f(R₀)). Hence dv/dr ≃ 2/f ≃ 1/(κ₋(r₋ − r)), which diverges logarithmically:
  the surface reaches r₋ at v = +∞ in finite proper time. So E lies on N1 = {v = +∞}.
* Every surface point before E has v_S < ∞.

*(iv) The exterior continuation.* The past ray toward the surface is, future-directed, an ingoing radial null ray. Its exterior
continuation from S is the ingoing null line v = v_S.
* If v_S < ∞: traced to the past this line runs out to the parent region. There t = v − r* with r* → ∞, so for large r it lies in
  the past of the slice Σ = {t = t₀}. By continuity it crosses Σ at finite r.
* If S lies after E (beyond N1): the line enters the region beyond the Cauchy horizon. In the maximal extension that region contains
  past-inextendible causal curves ending at the parent's timelike infinity without meeting Σ. Then p ∉ D⁺(Σ).

Both statements together are η_S < η_E, i.e. η + χ < η_E + χ_s.

*(v) Hidden.* The event horizon inside the ball is the outgoing ray η − χ = η_H − χ_s, where η_H is when the surface crosses r₊
inward. Post-bounce events have η > η_b > η_H, so they lie to its future and outside J⁻(I⁺). ∎

**Remarks.**
1. *Uniqueness.* Within spherical symmetry the reduced field equations are hyperbolic while h′(ψ) > 0 (Round 11, G11-05), and the dust
   equations are hyperbolic. So the solution in D⁺(Σ) is the unique spherically symmetric development of the data. Uniqueness against
   *non-spherical* perturbations needs a well-posed formulation of the full tower, which is open (P-R11-b). This is the main
   assumption behind "determined".
2. *Perturbations.* Mass inflation grows only as v → ∞ (Round 11). A point of 𝒟 at finite distance from the cone receives influx only
   from finite v, so the blueshift there is bounded.
3. *Stages.* The events of stage k (simultaneous in η) lie in 𝒟 exactly for χ < χ_k := η_E + χ_s − η_k. The centre completes stage k
   iff η_k < η_E + χ_s.

## 3. Corollary 13-2 (certified cycles)

**Method.** Stage times come from an ODE integration at three tolerances (10⁻⁸, 10⁻¹⁰, 10⁻¹²) and an independent quadrature. The
numerical error is the larger of the two spreads. "Certified" means the margin η_E + χ_s − η_k exceeds 10× that error.

| D, tower, μ/μ_ext, χ_s | S3 maximum: χ_k (volume fraction) | S5 second bounce: χ_k (volume fraction) | Numerical error |
|---|---|---|---|
| 10, Dymnikova, 1.05, 1.4 | 0.873 (5.4 %) | **0.0566** certified (2.5×10⁻¹²) | 1×10⁻¹² |
| 10, Dymnikova, 1.01, 1.5 | 0.959 (7.7 %) | **0.1595** certified (2×10⁻⁸) | 7×10⁻¹² |
| 10, Dymnikova, 1.10, 1.4 | 0.831 (3.7 %) | **0.0083** certified (8×10⁻²⁰) | 1×10⁻¹² |
| 10, Hayward, 1.20, 1.0 | 0.384 (0.05 %) certified | not reached | 1×10⁻¹² |
| 10, Hayward, 1.05, 1.4 | 0.636 (0.46 %) certified | not reached | 7×10⁻¹² |
| 5, Hayward, 20, 0.5 | not reached (ends while expanding) | — | 6×10⁻¹¹ |

**What the certification covers.** The second bounce of the Dymnikova cases is certified at the centre. In D = 10 the proper volume is
weighted by sin⁸χ, so the region that completes the full cycle is a tiny fraction of the ball: 2.5×10⁻¹² of the volume for the first
case. The maximum is reached by a larger fraction: 5.4 % for that case, 0.05 % for Hayward 1.2/1.0.

The fluid code reproduces the exact solution: centre cone arrival 2.6565 against 2.658, second bounce 2.6018 against 2.6008. So the
cone tracer is validated.

**Full-cycle window** (Dymnikova, D = 10; `full_cycle_boundary`): for μ/μ_ext ≤ 1.10, χ_s runs from 1.26–1.39 (rising with μ) up to
about 1.55. The best centre margin is 0.24 at μ = 1.002 μ_ext. The window closes at μ ≈ 1.13 μ_ext. In D = 9 it survives only at
μ ≈ 1.005 μ_ext.
