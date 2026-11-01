# Round 13 (GENESIS-13) — Robust finite-daughter formation

**Answers to the four priorities**

| # | Priority | Result |
|---|---|---|
| 1 | Prove the causal-domain result | **Done.** Theorem 13-1: an interior event lies in the domain of dependence iff η + χ < η_E + χ_s. Every cycle event is certified against it, including the Dymnikova second bounce (margins 0.008–0.16 against errors ≤ 7×10⁻¹²). Definition 13-0 makes "daughter region" and its stages precise |
| 2 | Realistic matter in D = 10 | **The cycle survives pressure and smooth inhomogeneity, but needs a nearly homogeneous core.** 240 fluid runs, all forming black holes. Uniform core: 48/48 reach maximum and recollapse inside the domain, 12 complete the full cycle. Centrally peaked profiles: 75/86 end while still expanding |
| 3 | Collider formation | **The same internal evolution does not follow.** Collider holes are near threshold (median μ/μ_ext = 1.044), the right window. But their content is two massless partons: unbound, so it never bounces in the core (Lemma 13-3) |
| 4 | Paper C and external review | **Paper C v2** (6 pages): established vs new claims, literature comparison, reproducibility section. A **review brief** doc is ready for reviewers (private until the project lead shares it) |

Detector, funding and collider phenomenology moved in parallel: RCF deadline re-checked (11 Nov, provisional); collider-hole mass and
spin distributions computed.

**Files**
* Data: `DATA/round13.json`. Figure: `FIGURES/round13.png`.
* Code: `babyuniverse/robust13.py`. Script: `scripts/round13.py`. Tests: `tests/test_round13.py`.
* Proof: `ROUND13/CAUSAL_DOMAIN_PROOF.md`. Collider analysis: `ROUND13/COLLIDER_FORMATION.md`.
* Paper: `PAPER/paperC_daughter_cycle.pdf`.

---

## 1. Causal domain (`ROUND13/CAUSAL_DOMAIN_PROOF.md`)

* **Definition 13-0.** The daughter region is 𝒟 = {p ∈ D⁺(Σ) \ J⁻(I⁺) on a fluid worldline after its bounce}. The stages are bounce,
  expansion, maximum, recollapse and second bounce, defined through the fluid expansion θ. A stage is completed when a neighbourhood of
  the centre reaches it inside 𝒟; its size is the fraction of the ball's volume that does.
* **Theorem 13-1.** Interior event (η, χ) ∈ D⁺(Σ) ⟺ η + χ < η_E + χ_s.
  * *Proof:* reduce to the 2D quotient. The surface has no turning point in the core and reaches r₋ only at v = ∞, at E. Past
    ingoing null lines from surface points before E cross Σ; those from points after E do not.
  * *Hidden:* post-bounce events lie inside the event horizon.
  * *Assumption:* uniqueness is within spherical symmetry. Off it, uniqueness needs a well-posed tower, which is open.
* **Corollary 13-2 (certified).** "Certified" means margin > 10× the numerical error (three ODE tolerances plus an independent
  quadrature).
  * Dymnikova, D = 10: the second bounce is certified at the centre for μ/μ_ext = 1.01, 1.05 and 1.10, with margins 0.160, 0.057 and
    0.008. The maximum is reached by up to 7.7 % of the volume.
  * Hayward, D = 10: the maximum is certified; the second bounce is not reached.
* **Full-cycle window** (Dymnikova, D = 10): μ ≤ 1.10 μ_ext and χ_s ≈ 1.26–1.55, closing at μ ≈ 1.13 μ_ext. In D = 9 it is only a
  sliver at 1.005 μ_ext.
* **Honest size statement.** The D = 10 volume element is sin⁸χ, so the full cycle is completed only by a tiny central volume fraction
  (2.5×10⁻¹² in the reference case).

## 2. Realistic matter (Proposition-N 13-5)

**Method.** The Round 6n Lagrangian code with four additions:
1. the exit event E;
2. the Cauchy-horizon cone traced through the fluid;
3. stage detection for tracked elements, with a 2 % re-expansion filter;
4. hiddenness: only events after the surface crosses r₊ count.

**Validation.**
* Exact dust: cone arrival at the centre 2.6565 against 2.658; second bounce 2.6018 against 2.6008.
* Verdicts identical at N = 60, 120 and 200 zones.

**Survey.**
* Configurations: D = 10, four (tower, μ/μ_ext, χ_s) cases, five profiles.
* Pressure: ε₀ = 0, 0.1, 0.5, 2 (central p/ρ ≤ 0.074).
* Inhomogeneity: smooth 5 % density modes (two seeds), or none.
* Outcome: 240 runs, all form a black hole; 230 complete; 10 break down (shell crossing, mostly peaked profiles at high pressure) and
  are reported separately.

| Profile (completed) | Full cycle | Maximum + recollapse | Expands, ends | No bounce in 𝒟 |
|---|---|---|---|---|
| uniform (48) | 12 | 36 | 0 | 0 |
| cored, flat to O(x⁴) (48) | 0 | 25 | 21 | 2 |
| parabolic (48) | 0 | 25 | 23 | 0 |
| peaked (1−x²)² (40) | 1* | 0 | 36 | 3 |
| Gaussian (46) | 0 | 6 | 39 | 1 |

\*Single highest-pressure run whose other two seeds broke down; not counted.

**Findings.**
* Pressure (up to p/ρ ≈ 0.07) and smooth 5 % inhomogeneity do not destroy the cycle. With a homogeneous core it always reaches its
  maximum and recollapses inside the domain; pressure slightly helps.
* The central density gradient decides. An order-one quadratic central gradient makes the central region re-expand later than the cone
  arrives, so the daughter ends while still expanding. **The daughter needs a nearly homogeneous core.**

## 3. Collider formation (`ROUND13/COLLIDER_FORMATION.md`)

* **Kinematics** (100 TeV, n = 6, M_D = 8 TeV, c_l = 0.70): μ/μ_ext median 1.044; 80 % below 1.1; 97 % below 1.2. Spin a* median ≈ 3.5.
* **Lemma 13-3.** Γ > 1 content has U² = Γ² − f > 0 in the core: no turning point. An open-FRW interior never bounces and crunches to a
  density caustic at bounded curvature (D = 10: ρ → 5×10⁵⁵ at τ = 4.5 l).
* **Proposition 13-4.** Collider holes contain unbound, massless content and do not host the bound daughter. Their interior holds only
  null matter streaming to the Cauchy horizon. A collider discovery would test the tower (Paper B), not the daughter (Paper C).
* **Anisotropy (effective GR scaling).** Amplification over a near-threshold cycle is (a_max/a_b)^{D−1}: 2.5 for the full-cycle case,
  106 for Hayward 1.2/1.0. Small departures from sphericity survive; collision data are O(1) anisotropic.
* **What would be needed.** Bound (Γ < 1), compact, nearly homogeneous collapse of mass just above M_ext. No known channel supplies it.

## 4. Paper C and external review

* `PAPER/paperC_daughter_cycle.pdf` (v2, 6 pages) contains: the definitions, Theorem 13-1 with proof, the timing law, the certified
  table, the fluid survey, the formation-channel section, the literature comparison and the limitations.
* **Newly verified references:**
  * Frolov–Markov–Mukhanov, PRD 41, 383;
  * Chakrabarty et al., arXiv:1909.07129;
  * Kelly–Santacruz–Wilson-Ewing, arXiv:2006.09325;
  * Idrissov–Quevedo, arXiv:2609.03278;
  * arXiv:2602.04956.
* **Established vs new:** the originality audit (`ROUND11/publications/ORIGINALITY_AUDIT.md`, C1–C6) and the review brief.
* **Review brief:** a doc listing the three claims, specific questions for reviewers, the established/new table, results and
  reproduction commands. It is private until you share it.
* **Publishing:** arXiv and Zenodo are unreachable from this environment. The evidence (code, data, tests) is on the branch.
  Submission and contacting reviewers are for the project lead.

## 5. Open problems

* **P-R13-a.** A characteristic code with fluid and a regular centre, to follow throats and multi-component interception.
* **P-R13-b.** Why do peaked profiles fail? Look for an analytic criterion on the central gradient.
* **P-R13-c.** A formation channel for bound, compact, near-threshold, nearly homogeneous collapse.
* **P-R13-d.** Uniqueness off spherical symmetry (well-posedness of the tower).
* **P-R13-e.** Rotation; the anisotropic sector of the tower.
