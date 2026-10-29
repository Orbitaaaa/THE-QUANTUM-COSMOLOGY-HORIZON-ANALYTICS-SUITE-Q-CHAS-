# Round 6n — model U beyond idealised collapse

**Task.** Test whether model U's daughter mechanism survives three things that idealised collapse leaves out: pressure, inhomogeneity
and shocks.

**Answer: yes, it passes every test we can run, with one stage still unresolved.**
* **What survives everything we can run.** We tested collapse with any fluid equation of state, uneven density profiles, pressure and
  shocks. In every case:
  * a horizon forms;
  * the matter bounces inside an untrapped core at the limiting curvature;
  * curvature stays below the bound;
  * total energy is conserved to round-off.
* **No pathology appeared.** There was no singularity, no curvature blow-up and no energy violation.
* **The bounced core does not escape cleanly.** In realistic, uneven collapse it meets the still-infalling envelope at the inner horizon.
  In some runs its outer layers then enter the white-hole region, opening a throat towards the trapped envelope; that is the
  daughter starting to separate.
* **What remains unresolved.** The evolution beyond that inner-horizon collision is not resolved by our method (P-R6n-a).

Code `hydro6n.py`; script `scripts/round6n.py`; tests `tests/test_round6n.py` (4); data `DATA/round6n.json`; figure `FIGURES/round6n.png`.

## 1. The equations [Proposition 6n-0]

We use comoving coordinates, ds² = −e^{2Φ}dt² + e^{2λ}dr² + R²dΩ²_{D−2}, with U = D_τR and Γ = R'e^{−λ}. In any theory whose spherical
sector has the Lovelock structure, which includes the quasi-topological towers of model U, four relations hold:
1. Generalised Misner–Sharp relation: Γ² − U² = (∇R)² = 1 − R²ψ(μ/R^{D−1}).
2. Unified first law: μ' = (D−1)ρR^{D−2}R' and D_τμ = −(D−1)pR^{D−2}U.
3. The tr equation is h'(ψ) times its GR form, so D_τΓ = −UΓp'/((ρ+p)R').
4. Euler equation: Φ' = −p'/(ρ+p).

Together these give the acceleration equation

  D_τU = −Γ²p'/((ρ+p)R') + Rψ(y) − ((D−1)/2)Rψ'(y)(y + p).

It reduces correctly in both limits:
* to GR's −M/R² − 4πpR in D = 4;
* to the generalised Friedmann equation of Round 6l for homogeneous dust.

**Numerical scheme.** The code is Lagrangian, of May–White type.
* Energy is conserved by construction: the boundary masses are evolved by the first law, so the constraint holds exactly and the total
  energy error is about 10⁻¹⁶.
* Tensor-type artificial viscosity acts only on non-homologous compression, so the dust Oppenheimer–Snyder solution is exact in the
  scheme.
* Thin zones are merged inelastically.

## 2. Results

| Test | Result | Label |
|---|---|---|
| Exact OS bounce, D = 5 and 10 | bounce radius matched to 6×10⁻⁶ and 4×10⁻⁶; bounce time to 4×10⁻⁵ and 7×10⁻⁴ | numerical validation |
| Any perfect fluid, closed FRW | **Proposition 6n-1**: because ψ ≤ 1/α, ȧ² = a²ψ − 1 < 0 for a < l, so every closed universe bounces at a_b ≥ l. Checked for w = −0.3 to 1 in D = 5 and 10 (a_b = 1.00005 l in every case) | Proposition |
| Inhomogeneous dust (LTB) | every shell bounces at R_b = √(kα) (within 0.3–1.8 %); inner shells bounce first (bounce times spread ×1.6–2.1), so the bounced core always meets the envelope | Proposition-N |
| Pressure, inhomogeneity, shocks (6 runs) | every run: horizon, then bounce inside an untrapped core; ψ_max ≤ 0.994/α; energy error 0 | Proposition-N 6n-2 |
| Inner-horizon collision | the bounced core expands to R ≈ l (the inner horizon) and collides with the infalling envelope (density spike 10³–10⁶, curvature still bounded). In 2 of 6 runs, shells then become anti-trapped (white-hole region) next to trapped envelope shells, with Γ → 0 between them: a **throat**. In the other runs a strong shock ends the comoving slicing first; which happens depends on resolution and viscosity | Proposition-N 6n-3 (with limits) |
| Pressure threshold | for central p/ρ ≳ 0.2 (D = 5, γ = 5/4) the star never forms a trapped surface: no black hole, so no daughter is needed | numerical |

**Interpretation.**
* **The kinematic reason for the throat.** Γ² = U² + (∇R)², so a trapped shell (U < 0) next to an anti-trapped shell (U > 0) forces
  Γ → 0 between them. The slice then passes through a neck joining the parent's black-hole region to the daughter's white-hole region.
  That is the semi-closed ("bag of gold") topology: the daughter connected to the parent through a throat.
* **The comoving slicing cannot cross a neck.** It requires Γ > 0, so continuing past this point needs a characteristic (double-null)
  code. In double-null coordinates trapping only changes the sign of r_u r_v and causes no breakdown.

**What this round corrects.** Di Filippo, Kubizňák & Srinivasan (PRD 114, 044061 (2026)) show that the junction conditions of pure-gravity
quasi-topological theories forbid null thin shells in these regular black holes.
* Our Round 6l Proposition-N 6l-6b assumed such shells: it used the Dray–'t Hooft–Redmount relation for crossing null shells. Its
  derivation is therefore not valid in model U (F63).
* The question it addressed, what happens when later matter meets the bounced core at the inner horizon, is answered here with
  continuous timelike matter instead.
* The answer: the curvature stays bounded, matter piles up at R ≈ l, and the daughter's white-hole region begins to form.
* This agrees in spirit with Bonanno, Khosravi & Saueressig (2021) and Bonanno, Panassiti & Saueressig (2025): some regular cores,
  Hayward-type included, are tamed by attractors rather than destroyed by mass inflation.

## 3. Status

| Requirement | Status |
|---|---|
| Bounce survives any equation of state | yes (proved) |
| Bounce survives inhomogeneity | yes (exact LTB, every shell) |
| Survives pressure and shocks with bounded curvature and exact energy conservation | yes (6 hydrodynamic runs) |
| Matter enters the daughter's white-hole region in realistic collapse | yes in 2 of 6 runs, with a throat forming; resolution-dependent |
| Evolution after the inner-horizon collision, to the new universe | not resolved by the comoving scheme; needs a characteristic code (P-R6n-a) |
| Later matter forms a second daughter by quiet addition (6l) | that derivation is withdrawn (F63); collision dynamics as above |

## Literature (Round 6n)
* Di Filippo, Kubizňák & Srinivasan, PRD 114, 044061 (2026), arXiv:2604.27980: no null thin shells in quasi-topological regular black holes.
* Frolov & Zelnikov, PRD 113, 084007 (2026), arXiv:2601.01861: colliding null shells in quasi-topological regular black holes.
* Bonanno, Khosravi & Saueressig, PRD (2021), arXiv:2010.04226: regular black holes with stable cores.
* Bonanno, Panassiti & Saueressig, arXiv:2507.03581: Cauchy-horizon (in)stability of regular black holes.
* Carballo-Rubio et al., arXiv:2402.14913: mass inflation without Cauchy horizons.
* Di Filippo, Kolář & Kubizňák, arXiv:2404.07058: inner extremality.
* May & White (1966): Lagrangian relativistic hydrodynamics.
* Misner & Sharp (1964).
