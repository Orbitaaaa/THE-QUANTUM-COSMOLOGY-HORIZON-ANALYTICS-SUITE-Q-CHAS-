# Round 6f — formation, crunching observers, the condensate beyond its fluid limit, one history, and problems 4 and 5

Problems posed (user, after Round 6e):

1. Complete the formation of the false-vacuum ball. Run the thermal transition, radiation escape and hand-off to gravitational collapse
   in one continuous evolution.
2. Resolve the remaining generalised-entropy question, especially the entropy accounting for observers whose histories end when the
   daughter recollapses.
3. Extend the field dynamics. Go beyond the simplest condensate approximation, study non-spherical perturbations, and establish the
   full stability range.
4. Combine and cross-check everything. Show one continuous formation history in which all necessary matter and gravitational fields
   evolve self-consistently.

Then the original problems:

* **Problem 4:** establish that nature supports the proposed fields and mechanism, not merely that the equations admit them.
* **Problem 5:** make the process technologically accessible at achievable energies, masses and physical conditions.

New code:
* `formation.py`. A perfect-fluid plasma with a field-dependent free energy, the scalar, the ghost condensate, a null shell and gravity,
  in one double-null solver. It has a corner stage, a regular-centre stage and an exact hand-off between them. Helpers cover the
  restoration threshold, cone normality and centre histories.
* `condensate_eft.py`. The nonlinear 1+1 ghost-condensate EFT with its exact linearisation, and the growth exponent on measured
  histories (k and l modes). It also has the l-mode probe (scalar and Regge–Wheeler proxies).
* `accessibility.py`. Problems 4 and 5 in numbers.
* `gsl_solution.ending_observer_deficits`.
* `route_condensate` gained `clip`/`dphi0`, for the seeded exit.
* `scripts/round6f.py`, `tests/test_round6f.py` (11 tests), `DATA/round6f.json`, `FIGURES/round6f.png`.

Labels as before: **Theorem**, **Proposition**, **Proposition-N** (numerical, with convergence evidence), **Literature**, **Assessment**,
**Conjecture**. Units G = c = ħ = 1. The simulations use l_H = 30 (r_N = 17.3, S_N = 942).

> **Status (paused at the user's request).**
>
> All results below are computed. Their raw outputs and scratch scripts are in `DATA/round6f_scratch/`, and the documents are written.
>
> Still to do when the round resumes:
> 1. Run `CODE/scripts/round6f.py`, which writes `DATA/round6f.json` and `FIGURES/round6f.png`. It now saves each finished item to
>    `DATA/round6f_partial.json`; the first attempt stopped at a hand-off check, which is now guarded. Then check the numbers in
>    Section 1.3 (centre runs) and Section 5.2 against the JSON.
> 2. Run `tests/test_round6f.py` in full. Four tests passed. The EFT test failed before the linearisation fix and now passes the same
>    check by hand. The others have not been run yet.
> 3. Run the whole suite (185 + 11 tests).
> 4. Update README, and add Addendum 8 to the paper.
> 5. Send the final ZIP.

---

## Summary

| Item | Answer | Key evidence |
|---|---|---|
| 1. Thermal formation, one GR evolution | **Done; the result is negative.** Heating cannot make the persistent, normal false-vacuum ball the route starts from. Proposition-N 6e-10 is **withdrawn.** | Restoration needs E ≥ 1.74–4.7 ρ_f (6f-1). The largest hot ball that is still normal is about the target size. In full GR, with plasma, field and gravity, the false-vacuum region shrinks or is trapped before the radiation leaves (6f-2). |
| 2. Observers who end at the crunch | **Resolved (theorem).** No GSL can hold for observers whose histories end at a singularity, with or without a daughter. | Their causal-horizon generators run into the endpoint, so the area drops to zero. Measured: the whole area, already for a plain black hole (6f-3). The transient daughter is GSL-consistent for every complete observer (6e-7). |
| 3. Beyond the IR fluid; non-spherical modes; full stability range | **Done, within the EFT and test-field proxies.** | Nonlinear 1+1 EFT = its exact linearisation; caustics only at seeds 10³ above quantum (6f-4). Non-spherical proxies show no growth beyond the controls (6f-5). Stability window universal in r_N, ∝ α^(−1/2): ×83 (initial charge), ×12 (self-generated) at α = 1 (6f-6). |
| 4. One continuous history, cross-checked | **From the ball to the hidden daughter: yes, in one evolution, two solvers agree. From nothing to the ball: no** (item 1). Exit and recollapse in the route geometry are beyond the grid. | Formation solver = condensate solver on the route to 10⁻⁵ (6f-8). Seeded exit needs the spinodal after throat formation (6f-9). |
| Problem 4 (nature) | **Not established, and the laboratory corner is excluded.** | The condensate's observational bound M ≲ 100 GeV forces throats of ≳ 10⁹ M_☉ (6f-10). No ingredient is observed; none is excluded at astrophysical scales. |
| Problem 5 (technology) | **Not accessible, by about 39 orders of magnitude, with no creation route at any energy.** | Any near-Nariai throat needs E·√ρ ≥ 8.9×10⁶⁴ J^{3/2} m^{−3/2}, independent of the vacuum scale (Theorem 6f-7). The best human figure is ~10²⁶; neutron-star mergers exceed it. |

---

## 1. Thermal formation of the ball in one GR evolution

### 1.1 The solver

`formation.py` evolves, together with gravity (metric −e^{2s} du dv + r²dΩ²):
* the vacuum field ϕ, with the tilted quartic of Round 5 (false vacuum −v, true vacuum +v);
* a perfect-fluid plasma whose free energy depends on ϕ,
  F(T, ϕ) = −a_L T⁴ − a_H T⁴ exp(−kg²(ϕ+v)²/T²);
  the light species are massless in both vacua; the heavy species are massless at the false vacuum and Boltzmann-suppressed in the true
  vacuum;
* the plasma–field friction of Ignatius–Kajantie–Kurki-Suonio–Laine, with the matching source so that the total stress tensor is
  conserved;
* the ghost condensate (the IR fluid of Round 6e);
* an ingoing null shell.

The fluid's two conservation laws use a θ-weighted (Preissmann, θ = 0.6) box scheme, solved by a 2×2 Newton step for (ln T, η_f). The
scalar's potential force is implicit.

There are two stages:
* **Regular-centre stage.** α = 1 cone gauge, with ψ = rϕ and a near-axis regular form of the r-equation.
* **Corner stage.** It takes over by an exact domain decomposition at a grid point, so the centre stage never has to follow the
  formation of a black hole.

**Verification (Proposition-N).**

* **Homogeneous limit (G → 0), against the ODE for field + plasma + friction.** Second order: ϕ error 1.3×10⁻³ → 3.0×10⁻⁴ → 6.0×10⁻⁵ of
  v; T error 6×10⁻⁴ → 2×10⁻⁴ → 8×10⁻⁵ (h = 0.02, 0.01, 0.005).
* **Constraints.**
  * Corner stage: first order in the plasma sector (6.3×10⁻³ → 2.1×10⁻³ → 8.4×10⁻⁴), because the θ-weighting is first order.
  * Centre stage with gravity: relative residual 8.2×10⁻⁵ → 2.5×10⁻⁵ → 1.3×10⁻⁵; E_uu 4.3×10⁻⁴ → 1.9×10⁻⁴ → 9.9×10⁻⁵.
* **Vacuum on the centre stage.** |s_axis| ≈ 3×10⁻¹¹ over the whole run, after the axis clock is set before the correction step (F44).
* **Hand-off overlap.** Where both stages are computed, the median differences are 10⁻¹⁰–10⁻⁷. The 95th percentile of ϕ is ~10⁻³ v; the
  largest differences sit only in late, near-spinodal regions.
* **Near-axis temperature.** Errors are 2.4% within 1–2 cells of the axis, 0.4% at 5 cells and < 0.1% at 20, from the r²-weighted
  fluxes. Accepted.

### 1.2 Proposition-N 6f-1 (restoration threshold)

Heating the true vacuum moves ϕ back to the false vacuum only above a total energy density E_min. Below it, field and plasma settle into
a stall equilibrium; the field must also be pushed over the steepest part of the potential.

In the homogeneous problem (f_H = 0.999, x_c = 30, Γ = 0.3 m):

| ε (potential tilt) | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| E_min/ρ_f | 4.70 | 2.88 | 2.29 | 2.00 | 1.81 |

The friction dependence is weak: E_min/ρ_f = 1.90, 2.19 and 2.25 for Γ/m = 0.03, 1 and 3 (ε = 4). The lowest value found is 1.74 (ε = 5.2,
Γ = 0.03 m).

**Correction.** The high-temperature free energy (c_T/2)T²(ϕ+v)² used for Proposition-N 6e-10 is invalid where it was applied
(m/T ≈ 15). With the free energy above, restoration always costs well over the false-vacuum energy itself.

### 1.3 Proposition-N 6f-2 (no persistent normal ball from heating)

Two things work against the hot ball:

* **Normality.** For E just above E_min, the largest hot ball whose outgoing cone from a regular centre is still normal has radius
  R ≈ 17.5 (ε = 4) to 18.25 (ε = 5.2). The route needs a ball of radius 0.595 l_H = 17.85. Normality of null-cone data is stricter than
  the spacelike estimate 2m/r < 1.
* **Dynamics.** Continuous centre-stage evolutions (plasma, field, gravity; radiation free to leave) of the largest normal balls:
  * E = 2.05 ρ_f, R = 17. The false-vacuum radius peaks near 11, the region is trapped by t ≈ 30, and it collapses to a black hole.
  * E = 2.3 ρ_f, R = 16. The false vacuum reaches ≈ 15, but the region is trapped.
  * ε = 5.2, E = 1.79 ρ_f. The false vacuum is intermittent (radius 3.6, 0, 13.1, 0 at successive times) and trapped by t = 40.

  The runs of `scripts/round6f.py` are recorded in `DATA/round6f.json` (item1.centre_runs); see Fig. (a).

The reason is the time ordering. Radiation needs about one light-crossing time ~R to leave. The region, with total density
ρ_f + E_thermal ≳ 2.7 ρ_f, collapses on (Gρ)^{−1/2} ≲ R. The margin is O(1), not many orders, but it is the wrong sign across
ε = 1–5.2, Γ/m = 0.03–3, f_H = 0.95–0.999 and x_c = 5–30.

**Withdrawal.** Proposition-N 6e-10 ("creation by supercooling") is withdrawn. Its free energy was invalid, and its "free-streaming
emptying" was not a GR evolution. Fig. 6e(d) shows the withdrawn model. The ball the route needs has **no known creation mechanism**
(P-R6f-a).

---

## 2. Observers whose histories end at a singularity

**Theorem 6f-3.**
* **Setting.** Wall's GSL and Hawking's area theorem apply to causal horizons whose generators are future-complete: the boundary of the
  past of a future-inextendible observer.
* **Statement.** If an observer's worldline ends at a singularity at a point p (in the spherical models here, where r → 0), the boundary
  of its causal past is the past light cone of p. Every generator ends at p, and the cross-sectional area falls monotonically to zero
  along them. No generalised entropy that includes A/4G can be non-decreasing there.
* **Proof.** Raychaudhuri focusing to the vertex, with the matter entropy bounded on the cone.
* **Generality.** This happens in every gravitational collapse, with or without a daughter.

**Measured (Proposition-N).** Largest decrease of A/4G along the past boundaries of observers that end next to r → 0:
* 1001 nats for a plain black hole (the route's data with no push);
* 966 for near-Nariai collapse;
* 1245 for the transient daughter.

Compare S_N = 942: in each case the loss is the full area.

**Consequence.** The entropy accounting of the transient daughter is complete:
* The GSL holds on every causal horizon of every complete observer. The parent's event horizon has non-decreasing area (Proposition-N
  6e-7).
* Observers who end at the daughter's crunch are in the same situation as observers who fall into any black hole. Their causal "horizon"
  is not a GSL horizon.
* There is nothing left to resolve; P-R6e-b is closed.

---

## 3. Beyond the IR fluid; non-spherical perturbations; the full stability range

### 3.1 Proposition-N 6f-4 (nonlinear EFT)

The full ghost-condensate EFT along the expanding direction of the Nariai band (dS₂ × S², flat slicing, a = e^{t/r_N}) has these
ingredients:
* L = M⁴P(X) − (α/2M²)(□ϕ)², with P = (X − 1)²/2;
* the Jeans term J = 8πα(μ/r_N)²;
* an optional steady tilt source σ.

It is evolved nonlinearly, pseudo-spectrally (`evolve_1p1`).

**Exact linearisation.**

    π_tt + (1 + A′/A)π_t + [εq² + (α/2M²)q⁴ − Jq²]/A π = 0,   A = P′ + 2XP″ = 2 + 3ε.

It agrees with the nonlinear solver to 3×10⁻⁴: 17.740 against 17.745 at ε₀ = −0.08, and 1.2976 against 1.2976 at ε₀ = −0.05.

**Correction.** The mode equation of Round 6d (`ghost_condensate.mode_growth`) uses A = 2 + ε and omits A′/A. This is an O(ε) difference,
< 1% for |ε| ≤ 3×10⁻³ (the window), so R101–R111 stand.

**Large seeds.** At the window edges (μ = 724 and 3175, r_N = 10⁶):
* Up to seed 10⁻⁴ (δX up to 7|ε|), the nonlinear amplification stays within 1–5% of linear, with no caustic.
* Caustics (X → 0; the known breakdown of the ghost condensate's single-valued phase) appear only at seed 10⁻³.
* The quantum seed is ~1/(2πμ²) ≈ 3×10⁻⁷.

### 3.2 Proposition-N 6f-5 (non-spherical perturbations, test-field proxies)

The l-mode master equation of a test field is integrated on each simulated background:

    ψ_uv = (r_uv/r)ψ − (e^{2s}/4)[l(l+1)/r² − κ 6M/r³]ψ.

* κ = 0 is the scalar proxy.
* κ = 1 is the Regge–Wheeler-type odd-parity proxy, with M = m − (4π/3)U r³. This is the mass parameter without the local vacuum energy,
  as in the Schwarzschild–de Sitter Regge–Wheeler potential.

Amplitudes are measured in two ways:
* for the scalar, as ϕ = ψ/r;
* for κ = 1, as the curvature ψ/r³.

**Result.** Away from the crushing region (r > r_N/2), across l = 0–20 and all five backgrounds:
* the largest amplification is ≤ 3.6;
* tipping backgrounds and their no-tip controls agree to 0.3% (route: 3.33 against 3.32);
* the band and transient backgrounds give 1.0.

The large raw amplifications (up to 10⁷) all sit at r → 0, at the singularities, and occur equally in the controls (F46).

**Assessment.** These are test-field proxies, not the coupled even- and odd-parity system with matter (P-R6f-c).

### 3.3 Proposition-N 6f-6 (full linear stability range)

The growth exponent G = ∫ max(0, −g + √(g² − ω²)) dt (g = (1 + A′/A)/2) is evaluated on the measured charge histories of Round 6e. It
covers both kinds of mode at once:
* comoving k along the expanding dS₂ direction, which redshift;
* l along the S², which do not redshift while r ≈ r_N.

Stability is G < 1 (one e-fold). Cross-check: e^G tracks the full ODE amplification to within a factor of about 2.

| history | α = 0.1 | α = 1 | α = 10 |
|---|---|---|---|
| initial charge | 1.2×10⁻³ – 0.10 (×83) | 3.8×10⁻⁴ – 3.2×10⁻² (×83) | 1.2×10⁻⁴ – 1.0×10⁻² (×83) |
| self-generated | 8.3×10⁻³ – 0.10 (×12) | 2.6×10⁻³ – 3.2×10⁻² (×12) | 8.3×10⁻⁴ – 1.0×10⁻² (×12) |

The table gives the window in M³r_N.

* **Universality.** The windows are identical for r_N = 10⁶, 10²⁰ and 10⁴⁰. Every rate depends on M³r_N alone.
* **Scaling.** Both edges scale as α^(−1/2).
* **Comparison with Round 6e (R111).**
  * The upper edge is about twice R111's analytic 1.7×10⁻², because Hubble friction in dS₂ limits the Jeans growth.
  * The self-generated window is ×12 wide, not the "factor 3" of Round 6e.
* **Corrections.**
  * With the code's normalisation the Jeans rate is 4π√α M³, not 2π√α M³ (Proposition 6d-2). The numerical window supersedes both.

---

## 4. One history, cross-checked

### 4.1 Proposition-N 6f-8 (two independent solvers, one route)

The route of Round 6e is: a false-vacuum ball, an imploding shell, a wall, no black hole and no condensate charge on the data, with the
tilt λ_f = 10⁻⁶ making the charge. It was evolved twice:
* by `condensate.evolve`;
* by the independent formation solver (Preissmann fluid box, implicit scalar), with a negligible plasma.

The results agree:
* both produce the daughter (n_big = 19997 against 19965);
* for r_N/2 < r < 3r_N: r agrees to 4.6×10⁻⁶ (median) and 3.6×10⁻⁴ (95th percentile); ϕ to 6×10⁻³ v (95th percentile); the charge to
  6×10⁻⁸.

Including the cold plasma remnant at its floor temperature turns the far end of the data ray marginally trapped (r_v = −1.3×10⁻⁶). The
route's data are tuned to 6×10⁻⁹ there, so any extra energy on the outer ray must be retuned.

### 4.2 The continuous route from the ball to the hidden daughter

* **One evolution.** The run starts from the ball: zero charge, normal data with no black hole and no bridge. It then passes through the
  self-generated NEC phase, tipping, and an inflating daughter behind the parent's event horizon (u_EH = 0; no daughter row reaches the
  outer edge). All fields are self-consistent: the scalar, the condensate with its tilt, the shell, gravity, and the plasma remnant in
  the formation solver.
* **Proposition-N 6f-9 (exit after throat formation).** The tilt that makes the charge also brings the false vacuum to its spinodal at
  T_exit = ρ_f(16/(3ε) − 1)/λ_f. In the route, the throat forms only at clock T ≈ 125. The exit must come later, which needs
  ρ_f(16/(3ε) − 1) > λ_f T_throat with λ_f ≥ λ_th.
  * For ε = 4 this is exactly marginal: λ_th T_throat ≈ 4.5×10⁻⁵ against ρ_f/3 = 4.4×10⁻⁵. A seeded run at λ_f = 10⁻⁶ exits at T ≈ 50
    near the centre, before any throat forms, and produces no daughter.
  * For ε = 2 the exit is delayed to T = 221. The route then works: normal data (shell energy retuned to E = 7.2330 so the data ray is
    again marginal at 6×10⁻⁹), a daughter, hidden. The control with λ_f = 0 makes no daughter.
* **What the grid cannot do.** Once the daughter's false vacuum passes the spinodal (T > 221), the corner-gauge cells are far larger in
  proper time than the field's roll-off time. The stiffness parameter h²e^{2s}|U″|/16 reaches ~4×10³.
  * Rows break off a little after the spinodal (T ≈ 250–270), before the roll to the true vacuum (F45).
  * A homogeneous-FRW hand-off was attempted. It fails because surfaces of constant condensate clock are not homogeneous slices: the
    charge spread is ~100%.
  * The exit and recollapse are therefore established only in the band geometry (Proposition-N 6e-7: exit at T = 232 against 221
    predicted, then recollapse). Following them in the route geometry needs a daughter-adapted solver (P-R6f-b).

**What item 4 delivers.**
* **From the ball onward:** one self-consistent evolution to a hidden inflating daughter, with every matter field and gravity, checked by
  two solvers.
* **Exit and recollapse:** computed in the band geometry and consistent with the route's clock.
* **Before the ball:** no history exists. The ball cannot be made thermally (item 1), and no other creation mechanism is known.

---

## 5. Problem 4: does nature support the fields and the mechanism?

Computation cannot establish what nature contains. What can be done is to confront every ingredient with what is observed, using the
route's computed requirements scaled to physical units (`accessibility.py`).

### 5.1 The ingredients

| ingredient | what the route needs | status in nature |
|---|---|---|
| Ghost condensate (a Lorentz-violating, healthy NEC-violating sector) | scale M with M³r_N in the window of 6f-6 | **Not observed.** Allowed only for M ≲ 100 GeV (nonlinear analysis, Arkani-Hamed–Cheng–Luty–Mukohyama–Wiseman 2005/07), or M ≲ 10 MeV in linear theory. Its vacuum is the condensed phase (X = 0 is a ghost), so it must fill all space today, and these bounds apply. The GSL objection (Dubovsky–Sibiryakov) is answered by accretion (Mukohyama 2009). |
| Stable NEC violation at all | ε M⁴ ~ 10⁻³/r_N² inside the throat | **Hint only.** DESI DR2 BAO with CMB and supernovae prefers evolving dark energy that crosses w = −1 (2.8–4.2σ, depending on the supernova sample). That would mean effective NEC violation in the dark sector at ρ ~ (meV)⁴ on Hubble scales. It is not established and is not the throat-scale violation needed. |
| A metastable false vacuum above ours | ρ_f = 1/(8π r_N²) | **None known.** The Standard Model electroweak and QCD transitions are crossovers (Kajantie et al. 1996–98: no transition for m_H ≳ 72 GeV). A false vacuum needs new physics. Its existence is unconstrained in general; an inflaton vacuum is limited to V^{1/4} ≲ 1.5×10¹⁶ GeV by r < 0.036 (BICEP/Keck 2021). |
| The tilt λ_f (soft shift-symmetry breaking confined to the false vacuum) | λ_f r_N³ ≈ 4×10⁻⁴ | **None known** (the ghost-inflation tilt is the closest analogue). |
| A strongly first-order thermal sector | — | Not needed now: item 1 shows that thermal creation fails anyway. |

### 5.2 Proposition 6f-10 (where the route can operate if the condensate exists)

The window M³r_N ≥ K_lo (6f-6) combined with an upper bound M ≤ M_max gives r_N ≥ K_lo/M_max³. For the initial-charge window at α = 1:

| bound on M | smallest r_N | Nariai mass | vacuum scale ρ_f^{1/4} | light-crossing time |
|---|---|---|---|---|
| 100 GeV (nonlinear) | 1.1×10¹³ m (75 AU) | 2.5×10⁹ M_☉ | 6.6 keV | 10 h |
| 10 MeV (linear) | 1.1×10²⁵ m (360 Mpc) | 2.5×10²¹ M_☉ | 6.6 meV | 1.2 Gyr |

The self-generated window raises these by 2.6×10⁻³/3.8×10⁻⁴ ≈ 7.

Conversely, a throat of laboratory mass needs a much larger condensate scale:
* 1 kg needs M ≳ 2×10¹⁵ GeV;
* 10⁶ kg needs M ≳ 2×10¹³ GeV.

Both are excluded by 11–13 orders of magnitude in M.

**Assessment.**
* **The laboratory corner is excluded.** If nature has a ghost condensate at all, the route can tip only supermassive-black-hole-scale
  throats or larger.
* **A natural venue would be cosmological.** The smallest allowed throat (≈ 10⁹ M_☉) equals the horizon mass at T ≈ ρ_f^{1/4} ≈ 7 keV.
  This is not a coincidence: both are ~1/√(Gρ). The candidate is a first-order transition at keV temperatures, during the end of BBN.
* **Such a transition is constrained but not excluded.** It is limited by N_eff (ΔN_eff ≲ 0.3, Planck 2018) and BBN unless its energy
  stays in a cold dark sector. Primordial black holes of ~10⁹ M_☉ are constrained to small abundance (Carr et al. 2021).
* **No observation could tell the result apart.** From outside, a daughter-hosting throat is a black hole (Phase VI).

**Verdict on problem 4.** Not established, and **cannot be established by computation**. The route needs at least three unobserved
ingredients: a condensate, a false vacuum and a tilt. The one partial hint (DESI's phantom crossing) is at the wrong scale. The parameter
region that would be lab-accessible is excluded by the condensate bounds. What remains possible is a cosmological, supermassive-scale
version that is unobservable event by event. The concrete tests are:
* gravity in a Higgs phase (modified Newtonian potential at L ~ M_Pl/M², a Jeans-like instability on a timescale ~M_Pl²/M³);
* hidden-sector keV-scale transitions (N_eff, BBN, gravitational waves);
* the ~10⁹ M_☉ primordial black-hole population.

---

## 6. Problem 5: is the process technologically accessible?

### 6.1 Theorem 6f-7 (the invariant)

Any route that ends in a near-Nariai throat of a false vacuum of density ρ_f must assemble at least the Nariai mass M_N = r_N/3 at energy
density ρ_f = 1/(8π r_N²), inside a region of size ~r_N. Therefore

    E_req √ρ_req ≥ (1/3)(8π)^{−1/2} E_P √ρ_P ≈ 8.9×10⁶⁴ J^{3/2} m^{−3/2},

whatever the vacuum scale. Lowering ρ_f raises the mass needed as ρ_f^{−1/2}. This is the black-hole formation invariant
M√ρ ~ (c²/G)^{3/2}. It is a statement about gravity, independent of the condensate and of any creation mechanism.

* **Proof.** Nariai: Λ = 8πρ_f = 1/r_N², M_N = r_N/3.

### 6.2 Proposition 6f-11 (the gap)

| capability | energy (J) | energy density (J/m³) | shortfall against the invariant |
|---|---|---|---|
| LHC Pb–Pb central collision | 1.8×10⁻⁴ | 2×10³⁶ (12 GeV/fm³) | 10^50.5 |
| record laser focus (1.1×10²³ W/cm²) | ~10² | 3.7×10¹⁸ | 10^53.7 |
| NIF ignition shot (8.6 MJ, 2025) | 8.6×10⁶ | ~5×10¹⁶ | 10^49.7 |
| largest thermonuclear explosion (50 Mt) | 2.1×10¹⁷ | ~10¹⁷ | **10^39.1** |
| core-collapse supernova (nature) | 3×10⁴⁶ | 2.4×10³⁴ | 10^1.3 |
| binary neutron-star merger (nature) | 4.8×10⁴⁷ | 9×10³⁴ | exceeds it |

The capability rows are orders of magnitude; sources are in LITERATURE.md (Round 6f).

**Further facts.**
* **The lightest throat.** At the largest vacuum scale considered (1.5×10¹⁶ GeV, the inflation bound), the throat is light: 1 g of
  mass-energy (8.6×10¹³ J, ~20 kt TNT). It must, however, sit in a region of 2×10⁻³⁰ m at 10¹⁰² J/m³, 66 orders above the LHC. It would
  also need a condensate at M ~ 10¹⁶ GeV, excluded by problem 4.
* **With an allowed condensate.** M ≤ 100 GeV demands ≥ 2.5×10⁹ M_☉ ≈ 4.5×10⁵⁶ J, about 10³⁶ years of today's world energy use.
* **No creation route at any energy.**
  * Classically, heating fails (item 1).
  * Farhi–Guth forbids NEC-respecting creation of an inflating region from regular data. The route needs only a normal ball, but no
    mechanism making one is known.
  * The quantum route is bounded by the project's frontier P ≤ 10⁻²³ (R42) and in practice is exponentially smaller.

**Verdict on problem 5.** Not accessible. The requirement is of compact-object class. Neutron-star mergers and supernovae reach it;
human technology falls short by about 39 orders of magnitude in the only combination that matters (energy × √energy density). Even
unlimited energy would not help, because no known process builds the required ball. This is not a statement about present technology
that a better machine could change. The invariant is fixed by G and c, and the creation step is missing in principle, not in
engineering.

---

## 7. What this does and does not establish

* **Established.**
  * Within classical GR plus the ghost-condensate EFT, from a normal false-vacuum ball onward: one self-consistent evolution to a hidden
    daughter, cross-checked by two solvers.
  * A full linear stability window that is universal in r_N, with the nonlinear EFT and non-spherical test-field proxies showing no new
    instability.
  * A complete GSL accounting: complete observers satisfy it, and crunching observers cannot satisfy any GSL, as in every collapse.
* **Not established, with the reason.**
  * The creation of the ball: thermal creation fails (6f-2), and nothing else is known.
  * The exit and recollapse in the route geometry: grid stiffness (F45).
  * The coupled non-spherical perturbation system.
  * That nature contains the ingredients: unobserved; the lab corner is excluded (6f-10).
  * Technological accessibility: excluded by the invariant (6f-7, 6f-11).
* **Problem 6** (verify formation experimentally) inherits problem 5. Verification also faces Phase VI: daughter formation is invisible
  event by event.

## Literature used (Round 6f)

All of these were checked by web search, through secondary sources or abstracts where marked in LITERATURE.md. They cover:
* ghost-condensate bounds;
* NEC and superluminality;
* the second law with Lorentz violation;
* DESI DR2;
* the electroweak crossover;
* BICEP/Keck and Planck N_eff;
* PBH constraints;
* NIF, the laser record and ALICE energy densities.
