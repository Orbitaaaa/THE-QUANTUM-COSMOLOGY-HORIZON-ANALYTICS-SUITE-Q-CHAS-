# Round 6e — the three immediate problems

Problems posed (user, after Round 6d):

1. Make the ghost condensate generate the required negative-energy pulse itself, in the right location, with full nonlinear stability.
2. Resolve the generalised-entropy problem and establish that the complete transition is physically consistent.
3. Explain how to create the initial false-vacuum ball, rather than assuming it already exists.

(Filed as "Round 6e" because the project already has a "Round 7" roadmap stage, `ROADMAP_ROUNDS5_7_12.md`.)

New code: `condensate.py` (the ghost condensate as a dynamical, back-reacting source in the double-null solver: geodesic flow, shift
charge, condensate clock, optional false-vacuum tilt), `fireball.py` (supercooling of a hot ball from the true vacuum),
`scripts/round6e.py`, `tests/test_round6e.py` (7 tests), `DATA/round6e.json`, `FIGURES/round6e.png`. Small additions:
`ghost_condensate.mode_growth` accepts a measured NEC history.

Labels as before: **Theorem** (proved), **Proposition** (proved under stated assumptions), **Proposition-N** (numerical, with convergence
evidence), **Literature**, **Assessment**, **Conjecture**. Units G = c = ħ = 1; the simulations use l_H = 30 (r_N = 17.3, S_N = πr_N² = 942)
and scale classically with r_N.

---

## Summary

| Problem | Answer | What was done |
|---|---|---|
| 1. Condensate makes the pulse itself, located, stable | **Yes, with one stated input (the self-generated window is tight)** | The condensate's own stress tensor tips the throat (no prescribed flux). Its charge sits where it is needed, and it is stable on the measured history. Shift symmetry forbids making the charge from nothing (Theorem 6e-1). A soft false-vacuum tilt makes it dynamically, from zero charge (Proposition-N 6e-3). The complete route closes with this self-generated charge. |
| 2. Generalised entropy and physical consistency | **Resolved for transient daughters. Impossible for eternal ones (theorem)** | Every eternal daughter grown from normal data violates the local GSL, whatever supplies the NEC violation (Theorem 6e-6, measured 7–164 nats). The same tilt that makes the charge also ends the daughter's inflation (Proposition-N 6e-7). The daughter then exits to Λ = 0 and recollapses inside the parent hole. The GSL then holds on every causal horizon of a complete observer. |
| 3. Create the false-vacuum ball | **Explained and computed, one link open** | No classical obstruction applies (Theorem 6e-8). Heating restores the false vacuum from the true vacuum. Free-streaming cooling then leaves an empty, momentarily static ball of radius 0.6 l_H, normal throughout for strong transitions (Proposition-N 6e-10). Confined radiation cannot do it (Proposition 6e-9). Open: one GR simulation of the emptying, handed off to the route. |

---

## Problem 1 — the condensate as its own source

### 1.1 The dynamical condensate (new solver capability)
Infrared limit of the ghost-condensate EFT L = M⁴P(X) − (α/2)M²(□φ)² − V(φ), around X₀ = ½ with φ normalised to proper time.
* **Stress tensor.** T_ab = n u_a u_b on top of a cosmological constant.
  * n = M⁴P′√(2X) is the shift charge. Its sign is free; n < 0 is the NEC-violating phase.
  * u_a = −∂_aφ/√(2X) is geodesic.
  * Pressure enters only at O(n²); α enters through the gradient terms.
* **Double-null form.** u_u = −(e^s/2)e^η and u_v = −(e^s/2)e^{−η}, with D = 4πr²n e^{s−η}:
  * geodesic flow: (η − s)_u + e^{2η}(η + s)_v = 0;
  * charge: D_u + (De^{2η})_v = 4πe^{2s}r²S;
  * clock: T_u + e^{2η}T_v = e^{s+η}.
* **Coupling to gravity.** r_uv gains G·T2_uv/r and s_uv gains −G·T2_uv/r², with T2_uv = De^{s+η}/4.
* **Scheme.** Box (diamond) transport, inside the existing three-pass diamond.

**Verification.**
* **Zero condensate.** Reproduces the diamond solver bit for bit.
* **Constraints with a back-reacting condensate.** Converge at second order: 5.9e-4 → 1.6e-4 → 4.2e-5.
* **Clock and source.** The clock is unit-normalised to 1e-7. The tilt source gives n = −λτ to 3e-4.
* **Maximum residual.** In the Nariai runs the maximum constraint residual sits at the boosted thin wall. It is just as large with no
  condensate at all, so it is inherited from the base solver, not introduced here.

### 1.2 Theorem 6e-1 (shift-charge conservation) [Theorem in the EFT]
The field equation is ∇_μ(n u^μ) = αM²□θ − V′(φ), where θ ≈ −□φ is the expansion of the condensate flow.
1. **(α → 0, V = 0)** n × (comoving volume) is conserved along every flow line. A NEC-phase element was NEC-phase at all earlier times.
   No local dynamics turns a charge-free condensate into a NEC-violating one.
2. **(α ≠ 0)** In a homogeneous patch, n = −αM²θ̇ + C/a³ = −3αM²Ḣ + C/a³.
   * This renormalises Newton's constant in the Ḣ equation, G_eff = G/(1 − 12πGαM²).
   * It is sign-locked to Ḣ, so it amplifies existing motion and cannot start a tipping.
   * Inside the stability window 12πGαM² ≲ 10⁻⁴α^{2/3}(r_N/10⁶)^{−2/3}.
3. **(V′ ≠ 0, soft shift breaking)** dn/dτ + θn = −V′: charge is made locally, with the sign of −V′. This is the mechanism of ghost
   inflation and of the stable NEC-violating phase of Creminelli–Luty–Nicolis–Senatore.

So "generate the pulse itself" has two honest readings. Either the condensate's own NEC-phase charge does the tipping, given that the
charge is part of the state (Proposition-N 6e-2). Or the condensate makes that charge dynamically, which requires V′ ≠ 0
(Proposition-N 6e-3).

### 1.3 Proposition-N 6e-2 (the condensate tips the throat by itself)
Setup:
* the near-Nariai light-wall band of Round 6d (δ = 10⁻⁴), with **no prescribed flux**;
* condensate at rest in the corner's geodesic frame;
* proper NEC-phase density n₀ throughout the static false-vacuum region on the data.

Results:
* **Threshold.** |n_th| r_N² = 1.94, 1.65 and ≈1.47 (×10⁻³) at h = 0.08, 0.04 and 0.02. The last is bracketed between 4.7 and
  5.1 ×10⁻⁶ in simulation units. Extrapolated: |n_th| r_N² ≈ (1.2–1.5)×10⁻³.
* **Same law as a prescribed pulse.** The weighted deficit at threshold, W = ∫πr²|n|e^{−u/τ}du ≈ 0.09–0.12, matches the Round 6d
  prescribed-pulse threshold for light walls (0.071–0.106).
* **Hidden daughter.** The daughter is hidden: there is no untrapped true vacuum at large r. The parent's event-horizon area never
  decreases.
* **Location** (all at h = 0.08):
  * Compact lumps at the corner obey charge × width = const: n_th = 3.3, 1.75 and 0.91 ×10⁻⁵ for widths 4, 8 and 16.
  * Removing the charge beyond v = 20 on the outgoing ray leaves the threshold unchanged (6.3–7.0 versus 6.5 ×10⁻⁶).
  * Removing it beyond 0.3 u_hole on the ingoing ray doubles the threshold.
  * So what matters is the charge within about one response time τ of the corner, on the ingoing side. That is exactly where the
    condensate sits: the corner's static frame is geodesic, so the condensate stays put and is diluted only at the dS₂ rate.
* **CLNS margin.** The CLNS stability ceiling |Ḣ| ≲ H² allows |n| ≲ 1/(4πr_N²) = 0.080/r_N². The required phase is about 50× below it.

### 1.4 Proposition-N 6e-3 (self-generated NEC phase)
Setup:
* the tilt is confined to the false vacuum: V = λ_f χ(ϕ) φ, where χ = 1 in the false vacuum and 0 in the true vacuum;
* **zero charge on the data**.

The condensate makes its own NEC phase wherever the false vacuum is.
* **Energy bookkeeping.** By energy conservation the false-vacuum energy rises with the condensate clock, ρ_f → ρ_f + λ_f T. The two
  pieces together have ρ + p = n < 0.
* **Neither piece alone is physical.**
  * The relative constraint residual is 2.1e-6 → 5.3e-7 → 2.3e-7 (h = 0.08, 0.04, 0.02) for the consistent pair.
  * It is 0.04–0.37, and does not converge, with either half alone.
  * The source alone does not tip; the drift alone does.
  * Only the conserved combination is reported as physics. Its strongly negative pressure also helps the tipping.
* **Threshold.** λ_th r_N³ = 4.7×10⁻⁴ (h = 0.08) and 3.8×10⁻⁴ (h = 0.04). These are final values with the stable (implicit) scalar step;
  the explicit step gave 5.5 and 4.5 ×10⁻⁴. The daughter is hidden and the parent area never decreases.

### 1.5 Proposition-N 6e-4 (stability on the measured history)
* **Method.** The NEC history along the corner worldline is taken from the simulations. Beyond the simulated range (≈1.1τ) it is
  continued by dS₂ dilution (initial charge) or held at its steady state with daughter expansion (self-generated). The full mode
  equation (gradient, k⁴ and Jeans terms) is then integrated over all (k, l) up to 3τ at r_N = 10⁶.
* **Initial charge.**
  * The integrated phase is I = ∫|n r_N²|dx = 0.0032.
  * Growth is 3.96×, 2.67× and 1.04× at 0.5, 1.5 and 5 times the analytic lower edge (M³r_N = I/4 ≈ 8×10⁻⁴).
  * The window is roughly 1×10⁻³ ≤ M³r_N ≤ 1.7×10⁻², about 17 wide.
* **Self-generated.**
  * I = 0.042: the steady NEC phase lasts, so the integrated phase is larger.
  * Growth is 3.1× at M³r_N = 5.2×10⁻³ and 1.97× at 1.6×10⁻². It is 11.6× at 5.2×10⁻², beyond the Jeans edge (1.7×10⁻²).
  * The window is narrow, roughly 6×10⁻³ ≤ M³r_N ≤ 1.7×10⁻², about 3 wide. It exists, but it is tight.
* **What "nonlinear stability" covers.** The radial (l = 0) sector is fully nonlinear in the simulation. All l ≥ 1 sectors are linear
  on the nonlinear background, with growth < e. Spherical symmetry prevents a nonlinear l ≥ 1 test (P-D2′).

### 1.6 Proposition-N 6e-5 (the complete route with a self-made push)
Setup:
* the Round 6d route data: no black hole, false-vacuum ball, imploding ordinary shell, light wall, normal true vacuum;
* the in-ball push replaced by the tilted condensate with zero initial charge.

Results:
* A daughter grows above λ_th = (3.5–3.8)×10⁻⁷ (bisected) and not below.
* **Hidden, by the event-horizon test.**
  * The parent's event horizon sits at u_EH = 0 to grid resolution.
  * None of the daughter's outgoing rows reaches the outer edge: 0 of 79 at λ_f = 10⁻⁶ and 0 of 215 at 3×10⁻⁶.
* **Transiently untrapped far exterior.** With the tilt, part of the far exterior (v ≈ 250–340, u ≈ 3.5–4.8) stays locally untrapped
  true vacuum for a while. It lies inside the parent's event horizon, and every one of its rows ends at the singularity. (The cruder
  "untrapped true vacuum at large r" count used in Round 6d therefore misreports this case; `condensate.hiding_report` gives the
  event-horizon test.)
* The parent's event-horizon deficit is 0.
* λ_th r_N³ ≈ 1.9×10⁻³ at the route's δ ≈ 0.009, about 5 times the δ = 10⁻⁴ value. That is consistent with the √δ scaling of the
  invariant cost.
* The push is now made in the band itself. The Round 6d penalty F36 (an in-ball push costing 1.5 M_N) no longer applies.

**Problem 1 status: solved in the sense above.**
* The pulse is the condensate's own stress.
* It sits in the band.
* It is stable (radial dynamics nonlinear, l ≥ 1 linear).
* Its charge is either part of the state or made by a soft, technically natural tilt.

Limits:
* the condensate is in its IR (dust) limit, with the gradient terms only in the linear analysis (P-R6e-a);
* the tilt is an explicit model choice (P-R6e-d).

---

## Problem 2 — generalised entropy and physical consistency

### 2.1 Theorem 6e-6 (eternal daughters violate the local GSL) [Theorem, classical sources, spherical symmetry]
Assume:
* the data are normal (r_u < 0 < r_v on the whole data surface);
* the sources are classical (S_out = 0);
* the daughter is eternal (it ends in de Sitter), with a complete observer γ in it.

Proof sketch:
1. By Theorem GSL-3 (Round 6b), the past of γ is compact in v: there is a last ingoing ray v_∞ < ∞ whose events reach γ.
2. The causal horizon ∂I⁻(γ) contains the ingoing generator v = v_∞, which starts on the data surface.
3. There r_u < 0, so its area decreases from the start, and the local GSL fails.

What follows:
* The argument does not care what supplies the NEC violation.
* **No choice of condensate, or of any other classical source, can make an eternal daughter GSL-consistent.**
* The size of the deficit is bounded by π(r_c² − r_b²) ≈ 3.3√δ S_N.

Measured largest decreases of A/4G on daughter-observer horizons:

| configuration | daughter deficit (nats) | parent event horizon |
|---|---|---|
| condensate, initial charge n₀ = −6×10⁻⁶ | 13.2 | 0 |
| self-generated, λ_f = 1.3×10⁻⁷ | 16.1 | 0 |
| self-generated, λ_f = 4×10⁻⁷ | 7.1 | 0 |
| route, self-generated, λ_f = 10⁻⁶ | 123 | 0 |
| route, self-generated, λ_f = 3×10⁻⁶ | 164 | 0 |
| route, in-ball push (Round 6d) | 50 | 0 |

These are 0.8–17% of S_N, all within the bound. They scale with S_N at fixed δ, so they are astronomical for macroscopic throats.
(The 0.6–2.4 nats quoted in Round 6d came from a single sampling of daughter points; the systematic sampling here gives larger numbers.
Treat this as a correction of scale, not of sign.)

One caution: when the charge also fills the region up to the wall and is strong (n₀ = −2×10⁻⁵), part of it crosses the parent's horizon,
which loses 0.32 nats. Charge confined to the band, or produced by the tilt, gives 0.

### 2.2 Proposition-N 6e-7 (the transient daughter: the same tilt ends inflation)
Mechanism:
* The tilt that makes the NEC phase also raises the false-vacuum energy: ρ_f + λ_f T.
* When the effective barrier ratio 4(ρ_f + λ_f T)/(λ_q v⁴) passes 16/3, the false vacuum becomes a hilltop.
* Any seed then sends the field to the true vacuum, with Λ_final = 0 = Λ_parent.
* The exit happens at clock time T_exit = ρ_f(16/(3ε) − 1)/λ_f, which the barrier ε moves independently of the tipping.
* The tilt profile is held flat beyond the false vacuum (C¹), so the roll can only go to the true vacuum. Without that, the tilted
  quartic has a spurious second false minimum.

Simulation (the route's ε = 4, λ_f = 2×10⁻⁷, so T_exit = ρ_f/(3λ_f) = 221; domain 520):
* **Seeded run** (the field 10⁻⁶ v off the false vacuum):
  * Identical to the control up to u ≈ 320, with a daughter of r_max ≈ 140–210.
  * At clock time T ≈ 225–250 the field rolls to the true vacuum, against 221 predicted.
  * The region then recollapses inside the parent black hole: r_max 139 → 66 → 42 → 30 → 16 by u = 520.
* **Control** (no seed, so it cannot exit): it keeps r_max at 76–427 over the same range.
* Final driver numbers:
  * First exit inside the band at clock time T = 232, against 221 predicted.
  * Largest radius by u-window: seeded 9,905 → 9,686 → 1,294 → 34; control 9,905 → 9,766 → 9,982 → 9,985 (the tracked limit is 10⁴).
* Both daughters are hidden (no outgoing row reaches the outer edge), and the parent's event-horizon area never decreases in either run.

(An earlier version of this run used ε = 2.5 with the explicit scalar scheme. It was contaminated by the instability recorded as F42, so it
is withdrawn.)

Consequences for the GSL:
* With Λ_final ≤ Λ_parent and a crunch, there are no complete daughter observers.
* The only causal horizons of complete observers are then the parent's, whose area never decreases (deficit 0).
* So the GSL in the form used by Wall and by Theorem GSL-3 (causal horizons of complete observers) holds for the complete transition.

**Caveat.** If the GSL is imposed on the past boundaries of every inextendible worldline, including those that end at the daughter's
crunch (the Jacobson–Parentani definition), those boundaries still carry a few-nat deficit (≈4 nats measured). The consistency claim
therefore rests on the complete-observer formulation (P-R6e-b).

A related point about lifetime. A short exit gives a short-lived daughter that crunches. Delaying the exit (smaller λ_f relative to ρ_f,
or a deeper barrier) gives more e-folds and an exponentially later crunch. It does not give an eternal daughter: that is exactly what the
GSL forbids from normal data.

### 2.3 Physical-consistency checklist for the complete transition
* **Energy–momentum conservation.** Constraint residuals converge at second order with the condensate, the tilt and the drift included.
* **Well-posedness and stability.** Proposition-N 6e-4 covers this (radial nonlinear, all l linear). The NEC phase is about 50× below
  the CLNS ceiling.
* **EFT validity.**
  * |ε| = |n|/(2M⁴) ≪ 1 inside the window.
  * Gradients ≪ M, since M r_N ≫ 1 for r_N ≳ 10³ ℓ_P.
  * Below the cutoff the π-phonons have group velocity < 1.
* **Naturalness.** The tilt λ_f breaks the shift symmetry softly, so it is technically natural: setting λ_f = 0 restores the symmetry.
* **Exterior.** The parent's event-horizon area never decreases, the daughter is hidden, and by Birkhoff the exterior cannot tell
  whether a daughter formed.
* **Generalised entropy.** Consistent for the transient daughter (2.2). Inconsistent for any eternal daughter (2.1).

**Problem 2 status.**
* For eternal daughters, the generalised-entropy problem cannot be resolved by any classical mechanism (Theorem 6e-6).
* It is resolved for transient daughters, and the condensate tilt supplies the exit automatically (Proposition-N 6e-7).
* Every other consistency check passes.

---

## Problem 3 — creating the false-vacuum ball

### 3.1 Theorem 6e-8 (no classical obstruction) [Theorem]
The route's data surface (Proposition 6d-5) is:
* normal everywhere: no trapped and no anti-trapped sphere;
* NEC-satisfying;
* asymptotically flat and bridge-free.

Its ball is sub-horizon, with the wall at r ≈ 0.6 l_H (exact value in the JSON). The obstructions do not apply:
* The Farhi–Guth / Penrose obstruction to making an inflating region needs an anti-trapped surface.
* The black-hole singularity theorem needs a trapped one.
* Neither is present, so NEC-respecting dynamics is not forbidden from producing these data.

This is the decisive difference from the classic "universe in the laboratory" problem: the ball does not need to inflate when it is
made. Inflation starts only after the band tips, and it is hidden behind the parent's horizon.

### 3.2 Proposition 6e-9 (confined radiation cannot empty the ball) [Theorem in the homogeneous model]
* **Statement.** A closed homogeneous ball of radiation plus vacuum energy that turns around and recollapses does so at a < l_H/√2, with
  ρ_rad ≥ ρ_f.
* **Proof.** Turnaround requires ρ_f a² + C/a² = 3/(8π). The recollapsing root lies to the left of the minimum of the left side, which
  gives C/a⁴ ≥ ρ_f.
* **Numerically** (closed FRW, ϕ starting in the true vacuum):
  * Every normal ball turns around with 1.7–19 ρ_f of radiation left and radius ≤ 0.42 l_H.
  * The others inflate, which is the anti-trapped case that Farhi–Guth forbids.
* **Consequence.** A ball kept hot by confined radiation cannot become the route's empty ball of radius 0.6 l_H. The radiation must
  leave.

### 3.3 Proposition-N 6e-10 (supercooling creates the ball from the true vacuum)
**Model.** The particles in the hot matter get their mass from ϕ, m_i² = g_i²(ϕ + v)², so they are massless at the false vacuum. The
leading (high-T) thermal correction then favours ϕ = −v:
F(ϕ, T) = −(π²/30)g T⁴ + (c_T/2)T²(ϕ + v)² + U(ϕ),
with dissipation Γ of the field's oscillations into the bath.

**Spinodal.** Above T₁ = 1.24 T_c (for the route's ε = 4; 1.375 T_c for ε = 2), the field cannot stay in the true vacuum.

**Stage A — restoration.** A momentarily static hot ball of radius 0.6 l_H at T = 1.4 T_c, with ϕ in the true vacuum:
* It is normal iff the transition is strong: 2m/r = 0.95 at α_PT = 2.5 and 0.68 at α_PT = 3, but 1.37 at α_PT = 2.
* The field moves to −v within t ≈ 6, against a light-crossing time of 18.

**Stage B — emptying.**
* The radiation free-streams out on the light-crossing time.
* ϕ = −v minimises F at every temperature (it is a stationary point of U and the minimum of the thermal term), so slow cooling leaves the
  field exactly there. The adiabaticity parameter is 0.004.
* The result is an empty, momentarily static false-vacuum ball of radius 0.6 l_H, with zero kinetic energy: the route's interior.
* The same holds with the route's ε = 4 potential (restored to −0.995 v, ending at −v).

**Energy.**
* The released radiation is E_rad = 2.9 (α_PT = 3) to 5.3 (α_PT = 2.5), against the shell the route needs, M_N(1 − δ) ≈ 5.7.
* Reflected back, the same radiation can supply 50–90% of the imploding shell.

**Assembled chain** (every stage NEC-respecting and normal):
1. a hot ball, with ϕ in the true vacuum;
2. heating restores the false vacuum;
3. free-streaming radiation leaves, and the field stays in the false vacuum;
4. the radiation (reflected) and/or an extra ordinary shell implodes onto the ball;
5. near-Nariai band;
6. the condensate's own NEC phase tips it;
7. a hidden, transient daughter.

**Problem 3 status.**
* The mechanism is explained and computed stage by stage, and no obstruction applies.
* **Open (P-R6e-c):** a single general-relativistic simulation of the emptying (free-streaming radiation and the moving wall, with
  reflection and timing) that hands off to the route's data. Only spherically symmetric radial null dust is available in the
  characteristic code, and isotropic free streaming needs more.

---

## Failed approaches (F37–F41)
* **F37.** Condensate data at rest in the static frame. Near both horizons the static observers are boosted without limit, so W_v
  blows up (W_v = 10³–10⁶). The physical frame is the corner's geodesic frame.
* **F38.** Splitting the tilt into "charge source only" or "vacuum drift only". Each violates energy conservation (relative constraint
  residual 0.04–0.37). Only the combination is physical.
* **F39.** Exits that do not work:
  * A fixed hilltop false vacuum drains the band before it tips.
  * An exactly stationary false vacuum never exits without a seed.
  * The unclipped tilted quartic can roll into a second positive-energy minimum.
  * The working exit is the clock-driven one (2.2).
* **F40.** A coherent classical lift (an ingoing massive-field packet, no dissipation).
  * In a random search the field crosses the barrier only transiently. Late false-basin points cover ≤ 10% of the interior.
  * With full gravity the packets collapse to black holes (2m/r up to 0.9).
  * Capture into a metastable basin needs a sink for the kinetic energy.
* **F41.** A confined-radiation fireball. It turns around before it cools (Proposition 6e-9).
* **F42 (numerical, found and fixed this round).** The explicit treatment of the scalar's potential force in the diamond scheme is
  unstable wherever m h e^s/2 > 1. There the scalar's proper Compton length is not resolved (s ≳ 2 inside the parent hole).
  * A field 10⁻⁶ v off the false vacuum then grows at the grid scale into a spurious true-vacuum intrusion. That removed the daughter in
    every seeded run.
  * A field sitting exactly at the minimum is an exact solution and is unaffected. So the unseeded results (in this round and in the
    Round 6c/6d solver) stand.
  * Fix in `condensate.py`: the potential force is implicit in ϕ_N (Newton step, trapezoid-like, stable for U″ > 0). Seeded and unseeded
    initial-charge runs now give identical daughters (test). The legacy `throat_evolution` kernel keeps the explicit form, so its results
    should not be used with off-minimum seeds in large-s regions.
  * Residual effect: the implicit step is stable but does not damp grid-scale oscillations. Where the light wall is crushed against the
    forming parent singularity (s ≈ 4–6), ϕ shows a checkerboard (≈ +3v alternating with −v) that does not resolve as h → 0 (it would
    need h ≲ 10⁻³). Hiding is therefore judged by points where the field actually sits at the true vacuum (|ϕ/v − 1| < 0.5): there are
    none at h = 0.08 and 0.04. The raw ϕ > 0 count is kept in the data as `phi_positive_flags`.

## Open problems
* **P-R6e-a.** The full ghost-condensate EFT (k² and k⁴ terms) simulated nonlinearly, beyond the IR fluid, together with the non-spherical
  modes (P-D2′).
* **P-R6e-b.** The GSL for observers who end at the daughter's crunch (Jacobson–Parentani horizons). Also long-lived transient daughters
  (many e-folds), and a persistent-neck (charged) variant.
* **P-R6e-c.** A single GR simulation of the fireball emptying, the reflection and the hand-off to the route.
* **P-R6e-d.** Microphysics of the tilt and the thermal sector. In the toy units the quartic is strongly coupled (λ ≈ 5×10⁷, because
  v ≪ ρ_f^{1/4}). The realistic scaling with r_N and the UV completion are open.
* **Closed or advanced:**
  * P-R6k: the condensate is now a dynamical source (IR limit).
  * P-R6l: a band-localised bulk source.
  * P-R6m: the ball's creation is explained.
  * P-R6n: the GSL accounting is replaced by the eternal/transient theorem.

## Literature used
* **Arkani-Hamed, Cheng, Luty, Mukohyama** (2004) hep-th/0312099: ghost condensation.
* **Arkani-Hamed, Creminelli, Mukohyama, Zaldarriaga** (2004) hep-th/0312100: ghost inflation, a soft tilt V(φ).
* **Creminelli, Luty, Nicolis, Senatore** (2006) hep-th/0606090: stable NEC violation with |Ḣ| ≲ H².
* **Farhi, Guth** (1987): the obstacle needs anti-trapped surfaces.
* **Wall** arXiv:1010.5513, with Theorem GSL-3 (Round 6b): the complete-observer GSL.
* **Guth** (1981): supercooling into a false vacuum.
* **Strong first-order transitions** (α_PT ≳ 1): standard; see e.g. Caprini et al. arXiv:1512.06239.

## What this does and does not establish
* **Established.** Within classical GR plus the ghost-condensate EFT (IR limit, with a soft false-vacuum tilt), plus standard
  finite-temperature physics, there is a chain from a hot ball in the true vacuum to a hidden, transient daughter. Each link is computed,
  and the links are energy-consistent and GSL-consistent for complete observers.
* **Not established.**
  * One continuous simulation of the whole chain (the fireball-to-route hand-off is missing).
  * Nonlinear non-spherical stability.
  * The full EFT beyond its IR fluid.
  * That nature contains a ghost condensate, a false vacuum or such a tilt (problem 4).
  * Any technological accessibility (problems 5 and 6).
