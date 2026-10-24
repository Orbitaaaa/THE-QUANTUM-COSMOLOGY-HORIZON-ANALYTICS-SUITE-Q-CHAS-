# Round 6d: one complete route, audited

**Brief.** Find one complete, physically consistent route from ordinary initial conditions to a growing daughter. Specifically:

* independently verify the Nariai-throat daughter;
* settle whether a localised ghost condensate stays stable in the near-Nariai geometry;
* derive the generalised entropy of our own solution, rather than assuming the ghost condensate violates the GSL;
* turn the two remaining "no" answers (ordinary quantum fields supply the push; preparation of the initial geometry) into "yes" if
  possible.

Theory only. **No apparatus is designed and no experimental procedure is proposed.** Labels as before.

**Files.**

* Code: `CODE/babyuniverse/throat_mol.py` (independent solver), `preparation.py` (complete route), `collapse.py` (regular-centre collapse with
  a second field), `ghost_condensate.py`, `gsl_solution.py`; additions to `throat_evolution.py` (light walls, ingoing classical flux,
  charge, corner re-tuning) and `quantum_throat.py`.
* Script: `CODE/scripts/round6d.py`.
* Data: `DATA/round6d.json`.
* Figure: `FIGURES/round6d.png`.
* Tests: `CODE/tests/test_round6d.py` (6).

---

## Where Genesis stands after Round 6d

| problem | Round 6c | Round 6d |
|---|---|---|
| Growing daughter | yes (one code) | **Yes, independently verified.** A second solver with a different numerical method and its own data integrator reproduces growth and hiding. Its threshold is higher (W between 0.053 and 0.071, against 0.031–0.035): the existence of a threshold is robust, its exact value to a factor of about 2. |
| Perturbation threshold | 270× below a flat-space flip (local energy) | **Preserved in the bridge-free configuration** when the NEC source sits in the band. In the single-evolution route it is not, because the push is injected inside the ball (Part III). |
| Hidden behind a trapped exterior | yes, but the exterior was already trapped in the data | **Yes, with an asymptotically connected exterior.** A light wall keeps the data normal and bridge-free. The exterior collapses into a parent black hole around the daughter, and the parent's event-horizon area never decreases. |
| Quantum fields supply the push | no (c ≥ 0.89 r_N² inside the EFT) | **No, for any c (theorem plus numerics).** The strong QNEC alone bounds every state's weighted deficit by c/48πτ. The same fields' vacuum polarisation raises the threshold faster than that bound grows (checked to c = 5 r_N²). Using lowest-Landau-level modes to get c_eff without extra species does not help. **Not on the critical path.** |
| Classical field supplies the push | conditional yes | **Ghost condensate: stable in the Nariai geometry.** Mode integration on dS₂ × S² over three response times gives growth ≤ 3.3× inside a window 86× wide in M³r_N. Its NEC phase switches itself off on the dS₂ Hubble time. The GSL cost is quantified below. |
| Preparation | no known method | **In theory, yes: a complete route from data containing no black hole is simulated.** A false-vacuum ball with a light wall, an imploding shell of ordinary radiation, and a classical NEC push give a hidden, inflating daughter. The data are normal, bridge-free and asymptotically connected, and satisfy the NEC except for the push. **Not a laboratory method.** |

---

## Part I — independent verification

### I.1 A second solver
`throat_mol` integrates the same equations by a different method:

* **Method of lines along the characteristics.** At each u it integrates the linear v-ODEs for (r_u, s_u, φ_u) with Heun steps.
* **Heun steps in u**, with substeps (the explicit u-update applies a Volterra operator of norm ~m²L_v/4).
* **Its own data integrator** (fixed-step RK4 on a 4× finer grid).

**Checks [Proposition-N].**

* Reissner–Nordström: second-order convergence; Misner–Sharp mass held to 1.3×10⁻⁶ at h = 0.1.
* Where both solvers have data, r agrees to all printed digits.

### I.2 The Nariai daughter, re-run

| pulse W (local) | diamond (Round 6c code) | method of lines (new) |
|---|---|---|
| 0.031 | no | no |
| 0.035 | grows (marginal), hidden | no |
| 0.044 | grows, hidden | no |
| 0.053 | grows, hidden | no |
| 0.071 | grows, hidden | grows, hidden |

The daughter, its hiding and the existence of a threshold are reproduced. The thresholds differ by a factor of 1.5–2: the diamond code gives
0.031–0.035 and the method of lines 0.053–0.071. The near-threshold growth criterion is sensitive to discretisation. Quantitative statements
that depend on W_th (the tipping law's constant, the ghost-condensate requirement, the quantum-field gap) inherit this factor. None of them
changes qualitatively under it: the ghost-condensate margins are about 100, and the quantum-field gap grows with c in either normalisation.

### I.3 A defect found and removed: the exterior was beyond a bridge
In every Round 6c run the false-vacuum wall's tension focused the outgoing data ray so strongly (Δr_v = −0.017 against r_v = +4×10⁻⁵) that
the true-vacuum exterior was already trapped on the data surface, with r decreasing outward into the parent's singularity. The parent's
asymptotically flat region was not on the data, and a configuration of that type is the Blau–Guendelman–Guth one beyond the Einstein–Rosen
bridge.

**Fix [Proposition + Proposition-N].**

* **Light wall.** A light wall (ε = 4, vphi = 1.5×10⁻³, δ = 10⁻⁴) has Δr_v = −3.6×10⁻³, smaller than r_v = 4.1×10⁻³ at the corner. The data
  are then normal everywhere outside the inner hole, with r increasing outward and no bridge.
* **Threshold.** W between 0.071 (no) and 0.106 (grows) at δ = 10⁻⁴.
* **Dynamics.** Because the configuration sits at its own Schwarzschild radius (m ≈ r/2 at the wall), the exterior collapses into a parent
  black hole by u ≈ 2.
* **With a push.** With W = 0.142, the daughter grows inside that black hole (n_big = 14,652; r up to 600 within the box). Without a push the
  whole configuration collapses.
* **Requirement unchanged.** False vacuum must continue along the ingoing ray to u ≳ 3u_hole (inner wall at u = 100 grows; at 10–60 it
  does not), as in Round 6c.

---

## Part II — the ghost condensate in the near-Nariai geometry

### Proposition 6d-1 (the NEC phase switches itself off) — [Theorem in the EFT]
Near the Nariai radius the geometry is dS₂(r_N) × S²(r_N), and f′(r_N) = 0 makes the static observer at r_N geodesic. In flat dS₂ slicing the
condensate φ = M²t has X = 1. Shift symmetry gives d/dt(aP′φ̇) = 0 with a = e^{t/r_N}. Hence a displaced condensate relaxes as
P′ = ε₀e^{−t/r_N}, and

* ρ + p = 2ε₀M⁴e^{−t/r_N};
* the null fluxes are N_u = N_v = πr²(ρ + p);
* the weighted deficit is W = 2πr_N²|ε₀|M⁴/(1/r_N + 1/τ).

### Proposition 6d-2 (stability window) — [Proposition, conditional on the EFT; checked numerically]
* **Perturbations.** L₂ = Aπ̇² − ε(∇π)² − (α/2M²)(∇²π)², plus the Jeans term −(αM²/2M_Pl²)k² (Arkani-Hamed–Cheng–Luty–Mukohyama).
* **Modes.** On dS₂ × S² a mode (k, l) has k_ph² = k²e^{−2t/r_N} + l(l+1)/r_N².
* **Window.** Tipping with W = W_th, while every instability stays within one e-fold over 3τ, requires
  $$1.9\times10^{-4}/\sqrt{\alpha P''}\ \le\ M^3r_N\ \le\ 1.7\times10^{-2}/\sqrt\alpha\qquad(\text{Planck units}),$$
  about 86× wide. It lies inside the EFT (Mr_N ≫ 1, |ε₀| ≪ 1) for r_N ≳ 10³ ℓ_P.
* **Numerical check.** Integrating every mode equation over 9.6 r_N, for r_N = 10⁴ and 10⁶ ℓ_P across the window, gives a largest
  amplification of 1.8–3.3×. A NEC phase 100× deeper grows by more than 100×.

**What this does and does not establish.**

* *Established:* a localised ghost-condensate NEC phase of the needed strength is stable, in the EFT sense, for the time it takes to tip a
  near-Nariai throat.
* *Not established:* a nonlinear simulation of the EFT itself (P-R6k), and the ghost condensate's UV completion and observational status.

---

## Part III — the complete route from data without a black hole (P-R6j)

### Proposition 6d-5 [Proposition-N]
Data (`preparation.build`), at a stage with no black hole:

* the corner inside a false-vacuum ball (pure de Sitter);
* along the outgoing ray, an **imploding shell of ordinary radiation** (null fluid, T_vv > 0) of local energy E = 7.26, then a light wall,
  then normal true vacuum;
* along the ingoing ray, de Sitter towards the centre.

Every sphere on the data surface is normal. r increases outward, so there is no bridge and no white hole. The NEC holds except for an optional
push N_u inside the ball. Outside the shell the mass parameter is M = M_N(1 − 0.009), the closest to Nariai that keeps the data normal.

Evolution:

* **No push.** The shell implodes and makes a black hole inside the false vacuum (more than 8×10⁵ trapped false-vacuum points). The exterior
  collapses into a parent black hole. No daughter.
* **Push N = 3 or 6.** No daughter.
* **Push N = 10.** The whole false-vacuum region becomes anti-trapped and inflates to r ≈ 4,500. Every row of the daughter era crosses only a
  thin sliver beyond the wall and then trapped true vacuum ending at r → 0. **The daughter is hidden inside the parent's black hole**, and the
  parent's event-horizon area never decreases.

**Caveats.**

* **The push is injected inside the ball.** The characteristic solver accepts sources only on the initial rays, and there f ≈ 0.75. The
  push's invariant energy there, about 1.5 M_N, is therefore not discounted. In the band's own frame the same flux carries
  W_band ≈ 10/λ ≈ 0.9, with λ ≈ 11 the boost between frames. A ghost-condensate phase located in the band would only need that.
* **The false-vacuum ball is part of the initial data.** Its creation (the "lift") is not simulated.
* **Regular centre.** With a regular centre, `collapse.py` shows that a shell of a second scalar field inside a false-vacuum bubble, from
  normal, asymptotically flat data with a regular centre, makes a black hole inside the false vacuum with the near-Nariai mass parameter
  (C/2 = 5.738–5.750 against M_N = 5.774). The near-horizon band is then too compressed in those coordinates to follow further.
* **Laboratory.** The ingredients are a black hole and a false vacuum tied by ρ_f = 1/(24πG²M²). This is a statement about what the laws
  permit, not a laboratory method.

---

## Part IV — the generalised entropy of our solution

With classical sources, S_out = 0 and S_gen = A/4G = πr² (Planck units).

### Proposition 6d-6 [Proposition-N + Theorem pieces]
* **Parent's event horizon.** Its area is non-decreasing in every run: bridge-free, isotropic (N_u + N_v) and route. The exterior GSL holds.
  The negative energy never crosses the parent's horizon outward, and the classical NEC violation stays inside.
* **Daughter observers' horizons** (boundaries of I⁻(p) for p deep in the daughter). These carry a deficit on their ingoing piece: the area
  decreases while the band's ingoing expansion is still negative, before the push arrives.
  * Numerically: 0.6 (push at u = 2) and 2.4 (push at u = 6) nats, against S_N = 942.
  * Analytically, the whole static-band crossing is bounded by ΔS ≤ π(r_c² − r_b²) ≈ 3.3√δ S_N.
* **Wall's theorem.** This is exactly the place where Wall's theorem (and Theorem GSL-3) bites. Its "weakly gravitating" hypothesis fails only
  at the level √δ, so the deficit is a fraction ~3√δ of the Nariai entropy rather than of order S_N.

### What this means for the ghost-condensate GSL dispute
* **Dubovsky–Sibiryakov and Feldstein.** Dubovsky–Sibiryakov (hep-th/0603158) argue that Lorentz-violating multi-speed theories violate the
  GSL. Feldstein (arXiv:0904.1212) builds ghost-condensate solutions with negative energy flowing into black holes and locates the pathology
  in NEC violation.
* **Mukohyama.** Mukohyama (arXiv:0901.3595) argues that the original ghost condensate respects the GSL, because its speed differences are
  suppressed by M²/M_Pl². arXiv:0908.4123 reaches a similar conclusion about the proposed violations.
* **Our solution** needs a generalised-entropy deficit of order √δ S_N, but only on the daughter observers' horizons, never on the parent's.
  * If the ghost condensate obeys the GSL in Mukohyama's sense, the deficit must be paid by S_out across those horizons. For a thermal 2D
    sector that requires c ~ S_N, the same scaling that defeats quantum pushes.
  * If it does not, as in Feldstein's mechanism, the route is consistent.
* **Bottom line.** The dispute is therefore decisive for the route, and the deficit it must accommodate is now quantified.

---

## Part V — quantum fields (the remaining "no")

* **Proposition 6d-3 [Theorem].** For every state of a 2D CFT obeying the strong QNEC with S′ → 0 before u₀,
  $$\int T_{uu}e^{-(u-u_0)/\tau}du\ \ge\ \frac1{2\pi}\int w\Big[\frac{S'}\tau+\frac6cS'^2\Big]\ \ge\ -\frac c{48\pi\tau},$$
  with equality for S′ = −c/12τ. This extends Proposition 6c-10 from conformal-map states to all states, and uses no QEI.
* **Proposition-N 6d-4 (threshold versus c).** In the Bunch–Davies-type background of the same fields the threshold rises:

  | x = c/r_N² | 0 | 1 | 2 | 3.3 | 5 |
  |---|---|---|---|---|---|
  | W_th/r_N | 1.85×10⁻³ | 5.1×10⁻³ | 7.8×10⁻³ | 1.14×10⁻² | 1.55×10⁻² |

  The capacity is 2.07×10⁻³x. The gap W_th − capacity grows (0.053, 0.064, 0.077, 0.088 at l_H = 30), so the two never cross.
* **Re-tuning the corner radius does not help.** At c = 300, η ∈ [−0.026, +0.02] leaves the threshold between 0.07 and 0.14.
* **Conclusion: no, and the obstacle is structural.** Fields that can supply negative energy back-react on the band by more than they can
  supply. The magnetic (lowest-Landau-level) loophole, which raises c_eff without extra species, avoids the species bound but not this.
* **The route does not need it:** the classical field of Part II supplies the push.

---

## Failed approaches and corrections
* **F32 (correction).** The Round 6c hiding results used data whose exterior was beyond a bridge (Part I.3). They are superseded by the
  bridge-free light-wall runs, which also grow and hide.
* **F33.** The semiclassical re-tuning of the corner radius against Bunch–Davies back-reaction: no effect on the threshold.
* **F34.** The method-of-lines solver without u-substeps: unstable for a massive field over a long v-line (Volterra norm m²L_v/4 ≈ 10³).
* **F35.** A regular-centre collapse in the r = v/2 cone gauge: the lapse jumps by e³ across a massive shell and the stiff scalar goes
  unstable. Fixed with the α = 1 gauge.
* **F36.** Emitting the push inside the false-vacuum ball: the route closes, but the invariant cost is about 1.5 M_N (no discount).

## Open problems
* **P-R6k.** A nonlinear ghost-condensate EFT simulation in the throat (the stability window is linear).
* **P-R6l.** A bulk NEC source in the band, so that the discount is kept in the single-evolution route.
* **P-R6m.** The "lift": simulate the creation of the false-vacuum ball by NEC-respecting dynamics.
* **P-R6n.** The GSL accounting for the ghost condensate: compute S_out of the condensate and its π-quanta across daughter-observer horizons,
  versus the required √δ S_N.
* **P-D2′.** Non-spherical stability of the tipping.
