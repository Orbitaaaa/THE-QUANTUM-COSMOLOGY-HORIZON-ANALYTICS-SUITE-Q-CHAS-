# Round 6j — what the invariant really bounds, and fifteen loopholes

The user asked for a programme organised around three targets:

| Target | Question |
|---|---|
| A. Gravitational compactness | Can a daughter-like region form without the classical trapped-surface condition? |
| B. High-energy initial state | Can an environment prepare the state with energy already present? |
| C. Classical spacetime | Can a quantum transition make the daughter with no classical microscopic throat? |

The user also pointed out a mathematical subtlety in Theorem 6i-11 and asked us to separate three things: a bound on the geometry, a bound on
the energy stored in a configuration, and a bound on the energy an apparatus must supply. Then came fifteen candidates in four groups.

All were computed or checked against the literature:
* Code: `accessibility6j.py` (Part I, Group A); `loopholes6j.py` (Groups B–D).
* `stable_shells.Shell` now allows a black hole inside the bubble.
* `scripts/round6j.py`, `tests/test_round6j.py` (11 tests), `DATA/round6j.json`, `FIGURES/round6j.png`.

Labels as before.

**Answer in one paragraph.**
* **The user is right about the distinction.** The invariant is a geometric identity at the first marginally trapped sphere, and with ρ ≥ 0 a
  bound on stored energy. It says nothing about the energy an apparatus must supply when nature already stores most of it.
* **Which density enters.** We derived it exactly: the areal-volume mean of the local energy density, in the frame with no radial energy
  flux. Non-uniformity and negative-energy components only make things harder.
* **Classical loopholes (Group A) and environment-prepared states (Group B):**
  * A field-dependent Planck mass costs more than it buys.
  * Higher curvature is a special case of that.
  * A finite pocket expands by 10⁻¹⁶.
  * Accretion and compact precursors cannot prepare a daughter-ready state.
  * A small supplied energy needs a natural precursor that is stable, old and within ~10⁻⁸ of threshold.
* **Quantum spacetime (Group C).** It gives a clean new law: a Euclidean channel that stays connected to us through a throat pays action
  ΔS = 3π³a₀² (8πG = 1), about 1.18 × (throat area)/4G. The throat-free limit is the disconnected no-boundary universe. That is the Q* lesson
  again, now exact.
* **The two doors that remain open are both in Group D.** Both are conditional on unestablished quantum gravity:
  * If black-hole interiors undergo a Frolov–Markov–Mukhanov transition, every black hole contains an inflating daughter, and making a
    daughter reduces to making a black hole.
  * If gravity is strong at 6–20 TeV (n ≥ 4 extra dimensions, still allowed), a 100 TeV collider makes black holes. But they sit at the
    fundamental scale (R_S M_D ≈ 1), so a semiclassical interior needs collisions ~10³–10⁴ times heavier.
* **Correction.** The Round 6c "partner vacuum" formula was off by a factor of 3 in density (F56).

---

## 1. Part I — geometry, stored energy, supplied energy

**Identity (the user's formula).**
* E√ρ̄ = √(3/32π) c⁶/G^{3/2} k^{3/2}, with k = 2GE/(c⁴R) (`compactness_identity`, verified to 10⁻¹²).
* Theorem 6i-11 is this identity at k = 1, the value it takes at a marginally trapped sphere. So it is a statement about geometry, not about
  an apparatus.

**Theorem (Misner–Sharp / Kodama identity; symbolic, `misner_sharp_identity`).**
* For ds² = −e^{2s}du dv + r²dΩ²: ∂_u m = r²e^{−2s}(G_uv r_u − G_uu r_v), and ∂_v m = r²e^{−2s}(G_uv r_v − G_vv r_u).
* In an orthonormal frame, along a radial slice: dm/dr = 4πr²(ρ + F r_t/r_n), with F the outward energy flux.
* So m(R) = ∫₀^R 4πr² ρ dr exactly in the frame with no radial flux (comoving, for fluids). On a time-symmetric slice it is the ordinary
  energy density.

**Proposition 6j-1 (which density; question A).**
* *The statement.* At the first marginally trapped sphere of regular data (radius R), the zero-flux energy density averaged over the areal
  ball is exactly ρ̄ = 3c⁴/(8πGR²).
* *Consequences.*
  * Some region has local density ≥ ρ̄.
  * Non-uniform profiles only raise the peak: 3.3× for a centrally peaked ball, 2.9× for a hollow shell, 3.0× for a negative core inside a
    positive shell (`which_density`).
  * Negative components raise the positive energy needed: ∫ρ₊ = Rc⁴/2G + ∫ρ₋, for example 1.04 × Rc⁴/2G for a negative core.
* *Why the formation argument holds.* Before the first marginal sphere, r is monotonic on the slice, so this integral is the energy inside.
* *A caveat.* A throat with no "inside" (a traversable wormhole) can have M_ADM = 0. But it cannot be formed from regular data by this
  route.

**Proposition 6j-2 (three different bounds; question B).** For a marginal sphere of radius R:
1. *Geometric (exact).* E_MS = Rc⁴/2G, ρ̄ = 3c⁴/(8πGR²).
2. *Stored (ρ ≥ 0).* The local matter energy in the ball is at least E_MS, because proper volume ≥ areal volume (binding lowers m).
3. *Supplied.* E_sup = E_stored − E_natural. The invariant does not bound it. Accretion-induced collapse of a neutron star is an
   astrophysical example: a small addition makes a black hole.

**Proposition 6j-3 (what a small supplied energy demands of nature).** E_sup ≪ E_stored needs a natural precursor within
δ = E_sup/(Mc²) of the threshold that satisfies all three conditions below. For E_sup = 10²⁰ J (about two months of world primary energy),
δ ≤ 6.4×10⁻⁹ at 1.7×10¹¹ kg and 5.6×10⁻²⁸ at a solar mass.
* *(a) Dynamically stable.* The near-Nariai critical state is not stable: it decides within τ = 3.2 r_N/c (4×10⁻²⁴ s for 1.7×10¹¹ kg).
* *(b) Old enough to be found.* A Hawking lifetime above 13.8 Gyr needs M > 1.7×10¹¹ kg.
* *(c) Not pushed over by natural perturbations.*

**Assessment (the trigger cannot be machine-made negative energy).**
* In a normal precursor of NEC matter, the trigger must violate the NEC (Theorem 6h-1).
* A quantum field holds |E₋| ≲ ħc/R in a region R over ~R/c (Ford–Roman scaling).
* The tipping law needs 1.85×10⁻³ Rc⁴/G, which is short by (R/l_P)²/540: 10^6.9 at 10⁻³⁰ m, 10^30.9 at 10⁻¹⁸ m.
* So a machine could at most switch on a natural NEC sector (P-R6h-d). It cannot supply the violation itself.

**Question C (no classical trapped sphere).** Theorem 6i-11 is classical. Quantum channels are bounded only by R10/R42 in their controlled
regimes, and by the new throat-area law (section 4).

**Correction (F56).**
* `nariai.partner_vacuum` used ρ_f = 1/(24πM²). The Nariai condition 9M²Λ = 1 gives ρ_f = 1/(72πM²), which is what `accessibility.py` has
  always used.
* Corrected values: ρ_f^{1/4} = 0.329 GeV (M_☉/M)^{1/2}, previously 0.433. A GUT vacuum pairs with a 2.2 g hole (previously 3.7 g), an
  electroweak vacuum with 3.6 Earth masses (previously 6).
* The Round 6c report, results and test are updated. No conclusion changes.

## 2. Group A — the classical gravitational assumptions

* **#1 Non-uniform matter.** Covered by Proposition 6j-1: no gain.
* **#2 Non-spherical data.**
  * Riemannian Penrose inequality (Huisken–Ilmanen 2001; Bray 2001): for time-symmetric data with ρ ≥ 0,
    M_ADM ≥ √(A/16π)c²/G for the outermost minimal surface.
  * That is k_A ≥ 1, so the same invariant holds with the area radius (`penrose_form`) [Theorem, cited].
  * The dynamical (non-time-symmetric) case is the open Penrose conjecture (P-R6j-a).
* **#3 Field-dependent gravitational coupling, Proposition 6j-4.**
  * *What it buys.* In S ⊃ F(φ)R/16πG, a core where F = (M_P,eff/M_P)² ≪ 1 lowers the requirement measured from outside by F^{3/4}. Closing
    10^39.5 needs F ≈ 10⁻⁵³, an effective Planck mass of about 40 eV.
  * *What it costs.* The core needs the canonical excursion Δχ = √(3/16π)|ln F| M_P. By the capacitance (Dirichlet) minimum its gradient
    energy is E ≥ 2πLΔχ², checked numerically (`dirichlet_check`).
  * *Consequence.* The core's own compactness is ≥ (3/4)ln²F: horizon-class for F < 0.32. At F = 10⁻⁵³ it is 1.1×10⁴.
  * *Verdict.* Making the strong-gravity region is itself black-hole formation, before it has bought anything.
  * *Scope.* f(R) gravity is the same theory (F = f′(R)) and gets the same result.
* **#4 Higher-curvature gravity.**
  * f(R): see #3.
  * Gauss–Bonnet: Round 6i (activation pushes toward collapse).
  * Weyl² (Stelle) gravity has a massive spin-2 ghost; we exclude it on unitarity grounds.
  * No gain.
* **#5 Finite-lived pocket, Proposition 6j-5.**
  * *Why it collapses quickly.* A classical pocket must stay sub-Hubble (R < l_f); otherwise Theorems 6h-1 and 6i-11 apply. It collapses
    within ~R/c.
  * *How much it expands.* Its interior expands by a fraction ~R/l_f before then: 3×10⁻¹⁶ for a TeV vacuum at 10⁻¹⁸ m, 1.2×10⁻¹⁴ at 10 TeV.
  * *Verdict.* There is no measurable spacetime creation as a milestone.

## 3. Group B — let an environment prepare the state

* **#6 Accretion-driven vacuum formation, Proposition 6j-6.**
  * *The precedent.* Marsden, Aurrekoetxea, Clough & Ferreira (PRD 111, L041501, 2025) show that accretion onto black holes classically makes
    symmetry-restored shells and triggers decay bubbles that grow outward at nearly c. Their bubbles are the parent-eating direction.
  * *What a daughter needs.* A vacuum whose Hubble length matches the hole (the near-Nariai partner) needs infalling density
    ρ_f = c⁶/(72πG³M²) at the horizon.
  * *What accretion supplies.*
    * Eddington flows are short of that by 10^18.8 (10 M_☉) and 10^19.7 (1.4 M_☉).
    * Even neutron-star merger debris (0.1 M_☉ within 3 r_s) supplies only 0.25–1.8% of it.
  * *Completion.* Theorem 6h-1 holds for collapse-formed holes, so a daughter still needs the NEC escape.
  * *Verdict.* State preparation is natural at best, and completion is not provided.
* **#7 Horizonless compact precursor, Proposition-N 6j-7.**
  * *The calculation.* The most compact self-bound false-vacuum star (MIT-bag equation of state, TOV, `bag_star`) has 2GM/c²R = 0.54,
    radius 0.28 l_f, and a vacuum share of 14% of its mass.
  * *Consequence.* That is sub-Hubble, so its collapse makes a black hole (Theorem 6h-1).
  * *Scope.* The smooth version removes the thin-shell pathology of Round 6i, not the obstruction.
  * *An observation.* The QCD bag constant B^{1/4} = 145 MeV is exactly the Nariai partner of a 5.2 M_☉ hole. Stellar black holes and the
    QCD scale are matched by the same compactness arithmetic. But the bag "vacuum" exists only inside quark matter, which is not
    vacuum-dominated.
* **#8 Evaporating-hole catalysis (Oshita–Yokoyama, complete), Proposition-N 6j-8.**
  * *Both steps, computed.* A bubble is made thermally around a hole of mass M_in (weight e^{−ΔM/T_H}). It then tunnels to the hidden
    inflating branch, using the thin shell with the hole inside.
  * *The exponent.* Σ = Σ̂(l_f/l_P)². Σ̂ = 0.013 (wall k_w = 0.5) near the Nariai bound, down to ~4×10⁻⁴ for light walls.
  * *GUT vacuum.* (l_f/l_P)² ≈ 5×10¹⁰ gives Σ ≳ 2×10⁷.
  * *When daughters actually occur.* Over a hole's evaporation the expected number is ~1.6×10⁴M_in² e^{−Σ}. That is O(1) only for
    l_f ≲ 10²–10³ l_P: a near-Planckian vacuum and a hole at the end of its life.
  * *Verdict.* Natural, Planck-scale.
* **#9 Cosmological critical state.**
  * Round 6h (R137): a trigger of ~10% of the patch vacuum energy.
  * Patches are decided within about their own Hubble time, in the early universe. None is waiting today.

## 4. Group C — quantum spacetime without a classical throat

* **#11 Wormhole and no-boundary nucleation, Proposition 6j-9 (analytic for small throats, numerical otherwise).**
  * *The model.* Euclidean axion wormholes with a cosmological constant (Giddings–Strominger + Λ; the "wineglass" family of Lavrelashvili &
    Lehners, arXiv:2603.11003; published as "Birth of inflationary universes via wineglass wormholes and their no-boundary relatives", PRD,
    August 2026), with a′² = 1 − (a²/3)(Λ + q²/a⁶).
  * *The no-boundary limit, reproduced.* The half-wormhole action S(q) tends to −12π²/Λ as the throat pinches off.
  * *The throat-area law.* The excess over no-boundary is ΔS = 3π³a₀² exactly (8πG = 1), independent of Λ to 10⁻⁵ over four decades. In
    ordinary units, ΔS = (3π²/8)(a₀/l_P)² ≈ 1.18 × 4πa₀²/4G.
  * *Analytic for a₀ ≪ Λ^{−1/2}.* With a₀⁴ = q²/3 the on-shell integrand is −6a√(1 − a₀⁴/a⁴). Since
    ∫₁^∞(1 − √(y²−1)/y) dy = π/2 − 1, ΔS = 12π²(a₀²/2)(π/2) = 3π³a₀² [Theorem in that limit].
  * *Consequences.*
    * A channel that stays connected to us through a throat pays action proportional to the throat's area.
    * The throat-free limit is the disconnected universe, created independently of us.
    * This agrees with Lavrelashvili–Lehners' statement that no-boundary instantons dominate.
* **#10 Quantum transition without a classical throat.**
  * The throat-area law and Q* (Round 6i) give the same answer from two directions: suppression vanishes exactly when the connection
    vanishes.
  * A well-defined Hilbert space for the shell exists (Rounds 3–5). The obstacle is not the existence of a quantum description. It is that
    the connected amplitude carries e^{−O(A_throat/4G)}.
  * No connected, throat-free channel was found.
* **#12 Holographic baby universes.**
  * Sasieta, Swingle & Vilar López (PRL 137, 041501, 2026; arXiv:2512.00149) construct a baby universe behind a black-hole horizon in JT
    gravity, dual to a thermal pure state of two SYK models with a heavy operator.
  * Their baby universe lives behind a horizon (a bridge), consistent with Theorem 6i-11.
  * Its supporting degrees of freedom are the black hole's microstates. Running the dual quantum circuit simulates the geometry; it does not
    create spacetime.
* **#13 More quantum degrees of freedom.**
  * Adding collective coordinates can only lower the minimal-action path.
  * But the barrier next to the inflating branch is fixed by the connection (throat-area law) and by "no free lunch" (R10).
  * Open (P-R6j-c).

## 5. Group D — a different gravitational theory

* **#15 Quantum-corrected interiors, Proposition 6j-10 (Frolov–Markov–Mukhanov limiting curvature).**
  * *The junction.* The Schwarzschild interior meets de Sitter (Hubble length l) across a spacelike surface at r₀ = (2Ml²)^{1/3}. There
    K_θθ is continuous (f_S = f_dS), the Kretschmann invariant is 12/l⁴ (half the de Sitter value), and M > l/2 is needed.
  * *The layer.* The jump of K_t^t gives a layer with tangential stress S_θ^θ = 0.06–0.18/l.
  * *The daughter.* The de Sitter region starts larger than its Hubble length (r₀/l = (2M/l)^{1/3}) and inflates.
  * *Consequence.* If quantum gravity behaves this way, every black hole with M > l/2 contains a daughter, made from the limiting-curvature
    vacuum. Making a daughter then reduces to making a black hole.
  * *Status.*
    * The NEC violation that Theorem 6h-1 demands sits in the quantum-gravity layer.
    * Competing models disagree. Regular black holes with static de Sitter cores have no daughter. Inner-horizon mass inflation questions
      long-lived regular interiors (Carballo-Rubio et al. 2022; PRL 2024).
* **#14 Dimensional crossover, Proposition 6j-11 (conditional).**
  * *What changes.* With n large extra dimensions and fundamental scale M_D, black holes form at M_BH ≳ few M_D. The 4D bridge invariant
    is replaced by a collider threshold.
  * *What is still allowed.*
    * LHC monojet: M_D ≳ 10.7–11.2 TeV (n = 2) to 5.2–5.9 TeV (n = 6–7).
    * CMS black-hole searches: semiclassical holes excluded below ~8.4–11.4 TeV.
    * Astrophysical bounds are strong for n = 2, 3 and weak for n ≥ 4.
  * *Collider reach.* A 100 TeV collider reaches M_D of order 15–20 TeV (scaling estimate). Geometric parton cross-sections are 15–40 pb at
    25–40 TeV.
  * *The catch.* R_S grows only as M^{1/(n+1)}. At M_BH ≈ 5 M_D, R_S M_D ≈ 1.0 (n = 4–7). Collider black holes are objects at the
    fundamental scale, with no semiclassical interior.
  * *What a semiclassical interior needs.* R_S = 3/M_D needs M_BH ≈ 960 M_D (n = 4) to 10⁴ M_D (n = 6), i.e. 10–100 PeV collisions.
    Cosmic rays reach ~0.75 PeV.

## 6. Synthesis: the doors that remain

| Door | What it needs | Testable? | Status |
|---|---|---|---|
| **Y1: FMM-type interiors** (#15) | Black-hole interiors that turn into inflating de Sitter at a limiting curvature | Only through quantum-gravity theory | Conjecture. In 4D, humanity still has to make a black hole (invariant). |
| **Y2: strong gravity at 6–20 TeV** (#14) | n ≥ 4 large extra dimensions, M_D ≲ 20 TeV | **Yes**: a 100 TeV collider probes it | Not excluded. Collider holes are fundamental-scale objects, and what they do inside is quantum gravity. |
| **Y1 + Y2** | Both, plus interiors at M_BH ≫ M_D | Partly | Would turn "make a daughter" into "make a big enough black hole": 10–100 PeV collisions |
| **Q*** (Round 6i) | A barely metastable vacuum; channels with sub-Planckian bridges | No | Outside every controlled regime; the throat-area law says connection costs action |
| Natural channels | Cosmological patches, evaporating near-Planckian holes | — | Nature, not humanity |

**What this round establishes.**
* Every remaining "yes" runs through one of two changes: gravity becomes strong at an accessible energy, or quantum gravity supplies the
  topology change and the NEC violation in black-hole interiors.
* The first is an experimental question that a 100 TeV collider can address.
* The second is a theoretical question about black-hole interiors.
* Within established physics (4D GR, quantum fields, semiclassical gravity), no route reaches humanity's scale. The invariant bounds stored
  energy, and the supplied energy can be small only with natural precursors that the computations above rule out.

## 7. Status

| Research objective | Status after Round 6j |
|---|---|
| Create the initial false-vacuum region | **Yes, by the universe** (Round 6g) |
| Establish a technologically accessible mechanism | **Not found within established physics.** The invariant is a bound on stored energy (Prop. 6j-2). Small supplied energies need natural precursors that are stable, old and within 10⁻⁸ of threshold: none exists (Prop. 6j-3). Classical loopholes cost more than they buy (Props. 6j-1, 6j-4, 6j-5). A connected quantum channel pays the throat-area action (Prop.-N 6j-9). **Two conditional doors remain**, both in quantum gravity, one of them testable at a 100 TeV collider (Props. 6j-10, 6j-11). |

## Literature used (Round 6j)

* *Accretion*
  * Marsden, Aurrekoetxea, Clough & Ferreira, PRD 111, L041501 (2025), arXiv:2403.17595: symmetry restoration and vacuum decay from accretion
    around black holes.
* *Wormholes and no-boundary*
  * Lavrelashvili & Lehners, arXiv:2603.11003 (March 2026); PRD (10 August 2026) "Birth of inflationary universes via wineglass wormholes and
    their no-boundary relatives". The Physics Letters B version cited in the request could not be found.
  * "Consistent evaluation of the no-boundary proposal", arXiv:2602.02682 (authors not confirmed here).
* *Holographic baby universes*
  * Sasieta, Swingle & Vilar López, PRL 137, 041501 (2026), arXiv:2512.00149: baby universes from thermal pure states in SYK.
* *Extra dimensions*
  * CMS, arXiv:2107.13021 (monojet; M_D > 10.7 TeV at n = 2 to 5.2 TeV at n = 7).
  * ATLAS, arXiv:2102.10874 (M_D > 11.2 TeV at n = 2, 5.9 TeV at n = 6).
  * CMS, arXiv:2604.10732 (black holes excluded below 8.4–11.4 TeV).
  * Hannestad & Raffelt, hep-ph/0304029 (astrophysical bounds).
  * Gamsızkan, arXiv:2404.10070 (black holes at 27 and 100 TeV).
  * Anchordoqui et al., hep-ph/0112247 (cosmic-ray black holes).
* *Penrose inequality*
  * Huisken & Ilmanen, J. Diff. Geom. 59, 353 (2001); Bray (2001); review: Bray & Chruściel, gr-qc/0312047; Bray & Khuri, arXiv:0905.2622.
* *Regular black holes*
  * Frolov, Markov & Mukhanov (1989/1990).
  * Carballo-Rubio, Di Filippo, Liberati, Pacilio & Visser, JHEP 09 (2022) 118, arXiv:2205.13556.
  * "Mass inflation without Cauchy horizons", PRL (2024).
