# Round 6g — the two remaining "no"s: creating the false-vacuum region, and a technologically accessible mechanism

Request (user, after Round 6f):

> Turn the last two rows of the status table into yes. Don't introduce new problems with them. You can invent a new machine, but it has to
> be doable within humanity.

| Research objective | Status before | Status after Round 6g |
|---|---|---|
| Create the initial false-vacuum region | Thermal route fails | **Yes, by the universe** (not by us). A supercooled first-order transition leaves false-vacuum patches. Simulated with gravity, patches larger than ≈ 0.58 l_H inflate into daughters hidden behind a black hole, with the GSL satisfied. The near-critical ones are exactly the route's ball. |
| Establish a technologically accessible mechanism | Current route exceeds accessible scales | **No, and now proved impossible for any machine.** Building any near-Nariai false-vacuum region needs ≳ 10⁵² W (or ≳ 10⁴³ N of force): 10³⁹ times all of humanity's power, above the conjectured maximum luminosity of any physical system. This is a property of gravity, not of present technology. |

The second row cannot honestly be made a yes. The rest of this report shows why, and audits every machine that could be proposed.

New code:
* `formation.frw_cone_data` (radiation-dominated FRW data on a light cone, in conformal or affine gauge) and `formation.cosmic_patch`
  (a supercooled false-vacuum patch in an expanding universe, evolved with gravity, with an outcome classifier). The centre stage now
  reports end codes.
* `accessibility.assembly_bounds`, `quantum_requirement`, `machine_audit`, `tev_gravity_conditional`.
* `scripts/round6g.py`, `tests/test_round6g.py` (5 tests), `DATA/round6g.json`, `FIGURES/round6g.png`.

Labels as before: **Theorem**, **Proposition**, **Proposition-N** (numerical), **Literature**, **Assessment**.

---

## A. The false-vacuum region, made by the universe

### A.1 Why the laboratory question had the wrong starting point

Every creation route tried so far (heating in Round 6f; the lift and the confined fireball in Round 6e) starts from flat, empty, true
vacuum and tries to put a false vacuum back. Two results block that:
* Farhi–Guth: a classical route from regular data needs anti-trapped surfaces.
* Theorem 6g-2 below: it needs Planck-class power or force.

The early universe does not start from flat vacuum.
* **It starts in the false vacuum.** At high temperature the symmetric phase is the vacuum everywhere.
* **It is already expanding.** Every region larger than the Hubble radius is anti-trapped.
* **Nothing has to be assembled.** When the transition proceeds by bubbles, regions where no bubble has nucleated yet stay in the false
  vacuum (supercooled), while the radiation around them is diluted by the expansion.

This "delayed decay" is studied as a source of primordial black holes (Kawana, Kim & Lu 2022; Liu, Bian, Cai, Guo & Wang 2022). Bubbles
nucleated during inflation and later left supercritical become baby universes behind black holes (Garriga, Vilenkin & Zhang 2016; Deng &
Vilenkin 2017, with numerical GR simulations).

### A.2 The simulation

* **Initial data.** A light cone from the centre of a patch in a radiation-dominated universe (`frw_cone_data`).
  * Total density is uniform, ρ₀, with the Hubble flow everywhere (H₀ = 1/η₀).
  * Inside comoving radius χ_p: false vacuum ρ_f, plus the rest as radiation.
  * Outside: true vacuum plus radiation (the latent heat).
  * Between them: the wall, at its own width.
  * f = ρ_f/ρ₀ is the vacuum fraction.
* **Matter content.** The same vacuum field and potential as the route (l_H = 30, v = 0.0015, ε = 4), a pure-radiation plasma, and
  gravity.
* **Evolution.** From the regular centre, in the affine cone gauge.

**Verification (Proposition-N).**
* Without a patch, the solver reproduces FRW. In conformal coordinates, r = aχ and s = ln a are matched to 6×10⁻³ → 4×10⁻³ → 2×10⁻³ and
  1.4×10⁻³ → 7×10⁻⁴ → 3×10⁻⁴ (h = 0.2, 0.1, 0.05): first order, as expected from the θ-weighted fluid scheme.
* In the affine gauge, the Misner–Sharp mass equals (4π/3)ρr³ to 9×10⁻³ → 4×10⁻³ (median; h = 0.1, 0.05).
* With patches, the median relative constraint residuals are 10⁻⁶–5×10⁻⁴.

### A.3 Proposition-N 6g-1 (critical size)

A patch either collapses into a black hole with no daughter (a trapped region forms and the false vacuum disappears), or its false vacuum
inflates (anti-trapped, with r growing by 10²–10⁴ inside the computed domain). The boundary between the two is:

| vacuum fraction f | collapses | inflates | critical radius / l_H |
|---|---|---|---|
| 0.9 (h = 0.05 and 0.1) | ≤ 0.569 | ≥ 0.588 | 0.57–0.59 |
| 0.5 | 0.42 | ≥ 0.566 | 0.42–0.57 |
| 1.0 | 0.55 | ≥ 0.65 | 0.55–0.65 |

* **Resolution.** The critical size is the same at h = 0.1 and h = 0.05.
* **The route's ball.** The ball that Rounds 6d–6f started from has radius 0.595 l_H, inside the critical band. The universe therefore
  produces that ball as its near-critical patches. Patch sizes are distributed around the Hubble radius at the transition (Kawana–Kim–Lu),
  so a population near any size exists.

### A.4 Proposition-N 6g-3 (the daughter is hidden and the GSL holds)

For f = 0.5 with a patch of 0.57 l_H, the false vacuum inflates (r up to 2.6×10³ in the domain).
* **Black hole outside the wall.** From u ≈ 30, a trapped true-vacuum region forms just outside the wall, at r ≈ 38. On 188 outgoing rays,
  the inflating false vacuum lies behind this trapped region: the daughter is behind a black hole as seen from the exterior universe.
  This is the Garriga–Vilenkin–Zhang structure, found here by the project's own solver.
* **Parent's event horizon.** It forms at u = 29.85. Its area never decreases (deficit 0).
* **Daughter observers.** Over the sampled daughter observers, the largest decrease of horizon area is 0.

**Comparison with Theorem 6e-6.** That theorem showed that an eternal daughter grown from normal data always violates the local GSL. Its
premise was normal data. Cosmological data are anti-trapped, and here the GSL holds even for the eternal daughter. The GSL problem of the
laboratory route does not arise in the cosmological route.

### A.5 Proposition-N 6g-4 (the condensate widens the window)

With the gentler potential ε = 2, a patch of 0.57 l_H collapses without the tilt. With the condensate tilt (λ_f = 10⁻⁶, which also
self-generates the NEC-phase charge), it inflates. The route's mechanism therefore also tips naturally made near-critical patches into
daughters.

The effect has limits:
* A smaller patch (0.53 l_H) still collapses with the same tilt.
* With ε = 4 the tilt does not tip a 0.57 l_H patch.

So the condensate lowers the critical size somewhat; it is not needed for supercritical patches.

### A.6 Does this introduce new problems?

| possible problem | status |
|---|---|
| Farhi–Guth / Penrose obstruction | Not violated. The cosmological data are anti-trapped. Their past singularity is the Big Bang (or inflation), as in standard cosmology. |
| GSL | Satisfied: parent event horizon non-decreasing, daughter observers' horizons non-decreasing (A.4). |
| Normal data, a black hole on the data | Not required: the patch is part of an expanding universe; the black hole forms during the evolution. |
| Hiding | The daughter is behind a trapped region (A.4). |
| The condensate's observational bound (problem 4) | Avoided for supercritical patches, which need no condensate. Only a first-order transition with strong supercooling is needed: standard beyond-Standard-Model physics, not established, testable. |
| Observable consequences | From outside, each daughter is a primordial black hole of about the horizon mass at the transition. Primordial black-hole abundances and the stochastic gravitational-wave background from the transition constrain the population. These are the same tests as for ordinary delayed-decay black holes (Phase VI: the daughter itself is invisible). |
| Numerical limits | The inflating interior becomes stiff on the grid (as in F45), so the computed domain ends after r grows by 10²–10⁴. In the f = 0.9 supercritical runs the exterior trapped region had not yet formed within the domain. For f = 0.5 it had. |

**Verdict for the creation row: yes, by nature.** The initial false-vacuum region is created by a cosmological first-order phase
transition. Within this project's own GR solver it then becomes a hidden daughter (or, near-critical, the route's ball). No new
problem is introduced. It is not created by humans: see part B.

---

## B. A technologically accessible mechanism: why no machine can do it

### B.1 Theorem 6g-2 (Planck power or Planck force)

* **Setting.** Let a region of radius R hold energy E with compactness k = 2GE/(Rc⁴). The route's near-Nariai ball, and any
  daughter-forming throat, has k ≈ 2/3. Such a region is gravitationally unstable, and the false vacuum is also pulled in by its own
  negative pressure. It evolves in t_dyn = αR/c, with α = O(1).
* **Statement.** A machine that builds it must do one of two things:
  1. Deliver E within t_dyn, which needs a power P ≥ E/t_dyn = (k/2α) c⁵/G ≈ 1.2×10⁵² W (k = 2/3, α = 1).
  2. Hold the region static while assembling it. The stress across a cross-section is then ~GE²/(c⁴R²) = (k²/4) c⁴/G ≈ 1.3×10⁴³ N.
* **Independence.** Both numbers are fixed fractions of the Planck power c⁵/G and the Planck force c⁴/G, whatever the vacuum energy: a
  denser vacuum gives a smaller region, but the same power and force.
* **Comparison.** Gibbons' maximum tension is c⁴/4G = 3×10⁴³ N; the conjectured maximum luminosity of any system is c⁵/4G = 9×10⁵¹ W.
* **Proof.** Energy conservation over t_dyn, and the static Tolman–Oppenheimer–Volkoff-type stress estimate for a region at
  compactness k. Order-of-magnitude in α, exact in its scaling.

| power | value |
|---|---|
| all of humanity's primary energy use | 1.9×10¹³ W |
| a hypothetical Dyson sphere around the Sun | 3.8×10²⁶ W |
| peak of the black-hole merger GW150914, the most powerful event ever measured | 3.6×10⁴⁹ W |
| conjectured maximum luminosity c⁵/4G | 9×10⁵¹ W |
| **needed to build the ball** | **1.2×10⁵² W** (10^38.8 × humanity) |

### B.2 Proposition 6g-5 (every machine, including invented ones)

Energy carried by quanta of energy E_q can be localised only to ħc/E_q. The invariant of Theorem 6f-7 then becomes
E_total · E_q ≥ ħc (E√ρ)_req^{2/3} = 6.3×10¹⁷ J².

| machine | E_total · E_q (J²) | shortfall |
|---|---|---|
| LHC: the whole stored beam, focused ideally | 4×10² | 10^15.2 |
| FCC-hh (100 TeV design): the whole beam, focused ideally | 7×10⁴ | 10^13.0 |
| NIF laser | 1×10⁻¹² | 10^29.7 |
| 50 Mt thermonuclear explosion | 3×10¹ | 10^16.3 |
| **Invented:** one year of all world energy in 100 TeV quanta, focused ideally | 1×10¹⁶ | 10^1.8 |

* **The invented machine.** It needs ~4×10²⁵ coherent 100 TeV quanta in a single 2×10⁻²¹ m spot within 7×10⁻³⁰ s. That means delivering
  9×10⁴⁹ W, more than the GW150914 peak. No focusing principle allows this. It is Theorem 6g-2 again, in quantum form.
* **A single collision.** It needs E_q ≥ 3.5×10¹⁸ GeV per particle. That is a 1.1-parsec-long linear collider at 100 GV/m, the
  plasma-wakefield record class, or a ring of 150 parsecs with 16 T magnets.

### B.3 Assessment 6g-6 (the loopholes, each checked)

* **Quantum tunnelling (Farhi–Guth–Guven).** It allows a small exterior mass, but the project's frontier is P ≤ 10⁻²³ even with Planck-scale
  energies in 100 l_P (R42, R50), and astronomically smaller at human energies.
* **Gravity strong at the TeV scale (large extra dimensions).** The Planck power and force would drop to M_D²c⁴/ħ-class values, which a
  collider reaches in one collision.
  * CMS excludes semiclassical black holes up to 10.1 TeV (JHEP 11 (2018) 042). A 100 TeV collider could probe further.
  * This would require unestablished physics, and the whole route rebuilt in higher dimensions, with a TeV-scale false vacuum and
    condensate. That is a set of new problems, so it does not qualify.
* **Analog systems** (false-vacuum decay in cold atoms and superfluids, analog horizons). These are doable and could test bubble and wall
  dynamics, but they do not create spacetime.
* **Steering nature.** Using an astrophysical collapse, or the universe's own phase transitions, is the only process that reaches the
  required power. It is not under human control.

**Verdict for the accessibility row: no, and it cannot be made a yes.** The obstacle is fixed by G and c, through the Planck power and
force, which no machine of any size or technology can reach. Only gravitational collapse and cosmology do. Calling it "a technology we
don't have yet" would be wrong.

---

## C. Updated status of all research objectives

| Research objective | Status |
|---|---|
| Model an expanding daughter spacetime | Achieved in simulations |
| Form a black hole dynamically | Achieved in simulations |
| Generate the necessary negative-energy pulse using a field | Achieved in the condensate model |
| Establish nonlinear condensate stability | Within the investigated regime (Round 6f) |
| Reproduce the principal result with two solvers | Agreement around 10⁻⁵ |
| Create the initial false-vacuum region | **Yes, by the universe** (cosmological first-order transition; simulated with gravity, hidden daughter, GSL holds) |
| Establish a technologically accessible mechanism | **No: proved impossible for any machine** (Planck power / force, Theorem 6g-2) |

## Literature used (Round 6g)

* Garriga, Vilenkin, Zhang (2016) JCAP 02, 064: black holes and the multiverse.
* Deng, Vilenkin (2017) JCAP 12, 044, arXiv:1710.02865: primordial black-hole formation by vacuum bubbles; supercritical bubbles become
  baby universes; numerical simulations.
* Kawana, Kim, Lu (2022) arXiv:2212.14037: PBH formation from overdensities in delayed vacuum transitions.
* Liu, Bian, Cai, Guo, Wang (2022) PRD 105, L021303: PBHs from first-order phase transitions (delayed decay).
* Gibbons (2002) Found. Phys. 32, 1891, hep-th/0210109: maximum tension c⁴/4G; the maximum luminosity c⁵/4G conjecture (Dyson).
* LIGO/Virgo GW150914: peak luminosity 3.6×10⁴⁹ W.
* CMS (2018) JHEP 11, 042, arXiv:1805.06013: semiclassical black holes excluded up to 10.1 TeV.
* Plasma-wakefield gradients: LBNL 7.8 GeV in 20 cm (≈ 400 GV/m average); AWAKE 0.2–1 GV/m.
