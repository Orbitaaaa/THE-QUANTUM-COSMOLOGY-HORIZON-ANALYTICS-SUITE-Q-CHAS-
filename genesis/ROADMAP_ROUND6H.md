# Round 6h — four routes toward a controllable mechanism, and what actually blocks one

The user reviewed Round 6g and pointed out three problems:
* its impossibility theorem rests on a dichotomy (rapid external delivery or static holding) that is not exhaustive;
* the maximum force and luminosity are conjectures with published counterexamples;
* the report's own numbers do not exceed those conjectures robustly.

The user proposed four routes to investigate:

| Route | Idea |
|---|---|
| A. Trigger a near-critical background | a small machine flips a configuration that is already almost critical |
| B. Black-hole-assisted formation | Oshita–Yokoyama: the black hole helps create the false vacuum and the wormhole |
| C. Dynamical assembly | the system becomes compact through its own evolution, with neither rapid injection nor static holding |
| D. Different gravitational physics | change the compactness law (e.g. TeV-scale gravity) |

Labels as before (Theorem, Proposition, Proposition-N, Literature, Assessment, Conjecture).

---

## 1. The correction (accepted)

* **The dichotomy is not exhaustive.** Proposition 6g-2 (formerly "Theorem") bounds two preparation strategies. It does not cover a
  dynamical history, or a small trigger on a pre-existing background.
* **The Planck limits are conjectures.**
  * Gibbons' maximum tension c⁴/4G: Jowsey & Visser (2021) give counterexamples in GR with fluid spheres near collapse. A later comment
    disputes them, and Faraoni argues that a black-hole-restricted version may survive.
  * The maximum luminosity c⁵/4G: Cardoso, Ikeda, Moore & Yoo (2018) find GR initial data with arbitrarily large luminosity. Requiring no
    past horizon seems to keep it below the Planck value, and critical collapse reaches about 0.2 of it.
* **The numbers do not exceed the conjectures robustly.**
  * The holding stress, 1.3×10⁴³ N, is below c⁴/4G = 3.0×10⁴³ N.
  * With α = 1.5 the power is 8.1×10⁵¹ W, below c⁵/4G = 9.1×10⁵¹ W.

`accessibility.assembly_bounds` and the Round 6g documents now state the bound as conditional, and the absolute claim is withdrawn (F51).

## 2. Theorem 6h-1 (the no-go that does hold for every strategy)

* **Statement (Penrose 1965; Farhi & Guth 1987; non-spherical through the time-reversed Penrose theorem).**
  * **Premises.** A spacetime that is globally hyperbolic, has non-compact Cauchy surfaces and obeys the null energy condition.
  * **The theorem.** If such a spacetime contains a closed anti-trapped surface, it is past null incomplete.
  * **The step to daughters.** An inflating daughter region contains closed anti-trapped surfaces: spheres larger than its Hubble radius.
  * **Conclusion.** A daughter can arise only with at least one of three ingredients:
    * (i) a violation of the null energy condition;
    * (ii) a past singularity, as in cosmology (the Big Bang);
    * (iii) a quantum transition.
* **Consequence.** A machine built of ordinary (NEC-respecting) matter, acting on a universe with a regular past, cannot create a daughter.
  This holds whatever its power, whatever its duration, and whatever the strategy: rapid, slow, dynamical, staged or triggered.
* **Wording sharpened in Round 6i (section 8).** "Regular past" is local: the argument closes once the past null generators of the anti-trapped
  surface reach their focal points, within an affine length ~R, which holds for any laboratory region in our (past-incomplete) universe.
  Round 6i adds the bridge invariant (Theorem 6i-11), which covers NEC-violating matter as well.
* **What it does not depend on.** It uses no assembly time, no maximum force and no maximum luminosity. That is why it, and not
  Proposition 6g-2, is the rigorous classical statement. Stage One already stated the obstruction precisely; this is its form for machines.
* **The project's escapes.**
  * (i) The ghost condensate (Rounds 6c–6f). Quantum fields cannot supply the classical violation, because of quantum energy
    inequalities (Round 6c).
  * (ii) The cosmological route (Round 6g).
  * (iii) The quantum frontier, P ≤ 10⁻²³ (R42).

## 3. Route A — trigger a near-critical background

* **In a normal (non-anti-trapped) background.**
  * By Theorem 6h-1, a trigger made of ordinary matter cannot turn collapse into a daughter.
  * The tipping law of Round 6c (R88) says the same: daughter iff W = −∫N e^{−u/τ}du > W_th > 0, which needs a net negative-energy flux.
  * In the route, the trigger is the condensate's own charge. That charge is made by the false-vacuum tilt, not by a machine
    (Theorem 6e-1: only explicit shift breaking makes charge). A machine could make it only if the condensate coupled to controllable
    matter, which is open.
* **In an anti-trapped (cosmological) background (Proposition-N 6h-2).**
  * An ordinary positive-energy trigger works. Here it enlarges a supercooled patch by converting a thin shell of true vacuum to false
    vacuum, which costs ρ_f per volume.
  * A patch of 0.60 η₀ clearly collapses into a black hole; one of 0.62 η₀ clearly inflates into a hidden daughter. The trigger is
    (0.62/0.60)³ − 1 ≈ 10% of the patch's vacuum energy.
  * The fine bisection puts the critical radius at 0.6175 η₀ (0.586 l_H) at h = 0.1. Close to it the outcome classifier is limited by the
    computed domain (0.615 η₀ is "undecided" at h = 0.05), so we can show the trigger shrinking only to about 1% in radius. Critical
    behaviour suggests it goes to zero at the threshold, but we do not claim that.
* **What humanity lacks.**
  * Our local spacetime is normal, not anti-trapped.
  * Supercooled patches existed only in the early universe, and the subcritical ones collapsed into black holes. In our runs, their false
    vacuum shrinks to nothing inside the horizon.
  * No near-critical background is known to be accessible today. A primordial black hole that still holds a near-critical false-vacuum
    interior could in principle be triggered by an infalling probe. Nothing suggests such interiors survive (P-R6h-a).

## 4. Route B — black-hole-assisted formation (Oshita–Yokoyama 2018)

* **The mechanism (Literature).**
  1. Hawking radiation thermally creates a thin-walled false-vacuum bubble around a black hole.
  2. The bubble wall triggers tunnelling to a wormhole-like configuration whose far side inflates.
  3. This happens once the black hole's mass is below a critical value of the order of the false vacuum's Hubble length. That is a
     black hole of the same mass class as the Nariai mass M_N.
  4. The rate is m e^{−B} with a Euclidean action B (Gregory–Moss–Withers).
* **Gaps.** The sources reachable here give no numerical probability, and the full paper could not be fetched from this session.
* **Assessment.**
  * It needs a black hole of about M_N in the final stage of evaporation.
  * Such holes are primordial (natural) if they exist. Making one is itself a horizon-class assembly: the same invariant, or TeV-scale
    gravity (route D).
  * The process is quantum and suppressed by e^{−B}.
  * So it is a possible natural channel, given primordial black holes and a false vacuum, but not a human mechanism.
* **Related work.** Oshita, Yamada & Yamaguchi (2019) show that horizonless compact objects catalyse downward vacuum decay efficiently.
  That is the opposite direction to the upward transition needed here.

## 5. Route C — dynamical assembly

* **Ordinary matter.** Theorem 6h-1 excludes it, whatever the dynamics.
* **With the condensate.** The route needs a near-Nariai false-vacuum ball, which is normal (not anti-trapped), so Theorem 6h-1 does not
  forbid it.
  * The ball necessarily holds mass M_N = r_N c²/(3G) within ~r_N. This is a geometric property of the near-Nariai throat (compactness
    ≈ 2/3), true for any history.
  * Gravity can supply the energy dynamically: collapse converts its own binding energy, with no external Planck power. But natural
    collapse compresses only stellar masses and above.
  * With the condensate observation allows (universal, M ≤ 100 GeV), M_N ≥ 2.5×10⁹ M_☉ (R128).

**What a smaller background would take** (`background_requirements`). The table assumes a condensate at a scale that observation forbids
for the universal condensate:

| false-vacuum scale | r_N | M_N | smallest localisable quantum ħc/r_N | quanta needed | time available r_N/c |
|---|---|---|---|---|---|
| 1.5×10¹⁶ GeV | 2×10⁻³⁰ m | 1 g | 9×10¹³ GeV | 6×10⁹ | 7×10⁻³⁹ s |
| 10¹⁴ GeV | 5×10⁻²⁶ m | 22 kg | 4×10⁹ GeV | 3×10¹⁸ | 2×10⁻³⁴ s |
| 10¹⁰ GeV | 5×10⁻¹⁸ m | 2×10⁹ kg | 41 GeV | 3×10³⁴ | 2×10⁻²⁶ s |
| 10³ GeV | 0.5 mm | 2×10²³ kg | 4×10⁻¹³ GeV | 3×10⁶² | 2×10⁻¹² s |

* **The trade-off.** Lower vacuum scales need gentler quanta but enormous masses. Higher scales need little energy, but every quantum must
  carry GUT-scale energy and the whole must arrive within ~10⁻³⁸ s.
* **Why this is hard.** For any machine that does not use gravity itself, that arrival time is what forces the conditional power bound.

**Conjecture-level escape (a new model, not a result).** Suppose the condensed phase exists only inside the false vacuum: a kinetic
function P(X, ϕ) whose condensed (healthy) phase is selected by the vacuum field. The observational bound constrains only the condensate in
our vacuum, so it would no longer apply, and the background could be grams to kilograms for high vacuum scales. Before this could count,
each of these must be solved:
* a ghost-free interpolation through the wall;
* the condensate's behaviour in the wall;
* the existence of a high-scale false vacuum;
* the creation of the ball itself (heating fails, Round 6f; the coherent lift fails, F40).

## 6. Route D — different gravitational physics

* With gravity strong at a scale M_D (large extra dimensions), horizon-class objects become collider-accessible above M_D.
* CMS excludes semiclassical black holes up to 10.1 TeV.
* The whole route would have to be rebuilt in higher dimensions, with a false vacuum and a condensate at that scale.
* It remains a conditional possibility, with new problems by construction.

---

## 7. What actually blocks a controllable mechanism

1. **Theorem 6h-1 (rigorous).** No machine of ordinary matter can create a daughter from a universe with a regular past, by any strategy.
   A mechanism needs NEC violation, a cosmological past singularity, or a quantum transition.
2. **With the condensate (the NEC escape).** A horizon-class background must exist: a near-Nariai false-vacuum ball of mass r_N c²/(3G).
   This is geometric and history-independent.
3. **With the condensate observation allows.** That background is ≥ 2.5×10⁹ M_☉.
4. **A smaller background** needs an unestablished false-vacuum-only condensate, plus a way to make the ball. No mechanism is known.
   External assembly within r_N/c costs ~c⁵/G-class power (a conditional bound), and natural collapse compresses only stellar masses and
   above.
5. **Ordinary-energy triggers** work only in anti-trapped (cosmological) backgrounds (Proposition-N 6h-2). Humanity has no access to one.
6. **Quantum and black-hole-assisted routes** need horizon-class black holes of the same mass class and are exponentially suppressed.

**Status.**

| Research objective | Status |
|---|---|
| Create the initial false-vacuum region | **Yes, by the universe** (Round 6g) |
| Establish a technologically accessible mechanism | **Not found.** Not proved impossible for every conceivable process (the Round 6g claim is withdrawn). Proved: no ordinary-matter machine can do it from a regular past (Theorem 6h-1). With the only known stable NEC violation at its allowed scale, the background is ≥ 10⁹ M_☉. Every known escape from that needs new, unestablished physics. |

## Literature used (Round 6h)

* Penrose (1965) PRL 14, 57; Farhi & Guth (1987) PLB 183, 149: closed trapped or anti-trapped surfaces and incompleteness; a growing
  false-vacuum bubble from regular data needs an initial singularity.
* Jowsey & Visser (2021) arXiv:2102.01831: counterexamples to the maximum force conjecture (disputed: arXiv:2109.07700; Faraoni,
  arXiv:2105.07929).
* Cardoso, Ikeda, Moore & Yoo (2018) PRD 97, 084013, arXiv:1803.03271: arbitrarily large luminosities from GR initial data; limited when
  past horizons are excluded; critical collapse ~0.2 Planck. Jowsey & Visser (2021) arXiv:2105.06650.
* Oshita & Yokoyama (2018) PLB 785, 197, arXiv:1601.03929: an inflationary universe out of an evaporating black hole, once its mass is
  below ~ the false vacuum's Hubble scale.
* Oshita, Yamada & Yamaguchi (2019) PLB 791, 149, arXiv:1808.01382: compact objects as catalysts for (downward) vacuum decay.
* CMS (2018) JHEP 11, 042: black holes excluded up to 10.1 TeV.
