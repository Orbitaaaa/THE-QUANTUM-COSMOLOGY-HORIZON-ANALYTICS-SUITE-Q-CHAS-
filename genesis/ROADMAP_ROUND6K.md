# Round 6k — the two doors, worked through: a conditional "yes"

The user asked us to treat the two remaining doors as hypotheses that change what a machine has to do, and to work toward the "yes":

| Door | Hypothesis | Missing discovery |
|---|---|---|
| Y1 | Black-hole interiors become expanding daughter spacetimes | A physically correct theory that produces such interiors |
| Y2 | Gravity is strong at particle-physics energies (extra dimensions), so microscopic black holes become accessible | The model realized in nature, and its black holes holding daughters |

The user also asked us to follow better leads wherever they appear.

Everything below is computed or checked against the literature:
* Code: `doors6k.py`, `scripts/round6k.py`, `tests/test_round6k.py` (6 tests).
* Data and figure: `DATA/round6k.json`, `FIGURES/round6k.png`.

**The answer.** Under two premises, both currently allowed by every experiment, the answer is yes. A machine humanity has already designed,
a 100 TeV proton collider (FCC-hh class), would make black holes whose interiors hold expanding daughter spacetimes. Each premise is stated
precisely below, together with its status. One is decided by that same collider; the other is a property of the correct theory of
gravity. The chain is not established, because its premises are not. It is also not excluded: no theorem, bound or experiment we found
rules it out.

## 1. Y1 — which definite theories put a daughter inside a black hole?

| Framework | What the interior becomes | Daughter | Standing |
|---|---|---|---|
| **Einstein–Cartan–Sciama–Kibble** (torsion sourced by fermion spin; Popławski 2010–2025) | collapse bounces; a new closed universe expands beyond the horizon; cycles grow with particle production until dark energy dominates | expanding | a classical theory, identical to GR below extreme densities; passes all tests |
| **Limiting curvature** (Frolov–Markov–Mukhanov 1989/90; Barrabès–Frolov 1996; Easson–Brandenberger 2001) | spacelike transition to de Sitter; possibly several disconnected de Sitter worlds | inflating | effective theory (higher-derivative constructions exist) |
| **Mimetic limiting curvature** (Chamseddine–Mukhanov 2017, 2019) | nested Schwarzschild "dolls" or a static de Sitter patch; stable remnants after evaporation | static core | effective theory |
| **Loop-quantum-gravity effective geometry** (Ashtekar–Olmedo–Singh 2018) | transition surface → white-hole region → new asymptotic regions | new regions, not inflating | effective quantum geometry |
| **Holography** (Cooper et al. 2018; Sasieta–Swingle–Vilar López 2026) | certain microstates hold a closed big-bang/big-crunch universe behind the horizon | crunching | controlled quantum gravity (AdS, JT/SYK) |

**What we computed.**
* **ECSK bounce [Proposition].**
  * For an unpolarised spin-½ fluid, ε_eff = ε − κs²/4 with s² = ħ²c²n²/8 and κ = 8πG/c⁴.
  * Matter of particle mass m therefore bounces at n_b = 4mc⁴/(πGħ²), i.e. ρ_b = 4m²c⁴/(πGħ²) = 3.9×10⁵⁸ kg/m³ for nucleons.
  * That is 10⁻³⁸ of the Planck density, so the bounce is semiclassical. Popławski quotes deviations from GR above ~10⁴⁵ kg/m³, a
    different matter content; we give both.
  * A horizon forms before the bounce for any collapsing mass above M_min = (3/4πρ_b)^{1/2}(c²/2G)^{3/2}. That is 4.3×10¹⁰ kg
    (2.7×10¹⁷ kg for 10⁴⁵ kg/m³).
  * Every astrophysical black hole is far above this. If ECSK is right, each of them already contains a bounce, and a new universe on its
    far side.
* **ECSK interior [Proposition-N 6k-1].**
  * The Oppenheimer–Snyder interior (closed dust FRW) with the spin term obeys ȧ² = A/a − S/a⁴ − 1 + a²/l_d².
  * It bounces at a_min ≈ (S/A)^{1/3} (numerics agree to 10⁻⁴) and re-expands beyond the horizon, into a new universe.
  * Without vacuum energy it returns to its original size and cycles.
  * A contracting region that carries vacuum energy with l_d < (3√3/2)A bounces and then inflates for ever. That criterion is analytic and
    checked numerically: l_d = 2.5 A inflates, 3.0 A cycles.
  * Consequence: a lasting cosmos from a stellar collapse needs growth from cycle to cycle (Popławski: particle production) or interior
    vacuum energy at the scale of the collapsing region.
* **Limiting curvature in D = 4 + n dimensions [Proposition 6k-2].**
  * The Schwarzschild–Tangherlini interior f = 1 − (r_s/r)^{n+1} meets de Sitter f = 1 − r²/l² at r₀ = (r_s^{n+1}l²)^{1/(n+3)}.
  * The junction exists iff r_s > l.
  * There r₀/l = (r_s/l)^{(n+1)/(n+3)} > 1: the de Sitter region starts larger than its Hubble length, so it inflates.

**What Y1 changes.** The machine no longer builds a false-vacuum bubble or forces a transition. It only has to make a black hole above the
theory's threshold:
* ECSK: M > 4.3×10¹⁰ kg.
* Limiting curvature: r_s > l.

In four dimensions that still meets the bridge invariant of Round 6i, unless a natural black hole is used (section 4).

## 2. Y2 — strong gravity at collider energies: what is allowed

* *LHC monojet.* M_D ≳ 10.7–11.2 TeV (n = 2) down to 5.2–5.9 TeV (n = 6–7) (CMS 2021, ATLAS 2021).
* *CMS black-hole search.* Semiclassical black holes and string balls are excluded below 8.4–11.4 TeV and 9.0–10.7 TeV (138 fb⁻¹; JHEP 08
  (2026) 098).
* *Astrophysics.* The bounds are strong for n = 2, 3 and weak for n ≥ 4.
* *What remains allowed.* n ≥ 4 with M_D ≈ 6–20 TeV.
* *Reach of a 100 TeV collider (scaling).* The LHC excludes black holes up to τ = M²/s ≈ 0.4–0.75. A 100 TeV collider with comparable
  luminosity reaches M_BH ≈ √0.5 × 100 TeV ≈ 71 TeV.

## 3. Y1 + Y2 — the conditional machine

**Proposition 6k-3 (conditional).** Assume P1 and P2:
* **P1:** n ≥ 4 large extra dimensions with M_D ≲ 15 TeV.
* **P2:** limiting-curvature interiors, with limiting length l = 1/M_D.

Then:
* A collider black hole holds an inflating de Sitter daughter once R_S > 1/M_D. With the Dimopoulos–Landsberg radius this gives
  M_BH > M_Y = 3.9 M_D (n = 4), 4.5 (n = 5), 4.7 (n = 6), 4.6 (n = 7).
* That is M_Y = 24–57 TeV for M_D = 6–12 TeV, and 68 TeV for M_D = 15 TeV (n = 5).
* All of these lie below the 71 TeV reach of a 100 TeV collider (τ at threshold 0.06–0.46).
* Geometric parton cross-sections at threshold are 8–34 pb.
* The interior reaches the junction within a proper time ~R_S/c. The hole's Hawking lifetime ~(1/M_D)(M/M_D)^{(n+3)/(n+1)} is about 7×
  longer, so the daughter forms before the hole evaporates.

**The premises, honestly.**

| Premise | Status |
|---|---|
| P1: gravity strong at 6–15 TeV (n ≥ 4) | Allowed. **Tested directly** by the 100 TeV collider: black-hole events, with high multiplicity and democratic decays. |
| P2: interiors become expanding daughters | True in ECSK and in limiting-curvature theories. Not true in effective LQG (white holes) or mimetic models (static cores). Decided by the correct theory of gravity. |
| P3: daughter forms before evaporation | Computed: lifetime/interior time ≈ 7 |
| P4: threshold within reach | Computed: 24–68 TeV < 71 TeV |
| P5: semiclassical control | **Not met at 100 TeV.** At threshold R_S ≈ l. Comfortable control (R_S ≥ 3l) needs M_BH ≈ 10³–10⁴ M_D (10–100 PeV collisions). |

**Partial empirical handles on P2 [Proposition-N 6k-4].**
* *The static-core class.* For a D-dimensional Hayward-type static de Sitter core, the temperature rises from zero at a remnant
  (μ_rem ≈ 10.8 l^{n+1}), peaks (T_max ≈ 0.19/l), then falls like GR (within 5% at large mass). A collider measuring the black-hole
  mass–temperature relation (Dimopoulos–Landsberg) and searching for stable remnants tests this class.
* *The expanding (FMM) class.* Its exterior evaporates as in GR. The daughter is invisible from outside, so P2 for this class rests on
  theory.

**Nature's version.** If P1 and P2 hold, cosmic-ray collisions with √s up to ~0.75 PeV already make such holes, and their daughters, in
the atmosphere every day. That is harmless: the daughters are hidden and the holes evaporate. The collider would be doing what nature
already does.

## 4. A better lead found this round: a natural bridge within reach

* **The hypothesis.** Scholtz & Unwin (PRL 125, 051103, 2020) show that the orbital anomalies attributed to Planet Nine, and a set of OGLE
  microlensing events, are consistent with a 5–10 Earth-mass primordial black hole at ~500 AU (horizon radius 4–9 cm).
* **Reaching it.** Witten (arXiv:2004.14192) shows that laser-launched probes at 0.001c could reach it in about a decade, for a launch
  system costing ~$0.5 B. Our flight time at 0.001c is 7.9 years.
* **If Y1 (ECSK or limiting curvature) holds:**
  * Such a hole (≫ 4.3×10¹⁰ kg) already holds a daughter, and humanity could reach a black hole and deliver matter into it within a
    human lifetime.
  * Whether delivering matter **creates a new** daughter, rather than feeding the existing one, is undecided. Barrabès–Frolov's
    multiple-world limiting curvature allows several disconnected de Sitter worlds; ECSK sends infalling matter through the existing
    bounce (P-R6k-c).
* **A coincidence worth recording.** With the corrected partner relation, a 5–10 Earth-mass hole is the Nariai partner of a 60–85 GeV
  vacuum. The Standard Model's electroweak vacuum energy is V(0)^{1/4} = 104 GeV, whose partner is 3.3 Earth masses. But the symmetric
  Higgs phase is not metastable at zero temperature, and restoring it thermally costs ~200× its vacuum energy (Round 6f). So this is a
  coincidence, not a route.

## 5. What would establish the "yes"

1. **P1 at a 100 TeV collider:** observe microscopic black holes (M_BH ≈ 25–70 TeV).
2. **P2 in theory:** a black-hole interior calculation in a theory that is otherwise consistent and recovers GR. ECSK already does this
   classically. Limiting curvature needs a complete action without pathologies. A quantum-gravity derivation (P-R6k-a) would settle the
   class.
3. **P2 empirically, partially:** the T(M) curve and remnant searches separate static cores from the rest.
4. **Planet Nine:** confirm or exclude the primordial-black-hole interpretation (probes, survey astrometry).
5. **P5:** 10–100 PeV collisions for semiclassically controlled daughters. That is 100–1000× the energy of the 100 TeV collider: a hard
   but definite engineering target, not 39 orders of magnitude.

## 6. Status

| Research objective | Status after Round 6k |
|---|---|
| Create the initial false-vacuum region | **Yes, by the universe** (Round 6g) |
| Establish a technologically accessible mechanism | **Conditional yes.** Under P1 (n ≥ 4 extra dimensions, M_D ≲ 15 TeV; allowed) and P2 (limiting-curvature or ECSK-type interiors; realized in definite theories, not yet shown to be nature's), a 100 TeV collider makes black holes above the daughter threshold (24–68 TeV < 71 TeV reach), whose daughters form before evaporation. P1 is tested by that collider; P2 is a property of the correct theory of gravity; semiclassical control needs 10–100 PeV. Without P1, the 4D bridge invariant stands (10^39.5 beyond technology) unless a natural black hole is used; a primordial black hole as Planet Nine would be reachable in ~10 years. |

## Literature used (Round 6k)

* *Einstein–Cartan–Sciama–Kibble*
  * Popławski, arXiv:1110.5019 (2011): spin–torsion averts singularities; each black hole's interior becomes a universe.
  * Popławski, arXiv:2008.02136, 2307.12190, 2509.11468 (IJMPA 40, 2544007, 2025): gravitational collapse with torsion into a universe in a
    black hole.
  * Popławski, arXiv:2007.11556: anisotropic (Kantowski–Sachs) version.
* *Limiting curvature*
  * Frolov, Markov & Mukhanov, PLB 216, 272 (1989); PRD 41, 383 (1990).
  * Barrabès & Frolov, PRD 53 (1996), hep-th/9511136; gr-qc/9607040.
  * Easson & Brandenberger, hep-th/0103019 (2001).
  * arXiv:2007.06664: an asymptotically de Sitter universe inside a Schwarzschild black hole.
* *Mimetic limiting curvature*
  * Chamseddine & Mukhanov, Eur. Phys. J. C (2017), "Nonsingular black hole"; "Black hole remnants", JHEP 10 (2019) 104, arXiv:1908.03498.
  * Ben Achour, Lamy, Liu & Noui, arXiv:1712.03876.
* *Loop quantum gravity*
  * Ashtekar, Olmedo & Singh, PRL 121, 241301 (2018); PRD 98, 126003 (2018).
* *Holography*
  * Cooper, Rozali, Swingle, Van Raamsdonk, Waddell & Wakeham, JHEP 07 (2019), arXiv:1810.10601: black-hole microstate cosmology.
  * Sasieta, Swingle & Vilar López, PRL 137, 041501 (2026).
* *Collider*
  * CMS, arXiv:2604.10732, JHEP 08 (2026) 098.
  * CMS arXiv:2107.13021; ATLAS arXiv:2102.10874.
  * Gamsızkan, arXiv:2404.10070 (black holes at 100 TeV).
  * Hossenfelder et al., hep-ph/0109085 (quasi-stable black holes and remnants).
* *Natural bridges*
  * Scholtz & Unwin, PRL 125, 051103 (2020), arXiv:1909.11090.
  * Witten, arXiv:2004.14192 (2020).
