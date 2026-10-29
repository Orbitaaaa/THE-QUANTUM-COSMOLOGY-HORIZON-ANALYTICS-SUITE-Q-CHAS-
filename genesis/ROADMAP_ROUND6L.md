# Round 6l — Unified quantum-gravity formation

The user asked for four things:
1. Unify P1 and P2: a consistent higher-dimensional model with singularity resolution, and whether its black-hole interiors actually
   inflate.
2. Collider production and evaporation across the allowed parameter space.
3. Multiple daughters from one black hole (P-R6k-c).
4. What experiments could distinguish.

They also asked for everything else needed to make the "yes" stronger.

Everything below is computed or checked against the literature:
* Code: `unified6l.py` (model U, interiors, daughters) and `collider6l.py` (parton distributions, rates, evaporation, discriminants).
* Script, tests, data, figure: `scripts/round6l.py`, `tests/test_round6l.py` (7 tests), `DATA/round6l.json`,
  `DATA/round6l_lumi*.json`, `FIGURES/round6l.png`.

**The answer in one paragraph.**
* **What is now established.** There is one consistent action that contains both doors. Call it model U: Einstein gravity in
  D = 4 + n dimensions plus the infinite tower of quasi-topological curvature terms of Bueno, Cano & Hennigar. It is pure gravity, with
  second-order spherical dynamics and a Birkhoff theorem. Its black holes are regular. Matter that falls in bounces at the limiting
  curvature and emerges through a white-hole region into a new universe, a daughter spacetime in the precise sense of a new asymptotic
  region.
* **At a collider.** These black holes form whenever the mass exceeds the horizon threshold M_ext = κ_n c_l^{n+1} M_D. A 100 TeV collider
  makes tens to ~10⁹ of them (conservative), provided M_ext lies in the window ≈ 9–50 TeV (LHC below, collider reach above). That window is a
  definite band of the theory's one free coefficient, the limiting length l = c_l/M_D. The interior forms 20–7000 times faster than the
  hole evaporates, and evaporation stops at a remnant.
* **What did not survive.** In model U the daughter does **not** inflate:
  * vacuum collapse leaves a static de Sitter core;
  * matter cores accelerate only for N ≈ ((D−3)/(D−1)) ln(R/l) e-folds, then recollapse;
  * lasting inflation needs a vacuum energy as dense as the hole itself.

  A second infall cannot simply add a second daughter. It joins the existing core, and after ~7–30 l of delay it drives the meeting region
  to the curvature bound. No experiment outside the horizon can see the daughter; experiments can test the class of interior that implies
  it.

## 1. Priority 1 — one model with both doors

### 1.1 Model U

* **The theory.** D = 4 + n dimensional gravity with fundamental scale M_D (P1; Dimopoulos–Landsberg convention G_D = 1/M_D^{n+2}). We add
  the quasi-topological tower of Bueno, Cano & Hennigar (PLB 861 (2025), arXiv:2403.04827).
* **Its known properties (Literature).**
  * In spherical symmetry the equations are second order.
  * A Birkhoff theorem holds, and Israel-type junction conditions exist (Bueno, Cano, Hennigar & Murcia, PRD 111, 104009 (2025); PRL 134,
    181401 (2025)).
  * The static solutions are unique and regular: Hayward-type, Bardeen-type and Dymnikova-type, as vacuum solutions with no matter fields.
* **How everything is fixed.** Everything is determined by one function h through the generalised Misner–Sharp relation:
  * **Static holes:** h(ψ) = μ/r^{D−1}, with ψ = (1 − (∇r)²)/r².
  * **Closed FRW regions:** h(H² + 1/a²) ∝ ρ. This is the form used by Bueno, Cano, Hennigar, Murcia & Vicente-Cano, arXiv:2505.09680.
* **The towers we use.** We implement two towers with a pole at ψ = 1/α, and limiting length l = √α = c_l/M_D:
  * Hayward: h = ψ/(1 − αψ).
  * Dymnikova: h = −ln(1 − αψ)/α.
* **ADD consistency [computed].** The extra-dimension radius is R = (M_P²/M_D^{n+2})^{1/n}/2π, from 4.5×10⁻¹⁶ m (n = 6, M_D = 8 TeV) to
  1.5×10⁻¹³ m (n = 4). That is 10³–10⁶ times larger than 1/M_D, so collider holes are D-dimensional and asymptotically flat.
* **Experimental consistency.** The tower acts only at curvatures ~1/l² ~ M_D². It changes nothing measured so far except black-hole
  formation itself.

### 1.2 Black holes of model U [Proposition, computed]

* **Structure.**
  * Near the centre, f(r) = 1 − r²ψ(r) → 1 − r²/α: a de Sitter core.
  * There are two horizons above the extremal mass, and none below.
  * Hayward tower: μ_ext = α^{(D−3)/2}((D−1)/(D−3))^{(D−1)/2}(D−3)/2, i.e. r_s,ext = 2.60 l (D = 4), 1.41 l (D = 10).
  * Dymnikova tower: r_s,ext = 1.96 l (D = 4), 1.20 l (D = 10).
* **Temperature.** T = ((D−3)r₊² − (D−1)α)/(4πr₊³). It vanishes at extremality, so evaporation ends at a remnant.
* **Vacuum collapse (colliding shock waves, thin shells).** Birkhoff fixes the interior to this static geometry: a trapped shell, an inner
  horizon, and a **static** de Sitter core of radius r₋ ≈ l (e.g. r₋ = 1.002 l at 50 μ_ext). The core does not inflate.

### 1.3 Do matter interiors inflate? [Proposition 6l-1]

* **The equation.** A closed dust interior obeys ȧ² = a²ψ(C/a^{D−1}) − 1.
* **Turning points.** They are exactly the horizon radii of the static hole of mass C (the same polynomial).
* **The bounce.** The interior bounces at a_b → l. This is universal: any tower with a pole at 1/α does it, and the Hayward and Dymnikova
  towers agree to 10⁻³.
* **Transient acceleration.** There follows an accelerating phase of N_acc = ln(a_infl/a_b) e-folds, with a_infl^{D−1} = 2αC/(D−3)
  (exact for Hayward):

  N_acc → ((D−3)/(D−1)) ln(a_max/l) + ln(2/(D−3))/(D−1).

  Numerically:
  * 0.88 e-folds for C = 10⁴ (D = 10);
  * **29.6 e-folds** inside a four-dimensional stellar hole (a_max/l = 2×10³⁸);
  * 36 e-folds for 10⁹ solar masses.
* **What follows.** The interior decelerates, re-expands to its original size a_max, and cycles. Even inside a stellar hole, about 30
  e-folds of acceleration end in a recollapse.
* **Verdict for P2 "inflating interiors" in the unified model:** no lasting inflation from matter alone.

### 1.4 What would make it inflate [Proposition 6l-2]

* **The criterion.** With an interior vacuum term Λ, the bounced region inflates for ever iff Λ > Λ_c(C).
  * **Low-curvature limit:** Λ_c = [2((D−3)/(D−1))^{(D−1)/2}/((D−3)C)]^{2/(D−3)}. In D = 4 this is l_Λ < (3√3/2)C, the same criterion as
    ECSK (Round 6k). It is also the Hayward extremality function with α → 1/Λ.
  * **Numerics:** they agree to 10⁻⁴ in D = 4. In D = 10 they lie 2–54 % above the limit (largest for the smallest holes).
* **What model U adds.** Because model U has no singularities, the Farhi–Guth past-singularity obstruction, which forced every earlier
  route through quantum tunnelling, is absent: a trapped false-vacuum region need not have come from a singularity.
* **What model U does not remove.** It does not remove the vacuum-energy requirement: the false vacuum must be as dense as the hole, which
  is the Nariai-partner relation of Rounds 6c–6j.

### 1.5 Einstein–Cartan, corrected [Propositions 6l-3 and 6l-4]

* **6l-3 (kinematic, any theory).** At a bounce of a radial timelike worldline, (∇R)² = (n·∇R)² ≥ 0: the bounce sphere is untrapped. A
  bounce inside a black hole therefore requires an untrapped core bounded by an inner apparent horizon. This is true in every
  singularity-resolving theory.
* **6l-4.** In ECSK the homogeneous ball's Misner–Sharp mass falls from (sin³χ/2)A to (sin³χ/2)a_b at the bounce (99 % lost for
  S = 10⁻⁶), because the spin pressure does work. So the ball cannot sit in a vacuum exterior of fixed mass.
  * The Round 6k picture of a homogeneous ball re-expanding beyond the horizon needs this correction (F61).
  * The realistic ECSK bounce is inhomogeneous: an untrapped core inside still-falling layers. That is the same causal structure as model U.

## 2. Priority 2 — production, rates, formation, evaporation

* **Parton distributions.** We use NNPDF2.3 QCD+QED LO (α_s = 0.130, the Pythia 8.3 default), read directly from the LHAPDF6 grid. The
  momentum sum is 0.998 and the valence numbers are 2.003 and 1.001. NNPDF3.1 LO serves as a cross-check: luminosities agree to ~1 % for
  τ ≤ 0.1 and to ~30 % at τ = 0.3.
* **Cross-section.** σ̂ = F_n π R_S², with the Yoshino–Nambu form factors F_n.
* **Trapped mass.**
  * **Conservative:** the head-on area bound y_n = 0.71, 0.67, 0.64, 0.62, 0.61, 0.60, 0.59, 0.58 for n = 0–7 (computed; n = 0 reproduces
    1/√2).
  * **Optimistic:** M = √ŝ.
* **Scale choice.** The factorisation scale is Q = M or 1/R_S; the spread is ≤ 20 %.

**Proposition-N 6l-6 (the window).** The model enters the rate only through M_ext.

| | requirement on M_ext | equivalent limiting length c_l = l M_D |
|---|---|---|
| LHC (13 TeV, 138 fb⁻¹, no black holes seen; simplified) | M_ext ≳ 8.7–9.7 TeV | Hayward ≳ 0.45–0.66, Dymnikova ≳ 0.54–0.76 |
| 100 TeV, 30 ab⁻¹, ≥ 10 events (conservative) | M_ext ≲ 46–51 TeV (optimistic: 79–84) | Hayward ≲ 0.63–0.80, Dymnikova ≲ 0.75–0.95 |
| 85 TeV (FCC-hh, 20 ab⁻¹) | M_ext ≲ 39–43 TeV (optimistic: 67–71) | Hayward ≲ 0.61–0.79, Dymnikova ≲ 0.73–0.92 |

(Ranges over n = 4–7 and M_D = 5–15 TeV.)

Inside the window the rates are large. Examples at 100 TeV, 30 ab⁻¹, conservative trapping:

| n | M_D | tower | c_l | M_ext | events |
|---|---|---|---|---|---|
| 6 | 8 TeV | Hayward | 0.70 | 34 TeV | 1×10⁵ |
| 6 | 8 TeV | Hayward | 0.60 | 11.5 TeV | 4×10⁸ |
| 6 | 8 TeV | Dymnikova | 0.70 | 11.4 TeV | 4×10⁸ |
| 4 | 8 TeV | Hayward | 0.70 | 43 TeV | 1.6×10³ |

With optimistic trapping the counts are 10–10³ times larger.

**Formation against evaporation [computed].**
* **Interior.** Matter crosses from the horizon to its bounce in 2.7 l of proper time, 0.7–2×10⁻²⁸ s.
* **Evaporation.** Brane-dominated Stefan–Boltzmann emission (g* = 106.75; Emparan–Horowitz–Myers capture radius) takes the hole to its
  remnant in 2×10⁻²⁷–10⁻²⁴ s. The ratio is 20–7×10³ over the scan: the daughter region forms long before the hole has evaporated. For
  comparison, a GR hole of the same mass is hotter and evaporates in ~10⁻²⁸ s, comparable to the interior time (the ratio of ~7 found in
  Round 6k).
* **Model U against GR, at equal mass.**
  * The remnant keeps the fraction M_ext/M of the mass (7–36 % for 120 TeV holes), seen as missing energy.
  * The hole emits 1.5–2.5× more quanta, which are softer: peak temperature 1.9–3× lower.
  * T_U stays more than 10 % below GR until r₊ ≳ 2.5–3.6 l, i.e. over essentially the whole collider range.

## 3. Priority 3 — multiple daughters from one black hole (P-R6k-c)

* **Second infall in model U [Proposition 6l-5, computed].**
  * In ingoing Eddington–Finkelstein coordinates, the bounced surface of the first collapse reaches the inner horizon only as v → +∞
    (r₋ − r ∝ e^{−κ₋v}).
  * So any later infall at finite advanced time meets it **inside the same core**:
    * after a delay of 2 l, it meets 0.16 l inside the inner horizon;
    * after 50 l, it meets 10⁻²¹ l inside;
    * after 500 l, 10⁻²¹² l inside.
  * The later matter enters and modifies the existing daughter. A theory with a Birkhoff theorem and a single regular centre leaves no
    room for a separate one. By Geroch's theorem, a disconnected daughter would need topology change, which a regular, causal spacetime
    forbids.
* **The meeting saturates the bound [Proposition-N 6l-6b].**
  * At the meeting, the Dray–'t Hooft–Redmount relation f_P f_F = f_L f_R holds. It is kinematic: it needs only continuity of r and a
    C⁰ metric, so it holds in model U.
  * It forces the region after the crossing to have f_F → −∞ (mass inflation; in GR the mass grows by 10¹³ after 50 l).
  * Model U bounds f_F ≥ 1 − r²/α, because ψ < 1/α.
  * The mass needed diverges, and the crossing region saturates at the limiting curvature (locally de Sitter, trapped), at a critical
    delay of 7–30 l (log-dependent on the mass ratio; D = 10: 9.3 l for 10⁻¹, 13.9 l for 10⁻³, 29.7 l for 10⁻¹²).
  * Beyond that delay no thin-shell solution exists.
  * This is exactly the regime in which Barrabès & Frolov constructed de Sitter worlds bounded by lightlike surfaces (PRD 53, 3215
    (1996); gr-qc/9607040).
  * **So:** additional matter cannot pass quietly into the existing daughter. It either seals or modifies it, or seeds a new de Sitter
    world. Which one requires thick-shell dynamics at the curvature bound (P-R6l-a).
* **Limiting-curvature (FMM) class [computed].**
  * The spacelike junction at r₀ = (r_s l²)^{1/3} has a continuous induced metric (match 10⁻¹¹), so Δt_dS = Δv.
  * On the contracting branch (continuous extrinsic curvature: the FMM choice), the daughter is one closed de Sitter world, and two
    entries have meeting futures.
  * On the expanding branch (which requires a spacelike shell), entries separated by Δv > l ln((r₀+l)/(r₀−l)) ≈ 2l(l/r_s)^{1/3} land in
    causally disjoint inflating patches: operationally new daughters.
* **Evaporation endpoint [Proposition-N 6l-7, closure].**
  * In spherical symmetry with a regular centre, a compact trapped region (complete evaporation, as in Hayward 2006 or Frolov 2014) sends
    everything back to our future null infinity: no daughter persists.
  * A remnant (T → 0 at M_ext, as in model U) means the trapped region never closes, so the continuation persists.

## 4. Priority 4 — what experiments can and cannot distinguish

| Observable | What it tests | Sensitive to the interior? |
|---|---|---|
| missing-energy monojet / monophoton (KK graviton emission) | P1: n and M_D | no |
| high-mass dilepton / diphoton tail (virtual gravitons) | P1: cut-off M_S | no |
| onset of high-multiplicity, high-ΣpT, democratic events | P1 and, through M_th/M_D, the class: GR-like few M_D, model U κ_n c_l^{n+1} M_D (13–56 c_l^{n+1}) | yes |
| T(M) from spectra against reconstructed mass | GR T ∝ M^{−1/(n+1)} against a regular core (T → 0 at M_ext; >10 % below GR up to r₊ ≈ 3 l) | yes |
| stable remnants (missing energy M_ext/M; heavy stable charged tracks if charged) | regular / mimetic cores against complete evaporation | yes |
| multiplicity at fixed mass | model U 1.5–2.5× GR | yes |
| final burst at M ~ M_D | not semiclassical | unreliable |
| the region beyond the core: new universe, inflating de Sitter, crunch, nothing | **nothing** | — |

**Exterior-equivalence lemma [Proposition-N 6l-8].**
* **Statement.** Two spacetimes formed by the same collapse that coincide on J⁻(I⁺) give identical classical data and an identical
  semiclassical out-state at I⁺ (Hawking quanta and all their exterior correlations). The data at I⁺ are propagated only from J⁻(I⁺),
  and the daughter lies outside it.
* **The logic that remains possible.**
  1. Measure M_D from monojets.
  2. Measure the black-hole onset M_th and the remnant with T → 0.
  3. These identify the regular-core class.
  4. In that class, the closure lemma makes the remnant a necessary condition for a persistent daughter.

  Experiments thus test the class that implies the daughter, never the daughter itself.

**Safety.** If stable remnants exist, those made by cosmic rays would long since have been captured in neutron stars and white dwarfs.
Giddings & Mangano (PRD 78, 035009 (2008)) show that for these scenarios the survival of observed neutron stars and white dwarfs excludes
any danger.

## 5. Status after Round 6l

| Question | Answer |
|---|---|
| One consistent model with P1 and P2 | **Yes**: model U, pure gravity, D = 4 + n (Literature for the theory, computed here for collider holes) |
| Do its holes contain a daughter? | **Yes, a new asymptotic region beyond a white hole** (bounce at the limiting curvature) |
| Do the interiors inflate? | **No lasting inflation**: static core for vacuum collapse; N ≈ ((D−3)/(D−1)) ln(R/l) transient e-folds for matter (29.6 for a stellar hole); lasting inflation needs a hole-density vacuum |
| Made at a 100 TeV collider? | **Yes, iff M_ext ≈ 9–50 TeV**, i.e. c_l in a window of width ~0.15–0.25 (e.g. 0.58–0.74 Hayward, 0.68–0.86 Dymnikova at n = 6, M_D = 8 TeV); up to ~10⁹ events (conservative) |
| Daughter forms before evaporation? | **Yes**, 20–7×10³ times faster; evaporation stops at a remnant, so the continuation persists |
| More than one daughter per hole? | **Not by quiet addition**: later matter joins the core; after 7–30 l of delay it saturates the curvature bound (the Barrabès–Frolov regime: open) |
| Can experiments see the daughter? | **No**; they can establish the class (M_D, onset, T(M), remnant) |

## 6. Literature used (Round 6l)

* *Model U*
  * Bueno, Cano & Hennigar, PLB 861 (2025), arXiv:2403.04827 (regular black holes from pure gravity, D ≥ 5).
  * Bueno, Cano, Hennigar & Murcia, PRD 111, 104009 (2025), arXiv:2412.02740 (thin-shell collapse, Birkhoff theorem, junction conditions,
    curvature bound); PRL 134, 181401 (2025), arXiv:2412.02742.
  * Bueno, Cano, Hennigar, Murcia & Vicente-Cano, arXiv:2505.09680 (Oppenheimer–Snyder collapse; the Friedmann equation h(Φ) ∝ ρ).
* *Related cosmologies*
  * Fernandes, PRD 112, 084028 (2025), arXiv:2504.07692 (singularity resolution and inflation from regularised Lovelock towers, 4D).
  * Ling & Yu, arXiv:2509.00137 (big bounce in quasi-topological gravity).
* *Multiple worlds and limiting curvature*
  * Barrabès & Frolov, PRD 53, 3215 (1996), hep-th/9511136; gr-qc/9607040.
  * Frolov, Markov & Mukhanov (1989/90).
* *Inner horizons*
  * Dray & 't Hooft (1985); Barrabès & Israel (1991); Poisson & Israel, PRD 41, 1796 (1990) (mass inflation).
  * Carballo-Rubio et al., arXiv:2101.05006; inner-extremal regular black holes, JHEP (2022).
* *Topology and evaporation*
  * Geroch, J. Math. Phys. 8, 782 (1967).
  * Hayward, PRL 96, 031103 (2006).
  * Frolov, JHEP 05 (2014) 049.
* *Collider*
  * Eardley & Giddings, PRD 66, 044011 (2002).
  * Yoshino & Nambu, PRD 67, 024009 (2003) (form factors 0.647–1.883).
  * Yoshino & Rychkov, PRD 71, 104028 (2005).
  * Emparan, Horowitz & Myers, PRL 85, 499 (2000).
  * NNPDF2.3, Nucl. Phys. B 867, 244 (2013); Pythia 8.3, SciPost Phys. Codebases 8 (2022).
  * FCC Feasibility Study Report (2025): FCC-hh ≈ 85 TeV.
  * Giddings & Mangano, PRD 78, 035009 (2008).
