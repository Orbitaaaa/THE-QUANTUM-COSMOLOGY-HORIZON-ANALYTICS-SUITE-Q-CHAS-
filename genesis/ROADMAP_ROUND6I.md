# Round 6i — seven routes that avoid manufacturing a false-vacuum ball, and the obstruction they share

The user proposed seven candidate routes, each with published precedent, that might replace "manufacture an enormous false-vacuum ball". The
user framed every route by three questions:

1. **Create the state.** Can something, natural or made, contain the high-energy vacuum configuration?
2. **Trigger.** Can a small intervention change its evolution?
3. **Complete the daughter.** Does the geometry actually inflate, within the applicable constraints?

The routes, in the user's priority order:

| # | Route |
|---|---|
| 1 | Curvature-activated field (scalar–Gauss–Bonnet) |
| 2 | Stable monopole to daughter (Sakai et al. 2006) |
| 3 | Breathing bubble (Guendelman–Sakai 2008) |
| 4 | Gauged / vacuum-dependent ghost condensate |
| 5 | Resonant tunnelling |
| 6 | Classical transitions from bubble collisions |
| 7 | Four-form membranes |

The user asked us to keep looking for further methods if these failed, and not to downplay results.

Every route was computed, not argued. Code: `gauss_bonnet.py`, `stable_shells.py`, `routes6i.py`, `scripts/round6i.py`, `tests/test_round6i.py`
(11 tests). Data: `DATA/round6i.json`. Figure: `FIGURES/round6i.png`. Labels as before (Theorem, Proposition, Proposition-N, Literature,
Assessment, Conjecture).

**Answer in one paragraph.**
* **None of the seven routes yields a mechanism humanity can operate.** Each fails at an identified, computed point, recorded in the ledger
  below.
* **One theorem explains why.** Every hidden daughter contains a sphere with 2Gm/(c²r) ≥ 1, a "bridge" (Theorem 6i-11). Making one
  classically from ordinary initial data is a black-hole formation event. That event obeys the scale-free invariant
  E√ρ̄ = √(3/32π) c⁶/G^{3/2} ≈ 2.3×10⁶⁵ J^{3/2} m^{−3/2}, whatever the vacuum scale, the route or the NEC-violating ingredient. The best
  technology is 10^39.5 short of it.
* **The quantum routes obey the project's "no free lunch" bound.** Creating a daughter is not exponentially easier than the decay of its
  own false vacuum (R10, R37); this round checks it for monopoles and membranes.
* **One door is not closed by any established result.** A barely metastable false vacuum with a light wall gives a semiclassically
  unsuppressed channel at TeV energies. But its bridge is about 10⁻¹⁶ Planck lengths across, which is far outside every regime where
  gravity is understood (Q*, section 8).
* **The largest real improvement of the round is in background mass.** A condensate whose scale depends on the vacuum lowers the required
  background from 2.5×10⁹ M_☉ to 9×10⁸ kg (conditional). The invariant is unchanged.

## 0. Ledger (state / trigger / completion)

| Route | Create the state | Trigger | Complete the daughter | Blocking point |
|---|---|---|---|---|
| 1 Gauss–Bonnet activation | Needs horizon-class curvature (lab short by 19 orders) | Activation raises F, so it pushes toward collapse (Prop. 6i-1) | Never reached: black hole or kinetic degeneracy in 152 nonlinear runs | Activation and daughter-tipping have opposite signs |
| 2 Monopole core | Natural if relic monopoles exist (none seen); stable below η_c ≈ 0.25–0.40 m_P | The channel at the monopole's own mass is disconnected from us (M ≪ Q); the connected channel needs M > Q ≈ 7 m_P, then B ≈ 3×10³–3×10⁴ | After tunnelling | Planck-density assembly, and e^{−B} |
| 3 Breathing bubble | Hidden case needs negative-energy dust (violates the NEC); every well crosses zero net surface energy at light speed | Tunnelling, B = B̂ (l_f/l_P)², B̂ ≈ 0.5–0.9 | After tunnelling | Exotic dust; consistency of the composite shell; e^{−B} |
| 4 Gauged ghost | Universal gauged condensate: window empty by 5 orders (Prop. 6i-8). Vacuum-dependent M(χ): background ≥ 9×10⁸ kg (Prop. 6i-9) | Condensate tipping (Rounds 6c–6f) | Simulated (Rounds 6e–6f) | Assembling the 10⁹ kg background within r_N/c |
| 5 Resonant tunnelling | Two-well shell potentials exist | At best turns B₁ + B₂ into max(B₂, B₁ − B₂) | As the underlying tunnelling | The barrier next to the inflating branch remains |
| 6 Classical transitions | Pockets of higher vacuum are kinematically allowed | Collision | Impossible in a normal parent (Prop. 6i-10); possible in an inflating parent | Needs an anti-trapped parent: the cosmological route |
| 7 Four-form membrane | Flux pocket is a Farhi–Guth wall (NEC saturated) | Quantum creation, B ≈ B_CdL | After tunnelling | Theorem 6h-1 classically; e^{−B} |

---

## 1. Route 1 — the curvature-activated field

Action (Doneva–Yazadjiev normalisation, G = c = 1), with the Gauss–Bonnet invariant
GB = R² − 4R_abR^ab + R_abcdR^abcd:

S = (1/16π)∫√−g [R − 2(∇φ)² − 2V(φ) + λ² f(φ) GB]

**The Gauss–Bonnet stress in spherical symmetry [Theorem, symbolic].**
* We computed Γ_ab for ds² = −e^{2s}du dv + r²dΩ² with F = λ²f(φ) (`gb_stress_double_null`).
* We verified the off-shell identity ∇^aΓ_ab = −(GB/2)∇_bF exactly with sympy.
* The null component is
  −Γ_uu = 4(2m/r³)∇_u∇_uF + (16e^{−2s}/r²)(∇_u∇_ur)(F_u r_v + F_v r_u),
  with m the Misner–Sharp mass.

**Proposition 6i-1 (activation pushes toward collapse).**
* *The flux.* On a vacuum-energy throat (r ≈ r_N, 2m/r³ = 1/r_N², ∇_u∇_ur = O(φ²)), the Round 6c tipping flux is
  N = 4πr²T_uu = r_N²φ_u² + 2F_uu (quadratic order, affine u).
* *A total derivative.* The Gauss–Bonnet part is a total second derivative, so its weighted deficit is
  W_GB = −2∫F_uu e^{−u/τ}du = −(2/τ²)∫(F − F₀)e^{−u/τ}du, plus boundary terms of the same sign when F grows.
* *Why activation raises F.* In a vacuum-energy throat GB = C² + 8Λ²/3 > 0. When the coupling is what destabilises φ, i.e. the mass
  m² = −(λ²/4)f″GB is negative, the field rolls so as to increase F.
* *Consequence.* So W_GB < 0: curvature activation can only push a throat toward the black hole.
* *The other arrangement.* The opposite case (f″ < 0, F falling during the roll) needs a potential that is tachyonic at low curvature. Such a
  field cannot be stored in the false-vacuum ball or in our vacuum.

**Proposition-N 6i-2 (the nonlinear homogeneous throat).**
* *The model.* Kantowski–Sachs metric ds² = −dt² + a²dx² + b²dΩ² with the coupling. We verified
  √−g GB = d/dt[8ȧ(1 + ḃ²)] symbolically, and derived the equations from the reduced Lagrangian (`ks_equations`). All runs start from exact
  Nariai data.
* *Scalarisation runs.* Coupling f = (1 − e^{−6φ²})/12 with λ² = 0.02–2 r_N². These end in a black hole (λ² = 0.02, 0.05, 0.1, 0.5), or in a
  degenerate kinetic matrix (λ² = 0.2, 0.3, 1, 2). The degeneracy is the homogeneous analogue of loss of hyperbolicity; it sets in when
  λ²|f′| ≈ b²/2 or through the HḞ terms.
* *Reversed-sign runs.* f″ < 0 with a hilltop potential that has a positive-energy minimum: all 144 runs end in a black hole (66) or the
  degeneracy (78). The control without Gauss–Bonnet collapses.
* *Discarded runs.* Runs with an unbounded potential that reached b = 3 r_N are not daughters: the potential turns negative before b
  doubles, and the solution then runs away. We discarded them.

**Assessment (where activation happens).**
* Activation needs λ²|GB| ~ 1/L² over the region size L.
* The bound is λ ≲ 0.28 km (GW230529). With it:
  * Earth's surface and dense laboratory masses (GB ≲ 10⁻⁴³ m⁻⁴) fall short by 19.0–19.5 orders of magnitude.
  * Neutron-star surfaces and stellar black-hole horizons fall short by 1.5 orders.

**Proposition 6i-3 (a curvature-made false vacuum).**
* *The bound.* V_eff = V − λ²fGB/2, so curvature can raise the vacuum by at most λ²|Δf|GB/2.
* *Consequence.* With GB = 12/l_curv⁴ at a horizon of radius l_curv, the lifted vacuum's Hubble length is l_f ≥ l_curv when |Δf| ≤ 1/2 and
  λ ≤ l_curv (the effective theory is valid).
* *Meaning.* A curvature-made false vacuum can fill a daughter-sized region only around an object that is already horizon-class at that
  scale.

**Decisive test (the user's question).**
* *Does it give NEC-violating stress?* Yes: the Gauss–Bonnet stress violates the NEC. Wormholes are known in dilatonic Einstein–Gauss–Bonnet
  theory, though they are unstable (Cuyubamba, Konoplya & Zhidenko 2018).
* *Does it give the right kind?* No. In the throat, the violation that curvature activation produces has the sign that makes a black hole.
* *Is it complete?* The nonlinear evolution never completes a daughter within the effective theory.

## 2. Route 2 — the monopole core

**Thin-shell model (after Sakai, Nakao, Ishihara & Kobayashi 2006).**
* Interior: de Sitter with density λη⁴/4.
* Wall: tension (√2/3)√λ η³.
* Exterior: magnetic Reissner–Nordström with Q² = 4π/e².

**Proposition-N 6i-5.** Static stable cores exist below a critical scale η_c (literature: O(0.3) m_P):

| λ | e | η_c / m_P |
|---|---|---|
| 1 | 0.5 | 0.25–0.26 |
| 1 | 1 | 0.32–0.33 |
| 0.1 | 0.5 | 0.39–0.40 |

GUT monopoles (η ≈ 10⁻³ m_P) are stable.

**The channels.**
* *At the monopole's own mass.*
  * Thin-shell exponent B ≈ 650/λ, independent of η for η ≤ 0.02 and falling to zero at η_c. That is 1.54× the flat thin-wall CdL action of
    the same wall (421/λ).
  * But here M ≈ 0.01–0.04 Q. The exterior has no horizon and no bridge, so the β_out < 0 branch is disconnected from our asymptotic region.
  * The far magnetic field (Gauss's law) forbids that outcome as a process in our universe.
* *Connected channel.* It needs M > Q ≈ 7 m_P, i.e. the core inside a magnetically charged black hole. Then
  B = 3.1×10³ (η = 0.01) to 6.5×10³ (η = 0.003) at M = 1.05 Q, rising with M to 2.8×10⁴.
* *Added energy.* With added mass up to 10× the static value, including past the horizon threshold, the barrier to the inflating branch
  persists. The charged Cauchy horizon gives no classical path.
* *NEC admixture.* Negative dust on the shell lowers B from 650 to 513, a reduction of 20% at most, before the stable core disappears. That
  happens where the net surface energy reaches zero. It opens no classical channel.

**Verdict.**
* *The state.* Nature might supply the state (relic monopoles), but none has been observed.
* *What remains.* The connected trigger needs Planck-density assembly (M > Q within ~10 l_P) and then e^{−3000} or less.

## 3. Route 3 — the breathing bubble

**The model.** Composite shell with surface energy m_s(R) = k_wR² + m_d (domain wall plus dust), de Sitter inside, Schwarzschild outside.

**Proposition-N 6i-4.**
* *Signs.* Bounded wells need k_w m_d < 0.
* *Which outcome.*
  * The hidden-daughter case (β_out < 0 on the large-R branch) needs k_w > 0, hence m_d < 0. That is negative-energy dust, which violates
    the NEC (σ + p = σ_d < 0).
  * k_w < 0 gives the bubble that eats the parent.
  * V is even in m_s and β_out is odd, which reproduces the Guendelman–Sakai dichotomy.
* *Every well crosses zero energy.* In 6000 random shells, all 713 wells contain the radius where m_s = 0. The shell crosses it at the speed
  of light, and that is where it changes side of the bridge (section 8).
* *Barrier.* B = B̂(l_f/l_P)², with B̂ ≈ 0.51–0.61 along the family (k_w, m_d) = (1, −0.1) and 0.84–0.86 along (2, −0.05). B̂ drops toward
  10⁻⁵ only at merger points of the well with neighbouring regions.
* *Mass.* The total mass can be about 10⁻³ l_f c²/G, two orders below the Nariai mass. The vacuum, wall and dust components are each of that
  order.

## 4. Route 4 — gauged and vacuum-dependent ghost condensation

**The NEC–charge identity.** In a ghost condensate ρ + p = 2XP′(X)M⁴ and the shift-charge density is J⁰ = 2P′(X)φ̇. So every NEC-violating
region (P′ < 0) is charged, with |ρ + p| = |J⁰|φ̇.

**Screening.** Gauged (Cheng, Luty, Mukohyama & Thaler 2006), the charge sources an electric field. For a uniform violating ball,
Coulomb/|ρ+p| = (P″|δ|/9)g²M²R². The violation is screened beyond R ~ 3/(gM√(P″|δ|)).

**Proposition 6i-8 (the gauged universal window is empty).**
* *The conditions.* Tipping a throat needs M³r_N ≥ 1.9×10⁻⁴ (Round 6f window) and an unscreened violation, r_N ≲ c/(gM). CLMT need
  M ≲ min(10¹² GeV, g²·10¹⁵ GeV).
* *Together.* M ≥ (K_lo²/(c²M_*))^{1/3}, which is 9.6×10¹⁶ GeV (c = 30) or 4.5×10¹⁷ GeV (c = 3).
* *Verdict.* The window is empty by 5.0–5.6 orders of magnitude. Gauging does not rescue the universal condensate.

**Proposition 6i-9 (vacuum-dependent scale, conditional).**
* *The model.* Let the condensate scale depend on the vacuum, M(χ): large inside the false vacuum, M_out ≤ 100 GeV in ours.
* *Why not an overall factor.* An overall factor F(χ) → 0 instead makes the strong-coupling scale M F² vanish in our vacuum: it is
  10⁻⁷⁵ GeV for M = 10¹³ GeV.
* *Wall thickness.* The effective theory stays valid across the wall only if the wall is thicker than 1/M_out. So r_N ≥ 1/M_out
  = 2.0×10⁻¹⁸ m.
* *The resulting background.*
  * M_N ≥ 8.9×10⁸ kg, with a 1.6×10¹⁰ GeV false vacuum and an interior condensate scale 1.4×10¹² GeV.
  * That is 30 orders of magnitude below the universal-condensate floor.
* *What remains.* Assembly: 5×10³³ quanta of 100 GeV within 6.6×10⁻²⁷ s inside 2×10⁻¹⁸ m.

## 5. Route 5 — resonant tunnelling

**Proposition 6i-6.**
* *Existence.* Two-well thin-shell potentials (well, barrier, well, barrier, hidden inflating branch, with m_s > 0 throughout) exist once
  charge and gas components are added: 13 clean cases in 20000 random shells.
* *Resonant masses.* Bohr–Sommerfeld levels of the two wells, with phases scaling as (l_f/l_P)², give discrete resonant masses
  (`resonant_masses`).
* *The double-barrier rate (WKB).*
  * Off resonance the exponent is B₁ + B₂.
  * On resonance it is max(B₂, B₁ − B₂) ≥ B₂ (coherent sharing versus incoherent hopping).
* *Verdict.* Resonance can remove the inner barrier at most, never the barrier next to the inflating branch. In the cases found the inner
  barrier is already tiny (B₁ ≲ 10⁻³ B̂), so the gain is negligible. This agrees with Copeland–Padilla–Saffin, who find no resonance for
  the homogeneous false vacuum.

## 6. Route 6 — classical transitions

**Kinematics.**
* The Dray–'t Hooft–Redmount relation f_A f_D = f_B f_C at a spherical collision allows a pocket of any Λ_C (`dtr_region_c`).
* Johnson & Yang (2010) show that in planar collisions with gravity the new vacuum cannot exceed the parent's.

**Proposition 6i-10 (completion is impossible in a normal parent).**
* *Sign constancy.* Where f_out > 0, β_out² = f_out + Ṙ² > 0, so β_out keeps its sign.
* *Consequence.* A pocket wall launched from a collision in a normal region has β_out > 0. It can never reach the hidden inflating branch,
  where β_out < 0. With positive tension it cannot grow into the parent either, so it collapses.
* *Check.* Every bound-branch turning point in `pocket_scan` has β_out > 0 and f_out > 0; every inflating branch has β_out < 0.
* *The exception.* In an anti-trapped (inflating) parent, collisions can make vacua up to the parent's. That is the cosmological route of
  Round 6g, which only nature runs.

## 7. Route 7 — four-form flux and membranes

**Proposition 6i-7.**
* *Same as a Farhi–Guth wall.* A Brown–Teitelboim membrane has charge q; the interior flux E + q raises the vacuum by q(E + q/2). It is a
  constant-tension wall with σ + p = 0 (NEC saturated), so its branches are those of Farhi–Guth.
* *The hidden branch.* The hidden inflating branch exists for every M below M_cr. Classically it starts at the white-hole singularity
  (Theorem 6h-1; Theorem 6i-11).
* *Quantum creation and "no free lunch".* B(M → 0) = 1.02 × the flat thin-wall CdL action for the same membrane, so the bound R10/R37 holds.
  B is small only when the high-flux pocket itself relaxes quickly.
* *What flux discreteness does and does not buy.* Discreteness makes the interior exactly vacuum-like and metastable. It does not change the
  cost.

## 8. What all seven share: the bridge invariant

**Theorem 6i-11 (bridge invariant; spherical symmetry).**
* *(1) Every hidden daughter contains a bridge.*
  * At the daughter's wall β_out < 0: the area radius decreases outward. It tends to infinity in our asymptotic region, so it has a minimum
    on every spacelike slice in between.
  * At that sphere n·∇r = 0, so ∇r is timelike or null and 1 − 2Gm/(c²r) = (∇r)² ≤ 0.
* *(2) Forming it is a black-hole formation event.*
  * Start from regular data (2m/r < 1 everywhere). The first such sphere is marginal: 2Gm = c²r.
  * There, with E = mc² and the mean density ρ̄ = 3mc²/(4πr³), we have E√ρ̄ = √(3/(32π)) c⁶/G^{3/2} ≈ 2.3×10⁶⁵ J^{3/2} m^{−3/2} exactly.
  * The value is scale-free. Negative energy lowers m, so it cannot help.
* *Assessment.* The matter energy that must be brought within r is at least m(r), since binding lowers m. The best technology
  (`accessibility_table`) is 10^39.5 short.
* *Natural events.* Core-collapse supernovae and neutron-star mergers reach the invariant (figure, panel e). They are natural black-hole
  formation events.

**Side-change lemma (thin-shell form).**
* β_out = (R(f_in − f_out) − m_s²/R)/(2m_s) is continuous where m_s ≠ 0, and β_out² = f_out + Ṙ².
* So on a classical trajectory, β_out changes sign only where m_s = 0 (crossed at the speed of light) or f_out ≤ 0.
* f_out ≤ 0 means inside a horizon. In a Schwarzschild region the future there is trapped, so the shell cannot re-expand.
* Check: 3351 of 3351 sign changes in 3000 random shells are explained.
* This is why the breathing bubble changes side only through m_s = 0, and why collision pockets and membranes cannot reach the hidden branch
  classically.

**Escapes from Theorem 6i-11.**
* *Quantum tunnelling.* Bounded by "no free lunch" (R10, R37) and, in the controlled regime M ≥ m_P, by P ≤ 10⁻²³ (R42).
* *A pre-existing bridge (a natural black hole).* Round 6c (Proposition 6c-15) found that the false vacuum must then extend inside the
  horizon over a horizon-scale range. Delivering that vacuum meets the same invariant.
* *A charged Cauchy horizon.* No help in the thin-shell model: the barrier persists (section 2).
* *An anti-trapped background.* The cosmological route, run by nature (Round 6g).

**Audit of Theorem 6h-1's wording.**
* Our universe has a past singularity, so "a universe with a regular past" must be read locally. The Penrose argument closes once the past
  null generators of the anti-trapped surface reach their focal points, within an affine length ≈ 2/|θ| ≈ R. A laboratory region today has a
  regular past over that length.
* For a closed universe, the boundary of the past of a small surface cannot be a Cauchy surface of the whole universe. So the conclusion
  stands; only the wording is sharpened.

## 9. New methods beyond the seven (examined this round)

* **Monopole plus Reissner–Nordström Cauchy horizon (escape iii above).** Computed: with added mass from 1× to 10× the static value,
  including past the horizon threshold, the forbidden band between the bound region and the inflating branch persists (`added_mass` in the
  data). No classical path.
* **Natural black holes as the bridge.** A stellar-mass black hole pairs with a QCD-scale false vacuum (ρ^{1/4} = 0.33 GeV (M_☉/M)^{1/2}, corrected in Round 6j,
  Round 6c). The region to be filled is horizon-scale, so its energy is about the hole's own mass at nuclear density: natural events only.
* **Q*: barely metastable vacuum with a light wall (outside controlled physics; Conjecture-level door).**
  * *The construction.* For a false vacuum ε above ours, choose the wall so that its flat thin-wall CdL exponent is B_CdL.
  * *An exact identity.* In the light-wall limit the hidden inflating branch at M → 0 starts at R₃ = 2κ/χ² = 3σ/ε = R_CdL.
  * *The exponent.* Semiclassically it is about B_CdL·Φ (R41).
  * *The numbers:*

    | False vacuum | B_CdL | Bubble radius | Bubble mass | Bridge size |
    |---|---|---|---|---|
    | 500 GeV | 1 | 3.5×10⁻¹⁹ m | 2.9 TeV | 4.7×10⁻¹⁶ l_P |
    | GUT | 1–10 | — | — | 0.01–0.05 l_P |

  * *The general relation.* Bridge/l_P ≈ 9.4 (vacuum scale/M_Pl) B_CdL^{3/4}. An unsuppressed channel with a classical bridge (≥ l_P) needs
    a near-Planckian vacuum.
  * *Status.* The channel is not excluded by any established result, but its semiclassical description is invalid by 16 orders of
    magnitude. It also needs a TeV-scale, barely metastable false vacuum, which no known physics provides.
  * *A consistency remark.* If such processes occurred, ultra-high-energy cosmic rays would have produced them for billions of years. The
    daughters would be hidden and harmless, so observation gives no constraint either way.
* **Vacuum-dependent condensate (route 4)** is the only route this round that lowers a requirement by many orders (background 2.5×10⁹ M_☉
  → 9×10⁸ kg). It leaves the invariant untouched.

## 10. Status

| Research objective | Status after Round 6i |
|---|---|
| Create the initial false-vacuum region | **Yes, by the universe** (Round 6g) |
| Establish a technologically accessible mechanism | **Not found.** All seven proposed routes are computed and fail at identified points (ledger). Every classical route must form a black hole at E√ρ̄ ≈ 2.3×10⁶⁵ (Theorem 6i-11), 10^39.5 beyond technology. Controlled quantum routes obey "no free lunch" and P ≤ 10⁻²³. The one door no established result closes (Q*) lies 16 orders of magnitude below the Planck length, where no calculation is valid. Not proved impossible for every conceivable process: Theorem 6i-11 is classical, and quantum gravity is unknown. |

## Literature used (Round 6i)

* *Spontaneous scalarisation and Gauss–Bonnet couplings*
  * Doneva & Yazadjiev, PRL 120, 131103 (2018); Silva et al., PRL 120, 131104 (2018): spontaneous scalarisation.
  * Ripley & Pretorius, PRD 99, 084014 (2019), arXiv:1902.01468: loss of hyperbolicity in spherical sGB collapse.
  * Thaalba, Bezares, Franchini & Sotiriou, arXiv:2306.01695: a Ricci coupling tames it.
  * Kanti, Kleihaus & Kunz, arXiv:1108.3003: wormholes without exotic matter in dilatonic EGB.
  * Cuyubamba, Konoplya & Zhidenko, PRD 98, 044040 (2018), arXiv:1804.11170: these wormholes are unstable.
  * GW230529 tests of GR, arXiv:2406.03568: √α_GB ≲ 0.28 km; arXiv:2405.13279: 0.26–0.30 km.
  * Nonlinear scalarisation dynamics: EPJC (January 2026); Herdeiro et al., arXiv:2603.24164 (phase structure). The user's "January 2026
    PRL" could not be identified.
* *Monopoles*
  * Sakai, Nakao, Ishihara & Kobayashi (2006), gr-qc/0602084: the universe out of a monopole in the laboratory.
  * Sakai et al., gr-qc/9808069: baby universes in monopole collisions.
  * Topological inflation threshold η ≳ m_P/4 for global defects (cited by Lepora & Martin).
* *Breathing bubbles*
  * Guendelman & Sakai (2008), arXiv:0803.0268: the universe out of a breathing bubble.
* *Gauged ghost condensation*
  * Cheng, Luty, Mukohyama & Thaler, JHEP 05 (2006) 076, hep-th/0603010: M < min(10¹² GeV, g²10¹⁵ GeV); the minimal theory M ≲ 100 GeV.
* *Resonant tunnelling*
  * Saffin, Padilla & Copeland, JHEP (2008), arXiv:0804.3801: resonant decay of an inhomogeneous (contracting-bubble) state.
  * Copeland, Padilla & Saffin: no resonant tunnelling in standard scalar QFT.
* *Classical transitions*
  * Easther, Giblin, Hui & Lim (2009), arXiv:0907.3234: classical transitions.
  * Johnson & Yang (2010), arXiv:1005.3506: gravity allows higher vacua than the colliding ones, but not above the parent.
  * Dray & 't Hooft (1985); Redmount (1985): the collision relation.
* *Four-form flux*
  * Brown & Teitelboim, PLB 195, 177 (1987); NPB 297, 787 (1988): membrane nucleation.
  * Bousso & Polchinski, JHEP 06 (2000) 006, hep-th/0004134.
* *Earlier rounds*
  * Penrose (1965); Farhi & Guth (1987); Farhi, Guth & Guven (1990); Fischler, Morgan & Polchinski (1990).
