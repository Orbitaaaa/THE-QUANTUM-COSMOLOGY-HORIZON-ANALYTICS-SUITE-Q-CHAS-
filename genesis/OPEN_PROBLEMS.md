# Open problems

Grouped by the phase that should attack them. Each item has a precise statement and a success criterion.

## Phase II — verification gaps

* **P-N1 (both ends proved, round 3 D2).** Prove NG4: in regime B, $\lvert\theta_D\rvert>\pi$. Both boundary values are now theorems: $F(t)>1$ at $M_{\rm cr}$, and $(\lvert J_0\rvert+\pi)/2\pi>1$ at $E_{AB}^-$. The interior is strictly monotone on 4800 points. Remaining (**P-N1′**): monotonicity of $EJ(E)$ on regime B.
* **P-T5 (mostly done, roadmap Part A).** A theorem for light walls (S3: $B=B_{\rm CdL}\Phi$), and in the double-flip class wherever the $\Delta t^W_E<0$ criterion holds (S5).
  Rigorous closed-form bounds S1 and S2. Round 3: the A/B end of the strip is a theorem (Prop. S5-b); certificates hold from $1.02E_{AB}$ (sampled);
  $J<0$ on 2100 points down to $10^{-7}$. Open (**P-T5″**): interval-verified certificates on $(E_{AB},E_\omega)$ and a quantitative remainder at the end.
* **P-S3 (new, round 3 D3).** Derive the correction coefficient $c(E)$ of S3 analytically: $-0.73<c<0$, $c(0^+)=-0.722221$ (conjecture $-13/18$), tail $-5E^{-3/2}$.
* **P-L1.** Compare T0–T3 with FGG's pseudomanifold rule and FMP's Hamiltonian formulae, and with
  de Alwis et al. / Céspedes et al., term by term (full texts required; arXiv is blocked in this environment).
* ~~P-G1~~ **done**: the atlas (roadmap Part G) and, in round 3, the Euclidean sections with conical defects, $I(\beta)=-(A_b+A_c)/4G$ for every $\beta$ (`sds_atlas.euclidean_section`). Original text: Explicit Penrose diagrams for the SdS exterior (patched Kruskal atlas at $r_b$ and $r_c$); Euclidean SdS
  pieces with two horizons (conical defects).

## Phase III — the quantum transition

* ~~P-Q2~~ **done (roadmap Part B):** exact integration of the momentum-space equation confirms $\sqrt2$ for three orderings, with corrections $h(A\ln h+B)$. Original text: The Stokes multiplier at the coalescence point is $\sqrt2$ for every ordering of the squared
  constraint, at leading local order (PHASE_III Prop. III-8). Confirm it with an exact solution of the $M=0$, $\kappa>\chi$ squared constraint
  (momentum-representation ODE), and compute prefactors.
* **P-Q9 (reformulated, roadmap Part C).** The orientation is configuration data, so it is not removed by squaring (Q9.1). The two squared families are two
  conventions for one daughter (Q9.2). The question is which Casoratian charge is probability (Q9.3/Q9.4), and it is moot outside the sliver $M_{\rm sd}<M$
  (Theorem S4). **Round 3:** under P2–P3 the question is settled in the reduced theory (Theorem A1: $e^{-B_{\rm an}}$ for double flip, no flux for single flip). Remaining **P-Q9′**: derive P3, i.e. show the Casoratian charge is the reduction of the full WDW Klein–Gordon current through the bridge. Original text: Is the exterior orientation of the wall superselected in a full theory? This decides between $B_{\rm an}$ and $B_{\rm sq}$ in the
  single-flip class (the double-flip class is unambiguous).
* **P-Q7.** Pseudo-differential quantization of the single sheet $H_-$ (real-angle convention): is the tunnelling amplitude governed by
  the branch locus, as with relativistic Schrödinger operators?
* **P-Q8 (partly done, roadmap C5).** The complex angles follow from analytic continuation of the WDW phase through the bulk forbidden zone. The measure does not. Original text: Derive the complex-angle convention (C2) and the Lebesgue measure (C3) from a covariant inner product
  (group averaging, or the crossed-product construction).
* **P-Q3.** Contour analysis (NEW_APPROACHES N3): show that Picard–Lefschetz contours of the fixed-energy minisuperspace path
  integral reduce to the real barrier class, or classify the Stokes phenomena that allow horizon-encircling contours.
* **P-T3 (partly done, roadmap D5).** There is one attempt per bubble (the bounded branch is not periodic in $t_\infty$), and the Kemble connection gives $P\simeq e^{-B}$. Open
  (P-T3′): the field-theoretic one-loop prefactor.
* **P-Q1.** Dirac (quantize-then-reduce) quantization with Kuchař variables; compare with Route A beyond WKB.
* **P-Q4 (revised, round 3 F17).** The roadmap-round claim (boundary-clock unitarity never produces the future-expanding daughter as a positive-norm subsystem) was a
  $t_\infty$-clock artefact. In the constraint's own charge the expanding double-flip daughter has positive flux. Open: a QFT-level statement of what the
  exterior can and cannot know (Wall, Freivogel et al.).

## Phase IV — fields instead of walls

* **P-T2.** Thick-wall versions of T2/T3 and of NG3 (Euclidean scalar-field FGG solutions).
* ~~P-D1~~ **done** (PHASE_IV): spherical Einstein–scalar evolution of laboratory bubbles and of both daughter branches; Prop. IV-1.
* **P-D2 (partly done).** Aspherical stability. Probe-wall amplification to $2GM$ is at most about 2. Open: self-gravitating aspherical perturbations
  of the bounded branch and of the inflating bounce branch; a 3+1 evolution.
* **P-D3.** Evolve the white-hole branch into its past (the past singularity required by Cor. IV-2); scan the potential family.

## Phase V — the full quantum system

* **P-Q5 (partly done, Phase V).** Reduced model: the daughter is not superselected, but it is invisible to the abelian exterior algebra
  (Theorem V-2) and needs its own clock after the edge (Theorem V-1). Open: the von Neumann-algebraic version with matter fields
  (crossed products, NEW_APPROACHES N2).
* ~~P-Q10~~ **done (roadmap Part D):** no boundary condition is needed. The edge is a coordinate singularity (Q10.1). The asymptotic clock folds back along every daughter (Q10.2), and the
  physically expanding daughter is a negative-flux channel of the boundary-clock evolution (Q10.3). Original text: The boundary condition at the daughter edge (interior horizon at finite $t_\infty$) and a daughter clock: derive them, rather than
  postulating an outgoing condition.
* **P-Q6.** Information: what the exterior observer's reduced state looks like after a transition (entropy bookkeeping:
  the Morse/endpoint areas of T1/T2 are suggestive).

## Phase VI — observability

* **P-O2 (scaling done, roadmap E4).** $f(M)\simeq\frac{4\pi}3\lambda(M_{\rm eq}/M)^{1/2}$. Open: GVZ formation efficiencies and the supercritical fraction. Original text: Derive the mass function of daughter-hosting PBHs from bubble nucleation during inflation (FRW embedding of the thin-wall and
  Phase IV dynamics); Phase VI supplies only the characteristic mass $M_*\simeq1/(2GL)$.
* **P-O1.** Rates of daughter production by vacuum transitions in the early universe (bubbles nucleated during inflation,
  Garriga–Vilenkin–Zhang; Deng–Garriga–Vilenkin) versus laboratory rates. Translate T5 into a statement about
  primordial-black-hole populations.

## Roadmap round (ROADMAP_STAGE1_2_4.md) — new items

* **P-Q9′.** Klein–Gordon current of the spherically symmetric Wheeler–DeWitt equation through the Einstein–Rosen bridge (imaginary part of the Kuchař phase).
* **P-T3′.** One-loop prefactor (wall and bulk fluctuations) of the laboratory transition.
* **P-T5′.** Analytic control of $\Delta t^W_E$ in the strip $E_{AB}<E<1.6E_{AB}$. (The bound $\Phi\ge\max(1,8/3E+2/3)$ is proved.)
* **P-A1.** The QEI species bound with curved-space, interacting QEIs and exact constants (currently a scaling argument with Ford–Roman constants).

## Round 3 additions

* **P-C4.** The field-theoretic burst at the transition (Bogoliubov coefficients of fields coupled to the wall), the only candidate signal present on the
  daughter branch rather than the collapse branch. Needs P-T3′.
* **P-C2.** Event-level discrimination (Prop. C1): a model of the wall's emission rate $\Gamma$ in a scalar-field wall, to turn $\Gamma\tau_{\rm collapse}$ into a number.
* **P-A2.** Extend Theorem A1 to the SdS exterior and to the non-flip window (E3), where $n_D$ can differ.

## Round 6 additions (ROADMAP_ROUND6_MECHANISM_ATLAS.md)

* **P-R6a.** The GSL and the area theorem in beyond-Horndeski genesis: is GSL violation unavoidable for self-amplifying NEC violation?
* **P-R6b.** A finite genesis seed whose boundary layer opens a throat with the G1 budget (otherwise Proposition G7 gives a parent-consuming
  visible bubble), and the hand-over to the daughter's vacuum energy. This is the first calculation Route G needs.
* **P-R6c.** A stabiliser compatible with a bridge, for Planck-density seeds (G5).
* **P-R6d.** The Kerr (rotating) route through the inner horizon in 3+1.
* **P-R6e.** Rates of flux-supported topology change (wineglass wormholes) seeded by laboratory fields.

## Rounds 5 and 7 additions (ROADMAP_ROUNDS5_7_12.md)

* **P-T2 (now central).** The thick-wall laboratory exponent. Perturbative theories reach small $B_{\rm CdL}$ only with thick walls, where S3 and S7 are not established.
* **P-C4 (sharpened).** The field's state after the transition on the daughter branch, from a full treatment. Round 7 shows the branch difference is set entirely by it
  (about 0.1 quanta per event under D1).
* **P-R7.** 3+1 greybody factors and higher multipoles for the Round 7 flux.
