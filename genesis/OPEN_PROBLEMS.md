# Open problems

Grouped by the phase that should attack them. Each item has a precise statement and a success criterion.

## Phase II — verification gaps

* **P-N1 (both ends proved, round 3 D2).** Prove NG4: in regime B, $\lvert\theta_D\rvert>\pi$. Both boundary values are now theorems: $F(t)>1$ at $M_{\rm cr}$, and $(\lvert J_0\rvert+\pi)/2\pi>1$ at $E_{AB}^-$. The interior is strictly monotone on 4800 points. Remaining (**P-N1′**): monotonicity of $EJ(E)$ on regime B.
* **P-T5 (mostly done, roadmap Part A).** A theorem for light walls (S3: $B=B_{\rm CdL}\Phi$), and in the double-flip class wherever the $\Delta t^W_E<0$ criterion holds (S5).
  Rigorous closed-form bounds S1 and S2. Round 3: the A/B end of the strip is a theorem (Prop. S5-b); certificates hold from $1.02E_{AB}$ (sampled);
  $J<0$ on 2100 points down to $10^{-7}$. Open (**P-T5″**): interval-verified certificates on $(E_{AB},E_\omega)$ and a quantitative remainder at the end.
* ~~P-S3~~ **done (Round 6c, Theorem 6c-S3):** $c(E)=2-2mI_0'/I_0+I_1/I_0$ in closed form, $c(0^+)=-13/18$ exactly, tail $-(16/\pi)E^{-3/2}$. Original text: Derive the correction coefficient $c(E)$ of S3 analytically: $-0.73<c<0$, $c(0^+)=-0.722221$ (conjecture $-13/18$), tail $-5E^{-3/2}$.
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
* **P-A1 (2D done, Round 6c):** exact weighted-ANEC bound $c/48\pi\tau$ for 2D CFTs (Prop. 6c-10) and its tipping corollary; 4D interacting QEIs remain open (no null QEI exists in 4D, Fewster–Roman). Original text: The QEI species bound with curved-space, interacting QEIs and exact constants (currently a scaling argument with Ford–Roman constants).

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
* ~~P-R7~~ **done (Round 6c, Proposition 6c-16):** Regge–Wheeler greybody factors for $l=0$–4; the branch difference is IR-finite in 3+1. Original text: 3+1 greybody factors and higher multipoles for the Round 7 flux.

## Round 6b additions (ROADMAP_ROUND6B_GSL.md)

* ~~P-R6f~~ **done (Round 6c):** the extremal throat cannot both hide and grow (Theorem 6c-1, Cor. 6c-2–6c-5); the Nariai throat does (Prop. 6c-9, 6c-12). Original text: 1+1 Einstein–Maxwell–scalar evolution of a bud in an extremal magnetic throat: its edges, hiding (GSL-10), persistence of the neck.
* **P-R6g (largely done, Round 6c):** LLL reduction and strong-QNEC saturation on curved 2D backgrounds proved (Prop. 6c-6/7); massive-mode QNEC on curved backgrounds remains [Literature only in Minkowski]. Original text: The QNEC and its strong 2D form for Landau-level fermions in a curved magnetic throat at large $q$, with massive-mode corrections.
* **P-R6h.** The sub-$T_{\rm brk}$ (Schwarzian) regime: anti-trapping and the bud as a quantum state in an SYK/JT model (escape E3; Route Q).
* **P-R6i (done modulo (S), Round 6c):** Theorem 6c-R6i proves GSL-3 without symmetry given the support form (S) of the local GSL; open: derive (S) from the GSL for non-smooth horizons. Original text: Theorem GSL-3 without spherical symmetry (needs a mean-convexity comparison for horizon cross-sections on normal data).
* **P-R6a (updated).** Beyond-Horndeski: the dynamical-entropy second law around time-dependent scalar backgrounds with $c_s\neq1$, and stability
  of a time-dependent inhomogeneous throat.

## Round 6c additions (ROADMAP_ROUND6C.md)

* **P-R6j (largely done, Round 6d):** the complete route from data without a black hole is simulated (Prop. 6d-5); open: the creation of the false-vacuum ball itself (P-R6m) and a band-localised push (P-R6l). Original text: Simulate the in-situ lift of a false vacuum around a black hole across its horizon (Proposition 6c-15) through to the tipping.
* **P-R6k.** The ghost condensate in the Nariai throat as an EFT simulation (gradient and Jeans dynamics included), not through its flux.
* **P-D2′.** Non-spherical stability of Nariai tipping (the $l\ge1$ modes on the Nariai $S^2$).
* **P-N1′ (attempted).** Reformulated as $-EJ=E^2\int sk(W-s^2)^{-1/2}ds$; the interval certificate is still open.

## Round 6d additions (ROADMAP_ROUND6D.md)

* **P-R6k.** Nonlinear ghost-condensate EFT simulation in the throat (Proposition 6d-2 is linear).
* **P-R6l.** A bulk NEC source located in the band, so the single-evolution route keeps the Nariai discount.
* **P-R6m.** The creation ("lift") of the false-vacuum ball by NEC-respecting dynamics.
* **P-R6n.** GSL accounting for the ghost condensate: S_out of the condensate and its quanta across daughter-observer horizons versus the
  required ~3√δ S_N (Proposition 6d-6); decides between the Mukohyama and Feldstein readings for this route.

## Round 6e additions (ROADMAP_ROUND6E.md)

* **Status of the Round 6d items.**
  * P-R6k advanced: the condensate is a dynamical, back-reacting source in its IR limit.
  * P-R6l done: the push is made by the condensate inside the band.
  * P-R6m largely done: creation by supercooling (Proposition-N 6e-10).
  * P-R6n replaced: Theorem 6e-6 (eternal daughters always violate the local GSL) and Proposition-N 6e-7 (transient daughter, GSL-consistent
    for complete observers).
* **P-R6e-a.** The full ghost-condensate EFT (k² and k⁴ terms, beyond the IR fluid) simulated nonlinearly, together with non-spherical modes.
* **P-R6e-b.** The GSL for observers who end at the daughter's crunch (Jacobson–Parentani horizons). Also long-lived transient daughters,
  and a persistent-neck (charged) variant.
* **P-R6e-c.** One GR simulation of the fireball emptying (free-streaming radiation, moving wall), with reflection and timing, handed off
  to the route's data.
* **P-R6e-d.** Microphysics of the false-vacuum tilt and of the thermal sector. That includes the strong quartic coupling of the toy
  units (λ ≈ 5×10⁷) and the realistic scaling with r_N.

## Round 6f additions (ROADMAP_ROUND6F.md)

* **Status of the Round 6e items.**
  * P-R6e-a done within the EFT and test-field proxies: nonlinear 1+1 EFT, l-mode proxies, full linear window (R122–R124).
  * P-R6e-b closed: Theorem 6f-3 (no GSL for observers ending at singularities, in any collapse).
  * P-R6e-c done, negative: heating cannot make the persistent normal ball (R119); Proposition-N 6e-10 withdrawn.
  * P-R6e-d partly answered by the problem-4 audit (R128); the microphysics of the tilt remains open.
* **P-R6f-a.** A creation mechanism for the route's ball. Thermal creation fails by an O(1) margin (R119). Candidates: NEC violation by
  the condensate during creation, a non-thermal (coherent plus dissipative) lift, a cosmological venue (late-blooming patches of a
  keV-scale first-order transition, R128).
* **P-R6f-b.** A daughter-adapted (cosmological-time) solver to follow the graceful exit and the recollapse in the route geometry (F45).
* **P-R6f-c.** The coupled non-spherical perturbation problem (even and odd parity, with the scalar, condensate and shell), beyond the
  test-field proxies of R123.
* **P-R6f-d.** Observational targets for problem 4: gravity in a Higgs phase at M ≲ 100 GeV (modified potential at L ~ M_Pl/M², growth
  on ~M_Pl²/M³), hidden keV-scale transitions (N_eff, BBN, gravitational waves), the ~10⁹ M_☉ primordial black-hole population.

## Round 6g additions (ROADMAP_ROUND6G.md)

* **Status.**
  * P-R6f-a (creating the ball) is answered by nature: a supercooled first-order transition leaves false-vacuum patches, and the
    near-critical ones are the route's ball; supercritical ones become hidden daughters (R131–R133).
  * There is no human route (Theorem 6g-2).
* **P-R6g-a.** Follow every supercritical patch until the exterior traps (f = 0.9 runs end first): the daughter-adapted solver of P-R6f-b.
* **P-R6g-b.** From nucleation statistics to the distribution of patch sizes, hence the abundance of daughters (as primordial black holes)
  and the gravitational-wave background, for a given transition scale.
* **P-R6g-c.** The cosmological route with a realistic plasma (heavy species, friction, latent heat in the walls) and statistically
  nucleated walls, beyond spherical symmetry.

## Round 6h additions (ROADMAP_ROUND6H.md)

* **P-R6h-a.** Near-critical false-vacuum interiors inside primordial black holes. Do any survive long enough for an infalling trigger to
  act? (Our subcritical runs: the false vacuum shrinks to nothing inside the horizon.)
* **P-R6h-b.** A condensate whose condensed phase exists only inside the false vacuum (P(X, ϕ)). Needs a ghost-free interpolation through
  the wall, the condensate's wall dynamics, and the observational bounds on its coupling to our vacuum. It would remove the 10⁹ M_☉ floor.
* **P-R6h-c.** The Oshita–Yokoyama rate in this project's potential: the Euclidean action, and the black-hole masses at which it is not
  negligible.
* **P-R6h-d.** A condensate coupled to controllable matter, which could let a machine source NEC-phase charge. Bounded by the
  ghost-condensate Lorentz-violation constraints.
* **P-R6h-e.** The near-critical scaling of the trigger in the cosmological background, which needs the daughter-adapted solver
  (P-R6f-b) to resolve outcomes close to the threshold.

