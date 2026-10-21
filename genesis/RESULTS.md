# Results (verified)

Only results that passed the checks in VERIFICATION.md are listed. Labels: **Theorem** (proof given),
**Proposition-N** (numerically established, proof open), **Reproduction** (known result re-derived independently).

## Reproductions (established theory, re-derived and checked)

* **Israel conditions, effective potential, mass formula, sign functions** for dS|SdS thin walls.
  JUNCTION_CONDITIONS.md; symbolic verification.
* **Blau–Guendelman–Guth global structure.** Inflating walls bounce beyond the Einstein–Rosen bridge, and visible
  bubbles collapse. GLOBAL_GEOMETRY.md; 15 branches, two independent constructions.
* **Farhi–Guth obstruction**, stated precisely with its exact hypotheses, plus an audit of 42 branches.
  NO_GO_THEOREMS.md Sec. 1.
* **Coleman–De Luccia thin-wall limit with gravity**, reached by thick-wall bounces with $O(w)$ convergence.

## Original results of Stage One

| # | result | label | file |
|---|---|---|---|
| R1 | **Lemma J1**: with σ > 0, squaring the Israel condition creates no spurious branch; the only spurious solutions of the squared equation are the negative-tension ones | Theorem (elementary) | JUNCTION_CONDITIONS.md §4 |
| R2 | **Chart-free itinerary** (Proposition G2): causal regions from (sgn f, ε, sgn Ṙ), valid for any number of horizons | Theorem + cross-validation | GLOBAL_GEOMETRY.md §3 |
| R3 | **Branch selection T0**: the canonical WKB phase is the analytic continuation along the barrier. The FGG slicing jump is a horizon monodromy ($P\to P\pm i\pi R/G$, $t\to t\pm i\beta/2$). *Refined in Phase III (R19–R23): the selection holds given the complex-angle convention; the sector normalisation remains a datum* | Proposition (analytic argument + two numerical criteria) | TUNNELLING.md §3 |
| R4 | **T1**: closed form at $M=0$, $B=S_E[\text{instanton}]+\frac1{4G}(\sum A_{\max}-\sum A_{\min})$; equals the CdL decay exponent of the false vacuum when the false vacuum is the larger cap | Theorem | TUNNELLING.md §4 |
| R5 | **T2 first law** $d\tilde B=\frac{{\rm Vol}^{\rm in}_E}{8\pi G}d\Lambda_{\rm in}+A^W_Ed\sigma-\Delta t^W_EdM$ | Theorem | TUNNELLING.md §4 |
| R6 | **T3 Smarr identity** $B=-\frac{\Lambda_{\rm in}}{8\pi G}{\rm Vol}^{\rm in}_E-\frac\sigma2A^W_E-\frac M2\Delta t^W_E+\frac{\mathcal A}{4G}$: a covariant expression for the FGG exponent that never uses the exterior pseudomanifold; integral representation from the $M=0$ instanton | Theorem | TUNNELLING.md §4 |
| R7 | **NG3 exterior embedding obstruction**: no manifold Euclidean interpolation for any transition from regular laboratory data to the inflating bounce | Theorem | NO_GO_THEOREMS.md §3 |
| R8 | **NG4**: regime-B transitions are also non-manifold (horizon-angle overshoot $\lvert\theta_D\rvert>\pi$) | Proposition-N | NO_GO_THEOREMS.md §3 |
| R9 | **Interior counter-example**: $\varepsilon_{\rm in}$ flips do not obstruct embedding; the obstruction needs a flip while the wall recedes from a Euclidean horizon | Proposition-N | NO_GO_THEOREMS.md §3 |
| R10 | **T5**: over laboratory transitions, $B\ge0.56\,B(0)$ for all γ, and $B\ge B(0)=B_{\rm CdL}$ for γ ≳ 0.2. Creating a daughter from laboratory data is not exponentially easier than vacuum decay | Proposition-N | TUNNELLING.md §4 |
| R11 | **Physical scale of the exponent**: $B(M\to0)$ follows the flat CdL scaling $\sigma^4/\rho^3$ at small κ/χ, not $(M_{\rm Pl}/H)^2$, but $B$ rises steeply with $M$. GUT example (σ^{1/3} = ρ^{1/4} = 10^{16} GeV): $B=133$ only for nanogram (sub-Planckian) masses, whose black holes are not semiclassical (Phase V); $B\gtrsim7\times10^4$ for $M\gtrsim10\,m_{\rm Pl}$, $B\sim10^5$–$10^7$ for milligram–gram masses | computation | TUNNELLING.md §4 T6 |
| R12 | **Unification**: the classical (FG), initial-data (NG2) and Euclidean (NG3) obstructions all stem from creating a marginally trapped surface (the bridge) in an asymptotically flat exterior | Theorem (each part) + Conjecture NG5 (general matter) | NO_GO_THEOREMS.md §4 |
| R13 | **Sector structure**: the four orientation sectors are sheets of one squared constraint, which communicate only under the barrier; the Hamiltonian generates imaginary shifts of the wall area by $8\pi\ell_P^2$ | derivation | QUANTIZATION.md §2–3 |

## Phase II results (PHASE_II_VERIFICATION.md)

| # | result | label |
|---|---|---|
| R14 | **Lemma P** (positivity, $0<\Theta_{\rm out}-\Theta_{\rm in}<\pi$) and the arccos formula $B=\frac2G\int R\arccos C\,dR$ | Theorem |
| R15 | **NG4 critical limit**: $\lvert\theta_D\rvert/\pi\to F(t)=\frac{3\sqrt3(t+1)^2}{4(t+2)\sqrt{t^2+2}}>1$ as $M\to M_{\rm cr}$ | Theorem |
| R16 | **T3′**: Smarr identity and first law for Schwarzschild–de Sitter exteriors ($\Lambda_{\rm out}$ term = excised exterior volume) | Theorem, verified to $10^{-21}$ |
| R17 | **Exact form of the FGG difference**: $B_{\rm an}-B_{\rm princ}$ = flip-radius area terms over $4G$ | numerical identity, 10 digits |

## Phase III results (PHASE_III_QUANTUM_TRANSITION.md)

| # | result | label |
|---|---|---|
| R18 | **III-1 horizon monodromy**: contour windings change $B$ by $\sum n_h\sigma_h(A(R_2)-A_h)/4G$; zero-net windings by $\pm(A_{h_2}-A_{h_1})/4G$ | Theorem + numerics ($10^{-4}$, $10^{-6}$) |
| R19 | **III-2 conventions**: analytic branch = complex (Sorkin) boost angles, satisfying a sign-free analytic constraint; principal branch = real angles, satisfying the sign-dressed constraint; related by a non-unitary anti-periodic intertwiner | Theorem |
| R20 | **III-4 branch point = flip = fold**: the true Hamiltonian is two-sheeted. The barrier trajectory touches its branch locus tangentially, exactly at the exterior flip radius, where the Euclidean wall folds in Killing time. Sector-diagonal true-Hamiltonian quantizations give zero transition amplitude | Theorem + numerics ($10^{-15}$) |
| R21 | **III-5 failure of the FGG principal construction**: the sign-dressed shift constraint has no non-zero analytic solution across a flip radius; $B_{\rm princ}<0$ in 33/36 models | Theorem + numerics |
| R22 | **III-6 no conserved current**: Casoratians of the area-shift constraint are $i\hbar G$-periodic; symmetric measures are functions of $e^{A/4\hbar G}$; Lebesgue is unique below that growth | Lemma |
| R23 | **Verdict**: no formulation removes the FGG ambiguity by its own rules; all place it in one datum (flipped-normal angle ⇔ sheet of $H$ ⇔ sector normalisation ⇔ fold action). The FGG principal splice fails (R21) | Theorem-level classification |
| R25 | **Lemma III-9 + Prop. III-8 (sector-blind quantization)**: $N=f_{\rm in}+f_{\rm out}-\kappa^2R^2$ is strictly concave, so $C=0$ occurs on the barrier iff the transition is single-flip. Every ordering of the squared constraint has Weber index $\nu=-\tfrac12$ there (Stokes multiplier $\sqrt2$), so its exponent is $B_{\rm sq}=\frac2G\int R\arccos\lvert C\rvert dR$. Hence **double flip: all consistent quantizations give $B_{\rm an}$ (ambiguity eliminated); single flip: $B_{\rm an}$ (orientation superselected) vs $B_{\rm sq}<B_{\rm an}$ (not)** | Theorem (lemma) + leading-order proposition + numerics ($10^{-15}$) |
| R24 | **Covariant action rule** for the laboratory pseudomanifold (GHY at ∞, background, bolt, folds as Hayward joints of angle π) equals $B_{\rm an}$; Lorentzian complex saddle ${\rm Im}\,T^*$ | Proposition, verified to $10^{-9}$ |

## Phase IV results (PHASE_IV_FIELD_DYNAMICS.md)

| # | result | label |
|---|---|---|
| R26 | Double-null Einstein–scalar code with a regular centre and singularity excision; exact in flat space, second order, validated against exact de Sitter and critical collapse | verified code |
| R27 | **Laboratory bubbles collapse to black holes** and converge to the Israel thin-wall dynamics as $w\to0$ ($R_{\max}\to R_1$, $r_{AH}\to2GM$, $A_{\rm in}-A_{\rm out}\to\kappa R$) | numerics |
| R28 | **Prop. IV-1**: in spherical symmetry with NEC matter, $r_u/\alpha^2$ is non-increasing along ingoing rays, so laboratory data can never produce an anti-trapped (inflating) region. Confirmed: 0 anti-trapped points in 10/10 laboratory runs. **Cor. IV-2**: every anti-trapped point lies to the future of initial data with $r_u>0$ (beyond a bridge, a white hole, or beyond a cosmological horizon) | Theorem + numerics |
| R29 | **Daughters exist as field-theory solutions**: the bounce branch ($\varepsilon_{\rm out}=-1$ throughout) from regular-centre data with a reversed-normal exterior, and the white-hole branch (regime C). Both inflate and converge to the thin-wall trajectories ($\dot R$ within 0.14–0.4%) | numerics |
| R30 | Thin-wall predictions accurate to about 1% for $w/R\lesssim0.05$; at $w/R=0.4$, $r_{AH}=0.79\cdot2GM$, but collapse with no anti-trapping persists | numerics |
| R31 | Aspherical growth before horizon formation is at most about 2 ($\ell=2$), in the probe Jacobi analysis (zero modes reproduced to $5\times10^{-6}$) | estimate (probe) |

## Phase V results (PHASE_V_QUANTUM_SYSTEM.md)

| # | result | label |
|---|---|---|
| R32 | **Theorem V-1**: the reduced theory with the asymptotic clock is a scattering problem. The parent channel ends at the black-hole horizon at $t_\infty=\infty$ ($t\simeq-2GM\ln(R-2GM)$, verified to $2\times10^{-4}$); the daughter channel ends at the interior de Sitter horizon at **finite** $t_\infty$ with $P\to-\infty$, an edge of the reduced phase space needing a boundary condition. The daughter needs its own clock | Theorem (WKB-level for self-adjointness) + numerics |
| R33 | **Theorem V-2**: the daughter is not superselected in the reduced theory, but it is invisible to every exterior observable (functions of $\hat H$) | Theorem |
| R34 | **Prop. V-3**: the transition preserves spatial topology ($\mathbb R^3\to\mathbb R^3$) and creates an Einstein–Rosen throat. Pinch-off needs a singularity (Geroch/Tipler) | Proposition |
| R35 | Semiclassical validity: WKB parameter $10^{-11}$ (GUT), $10^{-2}$ ($10^{18}$ GeV). Correction: Stage One's $B=133$ samples are sub-Planckian black holes. Total daughter probability, even with repeated attempts until evaporation, $\le10^{-184}$ | computation |

## Phase VI results (PHASE_VI_OBSERVABILITY.md)

| # | result | label |
|---|---|---|
| R36 | **Prop. VI-1**: a black hole hosting a daughter is indistinguishable from an ordinary one by any exterior measurement (spherical symmetry) | Proposition |
| R37 | **Prop. VI-2**: laboratory daughter probability $\le e^{-0.56B_{\rm CdL}}$. It is small only where the false vacuum is not metastable ("no free lunch") | Proposition (conditional on T5) |
| R39 | **Loophole search**: minimise $B$ over vacuum and wall scales under metastability, semiclassical black hole and WKB validity. Best laboratory case $B=107$, $P\le10^{-46.5}$ ($M=10\,m_{\rm Pl}$, $\rho^{1/4}\approx9\times10^{17}$ GeV), $10^{-656}$ for $M=100\,m_{\rm Pl}$. The optimum always sits at the edge of validity. Prop. IV-1 extends to any NEC matter | computation + Theorem |
| R38 | $M_{\rm cr}\to1/(2GL)$ as $\kappa/\chi\to0$ (analytic + $2\times10^{-6}$). Characteristic daughter-hosting black-hole mass $M_*\approx5.6\,{\rm g}\,(10^{16}{\rm GeV}/E_{\rm vac})^2$, which maps vacuum scales to PBH windows (asteroid window: $E_{\rm vac}\approx2\times10^5$–$7\times10^7$ GeV) | Theorem (limit) + scale estimate |

## Statement for a general audience

*We asked whether gravity lets a new, inflating universe grow out of a small bubble of false vacuum made in
an otherwise empty, flat universe.* Classically it cannot: a theorem of Farhi and Guth forbids it whenever
energy obeys the null energy condition. Quantum mechanically, the textbook calculation by Farhi, Guth and Guven ran
into a geometric object that is not a manifold. We proved that this object is unavoidable for every
"laboratory" starting point. We then found an exact identity, a Smarr relation, that gives the transition
exponent from well-defined geometric quantities without that object. The result: the probability is not
exponentially larger than the probability that the false vacuum simply decays. For grand-unification-scale
parameters and gram-scale bubbles the exponent is of order $10^6$ or much larger, i.e. a probability like
$e^{-1\,000\,000}$. Nothing in the analysis points to a practical route.
