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

## Roadmap round: Stages 1, 2 and 4 (ROADMAP_STAGE1_2_4.md)

| # | result | label |
|---|---|---|
| R40 | **Lemma S1 / Theorem S2**: $B_{\rm an}\ge B_{\rm sq}\ge\frac{2\kappa}G\int R^2\sqrt{V/f_{\rm in}f_{\rm out}}\ge\frac{2\pi\kappa^2R_1^2f_{\rm out}(R_1)}{G(L^2+4GM/R_1^3)\sqrt{f_{\rm in}(R_1)f_{\rm out}(R_2)}}$, asymptotically sharp for light walls | Theorem |
| R41 | **Theorem S3 (thin-barrier law)**: $B_{\rm lab}=B_{\rm CdL}\Phi((R_{\rm CdL}/R_M)^2)(1+O(\kappa^2/\chi^2))$, $\Phi\ge\max(1,8/3E+2/3)$, so $B_{\rm lab}\ge4\pi^2\sigma^2R_M^2/\rho$ (no $G$). Checked to $6.5\times10^{-6}$ | Theorem (fixed $E$, $\Phi$ bound proved); Prop-N (uniform to $M_{\rm cr}$) |
| R42 | **Corollary S7 (frontier)**: $B_{\rm lab}\ge\frac43\sqrt{\pi/2}\sqrt{B_{\rm CdL}(M/m_P)(R_M/\ell_P)}$. Controlled regime ($B_{\rm CdL}\ge10$, $M\ge m_P$, $R_M\ge100\ell_P$): $P\le10^{-23}$ | Theorem (given S3) |
| R43 | **Theorem S4**: double flip ⇔ $\kappa<\chi$ and $M<M_{\rm sd}=(\chi^2-\kappa^2)/2G\chi^3$ (boundary $R_2=1/\chi$). The FGG ambiguity lives only in the single-flip class | Theorem |
| R44 | **Theorem S5**: an explicit pointwise criterion for $\Delta t^W_E<0$ (sufficient: $E/E_{AB}\ge\max\omega(z_i)$, $\omega=z^2/12+7/6z$); hence $B\ge B_{\rm CdL}$ in the double-flip class where it holds | Theorem; thresholds Prop-N |
| R45 | **P-Q2 closed**: exact momentum-space integration gives Stokes multiplier $\to\sqrt2$ for three orderings, four parameter sets; controls $\nu=0,-1$ | Proposition-N (exact numerics) |
| R46 | **Q9.1–Q9.4**: orientation = boost-invariant configuration data; the squared families are two conventions for one daughter; Casoratian charges $q^{-n}$; sector-independent probability gives $2e^{-B_{\rm sq}}$, an entropy-weighted measure gives $e^{-B_{\rm an}}$ | Theorem (+ principle) |
| R47 | **Q10.1**: the daughter edge is a coordinate singularity ($\tilde P=P\mp\frac R{2G}\ln|f_{\rm in}|$); no boundary condition (corrects R32) | Theorem + numerics |
| R48 | **Q10.2–Q10.3**: $dt_\infty/d\tau<0$ on every daughter. The boundary clock folds back; unitary boundary evolution yields the time-reversed daughter; the expanding daughter is a negative-flux (Klein-type) channel, $|\beta|^2\simeq e^{-B}$ | Theorem (reduced model); interpretation flagged |
| R49 | **Stage 4 ledger**: seed BH (F15), QEI species bound $N\ge\frac83S_{\rm dS}$ (cutoff $0.35H$), dS parent $e^{-3\times10^{122}}$, GVZ cosmological route, monopoles, FMM. No laboratory mechanism beats $B\gtrsim B_{\rm CdL}$ | Propositions + literature |
| R51 | **SdS non-flip window**: near-side daughters ($\varepsilon_{\rm out}(R_2)=+1$, no NG3 obstruction) exist only for $\eta\gtrsim0.76$–$0.86\chi$; closed for our universe | Proposition-N |
| R52 | **SdS conformal atlas**: eight Kruskal-glued blocks; wall worldlines timelike, future-directed and continuous across all horizons (84/84 branches) | construction + numerical check |
| R53 | **Cosmological route (scaling)**: $f(M)\simeq\frac{4\pi}3\lambda(M_{\rm eq}/M)^{1/2}$; all dark matter in $10^{20}$ g daughter hosts needs $\lambda\approx10^{-16}$ ($B_{\rm inf}\approx37$) | Proposition (scaling) |
| R50 | **Five numbers** (controlled optimum): $E_{\rm req}\approx3\times10^{17}$ GeV, $L_{\rm req}\approx100\ell_P$, $\tau_{\rm req}\lesssim5\times10^{-42}$ s, $P\le10^{-23}$, $S_{\rm verify}=0$ | computation |

## Round 3: the three fundamental problems (ROADMAP_ROUND3_PROBLEMS_ABC.md)

| ID | result | status |
|---|---|---|
| R54 | **Theorem A1 (charge selection rule)**: under P1–P3 (reduced theory; geometric states analytic in the strip, each configuration once; probability = Casoratian charge) the net flip number is superselected in the flux. Double flip: $K$ constant, ordinary flux balance; single flip: $K\equiv0$, total reflection | Theorem (conditional on stated principles) |
| R55 | **Theorem A2 (exact toy)**: $\psi(u+ih)+\psi(u-ih)=-2\lambda u\psi$ is solved by $Z_{-iu/h}(i/\lambda h)$; $H^{(1)}$ is the unique solution geometric in both regions, $K\equiv0$, and $\lvert H^{(1)}(-u)\rvert=e^{\pi u/h}\lvert H^{(1)}(u)\rvert$: the decay is the analytic-angle exponent exactly, for every $h$ | Theorem + numerics |
| R56a | **Lemma A6**: every conserved nearest-shift current is $\Phi(q)K$ with $\Phi$ periodic, $q=e^{A/4\hbar G}$; P3 is equivalent to the absence of a current weight $\propto e^{\pm A/4\hbar G}$ | Theorem |
| R56 | **Corollary A3 (Problem A)**: $P=e^{-B_{\rm an}}(1+O(\hbar))$ for double flip (anchor-free); no flux into single-flip daughters. Supersedes $2e^{-B_{\rm sq}}$ (R46) | Theorem under P1–P3 |
| R57 | **Proposition A4/A5 (daughter clock)**: the expanding double-flip daughter carries positive conserved charge; both Kruskal null coordinates increase along parent and daughter on both sides of the wall; only $t_\infty$ folds back. Withdraws the physical reading of R48 (F17) | Theorem + numerics |
| R58 | **Proposition B1 (Problem B)**: dynamical assistance needs $\omega\gtrsim1/\Delta\tau_E=0.3$–$0.6\,E_{\rm vac}$ ($1.6\times10^{17}$ GeV at the frontier) | Proposition |
| R59 | **Proposition C1 (Problem C)**: event-level signature: collapse branch has a late tail $\propto e^{-u/4GM}$ with $\Gamma\tau_{\rm collapse}$ photons; daughter branch stops at the turning point (last redshift $\sqrt{f(R_1)}\approx0.74$). State-level no-go unchanged | Proposition + numerics |
| R60 | **Proposition S5-b (T5 strip, A/B end)**: $\lim_{E\to E_{AB}^+}J=J_0+\pi<0$ for all $\gamma$, via $(z^2+z+1)^2-z^2h'=(z-1)(z+1)(z^2-z+3)(z^2+z+1)/z^2$; interior of the strip: certificates from $1.02E_{AB}$, $J<0$ on 2100 points | Theorem (end); Prop-N (interior) |
| R61 | **NG4 at both ends**: $\lvert\theta_D\rvert/\pi\to(\lvert J_0\rvert+\pi)/2\pi>1$ at $E_{AB}^-$ (new) and $F(t)>1$ at $M_{\rm cr}$; strictly monotone on 4800 points | Theorem (ends); Prop-N (interior) |
| R62 | **S3, independent 25–35-digit check**: $B_{\rm an}=B(0)\Phi(E)[1+c(E)\kappa^2/\chi^2]$, $-0.73<c<0$ on all of regime A, $c(0^+)=-0.722221$ ($-13/18$?), $c\simeq-5E^{-3/2}$; S7 carries $1-0.73\kappa^2/\chi^2$ | Prop-N (+ conjecture) |
| R63 | **Euclidean SdS with conical defects**: $I(\beta)=-(A_b+A_c)/4G$ for every period $\beta$ (SdS Smarr cancels the $\beta$-terms) | Theorem + numerics |

## Round 6: the mechanism atlas (ROADMAP_ROUND6_MECHANISM_ATLAS.md)

| ID | result | status |
|---|---|---|
| R64 | **Theorem G1 (classical budget)**: anti-trapping from laboratory data needs $2Gm/r\vert_0+8\pi G\int r(-T_{uu})du/\alpha^2\ge1$; NEC-violating null force $c^4/8\pi G=4.8\times10^{42}$ N, independent of size; Casimir short by $\ge50$ orders | Theorem |
| R65 | **Theorem G2**: thin walls of any surface matter cannot cross the bridge classically (at $A_{\rm out}=0$, $V=f_{\rm out}$); the past of the flip is a singularity, a Cauchy horizon or a cosmological region | Theorem |
| R66 | **Proposition G3**: visible inflating bubbles need negative asymptotic tension; at $M=0$ their nucleation exponent equals the CdL exponent, so Minkowski would decay upward (survival needs $B\gtrsim560$); the parent is consumed | Theorem + Proposition |
| R67 | **Proposition G4**: Reissner–Nordström exteriors change the laboratory exponent by $<3\%$ but give a classical route through the Cauchy horizon above $M_{\rm cr}\approx0.16\,m_P^3/\sqrt\rho$ (continuous with FGG as $Q\to0$, rising with $Q$) | Proposition-N |
| R68 | **Proposition G5**: stabilised seeds at the frontier convert in $e^B/\omega$: $5\times10^{-19}$ s, 3 ns, 0.09 s at $R=100,200,300\,\ell_P$ | Proposition |
| R69 | **Proposition G6/G7 (Route G)**: self-amplifying NEC violation meets G1 automatically (budget ratio $6(L^2/h)^{1/3}$); a thin-edged seed in vacuum becomes a visible, parent-consuming bubble (interface energy turns negative at $HR=\sqrt{2GM/R}$); a hidden daughter needs a throat-opening NEC-violating boundary layer | Proposition (+ numerics) |
| R70 | **Problem B answer**: no more accessible route within GR + QFT + GSL (closure by G1, G2, S7, Wall, SCC); the only laboratory-scale class is self-amplifying NEC violation, which requires GSL and positivity violation; safety floor $\ge m_Pc^2\approx2$ GJ per hidden daughter | Assessment built on the theorems above |

## Rounds 5 and 7 (ROADMAP_ROUNDS5_7_12.md)

| ID | result | status |
|---|---|---|
| R71 | **Theory F1**: a tilted-quartic scalar derives $\rho,\sigma,w,m$ and the CdL exponent; the chain couplings → $B_{\rm lab}$ closes (gravitational bounce agrees with thin wall to 3–17%) | computation |
| R72 | **Proposition F1**: $B_{\rm CdL}\simeq6745/\epsilon^3\lambda$; the Round 2 frontier ($B_{\rm CdL}=10$) needs $\lambda\sim10^4$ (one loop ≫ tree); with $\lambda\le4\pi$, $\epsilon\le\frac12$: $B_{\rm lab}\ge1095$, $P\le10^{-476}$ | Proposition |
| R73 | **Round 7 radiation**: explicit 1+1 field, laboratory state: collapse branch has a negative precursor ($\approx-0.44F_H$), switch-on over $\approx3.5\times4GM$, deficit $\approx2F_H/\kappa$; late flux Hawking to 0.02% (check); the daughter branch (state D1, assumed) differs by $\approx0.1$ quanta per event | Proposition-N |
| R74 | **Inaccessibility**: the daughter's inflation and, at late times, which branch occurred cannot be established by exterior observers; only event-level statistics carry the branch | Theorem (Birkhoff, G1, G2) + numerics |

## Round 6b (ROADMAP_ROUND6B_GSL.md)

| ID | result | status |
|---|---|---|
| R75 | **Correction F21**: NEC violation does not imply GSL violation; "Route G must violate the GSL" withdrawn (a minimally coupled classical NEC-violating sector conflicts with black-hole thermodynamics as a theory, but genesis need not) | correction |
| R76 | **Lemma GSL-1**: in a globally hyperbolic spacetime, a relatively compact past slice of an observer forces relatively compact past slices on every earlier Cauchy surface | Theorem |
| R77 | **Proposition GSL-2 / Theorem GSL-3 (genesis dichotomy)**: with the GSL on causal horizons of complete observers and normal, weakly gravitating laboratory data, a complete daughter observer keeps a non-compact past (the neck never closes) and cannot end in a de Sitter phase; an eternal daughter and a laboratory-determined daughter exclude each other. de Sitter parent: $\Lambda_d\le\Lambda_p$ up to $O(\hbar)$ | Theorem (spherical symmetry); Corollary 3b Proposition |
| R78 | **Theorem GSL-4 (QNEC cost)**: flipping $\theta_-$ needs $-\int(r_0/r)S''_{\rm out}d\lambda\ge A_0\lvert\theta_0\rvert/4G\hbar$; $\sim S_{\rm BH}=10^{70}$ for a 1 m sphere in flat space; discounted by $\sqrt f$ near horizons; anti-trapped floor $m>Rc^2/2G$ | Theorem |
| R79 | **Proposition GSL-5**: in 2D with the strong QNEC, $\Theta'\le-\frac6cS'^2$ (exact quantum focusing) and an anti-trapped episode lasts $<c/6\lvert\Theta_0\rvert$; the 4D s-wave model focuses for $S_{\rm BH}\ge c/12$ | Theorem |
| R80 | **Propositions GSL-6/7**: in a near-extremal throat QNEC-saturating matter anti-traps spheres at depth $\rho$ iff $T\sinh(\rho/Q)<(c\eta/12)T_{\rm brk}$; verified analytically, by 4D s-wave integration and in the exact outgoing Bonnor–Vaidya–RN metric (flip iff $\delta m>rf/2$, the G1 budget) | Proposition + Proposition-N + exact metric |
| R81 | **Proposition GSL-8**: product solutions exist iff $\Lambda Q^2\le1/4$; a T-region grows iff $C>C_N$; extremal-throat bud $3/16$ (injected) or $0.2212$ (converted); a static throat in a fixed vacuum never grows; universal $M\simeq(0.07\text{–}0.16)m_P^3/\sqrt\rho$ | Theorem (criterion), exact values |
| R82 | **Proposition GSL-9**: the Cauchy-horizon route survives near extremality inside $[1/\sqrt L(1+s),1/\sqrt L(1-s)]$ (upper edge exact); width $\le r_+-r_-=2\kappa_-r_-^2$ | Proposition-N + exact edge |
| R83 | **Proposition GSL-10**: energy $\delta M$ kept in an extremal throat hides every region deeper than $\approx Q(1-\frac12\ln(2\delta M/Q))$; a growing bud deeper than $\approx2Q$ is hidden | Proposition (quasi-static) |
| R84 | **Remark GSL-11**: the Standard-Model electroweak corona stays below decompactification ($\Lambda Q^2\le0.075$; $\le0.226$ in an overestimate) | Assessment |
| R85 | **Problem B, GSL form**: the GSL forbids exactly eternal, data-determined (and disconnected) daughters; transient connected daughters and Cauchy-horizon daughters are GSL-compatible; the remaining obstacle is $\rho\gtrsim m_P^6/M^2$, common to all routes | Assessment built on R76–R83 |

## Round 6c (ROADMAP_ROUND6C.md)

| ID | result | status |
|---|---|---|
| R86 | **Theorem 6c-1 / Corollaries 6c-2 to 6c-5 (extremal throat)**: normal false vacuum is sub-horizon ($L_{\max}/2\ell_H\le0.78$); no throat discount for growth; the extremal throat can hide or grow, not both | Theorem + Proposition-N |
| R87 | **Proposition 6c-9 (Nariai throat)**: a near-Nariai SdS hole in a false-vacuum layer, tipped by outgoing negative energy, grows a hidden inflating daughter ($r>10^4r_N$, persistent trapped neck); converged; robust to wall position, tuning, Bunch–Davies back-reaction up to $c=150$ | Proposition-N |
| R88 | **Proposition 6c-12 (universal tipping law)**: all 25 runs obey "daughter iff $W=-\int Ne^{-u/\tau}du>W_{\rm th}$", $\tau=3.2r_N$, $W_{\rm th}=1.85\times10^{-3}r_N/G$ (local frame) | Proposition-N |
| R89 | **Invariant cost (correction F26)**: Misner–Sharp jump $\sqrt{f_0}W\approx(0.3\text{–}2)\times10^{-2}\sqrt\delta M_N$ ($(2.5\text{–}5)\times10^{-7}M_N$ at $\delta=10^{-8}$); $P_{\rm Hawking}=7.44\times10^{-5}/M^2$ reproduced; local energy $0.55\%M_N$, 270× below the flat-space flip | exact relation + Proposition-N |
| R90 | **Proposition 6c-10 / Corollary 6c-11 (P-A1)**: 2D CFT weighted ANEC $\int Te^{-u/\tau}\ge-c/48\pi\tau$ (sharp); tipping by a CFT needs $c\ge0.89(r_N/\ell_P)^2$ (species length = throat radius); vacuum polarisation is stabilising | Theorem (2D) + Proposition-N |
| R91 | **Proposition 6c-13**: with $T_{uu}\ge0$ the minimal anti-trapped radius never decreases, in any parent (de Sitter included) | Theorem (spherical) |
| R92 | **Proposition 6c-14 (ghost condensate)**: tipping needs $\lvert\rho+p\rvert\ge4.6\times10^{-3}\rho_f$ for $\gtrsim\tau$, 145–430× below the CLNS stability ceiling | Proposition (conditional on the EFT) |
| R93 | **Proposition 6c-15 (preparation)**: normal data with a regular past suffice; the hole must form inside the false vacuum (inner wall at $u\ge3u_{\rm hole}$); partner vacuum $\rho_f=1/72\pi G^3M^2$, $\rho_f^{1/4}=0.329$ GeV $(M_\odot/M)^{1/2}$ (corrected in Round 6j from $1/24\pi$, 0.433 GeV; F56) | Proposition + Proposition-N |
| R94 | **Propositions 6c-6/7/8 (P-R6g, P-R6h)**: LLL reduction $c=\lvert q\rvert N_f$; conformal vacua saturate the strong QNEC on curved 2D backgrounds (anomaly included); exact Schwarzian moments | Theorem |
| R95 | **Theorem 6c-R6i**: GSL-3 without spherical symmetry, by comparison at the touching point of the smallest enclosing normal leaf, given the support form of the local GSL | Theorem (conditional on (S)) |
| R96 | **Theorem 6c-S3 (P-S3 solved)**: $c(E)$ in closed form; $c(0^+)=-13/18$ exactly; $c\simeq-(16/\pi)E^{-3/2}$ | Theorem |
| R97 | **Proposition 6c-16 (P-R7)**: Regge–Wheeler greybody factors; the 3+1 branch difference is IR-finite (0.005–0.006 quanta per event per scalar) | Proposition-N |
| R98 | **Proposition 6c-17**: near threshold no exterior measurement distinguishes daughter from no daughter | Theorem + Proposition-N |

## Round 6d (ROADMAP_ROUND6D.md)

| ID | result | status |
|---|---|---|
| R99 | **Independent verification**: a method-of-lines solver with its own data integrator (second order; Reissner–Nordström mass to 1.3e-6) reproduces the Nariai daughter and its hiding; threshold W 0.053–0.071 vs 0.031–0.035 (agreement to a factor of 2) | Proposition-N |
| R100 | **Correction F32 / bridge-free data**: the Round 6c data had a trapped exterior (beyond a bridge); a light wall gives normal, asymptotically connected data, the exterior collapses into a parent black hole, and the daughter grows inside it (W = 0.142, δ = 1e-4) | Proposition-N |
| R101 | **Proposition 6d-1**: a displaced ghost condensate in the Nariai region relaxes as P′ = ε₀e^{−t/r_N}: its NEC phase switches itself off on the dS₂ Hubble time | Theorem (EFT) |
| R102 | **Proposition 6d-2**: ghost-condensate stability window 1.9e-4/√(αP″) ≤ M³r_N ≤ 1.7e-2/√α (86× wide); mode integration on dS₂×S² gives growth ≤ 3.3× over 3τ | Proposition + Proposition-N |
| R103 | **Proposition 6d-3**: the strong QNEC alone gives ∫T e^{−u/τ} ≥ −c/48πτ for every 2D CFT state | Theorem |
| R104 | **Proposition-N 6d-4**: the same fields' back-reaction raises the threshold faster than c/48πτ (no crossing to c = 5r_N²); quantum fields cannot push, with or without the species bound | Proposition-N |
| R105 | **Proposition 6d-5 (complete route)**: from data with no black hole (false-vacuum ball, light wall, imploding ordinary-radiation shell; normal, bridge-free, NEC except the push), a push gives a hidden daughter inflating to r ≈ 4,500 inside the parent's black hole; push injected inside the ball (no discount in this realisation) | Proposition-N |
| R106 | **Proposition 6d-6 (generalised entropy of the solution)**: the parent's event-horizon area never decreases; daughter-observer horizons carry a deficit ≤ 3.3√δ S_N (0.6–2.4 nats in the runs), the quantity the ghost-condensate GSL dispute must accommodate | Proposition-N + bound |

## Round 6e (ROADMAP_ROUND6E.md)

| ID | result | status |
|---|---|---|
| R107 | **Dynamical ghost condensate**: the IR condensate (geodesic flow, shift charge, clock) as a back-reacting source in the double-null solver. It reproduces the solver bit for bit at zero charge, constraints converge at second order (5.9e-4 → 1.6e-4 → 4.2e-5), and the potential force is implicit (F42) | verified code |
| R108 | **Theorem 6e-1 (shift charge)**: ∇(nu) = αM²□θ − V′. Without explicit breaking, each element's NEC-phase charge is conserved, and α only renormalises G (G_eff = G/(1 − 12πGαM²)). A charge-free condensate cannot make the pulse; a soft tilt V′ can | Theorem (EFT) |
| R109 | **Proposition-N 6e-2**: the condensate's own stress tensor tips the near-Nariai band, with no prescribed flux. Threshold \|n\|r_N² = 1.94, 1.65, ≈1.47 ×10⁻³ (h = 0.08, 0.04, 0.02); same weighted-deficit law (W ≈ 0.09–0.12); hidden; charge × width const for compact lumps; ≈50× below the CLNS ceiling | Proposition-N |
| R110 | **Proposition-N 6e-3 (self-generated)**: a tilt confined to the false vacuum makes the NEC phase from zero charge. λ_th r_N³ = 4.7 and 3.8 ×10⁻⁴ (h = 0.08, 0.04); only the energy-conserving pair (charge plus vacuum drift) satisfies the constraints (2e-6 → 2e-7, against 0.04–0.37 for either half) | Proposition-N |
| R111 | **Proposition-N 6e-4 (stability on the measured history)**: growth < e over 3τ for all (k, l). Window ≈ 1×10⁻³–1.7×10⁻² in M³r_N for an initial charge; a tight window, ≈ 6×10⁻³–1.7×10⁻², for the self-generated case | Proposition-N |
| R112 | **Proposition-N 6e-5 (complete route, self-made push)**: from data without a black hole, the tilted condensate (zero initial charge) gives a hidden daughter above λ_th = (3.5–3.8)×10⁻⁷ (λ r_N³ ≈ 1.9×10⁻³). u_EH = 0, no daughter row reaches the outer edge, parent area non-decreasing; the in-ball penalty F36 no longer applies | Proposition-N |
| R113 | **Theorem 6e-6**: with normal data and classical sources, every eternal daughter violates the local GSL on a complete observer's horizon (its ingoing generator starts on the data with decreasing area). Measured 7–164 nats (S_N = 942), whatever supplies the NEC violation | Theorem (spherical) |
| R114 | **Proposition-N 6e-7 (transient daughter)**: the tilt also destabilises the false vacuum at T_exit = ρ_f(16/(3ε) − 1)/λ_f. Simulated exit at T = 232 (predicted 221), then recollapse inside the parent (r_max 9,905 → 34); the unseeded control keeps growing. GSL-consistent for complete observers; observers who end at the crunch are open (P-R6e-b) | Proposition-N |
| R115 | **Theorem 6e-8**: the route's data are normal, sub-horizon (wall at 0.595 l_H), NEC-satisfying and bridge-free, so no Farhi–Guth or Penrose obstruction to making them applies | Theorem |
| R116 | **Proposition 6e-9**: a recollapsing closed ball of radiation plus vacuum turns around with ρ_rad ≥ ρ_f at a < l_H/√2, so confined radiation cannot empty the ball | Theorem (homogeneous model) |
| R117 | **Proposition-N 6e-10 (creation by supercooling)**: a momentarily static hot ball of radius 0.6 l_H with ϕ in the true vacuum is normal for α_PT ≥ 2.5. Heating (T > T₁ = 1.24 T_c) restores the false vacuum, and free-streaming cooling keeps it there with zero kinetic energy (adiabaticity 0.004). Released radiation 2.9–5.3, against a shell need of ≈5.7. Open: one GR simulation of the emptying, handed off to the route | **Withdrawn (Round 6f, R119):** the high-T free energy was invalid at m/T ≈ 15 and the emptying was not a GR evolution |
| R118 | **Proposition-N 6f-1 (restoration threshold)**: heating the true vacuum restores the false vacuum only above E_min = 4.70, 2.88, 2.29, 2.00, 1.81 ρ_f for ε = 1–5 (f_H = 0.999, x_c = 30); weakly friction dependent (1.90–2.25 for Γ/m = 0.03–3); lowest found 1.74 (ε = 5.2). Below it, field and plasma stall | Proposition-N (homogeneous) |
| R119 | **Proposition-N 6f-2 (no persistent normal ball from heating)**: the largest cone-normal hot ball at E ≳ E_min has R ≈ 17.5–18.25 against the route's 17.85; in continuous GR evolutions (plasma + field + gravity, regular centre) the false-vacuum region shrinks or is trapped before the radiation leaves. Proposition-N 6e-10 withdrawn; the ball has no known creation mechanism | Proposition-N |
| R120 | **Formation solver verified**: second order against the homogeneous ODE (ϕ 5.7e-3 → 1.5e-3 → 4.3e-4 v); constraints first order in the plasma sector (θ = 0.6 box), 8.2e-5 → 1.3e-5 on the centre stage; vacuum axis |s| ≈ 3e-11; centre-to-corner hand-off at h = 0.04: largest r difference 1.8e-5 (medians 1e-10–1e-7), not formable at h = 0.02 | Proposition-N (verification) |
| R121 | **Theorem 6f-3 (observers ending at singularities)**: their causal-past boundary is the past cone of the endpoint, whose generators all end there; its area falls to zero, so no generalised entropy containing A/4G is non-decreasing. Generic to every collapse. Measured: 1001 (plain black hole), 966 (near-Nariai collapse), 1245 (transient daughter) nats, S_N = 942. Closes P-R6e-b: the transient daughter is GSL-consistent for every complete observer | Theorem (+ Proposition-N) |
| R122 | **Proposition-N 6f-4 (nonlinear condensate EFT)**: the exact linearisation (A = P′ + 2XP″ = 2 + 3ε, friction 1 + A′/A) equals the nonlinear 1+1 solver to 3e-4; at the window edges no caustic up to seed 1e-4 (δX ≈ 7\|ε\|), caustics only at 1e-3 (quantum seed ≈ 3e-7). Correction: Round 6d's mode equation used A = 2 + ε without A′/A, an O(ε) (< 1%) difference in the window | Proposition-N |
| R123 | **Proposition-N 6f-5 (non-spherical, test-field proxies)**: scalar and Regge–Wheeler (mass without the local vacuum energy) l-modes, l = 0–20, on five backgrounds: away from r → 0 (r > r_N/2) amplification ≤ 3.6, equal between tipping backgrounds and no-tip controls (to 0.3%); the large raw values sit only at the singularities, as in the controls | Proposition-N (proxies) |
| R124 | **Proposition-N 6f-6 (full linear stability range)**: WKB growth exponent over redshifting k and non-redshifting l modes on the measured histories; universal in r_N (1e6, 1e20, 1e40: identical windows, G to 0.3%), edges ∝ α^(−1/2). α = 1: initial charge 3.8e-4 ≤ M³r_N ≤ 3.2e-2 (×83), self-generated 2.6e-3–3.2e-2 (×12). Supersedes R111's analytic upper edge and its "factor 3"; corrects the Jeans rate of Proposition 6d-2 to 4π√α M³ (code normalisation) | Proposition-N (EFT) |
| R125 | **Theorem 6f-7 (accessibility invariant)**: any near-Nariai throat of a false vacuum needs E√ρ ≥ (1/3)(8π)^(−1/2) E_P √ρ_P ≈ 8.9×10⁶⁴ J^{3/2} m^{−3/2}, independent of the vacuum scale (the black-hole formation invariant) | Theorem |
| R126 | **Proposition-N 6f-8 (two solvers, one route)**: the route evolved by `condensate.evolve` and by the independent formation solver agree (both form the daughter; r to 4.6e-6 median, ϕ to 6e-3 v, charge to 6e-8 at the 95th percentile); one evolution from the ball (zero charge) to the hidden inflating daughter with all fields | Proposition-N |
| R127 | **Proposition-N 6f-9 (exit after the throat)**: the spinodal T_exit = ρ_f(16/(3ε) − 1)/λ_f must come after throat formation (T ≈ 125), i.e. ρ_f(16/(3ε) − 1) > λ_f T_throat with λ_f ≥ λ_th: marginal at ε = 4 (seeded route exits early, no daughter), satisfied at ε = 2 (normal retuned data, hidden daughter). The post-spinodal roll is unresolved on the corner grid (stiffness ~4e3, F45) | Proposition-N |
| R128 | **Proposition 6f-10 (problem 4)**: with the condensate bound M ≲ 100 GeV (10 MeV linear), the route tips only throats with r_N ≳ 1.1×10¹³ m, M ≳ 2.5×10⁹ M_☉ (2.5×10²¹ M_☉), vacuum scale ≈ 6.6 keV (6.6 meV); laboratory masses (1–10⁶ kg) need M ≳ 2×10¹³–2×10¹⁵ GeV, excluded by 11–13 orders. No ingredient (condensate, false vacuum, tilt) is observed | Proposition (given the literature bounds) + Assessment |
| R129 | **Proposition 6f-11 (problem 5)**: the best human capability (a 50 Mt explosion) falls 10^39.1 short of the invariant (NIF 10^49.7, LHC 10^50.5); neutron-star mergers exceed it. Heating cannot create the ball at any energy (R119), and the quantum route is bounded by P ≤ 10⁻²³ (R42): not technologically accessible | Proposition + Assessment |
| R130 | **FRW on a light cone (verification)**: radiation-dominated FRW data in conformal null coordinates (s = ln a) and in the affine cone gauge (s = 0, comoving fluid at η_f = ln a) are reproduced by the centre-stage solver: r = aχ to 6e-3 → 2e-3 and s = ln a to 1.4e-3 → 3e-4 (h = 0.2 → 0.05); Misner–Sharp mass = (4π/3)ρr³ to 9e-3 → 4e-3 (affine) | Proposition-N (verification) |
| R131 | **Proposition-N 6g-1 (critical size of a supercooled patch)**: a false-vacuum patch left by a first-order transition in a radiation universe (uniform density, Hubble flow) collapses into a black hole with no daughter below, and inflates above, a critical radius of 0.57–0.59 l_H (f = 0.9; identical at h = 0.1 and 0.05), 0.42–0.57 l_H (f = 0.5), 0.55–0.65 l_H (f = 1). The route's ball (0.595 l_H) lies in the critical band: the universe makes it as its near-critical patches. Constraints 1e-6–3e-4 | Proposition-N |
| R132 | **Proposition-N 6g-3 (cosmic daughter hidden, GSL holds)**: f = 0.5, patch 0.57 l_H: the false vacuum inflates while a trapped true-vacuum region forms outside the wall (r ≈ 38) from u ≈ 30; on 188 outgoing rays the daughter lies behind it. Parent event horizon at u = 29.85 with non-decreasing area; daughter observers' horizons never lose area. Theorem 6e-6 does not apply (its premise, normal data, fails for cosmological data) | Proposition-N |
| R133 | **Proposition-N 6g-4 (the condensate widens the window)**: ε = 2, patch 0.57 l_H: collapses without the tilt, inflates with λ_f = 1e-6; 0.53 l_H still collapses; with ε = 4 the tilt does not tip a 0.57 l_H patch | Proposition-N |
| R134 | **Proposition 6g-2 (external assembly; corrected in Round 6h from "Theorem")**: building a region of compactness k ≈ 2/3 by delivering its energy from outside within αR/c needs P ≥ (k/2α)c⁵/G (1.2×10⁵² W at α = 1, 8.1×10⁵¹ W at α = 1.5); holding it static needs a stress ≈ (k²/4)c⁴/G = 1.3×10⁴³ N (below the conjectured c⁴/4G). A conditional engineering bound on two strategies, not a proof that no machine can exist; the maximum force and luminosity are conjectures with counterexamples in GR | Proposition (conditional) |
| R135 | **Proposition 6g-5 (every machine)**: in quantum form E_total·E_quantum ≥ ħc(E√ρ)^{2/3} = 6.3×10¹⁷ J². LHC beam (ideal focus) 10^15.2 short, FCC-hh 10^13.0, NIF 10^29.7, 50 Mt 10^16.3; an invented machine with one year of world energy in 100 TeV quanta, ideally focused, is 10^1.8 short and would have to deliver 9×10⁴⁹ W; a single-collision collider needs 3.5×10¹⁸ GeV per particle (a 1.1 pc linac at 100 GV/m). Loopholes (Assessment 6g-6): FGG tunnelling (P ≤ 1e-23 even at Planck energies), TeV-scale gravity (excluded below 10.1 TeV by CMS; would need the route rebuilt in higher dimensions), analog systems (no spacetime). Verdict (corrected in Round 6h): no accessible mechanism has been found; not proved impossible for every process | Proposition + Assessment |
| R136 | **Theorem 6h-1 (Penrose / Farhi–Guth for machines)**: in a globally hyperbolic spacetime with non-compact Cauchy surfaces obeying the NEC, a closed anti-trapped surface implies past null incompleteness; an inflating daughter contains one. So a daughter needs NEC violation, a past singularity (cosmology) or a quantum transition: no machine of NEC-respecting matter acting on a universe with a regular past can make one, by any strategy (rapid, slow, dynamical, triggered). Independent of assembly time and of any maximum force or luminosity | Theorem |
| R137 | **Proposition-N 6h-2 (a trigger in an anti-trapped background)**: in the cosmological background a positive-energy trigger works: a supercooled patch of 0.60 η₀ collapses, enlarged to 0.62 η₀ (trigger ≈ 10% of its vacuum energy) it inflates into a hidden daughter. Critical radius 0.6175 η₀ = 0.586 l_H (h = 0.1); closer to it the classifier is limited by the computed domain. In a normal background the trigger must violate the NEC (R136, R88) | Proposition-N |
| R138 | **Assessment 6h-3 (routes B–D)**: black-hole-assisted creation (Oshita–Yokoyama) needs a black hole of about the Nariai mass class at the end of its evaporation and is quantum-suppressed; dynamical assembly with the condensate needs the near-Nariai ball, mass r_N c²/3G (geometric, any history), ≥ 2.5×10⁹ M_☉ for the observationally allowed condensate; TeV gravity is excluded below 10.1 TeV. Smaller backgrounds need a false-vacuum-only condensate (conjecture, open) and GUT-scale quanta arriving within r_N/c (1 g, 6×10⁹ quanta of 9×10¹³ GeV in 7×10⁻³⁹ s at 1.5×10¹⁶ GeV) | Assessment |


### Round 6i

| # | result | label |
|---|---|---|
| R139 | **Gauss–Bonnet stress in spherical symmetry**: Γ_ab for ds² = −e^{2s}du dv + r²dΩ² with F = λ²f(φ), verified by the off-shell identity ∇^aΓ_ab = −(GB/2)∇_bF; −Γ_uu = 4(2m/r³)∇_u∇_uF + (16e^{−2s}/r²)(∇_u∇_ur)(F_u r_v + F_v r_u) | Theorem (symbolic) |
| R140 | **Proposition 6i-1 (curvature activation pushes toward collapse)**: on a vacuum-energy throat the tipping flux is N = r_N²φ_u² + 2F_uu; the Gauss–Bonnet part is a total derivative, W_GB = −(2/τ²)∫(F − F₀)e^{−u/τ}du + boundary terms of the same sign, negative whenever F grows. With GB = C² + 8Λ²/3 > 0, a field the coupling destabilises rolls to increase F: activation can only push toward the black hole. The opposite arrangement needs a field that is tachyonic at low curvature (cannot be stored) | Proposition |
| R141 | **Proposition-N 6i-2 (nonlinear homogeneous throat)**: Kantowski–Sachs with the coupling (√−g GB = d/dt[8ȧ(1 + ḃ²)] verified), exact Nariai data: scalarisation (λ² = 0.02–2 r_N²) ends in a black hole or a degenerate kinetic matrix (loss of predictivity at λ²|f′| ≈ b²/2 or through HḞ); the reversed sign with a positive-energy minimum: 0 daughters in 144 runs (66 black holes, 78 degenerate); unbounded-potential "daughters" discarded (V < 0 first) | Proposition-N |
| R142 | **Assessment 6i (activation) and Proposition 6i-3 (curvature-made vacuum)**: activation needs λ²|GB| ~ 1/L²; with λ ≲ 0.28 km laboratories fall short by 19.0–19.5 orders, neutron stars and stellar black holes by 1.5. A vacuum lifted by curvature has l_f ≥ l_curv (|Δf| ≤ 1/2, λ ≤ l_curv): only around an object already horizon-class at that scale | Assessment + Proposition |
| R143 | **Proposition-N 6i-5 (monopole core, thin shell)**: static stable cores below η_c = 0.25–0.26 (λ = 1, e = 0.5), 0.32–0.33 (e = 1), 0.39–0.40 m_P (λ = 0.1). At the monopole's own mass B ≈ 650/λ (η-independent, 1.54 × the flat CdL 421/λ), but M ≪ Q: no bridge, the β_out < 0 branch is disconnected from us (Gauss's law forbids it). Connected channel: M > Q ≈ 7 m_P (a magnetic black hole), B = 3.1×10³–2.8×10⁴. Added mass up to 10× and charged horizons leave the barrier in place; negative dust lowers B by ≤ 20% | Proposition-N + Assessment |
| R144 | **Proposition-N 6i-4 (breathing bubble)**: wells need k_w m_d < 0; the hidden-daughter case needs k_w > 0 and negative-energy (NEC-violating) dust, the other sign eats the parent (V even in m_s, β_out odd). All 713 wells in 6000 random shells contain the zero of the net surface energy, crossed at light speed. B = B̂(l_f/l_P)², B̂ ≈ 0.5–0.9 (→ 0 only at merger points); total mass ~10⁻³ l_f c²/G | Proposition-N |
| R145 | **Proposition 6i-8 (gauged ghost window empty)**: every NEC-violating condensate region is charged (ρ + p = 2XP′M⁴, J⁰ = 2P′φ̇); gauged, Coulomb/|ρ+p| = (P″|δ|/9)g²M²R². Tipping (M³r_N ≥ K_lo), no screening (r_N ≲ c/gM) and CLMT (M ≲ g²10¹⁵ GeV) need M ≥ 10^17.0–10^17.6 GeV > 10¹² GeV | Proposition |
| R146 | **Proposition 6i-9 (vacuum-dependent condensate, conditional)**: M(χ) large inside, M_out ≤ 100 GeV outside; EFT across the wall needs r_N ≥ 1/M_out = 2×10⁻¹⁸ m: M_N ≥ 8.9×10⁸ kg (vacuum 1.6×10¹⁰ GeV, interior scale 1.4×10¹² GeV), 30 orders below the universal floor. An overall factor F(χ) → 0 instead drives the strong-coupling scale MF² to zero. Assembly still needs 5×10³³ quanta of 100 GeV within 6.6×10⁻²⁷ s | Proposition (conditional) |
| R147 | **Proposition 6i-6 (resonant tunnelling)**: two-well shell potentials exist (charge + gas components; 13 clean in 20000); on resonance the exponent is max(B₂, B₁ − B₂) ≥ B₂ instead of B₁ + B₂: the barrier next to the inflating branch is never removed | Proposition |
| R148 | **Proposition 6i-10 (collision pockets)**: DTR allows a pocket of any Λ_C (planar, Johnson–Yang: ≤ parent); where f_out > 0 the sign of β_out is fixed, so a pocket launched in a normal region never reaches the hidden branch and collapses; in an inflating parent the cosmological route applies | Proposition |
| R149 | **Proposition 6i-7 (four-form membrane)**: a Brown–Teitelboim membrane is a Farhi–Guth wall (NEC saturated); hidden inflating branch for M < M_cr, classically from the white-hole singularity; B(M → 0) = 1.02 × flat CdL ("no free lunch" holds) | Proposition |
| R150 | **Theorem 6i-11 (bridge invariant)**: every hidden daughter contains a sphere with 2Gm/(c²r) ≥ 1 (the area radius has a minimum between the wall, β_out < 0, and our asymptotic region); from regular data its first appearance is a black-hole formation event with E√ρ̄ = √(3/32π)c⁶/G^{3/2} = 2.3×10⁶⁵ J^{3/2}m^{−3/2} (exact, scale-free, any route, any NEC-violating ingredient), 10^39.5 beyond technology. Side-change lemma: β_out changes sign only at m_s = 0 or f_out ≤ 0 (3351/3351). Escapes: quantum, a pre-existing bridge (natural black holes), charged Cauchy horizons (no help), anti-trapped backgrounds | Theorem (spherical) + Assessment |
| R151 | **Q* (Conjecture-level door)**: a barely metastable false vacuum with a light wall has R₃ = R_CdL exactly and a semiclassical exponent ~B_CdL; at 500 GeV, B_CdL = 1: 3.5×10⁻¹⁹ m, 2.9 TeV, but a bridge of 4.7×10⁻¹⁶ l_P (bridge/l_P ≈ 9.4 (scale/M_Pl) B^{3/4}): outside R42's controlled regime and semiclassical gravity; needs an unknown TeV false vacuum | Conjecture + Assessment |


### Round 6j

| # | result | label |
|---|---|---|
| R152 | **Compactness identity**: E√ρ̄ = √(3/32π)c⁶/G^{3/2}k^{3/2}, k = 2GE/c⁴R; Theorem 6i-11 is its value at k = 1 (a marginally trapped sphere): a statement about geometry | Theorem (identity) |
| R153 | **Misner–Sharp / Kodama identity**: ∂_u m = r²e^{−2s}(G_uv r_u − G_uu r_v) (and u↔v); orthonormal form dm/dr = 4πr²(ρ + F r_t/r_n). The density in Theorem 6i-11 is the areal-volume mean of the zero-flux-frame energy density | Theorem (symbolic) |
| R154 | **Proposition 6j-1 (which density)**: at the first marginal sphere some region has density ≥ ρ̄ = 3c⁴/8πGR²; non-uniform profiles raise the peak (3.3× Gaussian, 2.9× shell); negative components raise the positive energy needed (∫ρ₊ = Rc⁴/2G + ∫ρ₋) | Proposition |
| R155 | **Proposition 6j-2 (three bounds)**: geometric (exact), stored (≥ E_MS for ρ ≥ 0, proper ≥ areal volume), supplied (E_stored − E_natural; not bounded by the invariant) | Proposition |
| R156 | **Proposition 6j-3 (small supplied energy)**: needs a natural precursor within δ = E_sup/Mc² of threshold that is stable (near-Nariai states decide in 3.2 r_N/c), older than 13.8 Gyr (M > 1.7×10¹¹ kg) and not pushed over naturally: δ ≤ 6.4×10⁻⁹ for 10²⁰ J. A machine cannot supply the NEC-violating trigger: Ford–Roman |E₋| ≲ ħc/R is short of 1.85×10⁻³Rc⁴/G by (R/l_P)²/540 | Proposition + Assessment |
| R157 | **Penrose form (non-spherical)**: Riemannian Penrose inequality (Huisken–Ilmanen; Bray): time-symmetric, ρ ≥ 0 ⇒ k_A ≥ 1, the same invariant with the area radius; dynamical case open | Theorem (cited) |
| R158 | **Proposition 6j-4 (field-dependent G)**: a core with F = (M_P,eff/M_P)² gains F^{3/4} but needs Δχ = √(3/16π)\|ln F\|M_P, gradient energy ≥ 2πLΔχ² (capacitance), compactness ≥ (3/4)ln²F: horizon-class for F < 0.32 (F ≈ 10⁻⁵³ needed: compactness 1.1×10⁴). f(R) included | Proposition |
| R159 | **Proposition 6j-5 (finite pocket)**: a classical pocket is sub-Hubble and collapses within ~R/c, expanding by ~R/l_f (3×10⁻¹⁶ for a TeV vacuum at 10⁻¹⁸ m) | Proposition |
| R160 | **Proposition 6j-6 (accretion)**: a daughter-ready (near-Nariai partner) vacuum needs infalling density c⁶/72πG³M² at the horizon: Eddington flows are 10^18.8–10^19.7 short, merger debris supplies 0.25–1.8%; completion still needs the NEC escape | Proposition |
| R161 | **Proposition-N 6j-7 (compact precursor)**: the most compact bag star has 2GM/c²R = 0.54, R = 0.28 l_f, vacuum share 14%: sub-Hubble, collapses to a black hole. The QCD bag constant (145 MeV) is the Nariai partner of a 5.2 M_☉ hole | Proposition-N |
| R162 | **Proposition-N 6j-8 (Oshita–Yokoyama, both steps)**: Σ = Σ̂(l_f/l_P)², Σ̂ = 0.013 (k_w = 0.5) to ~4×10⁻⁴ (light walls); GUT vacuum Σ ≳ 2×10⁷; daughters from evaporating holes only for l_f ≲ 10²–10³ l_P | Proposition-N |
| R163 | **Proposition 6j-9 (throat-area law; analytic for a₀ ≪ Λ^{−1/2}: ∫₁^∞(1 − √(y²−1)/y)dy = π/2 − 1)**: Euclidean axion wormholes with Λ: S(q) → −12π²/Λ (no-boundary) as the throat pinches off, and ΔS = S(q) − S(0) = 3π³a₀² (8πG = 1) independent of Λ to 10⁻⁵, i.e. ≈ 1.18 × 4πa₀²/4G: a connected channel pays action in proportion to its throat area; throat-free = disconnected | Theorem (small throats) + Proposition-N |
| R164 | **Proposition 6j-10 (FMM interiors)**: Schwarzschild meets de Sitter (Hubble l) across a spacelike surface at r₀ = (2Ml²)^{1/3} (K_θθ continuous, M > l/2), Kretschmann 12/l⁴; the de Sitter region starts at r₀/l = (2M/l)^{1/3} > 1 and inflates. If quantum gravity does this, every black hole holds a daughter | Proposition (conditional on the FMM hypothesis) |
| R165 | **Proposition 6j-11 (extra dimensions)**: with n ≥ 4 and M_D ≈ 6–20 TeV (allowed), black holes form at a few M_D: a 100 TeV collider makes them. But R_S M_D ≈ 1.0 at 5 M_D (n = 4–7): fundamental-scale objects; a semiclassical interior (R_S = 3/M_D) needs 960 M_D (n = 4) to 10⁴ M_D (n = 6) | Proposition (conditional) |
| R166 | **Correction F56**: the Round 6c partner vacuum is ρ_f = 1/(72πM²) (Nariai 9M²Λ = 1), not 1/(24πM²): ρ_f^{1/4} = 0.329 GeV (M_☉/M)^{1/2}; GUT partner 2.2 g, electroweak 3.6 Earth masses | Correction |


### Round 6k

| # | result | label |
|---|---|---|
| R167 | **ECSK bounce**: an unpolarised spin-½ fluid (ε_eff = ε − κs²/4, s² = ħ²c²n²/8) of particle mass m bounces at ρ_b = 4m²c⁴/(πGħ²) = 3.9×10⁵⁸ kg/m³ for nucleons (10⁻³⁸ of Planck: semiclassical); a horizon forms first for M > (3/4πρ_b)^{1/2}(c²/2G)^{3/2} = 4.3×10¹⁰ kg. Every astrophysical black hole holds a bounce if ECSK holds | Proposition |
| R168 | **Proposition-N 6k-1 (ECSK interior)**: Oppenheimer–Snyder closed dust with spin term bounces at a_min ≈ (S/A)^{1/3} and re-expands beyond the horizon into a new universe; without vacuum energy it cycles; a contracting region with vacuum energy l_d < (3√3/2)A bounces and then inflates | Proposition-N |
| R169 | **Proposition 6k-2 (limiting curvature, D dimensions)**: Schwarzschild–Tangherlini meets de Sitter at r₀ = (r_s^{n+1}l²)^{1/(n+3)}, exists iff r_s > l, with r₀/l > 1 (inflates) | Proposition |
| R170 | **Proposition 6k-3 (the conditional machine)**: with P1 (n ≥ 4, M_D ≲ 15 TeV; allowed) and P2 (limiting-curvature interiors, l = 1/M_D), black holes hold daughters for M_BH > M_Y = 3.9–4.7 M_D = 24–68 TeV, below the ~71 TeV reach of a 100 TeV collider; parton σ 8–34 pb; daughters form before evaporation (lifetime/interior time ≈ 7); semiclassical control needs 10³–10⁴ M_D | Proposition (conditional) |
| R171 | **Proposition-N 6k-4 (collider test of the interior class)**: a D-dimensional static de Sitter core gives T(M) rising from zero at a remnant (μ ≈ 10.8 l^{n+1}), peaking (0.19/l) and approaching GR (5%); expanding (FMM) interiors look like GR outside | Proposition-N |
| R172 | **Natural bridge within reach**: a Planet-Nine primordial black hole (5–10 Earth masses, ~500 AU; Scholtz–Unwin 2020) is reachable by laser-launched probes in ~8–10 years (Witten 2020); under Y1 it already holds a daughter. Its Nariai partner vacuum is 60–85 GeV, near the electroweak V(0)^{1/4} = 104 GeV (partner 3.3 Earth masses): a coincidence, not a route | Assessment |


### Round 6l

| # | result | label |
|---|---|---|
| R173 | **Model U**: D = 4 + n Einstein gravity + the quasi-topological tower (Bueno–Cano–Hennigar) is one pure-gravity action containing P1 (strong gravity at M_D) and P2 (singularity resolution). Its holes are regular (de Sitter core), with two horizons above μ_ext = α^{(D−3)/2}((D−1)/(D−3))^{(D−1)/2}(D−3)/2 (Hayward), T = ((D−3)r₊² − (D−1)α)/(4πr₊³) → 0 there (remnant). ADD radius 10³–10⁶ × 1/M_D | Proposition (theory: Literature) |
| R174 | **Proposition 6l-1 (matter interiors)**: closed dust obeys ȧ² = a²ψ(C/a^{D−1}) − 1, bounces at a_b → l (any tower with a pole), accelerates for N_acc → ((D−3)/(D−1)) ln(a_max/l) + ln(2/(D−3))/(D−1) e-folds (29.6 inside a 4D stellar hole, 36.5 for 10⁹ M_sun), then recollapses and cycles. Vacuum collapse leaves a static core of radius ≈ l. No lasting inflation | Proposition |
| R175 | **Proposition 6l-2**: lasting inflation iff interior vacuum Λ > Λ_c(C) = [2((D−3)/(D−1))^{(D−1)/2}/((D−3)C)]^{2/(D−3)} in the low-curvature limit (= Hayward extremality with α → 1/Λ; numerics within 10⁻⁴ in 4D, 2–54 % higher in D = 10). Model U removes the Farhi–Guth past-singularity obstruction but not the hole-density vacuum requirement | Proposition |
| R176 | **Proposition 6l-3 (bounce lemma, any theory)**: at a bounce (∇R)² = (n·∇R)² ≥ 0, so a bounce inside a black hole needs an untrapped core inside an inner apparent horizon. **6l-4**: a homogeneous ECSK ball loses 99 % of its Misner–Sharp mass by the bounce (spin pressure work) and cannot sit in a vacuum exterior: the bounce is inhomogeneous (corrects R168's picture, F61) | Proposition |
| R177 | **Proposition-N 6l-6 (collider window)**: with NNPDF2.3 LO (sum rules 0.998, 2.003, 1.001), Yoshino–Nambu form factors and head-on trapped fractions y_n = 0.71–0.58, model-U daughter-bearing holes are produced at 100 TeV / 30 ab⁻¹ (≥ 10 events, conservative) iff M_ext ≲ 46–51 TeV, and are compatible with the LHC iff M_ext ≳ 8.7–9.7 TeV: a window c_l = l M_D ≈ 0.45–0.80 (Hayward) or 0.54–0.95 (Dymnikova) over n = 4–7, M_D = 5–15 TeV; up to ~10⁹ events | Proposition-N |
| R178 | **Formation vs evaporation**: bounce 2.7 l (0.7–2×10⁻²⁸ s) after horizon crossing; evaporation to the remnant 2×10⁻²⁷–10⁻²⁴ s (ratio 20–7×10³); remnant keeps M_ext/M; 1.5–2.5× more, softer quanta than GR; T within 10 % of GR only for r₊ ≳ 2.5–3.6 l | Proposition-N |
| R179 | **Proposition 6l-5 (second infall)**: the bounced first collapse reaches the inner horizon only at v → ∞, so later infall meets it inside the same core (10⁻²¹ l inside after a 50 l delay): it joins and modifies the existing daughter; no separate one with a Birkhoff theorem and one regular centre (Geroch) | Proposition |
| R180 | **Proposition-N 6l-6b (DTR saturation)**: at that meeting the Dray–'t Hooft–Redmount relation (kinematic) demands f → −∞ (GR mass inflation ×10¹³ after 50 l); model U bounds f ≥ 1 − r²/α, so after a delay of 7–30 l the crossing region saturates at the limiting curvature (locally de Sitter, trapped) and no thin-shell solution exists: the Barrabès–Frolov regime | Proposition-N |
| R181 | **FMM patches**: induced-metric continuity gives Δt_dS = Δv; on the expanding de Sitter branch, entries separated by Δv > l ln((r₀+l)/(r₀−l)) ≈ 2l(l/r_s)^{1/3} land in causally disjoint patches; on the contracting (FMM) branch, one closed de Sitter world | Proposition |
| R182 | **Proposition-N 6l-7/6l-8**: (closure) a compact trapped region with a regular centre returns everything to our future null infinity: no daughter persists without a remnant-like endpoint; (exterior equivalence) interiors that agree on J⁻(I⁺) give identical exterior data and Hawking out-states: experiments test the class (M_D from monojets; onset M_th/M_D; T(M); remnants), never the daughter | Proposition-N |


### Round 6m

| # | result | label |
|---|---|---|
| R183 | **Global daughter spacetime**: OS collapse in model U glued exactly to the static exterior (surface = geodesic with E = cos χ_b); maximal extension I → II → III → IV → I' computed in diamond patches matched through the finite null coordinate | Proposition |
| R184 | **Proposition 6m-1**: max Kretschmann × l⁴ is mass-independent over μ = 10⁻³–10⁸ (D = 4: 29–34; D = 5: 45–48; D = 10: 188 Hayward, 919 Dymnikova); the closed interior peaks at the de Sitter value 2D(D−1)/l⁴ at the bounce | Proposition-N |
| R185 | **Proposition 6m-2**: for μ > μ_ext the ball surface crosses r₊, r₋, bounces where f = E² > 0 and reaches the new universe I' in finite proper time (2.5–35 l); for μ < μ_ext it bounces back into the parent | Proposition |
| R186 | Both horizons regular (Kruskal metric function equal finite limits from both sides, ratio 1 ± 10⁻⁴); radial and non-radial timelike and null geodesics complete | Proposition-N |
| R187 | **Proposition 6m-3**: the daughter IV ∪ I' lies outside J⁻(I⁺) of the parent, has its own null infinity and ADM mass μ; the finite daughter (bounced closed FRW) grows in volume ×1.3×10⁵ (C = 10⁴, D = 10) to ×8×10¹¹⁴ (4D stellar) | Proposition |


### Round 6n

| # | result | label |
|---|---|---|
| R188 | **Proposition 6n-0 (model-U hydrodynamics)**: Γ² − U² = 1 − R²ψ(μ/R^{D−1}), first law, tr equation (h' × GR) and Euler give D_τU = −Γ²p'/((ρ+p)R') + Rψ − ((D−1)/2)Rψ'(y + p); GR and Friedmann limits reproduced; energy-conserving Lagrangian scheme matches the exact OS bounce to 10⁻⁵ | Proposition + numerical |
| R189 | **Proposition 6n-1**: every closed FRW universe with any perfect fluid bounces at a_b ≥ l (ψ ≤ 1/α); verified w = −0.3–1, D = 5, 10 | Proposition |
| R190 | Inhomogeneous dust: every shell bounces at √(kα) (≤ 1.8 %), inner shells first (spread ×1.6–2.1) | Proposition-N |
| R191 | **Proposition-N 6n-2/3**: with pressure, inhomogeneity and shocks, every run forms a horizon and bounces in an untrapped core with ψ ≤ 0.994/α and exact energy conservation; the bounced core meets the envelope at the inner horizon (R ≈ l); in 2/6 runs shells enter the white-hole region with a throat (Γ → 0) to the trapped envelope; beyond, the comoving slicing ends (resolution-dependent) | Proposition-N |
| R192 | Pressure threshold: central p/ρ ≳ 0.2 (D = 5, γ = 5/4) prevents trapping (no black hole) | numerical |


### Round 6o

| # | result | label |
|---|---|---|
| R193 | **One spacetime per event**: formation (trapped surface, ~10⁻²⁸ s) → partons through the regular centre (~10⁻²⁸ s) → model-U Vaidya evaporation (10⁻²⁶–10⁻²⁵ s) → zero-temperature remnant | Proposition-N |
| R194 | **Proposition-N 6o-1**: μ − μ_ext ∝ v^{−1.00±0.03}, r₊ − r₋ ∝ v^{−0.50±0.02} (24 cases): the trapped region never closes; the bounced partons stay below r₋ and reach the parent-chart boundary at finite affine parameter 2.0–2.1 l: the daughter persists through evaporation | Proposition-N |
| R195 | **Event generator** (100 TeV, n = 6, M_D = 8 TeV, M_ext = 11.5 TeV, σ ≈ 13.7 pb): model U emits 0.2–0.3 quanta per event (no quantum in 76–85 %), missing fraction 0.95–0.97; a GR hole at the same threshold emits ~2 and loses none. Near-threshold model-U production is largely invisible | Proposition-N |

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
