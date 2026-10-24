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
| R93 | **Proposition 6c-15 (preparation)**: normal data with a regular past suffice; the hole must form inside the false vacuum (inner wall at $u\ge3u_{\rm hole}$); partner vacuum $\rho_f=1/24\pi G^2M^2$, $\rho_f^{1/4}=0.433$ GeV $(M_\odot/M)^{1/2}$ | Proposition + Proposition-N |
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
