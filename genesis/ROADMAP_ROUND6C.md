# Round 6c: making the throat work, and everything that remained open

**Brief.** Six categories:

* I. the throat mechanism (P-R6f: boundary dynamics, persistence, growth, hiding, backreaction, stability);
* II. quantum throat physics (P-R6g, P-R6h, P-R6i, P-R6a, P-A1, a source for magnetic charge);
* III. the vacuum-energy obstacle (items 37–44);
* IV. foundational mathematics (P-Q9′/P-Q8, P-Q1, P-Q3, P-Q7, P-T3′, P-Q4/5/6, P-N1′, P-T5″, P-S3, P-L1, P-A2);
* V. matter dynamics and observables (P-T2, P-D2, P-D3, P-C4, P-C2, P-R7, P-O1/O2, signatures);
* VI. alternatives (P-R6b–e, Route Q, Route D).

A mid-round follow-up asked for definite answers to six questions about the throat (table below). Everything is theory. **No apparatus is
designed and no experimental procedure is proposed.** Labels as before: [Theorem], [Proposition], [Proposition-N] (numerical), [Literature],
[Assessment], [Conjecture].

**Files.**

* Code: `CODE/babyuniverse/throat_evolution.py`, `obstructions.py`, `quantum_throat.py`, `nariai.py`, `greybody.py`, and
  `crosscheck.s3_coefficient` / `crosscheck.s3_small_E_limit`.
* Scripts: `CODE/scripts/round6c_throat.py`, `round6c_nariai.py`, `round6c_answers.py`.
* Data: `DATA/round6c_throat.json`, `round6c_nariai.json`, `round6c_answers.json`.
* Figures: `FIGURES/round6c_throat.png`, `round6c_nariai.png`, `round6c_answers.png`.
* Tests: `CODE/tests/test_throat_evolution.py` (5) and `test_round6c.py` (8).

---

## The answers

| question | answer | status |
|---|---|---|
| Can a near-Nariai throat support a growing daughter? | **Yes.** It inflates to $r>10^4r_N$ behind a persistent trapped neck. The result is converged in $h$ and robust to the wall position, to the tuning, to Bunch–Davies back-reaction up to $c=150$, and to replacing the false vacuum inside the hole by true vacuum at depth $u\ge3u_{\rm hole}$. | Proposition-N |
| Is the required perturbation much smaller? | **Yes, with a correction to how it was stated.** One law governs all 25 runs: the daughter grows iff the exponentially weighted NEC deficit $W=-\int N e^{-u/\tau}du$ exceeds $W_{\rm th}=1.85\times10^{-3}r_N/G$, with $\tau=3.2r_N$. This is local energy in the corner's geodesic frame. The Round 6c checkpoint's "0.6–0.9% of $M_N$" was this local number, not a mass. Flat space needs $r/2G$ just to flip a sphere, 270 times more. The gauge-invariant (Misner–Sharp) cost is smaller still: $\sqrt{f_0}\,W\approx(0.3\text{–}2)\times10^{-2}\sqrt\delta\,M_N$, which is $(2.5\text{–}5)\times10^{-7}M_N$ at $\delta=10^{-8}$. It is bought with tuning $\delta$. | Proposition-N + exact relation |
| Is the daughter hidden behind a trapped exterior? | **Yes**, in every growing run: no untrapped true-vacuum point at large $r$, and a trapped neck that persists. | Proposition-N |
| Can QNEC-constrained quantum fields supply the perturbation? | **No, not inside the effective theory.** A 2D CFT can supply at most $W=c/48\pi\tau$ (Proposition 6c-10, sharp, from the Fewster–Hollands QEI). Tipping therefore needs $c\ge0.89(r_N/\ell_P)^2$, at which point the species length reaches the throat radius. The vacuum polarisation of the same fields is stabilising: it never tips the throat (tested at $c=30$–$1000$) and it raises the threshold. | Theorem (2D) + Proposition-N |
| Is there a physically consistent classical field that can? | **Yes, conditionally: the ghost condensate.** It needs $\lvert\rho+p\rvert\ge4.6\times10^{-3}\rho_f$ held for $\gtrsim\tau$. CLNS's stability ceiling for a softly broken ghost condensate, $\lvert\dot H\rvert\lesssim H^2$, lies 145–430 times above that. It is Lorentz-violating, and its GSL is not universal, as GSL-3 requires of any eternal daughter. No NEC-respecting field can do it in any parent (Proposition 6c-13). | Proposition (conditional on the EFT) + Theorem |
| Can the initial geometry be prepared? | **In theory, yes. No laboratory route exists.** The configuration is normal and can have a regular past: a hole formed inside false vacuum, or a vacuum lifted in situ across the horizon. Requirements: the outer wall inside the static band; $\delta\lesssim10^{-4}$; false vacuum continuing inside the horizon; and the hole tied to the vacuum by $\rho_f=1/(24\pi G^2M^2)$, i.e. $\rho_f^{1/4}=0.433\,{\rm GeV}\,(M_\odot/M)^{1/2}$. A GUT vacuum pairs with a 3.7 g black hole. The assembly process itself has not been simulated (P-R6j). | Proposition + Proposition-N |

Two further results with definite answers: **P-S3 is solved** (Theorem 6c-S3, $c(0^+)=-13/18$ exactly, tail $-(16/\pi)E^{-3/2}$), and **P-R7 is
done** (3+1 greybody factors; the branch difference is IR-finite in 3+1).

---

## Part I — the throat mechanism (P-R6f)

### I.1 The solver
`throat_evolution` integrates the spherically symmetric Einstein–Maxwell–scalar–CFT equations in double-null form:

* a magnetic charge;
* the tilted-quartic scalar of Round 5;
* $c$ conformal fields with the Polyakov stress tensor;
* optionally a classical null fluid of either sign;
* a second-order diamond scheme with three-pass predictor–corrector.

**Validation [Proposition-N].**

* Reissner–Nordström: Misner–Sharp mass conserved to $4\times10^{-6}$ at $h=0.1$ and $1\times10^{-6}$ at $h=0.05$ (second order); constraints at
  $10^{-6}$–$10^{-7}$.
* The back-reacted QNEC flip of Round 6b: flips at $c=50$, not at $c=30$, as predicted.

### I.2 What the extremal (AdS$_2$) throat cannot do
* **Theorem 6c-1 (normal false vacuum is sub-horizon).** In a region of vacuum energy with radial magnetic field, $f=f_C(r)=1-C/r+Q^2/r^2-\Lambda r^2/3$
  exactly. On normal data a false-vacuum segment starting at $f_a$ has proper length $\le L_{\max}=\int dr/\sqrt{f_C}$. The largest
  $L_{\max}/2\ell_H$ over all depths is 0.78, 0.53 and 0.41 for $\Lambda Q^2=0.3$, 0.5 and 1. Normal data never contain a super-Hubble segment.
* **Corollary 6c-2 (no throat discount for growth).** Growth needs anti-trapped false vacuum of extent $\gtrsim\ell_H$ (verified: segments with
  $L/\ell_H\ge1$–2 reach $\mathcal I^+$; shorter ones are eaten). Making such a segment from normal data costs the G1 budget at $f=O(1)$, which
  QNEC matter supplies only above the species bound.
* **Propositions 6c-3 to 6c-5.**
  * Creating the vacuum energy afterwards focuses rays by $xL/2Q$.
  * Field-restored coronae saturate at $\Lambda Q^2=1/4$ and never exceed it.
  * Growth needs the inner end outside the throat's mouth, where GSL-10 hides nothing.

The extremal throat can hide or it can grow. It cannot do both.

### I.3 The positive result: the Nariai (dS$_2$) throat
**Proposition 6c-9 [Proposition-N].** Configuration:

* a Schwarzschild–de Sitter hole with $M=M_N(1-\delta)$, $M_N=r_N/3G$, $r_N=1/\sqrt\Lambda$;
* inside a false-vacuum layer whose outer wall sits in the static band $(r_b,r_c)$;
* true vacuum outside;
* data normal.

An outgoing negative-energy flux on the ingoing ray tips it.

The mechanism:

* The flux makes the enclosed region heavier than Nariai, $m_{\rm in}=m_{\rm out}+\lvert E\rvert$.
* It also turns $r_u$ positive, so the super-Nariai Kantowski–Sachs region that results is anti-trapped and expands rather than collapses.

Outcome:

* The anti-trapped false vacuum reaches $r>10^4r_N$.
* The parent's true vacuum around it is trapped: the neck persists and the daughter is hidden.
* Positive energy instead produces the trapped (collapsing) branch, and never a daughter (tested at three amplitudes).

### I.4 The universal tipping law
**Proposition 6c-12 [Proposition-N].** The 25 runs comprise:

* Gaussian pulses at 9 amplitudes;
* 5 runs with late emission;
* negative/positive pairs at 5 separations;
* 6 top-hats of duration 10–100.

They are separated by a single criterion:
$$\text{daughter}\iff W=-\int N(u)\,e^{-u/\tau}du>W_{\rm th},\qquad\tau\in[2.8,3.9]\,r_N,\quad W_{\rm th}=(1.73\text{–}1.98)\times10^{-3}\,r_N/G .$$

* **Interpretation:** linear response of an instability. What counts is the exponentially weighted deficit, not the total. An equal positive
  pulse cancels the negative one if it follows within $\sim0.6\tau$ (dipoles at $D\le20$ fail, at $D\ge40$ succeed).
* **Scaling:** classical scale invariance makes $W_{\rm th}$ and $\tau$ proportional to $r_N$.
* **Use:** this law maps any NEC source onto the threshold.
* **Check with extended sources.** A source at rest in the corner frame gives $N=\pi r^2\lvert\rho+p\rvert$. The top-hats with
  $\lvert\rho+p\rvert=1.2\%\rho_f$ for $0.7\tau$ and $0.96\%\rho_f$ for $1.8\tau$ tip the throat. The one with $0.48\%\rho_f$ for $1.8\tau$
  does not, as the law predicts: it needs $0.46\%/(1-e^{-1.8})=0.55\%$.

### I.5 How big the perturbation is (correction F26)
* **Local energy.** $W$ is energy measured by the geodesic observer at the Nariai radius ($f'(r_N)=0$): $W_{\rm th}=0.55\%$ of $M_N$.
* **Invariant energy.** The Misner–Sharp jump follows from $m=\frac r2(1+4r_ur_ve^{-2s})$ and $\Delta r_u=(G/r)\int\lvert N\rvert du$:
  $\Delta m=2r_ve^{-2s}\int\lvert N\rvert du=\sqrt{f_0}\,W$, with $f_0=2\delta/3$. This is $4.7\times10^{-7}M_N$ at $\delta=10^{-8}$.
  Across the scan ($\delta=10^{-10}$–$10^{-3}$), $\Delta m_{\rm th}\approx(0.3\text{–}2)\times10^{-2}\sqrt\delta\,M_N$; the local threshold grows only like $\delta^{0.15}$.
* **Flat space.** Anti-trapping a sphere of radius $r$ needs local energy $\ge r/2G$, i.e. 270 times $W_{\rm th}$ at $r=r_N$.

The checkpoint's "about 400× cheaper" and "0.6–0.9% of $M_N$" referred to local energy with a coarser bracket. They are replaced by the
statements above.

### I.6 Robustness, back-reaction, stability (R6f-A to F)
* **A (boundary dynamics).** Light thin walls, natural width. The wall position is varied across the band (0.5, 0.9, 0.95, 0.995).
* **B (persistence).** The trapped neck persists to the end of every growing run.
* **C (growth).** Coherent anti-trapped false vacuum with area $\ge1$ and $r>5r_N$, converged between $h=0.1$ and $0.05$.
* **D (hiding).** No untrapped true vacuum at large $r$ in any growing run.
* **E (back-reaction).** Bunch–Davies Polyakov back-reaction:
  * growth survives at $c=20$, 60 and 150 with $W=0.071$, i.e. $k/r_N^2$ up to $1.3\times10^{-2}$;
  * at $c=300$ the threshold rises above 0.071, and 0.18 still grows;
  * vacuum polarisation alone never tips the throat ($c=30$–1000).
* **F (stability).**
  * Sub-threshold perturbations of either sign leave the hole.
  * The growing branch is stable to 10% changes of the segment's length, depth and wall width (Round 6c-I scan).
  * Non-spherical stability is open (P-D2).

---

## Part II — quantum throat physics

### II.1 P-R6g: Landau-level fermions and the strong QNEC
* **Proposition 6c-6 [Theorem].** On $M_2\times S^2_a$ with $q$ flux quanta, the Dirac operator has $\lvert q\rvert$ zero modes (index theorem).
  The tower sits at $\pm\sqrt{n(n+\lvert q\rvert)}/a$ with degeneracy $\lvert q\rvert+2n$. A 4D Dirac fermion is therefore $\lvert q\rvert$
  massless 2D fermions ($c=\lvert q\rvert$) plus massive modes, with gap $\ge\sqrt{1+\lvert q\rvert}/a$. A 4D mass passes to all modes.
  * Massive-mode corrections are suppressed for pulses of scale $a/\eta$ with $\eta\ll\sqrt{1+q}$.
  * For massive 2D fields the QNEC holds in Minkowski space by Ceyhan–Faulkner [Literature].
* **Proposition 6c-7 [Theorem].** Conformal vacua of any coordinate on any 2D metric $-e^{2\sigma}du\,dv$ saturate the strong QNEC
  $T_{ll}=\frac1{2\pi}(S_{ll}+\frac6cS_l^2)$ in the affine parameter, Weyl anomaly included. This is verified symbolically: the residual is 0.

### II.2 P-R6h: the Schwarzian regime
**Proposition 6c-8 [Theorem].** The Schwarzian partition function is $Z\propto\beta^{-3/2}e^{2\pi^2C/\beta}$ with $C=Q^3$. It gives
$\langle E\rangle=2\pi^2CT^2+\frac32T$ and ${\rm Var}\,E=4\pi^2CT^3+\frac32T^2$. Hence $\delta E/E=\sqrt{2T_{\rm brk}/T+1.5T_{\rm brk}^2/T^2}$,
checked against the density of states to $10^{-3}$.

Below $T_{\rm brk}$ the trapped/normal status of the deep throat is not classical. The bud must then be treated as a quantum state of JT/SYK
(Route Q), where "anti-trapped" is not an observable. For the Nariai throat the analogous statements are:

* **Fluctuations.** Radion fluctuations are $\ell_P/r_N$, below $W_{\rm th}/r_N=1.85\times10^{-3}$ once $r_N\gg540\,\ell_P$.
* **One loop.** The one-loop parameter is $Gc/12\pi r_N^2$.

### II.3 P-R6i: Theorem GSL-3 without spherical symmetry
**Theorem 6c-R6i [Theorem, conditional on (S)].** Hypotheses:

* **(H1)** Global hyperbolicity, an asymptotically flat Cauchy surface $\Sigma_0$, and a complete observer $\gamma$.
* **(H2)** The exterior of a compact core is foliated by closed leaves $S_R$, $R\ge\ell_*$, with ingoing expansion $\theta_-(S_R)\le-\theta_*<0$
  (normal data, no symmetry).
* **(S)** The support form of the local GSL: at each point $q$ of $H(\gamma)$ and for each $\epsilon>0$ there is a smooth past-support null
  hypersurface $W$ with $\theta_W(q)\ge-4G\hbar\,s(q)-\epsilon$, where $s$ is the outside-entropy variation density.
* **(H4)** Local weak gravitation: $4G\hbar\,s<\theta_*$.

Conclusion: if the past slice $P_0=I^-(\gamma)\cap\Sigma_0$ is relatively compact, it lies inside $S_{\ell_*}$. With Lemma GSL-1, which uses
no symmetry, parts (a)–(c) of Theorem GSL-3 hold without spherical symmetry.

*Proof.*

1. Let $R_*=\min\{R:\bar P_0\subset B_R\}$ and suppose $R_*\ge\ell_*$. The leaf $S_{R_*}$ touches $\bar P_0$ at a point $p\in\partial P_0\subset H(\gamma)$.
2. Every point of $I^-(\gamma)\cap J^+(\Sigma_0)$ lies in $D^+(\bar P_0)\subset D^+(\bar B_{R_*})$. Indeed, a past causal curve reaching $\Sigma_0$
   outside $\bar P_0$ would put its endpoint in $I^-(\gamma)$.
3. Hence $H(\gamma)$, and with it every past-support hypersurface $W$ at $p$, lies on the past side of the ingoing null hypersurface $N$ of
   $S_{R_*}$, and touches it at $p$.
4. The cross-sections $W\cap\Sigma_0$ and $S_{R_*}$ are tangent at $p$, with the first inside the second. So their mean curvatures satisfy
   $H_W\ge H_N$, while ${\rm tr}_\sigma K$ is the same. Therefore $\theta_W(p)=-H_W+{\rm tr}_\sigma K\le\theta_-(S_{R_*},p)\le-\theta_*$.
5. This contradicts (S) and (H4). $\square$

On (S):

* For a $C^2$ horizon it is exactly Wall's GSL for arbitrary cross-sections.
* For $\hbar\to0$ it is the Hawking–Galloway area theorem in support form (Galloway, Ann. Henri Poincaré 1 (2000) 543), which follows from
  the NEC and future-complete generators.
* Its semiclassical justification is the quantum-focusing analogue of that argument, which is the content of the conjecture.
* The de Sitter-parent corollary GSL-3b generalises the same way: $\partial P_0$ cannot lie inside any normal leaf.

### II.4 P-R6a: beyond Horndeski
**[Assessment + Literature.]**

* **$F(\phi)R$.** These theories reduce to the Einstein frame (Jacobson–Kang–Myers). There Farhi–Guth applies. A Jordan-frame daughter
  anti-trapped with $\phi\to$ const would be anti-trapped in the Einstein frame too. So a persistent daughter needs $\Omega^2=1-8\pi G\xi\phi^2\to0$,
  i.e. strong coupling.
* **Horndeski.** Flat non-singular cosmologies suffer gradient or tensor pathologies at some time (Kobayashi 2016, arXiv:1606.05831).
  Transient, localised NEC violation is not excluded.
* **Beyond Horndeski.** These admit stable wormholes (Franciolini et al., arXiv:1811.05481). Angular-gradient issues are under debate
  (Mironov–Rubakov–Volkova).
* **The GSL.** When species propagate at different speeds ($c_s\neq1$ around a time-dependent background) there is no universal GSL
  (Dubovsky–Sibiryakov, hep-th/0603158; Eling–Foster–Jacobson–Wall, hep-th/0702124).
* **Answer.** The dynamical-entropy second law is not available in that class. Theorem GSL-3 therefore does not bind it, and stability,
  not the second law, is the operative constraint.

### II.5 P-A1: quantum energy inequalities with exact constants
* **Proposition 6c-10 [Theorem, 2D CFT].** For a unitary 2D CFT in a state equal to the vacuum before $u_0$, the Fewster–Hollands null QEI
  $\int Tf^2\ge-\frac c{12\pi}\int f'^2$ with $f^2=e^{-(u-u_0)/\tau}$ gives
  $$\int T_{uu}e^{-(u-u_0)/\tau}du\ge-\frac c{48\pi\tau}.$$
  The conformal-map states with $A=1/\tau$ saturate it to $7\times10^{-4}$ (numerical check).
* **Corollary 6c-11 [Theorem + Proposition-N].** With the tipping law, $c\ge48\pi\tau W_{\rm th}=0.89(r_N/\ell_P)^2$.
* **In 4D.** There is no null QEI on a single generator even for free fields (Fewster–Roman, gr-qc/0209036). For the sphere-integrated s-wave
  flux, which is what tips the throat, each angular mode is a 2D field. Modes with $l\lesssim r_N/\tau$ act like massless ones, so 4D free
  fields give $c_{\rm eff}\approx N_{\rm species}\times O(1)$ and the conclusion stands.
* **Interacting theories.** The QNEC and the ANEC hold in Minkowski space (Ceyhan–Faulkner; Faulkner et al.). State-independent QEIs are not
  established.

### II.6 A source for large magnetic charge
**[Assessment.]**

* Dirac quantisation gives $q=2\sqrt\alpha\,Q/\ell_P$. A throat of radius $Q\sim r$ carries $q\approx0.17\,r/\ell_P$ quanta.
* It supports $c_{\rm LLL}=qN_f$ lowest-Landau-level 2D fermions (Maldacena, arXiv:2004.06084). This grows only linearly in $r$, far below
  what a quantum tip needs.
* Primordial monopoles, or magnetic black holes formed from them, are the only known source.
* **The Nariai route needs no charge at all.** The positive result above is neutral.

---

## Part III — the vacuum-energy obstacle (items 37–44)

* **37. Field theory for the false vacuum.**
  * Round 5's Theory F1 (tilted quartic, light thin walls, $\epsilon=2$) is the false vacuum used in every run.
  * The Nariai route adds one relation, the partner condition $\rho_f=1/(24\pi G^2M^2)$ [Theorem: $r_N=3GM$ is the photon sphere].
  * An eternal daughter additionally needs a metastable vacuum (decay rate per Hubble volume below $\sim H^4$).
* **38. Natural or dynamical configurations reaching the threshold.**
  * Dynamically: collapse inside a false-vacuum region, or a vacuum lifted in situ around a hole (Proposition 6c-15).
  * The threshold itself needs $W\ge W_{\rm th}$ of NEC deficit (Proposition 6c-13: unavoidable).
  * Nature supplies near-Nariai holes only in de Sitter, through Bousso–Hawking nucleation at a rate $\sim e^{-S_{\rm dS}/3}$ [Literature].
* **39. Energy–momentum constraints on preparing it.**
  * **Proposition 6c-13 [Theorem, spherical symmetry].** If $T_{uu}\ge0$, then along ingoing rays $r_u/\alpha^2$ is non-increasing. An
    anti-trapped point $p$ therefore lies on an ingoing ray that was anti-trapped all the way back to $\Sigma_0$, with $r(p)\ge r(q)$.
  * So the minimal radius of anti-trapped spheres never decreases. This holds in any parent, including de Sitter, where it is the Hubble
    radius. A daughter of Hubble radius $\ell_H$ smaller than the parent's minimal anti-trapped radius needs $T_{uu}<0$.
  * The prepared configuration itself is normal and NEC-compatible.
  * The tipping needs $\Delta m_{\rm MS}=\sqrt{f_0}W_{\rm th}$.
* **40. Stability.**
  * The near-Nariai band is an unstable equilibrium with response time $\tau=3.2r_N$.
  * NEC-respecting perturbations send it to the trapped branch (it remains a black hole).
  * Bunch–Davies one-loop effects are stabilising.
  * Hawking lifetime $\sim M^3\gg r_N$.
  * The growing daughter is robust (Part I.6).
* **41. Lower-energy alternatives.**
  * The invariant tipping cost $\propto\sqrt\delta$ can be made arbitrarily small by tuning. The hole's mass cannot.
  * $Mc^2=c^2/\sqrt{24\pi G^2\rho_f}$ is the universal law $M\sim m_P^3/\sqrt\rho$ again.
  * Lower vacuum energy means a heavier hole: $(0.43\,{\rm GeV})^4\leftrightarrow M_\odot$, $(100\,{\rm GeV})^4\leftrightarrow6\,M_\oplus$,
    $(10^{16}\,{\rm GeV})^4\leftrightarrow3.7$ g.
* **42. Non-spherical and quantum-gravity corrections.**
  * Theorem 6c-R6i removes the symmetry from GSL-3.
  * Quantum-gravity corrections: one loop $Gc/12\pi r_N^2\approx10^{-11}c$ for the 3.7 g case; fluctuations $\ell_P/r_N\approx2\times10^{-6}$.
  * Non-spherical stability of the tipping: open. The $l\ge1$ modes on the $S^2$ of Nariai are massive, $m^2\sim l(l+1)/r_N^2$ [Assessment].
* **43. Predictability and EFT validity.**
  * The runs are deterministic evolutions inside the domain of dependence of normal data.
  * The classical treatment holds for $r_N\gg540\ell_P$.
  * The ghost-condensate EFT holds inside its window (Proposition 6c-14).
  * Without GSL violation an eternal daughter could not be data-determined (GSL-3). With the ghost condensate, whose GSL is not universal,
    it is.
* **44. Formation rate or deterministic evolution.**
  * **Deterministic evolution is shown** (Proposition 6c-9 with a classical source).
  * A spontaneous quantum tip would need a radion fluctuation of $W_{\rm th}/r_N$ against a width $\ell_P/r_N$. Its probability is
    $\sim\exp[-\frac12(1.85\times10^{-3}r_N/\ell_P)^2]$, which is $e^{-4\times10^5}$ for the 3.7 g case [Assessment].

---

## Part IV — foundational mathematics

### Theorem 6c-S3 (P-S3 solved)
In units $R=\kappa\rho$, $2GM=\kappa^3m$ ($m=4E^{-3/2}$), the barrier is $V=1-X^2+\kappa^2V_1$ with $X=m/\rho^2-\rho/2$ and $V_1=-(\rho^2/2+m/\rho)$.
Expanding the junction angle and $B(0)=\pi\kappa^4(1-2\kappa^2)/8G$ gives
$$B_{\rm an}=B(0)\Phi(E)\,[1+c(E)\kappa^2/\chi^2],\qquad c(E)=2-2m\frac{I_0'}{I_0}+\frac{I_1}{I_0},$$
$$I_0=\int\rho^2\sqrt{1-X^2}\,d\rho,\quad I_1=\int\rho^2P\Big(\sqrt{1-X^2}-\frac1{2\sqrt{1-X^2}}\Big)d\rho,\quad I_0'=-\int\frac X{\sqrt{1-X^2}}d\rho .$$

**Asymptotics.**

* **As $E\to0$.** The barrier is a thin shell, $\rho=\rho_*(1+y)$ with $(1+y)^{-2}-(1+y)=2X/\rho_*$. Series inversion and Chebyshev integrals
  give $I_0=\pi\rho_*^2/3+\pi/3$ and $I_1=-25\pi\rho_*^2/54-2\pi/3$. Hence
  $$c(0^+)=\tfrac23-\tfrac{25}{18}=-\tfrac{13}{18}$$
  exactly, which proves the Round 3 conjecture.
* **As $E\to\infty$.** $I_0=\pi/2+2m$ and $I_1+2I_0=2m$, built from $\frac{20}3m$, $-\frac23m$ and $-4m$. So
  $c=-(16/\pi)E^{-3/2}+O(E^{-9/4})$. The Round 3 fit $-5E^{-3/2}$ was this law, since $16/\pi=5.093$.

**Checks.**

* The closed form matches the Round 3 25–35-digit values at $\kappa/\chi=0.01$ to $<6\times10^{-5}$, i.e. their $O(\kappa^2)$ remainder.
* $c\cdot E^{3/2}\to-5.0928$ at $E=10^6$.

### Status of the others
* **P-N1′.** Reformulated: $-EJ=E^2\int sk(W-s^2)^{-1/2}ds$ on regime B, both ends proved. The monotonicity proof (an interval-arithmetic
  certificate over $(\gamma,E)$ with a singular integrand) was not closed this round. **Open.**
* **P-T5″.** Open; unchanged.
* **P-Q9′/P-Q8, P-Q1, P-Q3, P-Q7, P-T3′, P-Q4/5/6, P-A2.** Not advanced this round. Their statements in OPEN_PROBLEMS stand. Nothing in
  Round 6c depends on them, because the Nariai route is classical.
* **P-L1.** The full texts are now reachable (web access works in this environment). The comparison itself remains to be done. **Open.**

---

## Part V — matter dynamics and observables

### Proposition 6c-16 (P-R7): 3+1 greybody factors
`greybody.transmission` solves Regge–Wheeler for spin 0. Checks:

* $\Gamma_0/16M^2\omega^2=1.034$, 1.16 at $M\omega=0.005$, 0.02, i.e. $\to1$ (absorption cross-section $\to A_H$);
* $\sum(2l+1)\Gamma_l/27M^2\omega^2=0.97$, 1.02 at $M\omega=1$, 1.5 (geometric optics).

Results:

* **Hawking power.** For one massless scalar, $P=7.44\times10^{-5}/M^2$ (l = 0: 90%, l = 1: 9.7%, l = 2: 0.24%). This reproduces the
  standard value (Page, Elster), which validates the solver. The number rate is $\dot N=6.64\times10^{-4}/M$.
* **The branch difference in 3+1.** Round 7's transient is a delay of the switch-on by $D/\kappa=4DM$ ($D=1.9$–2.2). In 3+1 it removes
  $\dot N\cdot4DM$ = 0.005–0.006 quanta per event per scalar degree of freedom.
* **Why 3+1 matters.** This is IR-finite, because $\Gamma_0\propto\omega^2$. The 1+1 number diverges logarithmically, so Round 7's
  "0.1 quanta" (energy over mean frequency) is replaced by this.

### Signatures
**Proposition 6c-17 [Theorem + Proposition-N] (threshold indistinguishability).**

* Outcomes switch at $W=W_{\rm th}$, while every exterior datum (mass, emitted flux) depends continuously on $W$.
* Near threshold, therefore, no exterior measurement of finite precision distinguishes "daughter formed" from "no daughter". Far from
  threshold, Birkhoff fixes the exterior to Schwarzschild in both cases.
* The only exterior record is the source's own flux, $\Delta m=\sqrt{f_0}W$, and it is the same on both sides of the threshold.

### Status
* **P-T2, P-D2, P-D3, P-C2, P-C4.** Not advanced. P-D2 now has a sharper target: non-spherical stability of Nariai tipping.
* **P-O1/O2.** Natural near-Nariai holes arise only in de Sitter phases (Bousso–Hawking). In the present asymptotically flat universe there is
  no population.

---

## Part VI — alternatives

* **Route D** is now concrete: the ghost condensate in the Nariai throat (Proposition 6c-14), at 145–430 times below its stability ceiling.
  Lorentz-invariant Galileon/beyond-Horndeski NEC violation is stable only with superluminal modes and anisotropy (Dubovsky–Grégoire–Nicolis–Rattazzi,
  hep-th/0512260).
* **Route Q.** Proposition 6c-8 marks where it starts.
* **P-R6b** (seed boundary layer) is superseded for the throat problem by the Nariai route, which needs no seed.
* **P-R6c** (stabiliser): the CLNS ghost condensate is the candidate stabiliser for NEC violation.
* **P-R6d, P-R6e.** Open.

---

## Failed approaches recorded this round
* **F26 (correction).** The pulse was quoted as "0.6–0.9% of $M_N$" and "about 400× cheaper". That was local energy in the corner frame, not
  a mass. The invariant cost is $\sqrt{f_0}$ times it (Part I.5).
* **F27.** A designed QNEC-saturating CFT pulse (step profile) leaves a permanent flux $(c/48\pi)Y^2$ that outweighs its negative lobe.
  Even optimal conformal-map pulses (rise, slow fall) fail up to $c=500$, as Corollary 6c-11 requires.
* **F28.** Vacuum polarisation as the source (Bousso–Hawking anti-evaporation in the Bunch–Davies state): stabilising; no tip at $c=30$–1000.
* **F29.** Positive-energy pulses: the trapped branch, never a daughter (Farhi–Guth, Proposition 6c-13).
* **F30.** A hole formed in true vacuum with only a thin false-vacuum layer (inner wall at $u\le60$, i.e. $\lesssim1.7u_{\rm hole}$): no growth
  even at 4× the pulse.
* **F31.** Non-minimal coupling $\xi\phi^2R$ as the classical source: conformally equivalent to GR with NEC matter, so no persistent daughter
  below $\Omega\to0$.

## Open problems added
* **P-R6j.** Simulate the in-situ lift of a false vacuum around a hole (preparation, Proposition 6c-15) through to the tipping.
* **P-R6k.** Ghost condensate in the Nariai throat: simulate the EFT itself, including its gradient and Jeans dynamics, rather than its flux.
* **P-D2′.** Non-spherical stability of Nariai tipping.
