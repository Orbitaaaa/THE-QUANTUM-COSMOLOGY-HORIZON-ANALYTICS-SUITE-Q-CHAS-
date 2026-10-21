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
| R3 | **Branch selection T0**: the canonical WKB phase is the analytic continuation along the barrier. The FGG slicing jump is a horizon monodromy ($P\to P\pm i\pi R/G$, $t\to t\pm i\beta/2$) | Proposition (analytic argument + two numerical criteria) | TUNNELLING.md §3 |
| R4 | **T1**: closed form at $M=0$, $B=S_E[\text{instanton}]+\frac1{4G}(\sum A_{\max}-\sum A_{\min})$; equals the CdL decay exponent of the false vacuum when the false vacuum is the larger cap | Theorem | TUNNELLING.md §4 |
| R5 | **T2 first law** $d\tilde B=\frac{{\rm Vol}^{\rm in}_E}{8\pi G}d\Lambda_{\rm in}+A^W_Ed\sigma-\Delta t^W_EdM$ | Theorem | TUNNELLING.md §4 |
| R6 | **T3 Smarr identity** $B=-\frac{\Lambda_{\rm in}}{8\pi G}{\rm Vol}^{\rm in}_E-\frac\sigma2A^W_E-\frac M2\Delta t^W_E+\frac{\mathcal A}{4G}$: a covariant expression for the FGG exponent that never uses the exterior pseudomanifold; integral representation from the $M=0$ instanton | Theorem | TUNNELLING.md §4 |
| R7 | **NG3 exterior embedding obstruction**: no manifold Euclidean interpolation for any transition from regular laboratory data to the inflating bounce | Theorem | NO_GO_THEOREMS.md §3 |
| R8 | **NG4**: regime-B transitions are also non-manifold (horizon-angle overshoot $\lvert\theta_D\rvert>\pi$) | Proposition-N | NO_GO_THEOREMS.md §3 |
| R9 | **Interior counter-example**: $\varepsilon_{\rm in}$ flips do not obstruct embedding; the obstruction needs a flip while the wall recedes from a Euclidean horizon | Proposition-N | NO_GO_THEOREMS.md §3 |
| R10 | **T5**: over laboratory transitions, $B\ge0.56\,B(0)$ for all γ, and $B\ge B(0)=B_{\rm CdL}$ for γ ≳ 0.2. Creating a daughter from laboratory data is not exponentially easier than vacuum decay | Proposition-N | TUNNELLING.md §4 |
| R11 | **Physical scale of the exponent**: $B(M\to0)$ follows the flat CdL scaling $\sigma^4/\rho^3$ at small κ/χ, not $(M_{\rm Pl}/H)^2$, but $B$ rises steeply with $M$. GUT example (σ^{1/3} = ρ^{1/4} = 10^{16} GeV): $B=133$ only for picogram masses, $B\sim10^5$–$10^7$ for milligram–gram masses | computation | TUNNELLING.md §4 T6 |
| R12 | **Unification**: the classical (FG), initial-data (NG2) and Euclidean (NG3) obstructions all stem from creating a marginally trapped surface (the bridge) in an asymptotically flat exterior | Theorem (each part) + Conjecture NG5 (general matter) | NO_GO_THEOREMS.md §4 |
| R13 | **Sector structure**: the four orientation sectors are sheets of one squared constraint, which communicate only under the barrier; the Hamiltonian generates imaginary shifts of the wall area by $8\pi\ell_P^2$ | derivation | QUANTIZATION.md §2–3 |

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
