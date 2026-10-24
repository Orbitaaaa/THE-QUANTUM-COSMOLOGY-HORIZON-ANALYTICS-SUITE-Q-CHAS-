# Round 6b: GENESIS — a serious attempt to get around the GSL

**Programme.** Round 6 concluded that the only laboratory-scale route, self-amplifying NEC violation, "must violate the GSL". The brief for
this round corrects that premise. Violating the null energy condition (NEC) does not by itself violate the generalised second law (GSL):
Gao–Jafferis–Wall traversable wormholes violate the averaged NEC and respect the GSL. The brief proposes five routes:

* **H**: a daughter that never disconnects (Visser's umbilical cord);
* **D**: beyond-Horndeski genesis with an inhomogeneous throat;
* **S**: keep $\Delta S_{\rm gen}\ge0$ with the right gravitational entropy and find which premise of Wall's theorem fails;
* **Q**: the daughter as a quantum spacetime state;
* a separately labelled theory that violates the GSL directly.

It asks for a **good** answer to Problem B in this form. The routes were treated as recommendations. Where the analysis showed a sharper
formulation, that formulation was followed.

Everything below is theory. **No apparatus is designed and no experimental procedure is proposed.** Labels as in Round 6:

* **[Theorem]**: proved under the stated hypotheses;
* **[Proposition]**: derived at the stated level of rigour;
* **[Proposition-N]**: numerical;
* **[Literature]**: a result of the cited work, not re-derived here;
* **[Assessment]**: a judgement;
* **[Conjecture]**.

**Files.**

* Code: `CODE/babyuniverse/gsl.py`.
* Script: `CODE/scripts/round6b_gsl.py`.
* Data: `DATA/round6b_gsl.json`.
* Figure: `FIGURES/round6b_gsl.png`.
* Tests: `CODE/tests/test_gsl.py` (5).

---

## The answer to Problem B (GSL form)

**Yes, the GSL can be kept. The GSL is not what makes laboratory genesis hard, and this round proves exactly what it does forbid.**

1. **What the GSL forbids (Theorem GSL-3, the genesis dichotomy).**
   * Assume the GSL on causal horizons of complete observers, global hyperbolicity, and laboratory data that are normal and weakly
     gravitating.
   * Then a complete daughter observer must keep a non-compact past: arbitrarily distant parts of the parent stay in its causal past, so
     the neck never closes.
   * Such an observer can never end in a de Sitter phase.
   * Consequently **an eternal daughter and a daughter determined by the laboratory data exclude each other.** This is the precise content of
     Wall's "cannot restart inflation", re-derived here with self-contained proofs in spherical symmetry.

2. **What it permits: two classes, each with an exact price.**
   * **(a) Transient, connected daughters (Route H made precise).**
     * The daughter inflates for a while, then leaves inflation with $\Lambda_{\rm final}\le\Lambda_{\rm parent}$ (or crunches).
     * The neck must persist forever. An extremal magnetic black hole never evaporates, so it supplies this automatically.
     * The GSL places no further constraint. The price is set by the QNEC instead (Theorem GSL-4). Flipping the ingoing expansion costs a
       change in the outside-entropy gradient of $S_{\rm BH}\,|\theta_-|$.
     * In weakly gravitating regions that cost is $\sim10^{70}$ for a 1 m sphere and cannot be paid. In a near-extremal throat it is
       discounted by $\sqrt f\to0$.
     * In a throat of radius $Q$ at temperature $T$, QNEC-saturating matter with central charge $c$ anti-traps spheres at proper distance
       $\rho$ from the horizon iff $T\sinh(\rho/Q)<(c/12)\,T_{\rm brk}$. Here $T_{\rm brk}=1/2\pi^2Q^3$ is the Schwarzian scale.
       This is shown three independent ways, including an exact Einstein–Maxwell–null-fluid metric (Propositions GSL-5 to GSL-7).
   * **(b) Eternal daughters beyond a Cauchy horizon (escape E2).**
     * The classical charged-bubble route of Proposition G4 is GSL-compatible, because global hyperbolicity fails.
     * The daughter is not determined by the laboratory data.
     * Near extremality the route survives, but its mass window has width $\le r_+-r_-=2\kappa_-r_-^2$ (Proposition GSL-9). Weakening the
       inner-horizon blueshift closes the window, so extremality buys a milder Cauchy horizon only with fine tuning.

3. **What remains hard is energy density, not entropy (Proposition GSL-8).**
   * A flipped region grows into a daughter only if its mass exceeds the charged-Nariai bound for the local vacuum energy.
   * For an extremal throat this means $\Lambda_{\rm loc}Q^2>3/16$ if the vacuum energy is injected, $0.2212$ if it is converted locally.
     A static throat in a fixed vacuum never grows.
   * Every route therefore meets the same law, $M\simeq(0.07\text{–}0.16)\,m_P^3/\sqrt\rho$:
     * Nariai floor: $0.0665$;
     * throat bud: $0.0864$;
     * thin-wall classical route: $0.155$.
   * **The throat solves the GSL problem. It does not solve the energy-density problem.** A throat of radius $10^{-13}$ m (mass
     $1.3\times10^{14}$ kg) still needs $\rho^{1/4}\approx5\times10^{7}$ GeV inside it.

4. **Route D** (classical NEC violation from non-minimal couplings) is not bounded by the QNEC, so Theorem GSL-4 does not apply to it. The
   dichotomy GSL-3 permits its transient daughters. Its obstacles are dynamical: stability of a time-dependent inhomogeneous throat,
   strong coupling, and unitarity below the cutoff. **Direct GSL violation** stays a separate, contested class (multi-speed Lorentz
   violation). It is not adopted.

**Correction F21** (Round 6, A3 and summary item 5): "the sector must violate the GSL" is withdrawn.

* What remains true: a minimally coupled *classical* NEC-violating sector conflicts with black-hole thermodynamics *as a theory*, because
  its stress can be dropped into a black hole and lower the area with no compensating entropy.
* But the genesis process itself does not need to violate the GSL. The GSL constrains the daughter's *fate* (Theorem GSL-3), not its
  formation.

---

## 1. Setting

* **Spacetime.** $(M,g)$ is four-dimensional and time-oriented. It solves the semiclassical equations
  $G_{ab}+\Lambda g_{ab}=8\pi G\langle T_{ab}\rangle$.
* **Observers.** An *observer* is a future-inextendible timelike curve $\gamma$. It is *complete* if it is future-complete.
* **Causal horizon** (Jacobson–Parentani, Wall): $H(\gamma)=\partial I^-(\gamma)$ for a complete observer.
* **GSL** (Wall's assumption, which we keep). On every causal horizon,
  $$S_{\rm gen}(\sigma)=\frac{A(\sigma)}{4G\hbar}+S_{\rm out}(\sigma)\quad(+\text{EFT corrections of Wall/HKR type})$$
  is non-decreasing along future-directed cross-sections $\sigma$. It is *not* assumed on arbitrary null surfaces; Wall's erratum notes
  counterexamples there.
* **QNEC.** $T_{kk}\ge\frac{\hbar}{2\pi}S''_{\rm out}$ per unit area. Integrated over a round sphere and combined with strong subadditivity,
  this gives $A\,T_{kk}\ge\frac{\hbar}{2\pi}\,d^2S_{\rm out}/d\lambda^2$ for uniform deformations. The off-diagonal second variations are
  $\le0$.
* **Strong 2D QNEC.** For a 2D CFT, $T_{kk}\ge\frac{\hbar}{2\pi}(S''+\frac6cS'^2)$; conformal vacua saturate it (Wall, arXiv:1105.3520).
* **Laboratory data.** $\Sigma_0$ is a Cauchy surface of the globally hyperbolic development, asymptotically flat (or asymptotically de Sitter
  for Corollary GSL-3b).
  * *Normal:* every round sphere has $\theta_+>0>\theta_-$.
  * *Weakly gravitating at scale $\ell_*$:* for every round sphere of radius $\ge\ell_*$ in $\Sigma_0$,
    $|dS_{\rm out}/d\lambda|<A|\theta_-|/4G\hbar$ on its ingoing congruence. This is the Bekenstein-type statement that matter entropy
    gradients are far below the black-hole entropy gradient of the sphere.
* **Past slices.** For a Cauchy surface $\Sigma$, $P_\Sigma:=I^-(\gamma)\cap\Sigma$. A daughter observer is *dependent* if no $P_\Sigma$ is
  relatively compact.

## 2. What the GSL forbids

### Lemma GSL-1 (compact pasts propagate backward) — [Theorem]
In a globally hyperbolic spacetime, let $\gamma$ be a future-inextendible timelike curve that meets a Cauchy surface $\Sigma_L$. If
$P_{\Sigma_L}$ is relatively compact, then $P_\Sigma$ is relatively compact for every Cauchy surface $\Sigma\subset J^-(\Sigma_L)$.

*Proof.*
1. Let $p\in P_\Sigma$, so $p\ll\gamma(s)$ for some $s$. Increasing $s$ keeps this true, so take $\gamma(s)\in I^+(\Sigma_L)$.
2. The timelike curve from $p\in J^-(\Sigma_L)$ to $\gamma(s)$ crosses $\Sigma_L$ at a point $q$ with $q\ll\gamma(s)$. Hence
   $q\in P_{\Sigma_L}\subset K:=\overline{P_{\Sigma_L}}$ and $p\in J^-(K)$.
3. Thus $P_\Sigma\subset J^-(K)\cap\Sigma$. For compact $K$ this set is compact in a globally hyperbolic spacetime, since
   $J^-(K)\cap J^+(\Sigma)$ is compact. $\square$

### Proposition GSL-2 (no compact past above normal data) — [Theorem, spherical symmetry]
Hypotheses:
* the spacetime is spherically symmetric, globally hyperbolic and semiclassical, and obeys the GSL;
* $\gamma$ is the regular centre and is complete;
* $\Sigma_0$ is normal and weakly gravitating at scale $\ell_*$.

If $P_0=P_{\Sigma_0}$ is relatively compact, its outer boundary is a round sphere of areal radius $r_0<\ell_*$. Otherwise the GSL fails on
$H(\gamma)$ at that sphere.

*Proof.*
1. $P_0$ is open, spherically symmetric and contains $\gamma\cap\Sigma_0$. Let $S_0$ be its outer boundary sphere, so $S_0\subset H(\gamma)$.
2. For $p\in H(\gamma)$, $I^+(p)\cap I^-(\gamma)=\emptyset$. Otherwise $p\ll q\ll\gamma$ would put $p$ in the open set $I^-(\gamma)$.
3. Suppose the ingoing null hypersurface $N_-$ from $S_0$ reached the centre at a point $c$. The centre points to the future of $c$ lie in
   $I^+(p)$ for every $p\in S_0$, by symmetry, and they include late points of $\gamma$. That contradicts step 2.
4. So $H(\gamma)$ contains $N_-$ near $S_0$. On normal data $\theta_-<0$, so $dA/d\lambda=\theta_-A<0$ there.
5. The GSL needs $dS_{\rm out}/d\lambda\ge A|\theta_-|/4G\hbar$. For $r_0\ge\ell_*$ this contradicts weak gravitation. $\square$

### Theorem GSL-3 (the genesis dichotomy) — [Theorem, spherical symmetry]
Under the hypotheses of GSL-2, with $\Sigma_0$ asymptotically flat and the centre observer of the daughter complete:

* **(a)** No Cauchy surface to the future of $\Sigma_0$ has a relatively compact past slice $P_\Sigma$. The only exception is the
  non-semiclassical case where $P_0$ is a ball of radius $<\ell_*$.
* **(b)** In particular the daughter cannot end in a de Sitter phase. Here "ending in a de Sitter phase" means that some late $P_\Sigma$ is
  bounded by a compact cosmological horizon.
* **(c)** For every $R$, some point of $\Sigma_0$ at areal radius $>R$ lies in the observer's past. The daughter stays causally downstream
  of the entire parent, so the neck never closes.

*Proof.* If some $P_\Sigma$ were relatively compact, Lemma GSL-1 would make $P_0$ relatively compact. Proposition GSL-2 then gives either
GSL violation or $r_0<\ell_*$. Parts (b) and (c) restate (a). $\square$

**Corollary GSL-3a (eternal vs determined).**
* Apply the theorem to the domain of dependence $D(\Sigma_0)$, which is globally hyperbolic. A daughter that ends in a de Sitter phase for a
  complete observer cannot lie in $D(\Sigma_0)$.
* An eternal daughter therefore lies, at least partly, beyond the Cauchy horizon of the laboratory data. Its existence is not determined by
  those data.

**Corollary GSL-3b (de Sitter parent) — [Proposition, spherical symmetry].**
* With a de Sitter parent ($\Lambda_p>0$) the past slice can be compact. The GSL then forces its boundary to have $\theta_-\ge0$, i.e. to
  lie at or beyond the parent's cosmological apparent horizon.
* Monotonicity of $S_{\rm gen}$ along $H(\gamma)$ to late times gives
  $A_d\ge A_c(\Sigma_0)-4G\hbar\,\Delta S_{\rm out}$. To leading order, $\Lambda_d\le\Lambda_p\,[1+O(GM\sqrt{\Lambda_p})+O(G\hbar\Lambda_p\Delta S_{\rm out})]$.
* A daughter may end with less vacuum energy than its parent, never more. This is Wall's "cannot restart inflation", made quantitative.

**Relation to the literature.**
* Wall (arXiv:1010.5513) classifies babies as dependent (permitted; the black hole cannot evaporate completely), disconnected (forbidden by
  the GSL) and transcendental (forbidden by global hyperbolicity) [Literature].
* GSL-2 is the disconnected case. GSL-1 together with GSL-2 gives the de Sitter exclusion.
* The new point used in this round is that **transient** dependent daughters are permitted. Nothing in the GSL forbids them.
* Bousso's robust singularity theorem (PRL 2025, arXiv:2501.17910) excludes semiclassical bounces inside black holes [Literature]. A
  persistent neck therefore cannot be a bouncing trapped neck. It must be untrapped or marginal, like the throat of an extremal hole, or it
  must pass a Cauchy horizon.

**The escapes.** Every daughter that is not excluded uses at least one of these:
* **E1:** transient daughters ($\Lambda_{\rm final}\le\Lambda_{\rm parent}$), or crunching ones (incomplete observer);
* **E2:** Cauchy horizons;
* **E3:** a Planckian past, i.e. a non-semiclassical origin;
* **E4:** a modified gravitational entropy.

E4 adds nothing for eternal daughters:
* the Wall/HKR higher-derivative corrections are EFT-suppressed;
* for scalar–tensor entropies $\int F(\phi)\,dA/4G$, the step "$\theta_-<0\Rightarrow S_{\rm gen}$ decreases" in GSL-2 survives wherever $F$
  is nearly constant on $\Sigma_0$, i.e. outside strong coupling.

## 3. The price inside the domain of dependence: Theorem GSL-4

Inside $D(\Sigma_0)$, any anti-trapped sphere needs $\int T_{uu}<0$ on its ingoing ray (Farhi–Guth; Theorem G1). Along that ray
$\alpha^2\partial_u(r_u/\alpha^2)=-4\pi Gr\,T_{uu}$ keeps the sign of $r_u$ if $T_{uu}\ge0$.

### Theorem GSL-4 (QNEC cost of genesis) — [Theorem]
Hypotheses: a spherically symmetric semiclassical spacetime; an ingoing radial null ray with affine parameter $\lambda$; areal radius
$r(\lambda)$; $\theta=2r'/r$; the integrated QNEC.

Then
$$r''\le-\frac{G\hbar}{2\pi r}\,S''_{\rm out},\qquad\text{so a flip }\theta(\lambda_0)<0\le\theta(\lambda_1)\text{ needs}\qquad
-\int_{\lambda_0}^{\lambda_1}\frac{r_0}{r}\,S''_{\rm out}\,d\lambda\ \ge\ \frac{A_0}{4G\hbar}\,|\theta_0| .$$

*Proof.*
1. Raychaudhuri for round spheres, with zero shear, gives $2r''/r=-8\pi GT_{kk}$.
2. The integrated QNEC gives $T_{kk}\ge\hbar S''/(2\pi\cdot4\pi r^2)$.
3. Integrate and use $r'_1-r'_0\ge|r'_0|$. $\square$

Both sides scale in the same way under $\lambda\to c\lambda$. `qnec_raychaudhuri` verifies equality for a QNEC-saturating state and strict
inequality with slack.

**What the cost means.**
* **Flat space.** In a static frame the ingoing normal $k=(u-n)/\sqrt2$ has $|\theta_-|=\sqrt{2f}/r$. For a 1 m sphere in flat space,
  $S_{\rm BH}=1.2\times10^{70}$. The entropy gradient must change by $3.9\times10^{35}$ per Planck length of affine parameter, i.e. by
  $\sim S_{\rm BH}$ across the sphere.
* **Localisation — [Proposition].** For weakly gravitating matter of energy $E$ inside radius $r$, the Bekenstein bound
  $S\le2\pi Er/\hbar=(2GE/r)\,S_{\rm BH}$ caps the available entropy far below $S_{\rm BH}$. The cost can be met only where $2Gm/r\to1$:
  in the near-horizon zone of a black hole. There $|\theta_-|\propto\sqrt f\to0$ (the throat discount).
* **Anti-trapped floor — [Theorem, elementary].** An anti-trapped sphere of radius $R$ has $1-2Gm/R<0$. Its Misner–Sharp mass therefore
  exceeds $Rc^2/2G$:
  * $6.7\times10^{26}$ kg at 1 m;
  * $6.7\times10^{13}$ kg at $10^{-13}$ m.

  A homogeneous region needs $\rho>3c^4/8\pi GR^2$, i.e. $\rho^{1/4}>29$ GeV at 1 m and $9\times10^7$ GeV at $10^{-13}$ m.

## 4. Two dimensions, exactly (the semiclassical limit of Route Q)

In the s-wave/JT description of a throat, the dilaton in entropy units is $\phi=(A-A_0)/4G\hbar$. It obeys
$\phi''=-\frac{2\pi}{\hbar}T^{(2)}_{kk}$, up to the classical term $2\pi r'^2/G\hbar$, which is negligible in a throat.

### Proposition GSL-5 (exact quantum focusing and a duration bound) — [Theorem]
Assume 2D CFT matter obeying the strong QNEC.
* **(a)** The quantum expansion $\Theta=\phi'+S'_{\rm out}$ obeys $\Theta'\le-\frac6c S'^2_{\rm out}\le0$.
* **(b)** An anti-trapped episode ($\phi'>0$) that starts quantum-normal ($\Theta_0<0$) lasts an affine length $<c/(6|\Theta_0|)$.
* **(c)** In the 4D s-wave model with a saturating state,
  $\Theta_{4D}'=-\frac12\Theta_{4D}^2+s^2\big(\frac12-\frac{6}{c}S_{\rm BH}\big)\le0$ whenever $S_{\rm BH}\ge c/12$. Here
  $s=4G\hbar S'/A$.

*Proof.*
* (a) $\Theta'=\phi''+S''=-\frac{2\pi}\hbar T+S''\le-\frac6cS'^2$.
* (b) Let $G(\lambda)=\frac6c\int_0^\lambda S'^2$. Then $\Theta\le\Theta_0-G$.
  * $\phi'>0$ needs $-S'>|\Theta_0|+G$, hence $G'>\frac6c(|\Theta_0|+G)^2$.
  * Integrating over the episode, $u=|\Theta_0|+G$ satisfies $1/u_s-1/u_e>6\Delta/c$, so $\Delta<c/6|\Theta_0|$.
* (c) Substitute the saturated strong QNEC into $\theta'=-\theta^2/2-8\pi GT_{kk}$ and complete the square. $\square$

In the conformal vacuum relative to a reparametrisation $X(\lambda)$, with $y=\ln X'$:
$S'=-\frac c{12}y'$ and $T=\frac{c\hbar}{24\pi}(\frac12y'^2-y'')$.

* A step of $y'$ by $Y$ kicks $\phi'$ by $\frac c{12}Y$.
* The positive tail $\frac{c\hbar}{48\pi}Y^2$, the quantum interest, refocuses within $\sim2/Y$.
* A bump leaves the ray *more* focused than before. The test checks this.

## 5. The throat

### Proposition GSL-6 (throat discount and the flip window) — [Proposition, verified three ways]
Setting: near-extremal Reissner–Nordström with charge $Q$ and temperature $T$, $QT\ll1$, $r_+r_-=Q^2$.

* At proper distance $\rho$ from the horizon, in the AdS$_2$ region,
  $$\sqrt f=2\pi QT\,\sinh(\rho/Q)\,[1+O(QT)].$$
  The exact profile agrees to $10^{-8}$.
* Take a QNEC-saturating CFT state with central charge $c$ whose reparametrisation steps by $\eta/Q$ over the static observer's proper time.
  The spheres at depth $\rho$ become anti-trapped iff
  $$\sqrt f<\frac{c\,\eta}{12\pi Q r}\quad\Longleftrightarrow\quad T\sinh(\rho/Q)<\frac{c\eta}{12}\,T_{\rm brk},\qquad T_{\rm brk}=\frac1{2\pi^2Q^3}.$$
* The episode lasts at most $c/(12\pi r\sqrt f)$ of local proper time. Inside the window that is $\ge Q$.
* Semiclassical control (small Schwarzian fluctuations) needs $T\gg T_{\rm brk}$. So the **flip window**
  $T_{\rm brk}\ll T<\frac{c\eta}{12}T_{\rm brk}/\sinh(\rho/Q)$ is non-empty iff $c\eta\gg12\sinh(\rho/Q)$. **A large number of light
  species is the essential ingredient.**

*Derivation.*
1. Along the ingoing ray $du/d\lambda=\sqrt2/\sqrt f$, and a local proper-time scale $Y_{\rm loc}$ has $Y_u=\sqrt f\,Y_{\rm loc}$.
2. The kick $\Delta r'=(G/r)\frac{c\hbar}{24\pi}Y_\lambda$ must exceed $|r'_0|=\sqrt{f/2}$.
3. With $Y_\lambda=\sqrt2\,Y_{\rm loc}$ this is $cY_{\rm loc}/(12\pi r)>\sqrt f$.

**Three independent checks** ($Q=10^3$, $T=1.97\,T_{\rm brk}$, $\rho=Q$, $\eta=1$; `round6b_gsl.json`):

| $c$ | analytic | 4D s-wave CFT integration | exact Vaidya | kick / need | $\delta m$ / budget | max $\Delta\Theta$ |
|---|---|---|---|---|---|---|
| 5 | no | no | no | 0.180 | 0.178 | $-1\times10^{-18}$ |
| 20 | no | no | no | 0.718 | 0.714 | $-3\times10^{-18}$ |
| 40 | yes | yes | yes | 1.437 | 1.428 | $-6\times10^{-18}$ |
| 100 | yes | yes | yes | 3.59 | 3.57 | $-1\times10^{-17}$ |

The anti-trapped durations (427 and 1019) respect the bound of GSL-5 (2025 and 5036). Quantum focusing holds along every ray.

**Where large $c$ comes from — [Literature + Assessment].**
* In a magnetic throat, each light charged fermion contributes $q$ two-dimensional modes: the lowest-Landau-level degeneracy
  (Maldacena–Milekhin–Popov, arXiv:1807.04726). So $c\approx N_fq$.
* The geometric charge is $Q=q\,\ell_P/2\sqrt\alpha$.
* The modes are light only if $m_fQ\lesssim1$. For electrons that means $Q\lesssim\lambda_C\approx3.9\times10^{-13}$ m.
* **Example:** $Q=10^{-13}$ m, i.e. $q\approx1.1\times10^{21}$, $c\sim10^{21}$, extremal mass $1.3\times10^{14}$ kg.
  * $T_{\rm brk}=3\times10^{-35}$ K;
  * $T_{\rm flip}(\rho=Q)=2\times10^{-15}$ K;
  * so the window spans 20 decades.

  Magnetic black holes are hypothetical (they need monopoles). These numbers are theoretical scales, not a device.

### Proposition GSL-7 (the exact metric) — [Theorem for the metric; Proposition for the CFT flux]
The outgoing Bonnor–Vaidya–Reissner–Nordström metric
$$ds^2=-\Big(1-\frac{2m(u)}r+\frac{Q^2}{r^2}\Big)du^2-2\,du\,dr+r^2d\Omega^2,\qquad T_{uu}=-\frac{m'(u)}{4\pi r^2},$$
is an exact solution of Einstein–Maxwell plus a null fluid.

* With $k_+=\partial_r$ and $k_-=\partial_u-\frac f2\partial_r$ ($k_+\cdot k_-=-1$), the expansions are $\theta_+=2/r$ and
  $\theta_-=-f(u,r)/r$. A sphere is anti-trapped iff
  $$m(u)-M>\tfrac12 r f(r).$$
  This is exactly the budget of Theorem G1 in mass units, and **it vanishes at an extremal horizon.**
* Taking the CFT flux $m'(u)=-\frac c{24\pi}(\frac12y'^2-y'')$ gives $\delta m\le cY_u/24\pi$. That reproduces GSL-6 identically (table
  above).
* The flipped spheres lie between the source depth and the threshold depth. The flux is outgoing, so it never crosses the future event
  horizon, and the GSL on that horizon is untouched.

### Proposition GSL-8 (growth: the charged-Nariai criterion) — [Theorem for the criterion; values exact]
* **Product solutions.** Einstein–Maxwell–$\Lambda$ admits $M_2\times S^2$ with $1/a^2=\Lambda+Q^2/a^4$. The two branches are AdS$_2$ (the
  throat) and dS$_2$ (charged Nariai).
  * They exist iff $\Lambda Q^2\le\frac14$.
  * At $\frac14$ they merge into flat$\times S^2$ with $a^2=2Q^2$.
  * Beyond it no static throat exists. These identities are checked with sympy.
* **Growth.** A T-region with $\dot r^2=-f_C$, where $f_C=1-C/r+Q^2/r^2-\Lambda r^2/3$, grows without bound from $r_0$ iff
  $$C>C_N:=\max_{r\ge r_0}\Big(r+\frac{Q^2}r-\frac{\Lambda r^3}3\Big).$$
  * For $Q=0$ this is $M>1/(3\sqrt\Lambda)$, the Nariai mass.
  * Late expansion is accelerating ($\ddot r>0$): inflation, distinct from the decelerating, recollapsing Λ = 0 flip of GSL-5.
* **Extremal throat ($r_0=Q$).**

  | how the vacuum energy appears | growth iff |
  |---|---|
  | injected from outside ($C$ stays $2Q$) | $\Lambda Q^2>3/16$ exactly |
  | converted locally at fixed Misner–Sharp mass ($C=2Q-\Lambda Q^3/3$) | $\Lambda Q^2>0.22121$ |
  | static throat in a fixed vacuum | never: $C=g(a)<C_N$, because $g$ has a local minimum at the AdS$_2$ radius and its maximum at the Nariai radius |

* **The universal law.** For a given vacuum energy, the masses are
  * $M_N=0.0665\,m_P^3/\sqrt\rho$ (Nariai);
  * $0.0864\,m_P^3/\sqrt\rho$ (throat bud, $3/16$);
  * $0.155\,m_P^3/\sqrt\rho$ (thin-wall classical route, G4).

  The required density at fixed mass changes by at most a factor 5.4 between routes.

Failed attempt (F22): "switching on" Λ in a throat at rest is not valid initial data. For a homogeneous region the radial null constraint
reads $\ddot r=(\dot a/a)\dot r$, so marginal data ($\dot r=0$) cannot accelerate. The correct analysis is the RN–dS T-region above.

### Proposition GSL-9 (Cauchy-horizon route near extremality) — [Proposition-N + exact edge]
For the charged bubble of Proposition G4:
* The outer flip of $A_{\rm out}$ solves $(\chi^2+\kappa^2)R^4=2MR-Q^2$. It must lie between the horizons.
* At $Q=qM$ this confines $M$ to $[1/\sqrt L(1+s),\,1/\sqrt L(1-s)]$, with $L=\chi^2+\kappa^2$ and $s=\sqrt{1-q^2}$. The upper edge is
  attained to $10^{-5}$.
* The lower edge is raised slightly by the barrier condition. So the window width satisfies $0.87\le\text{width}/(r_+-r_-)\le1$, with
  $r_+-r_-=2\kappa_-r_-^2$.
* As $q\to1$ the window collapses onto $M=1/\sqrt{\chi^2+\kappa^2}$.

| $1-q$ | window $M\chi$ ($\kappa=0.3\chi$) | $\kappa_-r_-$ | window $M\chi$ ($\kappa=\chi$) |
|---|---|---|---|
| $10^{-1}$ | 0.674–1.698 | 0.773 | 0.523–1.254 |
| $10^{-2}$ | 0.843–1.115 | 0.164 | 0.638–0.823 |
| $10^{-3}$ | 0.918–1.003 | 0.047 | 0.685–0.740 |

**Consequence.**
* Mass inflation grows like $e^{\kappa_-v}$. A route with a weak inner-horizon blueshift needs $M$ tuned to fractional precision
  $\sim\kappa_-r_-$.
* Even at exact extremality the Cauchy horizon is extendible only weakly, in $C^{0,1/2}\cap W^{1,2}_{\rm loc}$, and **non-uniquely**
  (Gajic–Luk, arXiv:1709.09137) [Literature].
* The quantum stress tensor diverges as $C/V^2$ at Reissner–Nordström–de Sitter Cauchy horizons, generically (Hollands–Wald–Zahn,
  arXiv:1912.06047) [Literature].
* E2 is GSL-compatible, but it gives up determinism.

### Proposition GSL-10 (deep throats hide what happens in them) — [Proposition, quasi-static]
* In extremal Reissner–Nordström, depth $\rho$ below the throat mouth ($r=2Q$) corresponds to $r-Q\simeq Qe^{1-\rho/Q}$.
* Energy $\delta M$ that stays inside moves the event horizon to $r_+-Q=\delta M+\sqrt{2Q\delta M+\delta M^2}$.
* So a region at depth $\rho$ ends up inside the event horizon once $\delta M>\frac Q2e^{2-2\rho/Q}$.
* A bud that meets GSL-8 carries $\gtrsim Q/32$, the injected vacuum energy. **Any growing bud deeper than $\rho\approx2\,Q$ is hidden** (`hiding_depth`: 1.98 $Q$ exactly; the logarithmic form, valid for $\delta M\ll Q$, gives 2.39 $Q$).
  Visible bubbles (G3, G7) are avoided automatically.

### Remark GSL-11 (the electroweak corona: a near miss) — [Assessment with explicit approximations]
* Magnetic black holes with $eB>m^2$ at the horizon carry a corona of symmetric-phase vacuum (Maldacena, arXiv:2004.06084). This is a
  *natural* false vacuum inside a throat.
* Self-consistency with the AdS$_2$ radius gives $\Lambda Q^2\le K\,y(\Lambda Q^2)^2$, with $y=Q^2/a^2$ and
  $K=8\pi\alpha\rho_{\rm EW}/m_{\rm th}^4$.

| threshold $m_{\rm th}$ | $K$ | $\Lambda Q^2$ max |
|---|---|---|
| $m_H$ (full restoration) | $\pi\alpha/2\lambda=0.089$ | 0.075 |
| $m_W$, with the full $\rho_{\rm EW}$ (an overestimate) | 0.52 | 0.226 |

* Both stay below $1/4$, so the Standard Model throat always exists.
* A corona that is part of the static solution never buds (GSL-8). The near miss is a coincidence of order one: corona size and
  decompactification size both scale as $1/\sqrt{G\rho_{\rm EW}}$.

## 6. The routes of the brief, answered

**Route H** (research question: can regular, asymptotically flat data evolve into an expanding interior with a persistent, nontraversable
neck without violating the GSL?).
* **Answer: the GSL permits it exactly when the interior is transient for complete observers. It then *requires* the neck to persist**
  (GSL-3).
* Formation needs NEC violation on ingoing rays (G1), priced by GSL-4: impossible in weakly gravitating regions, possible in near-extremal
  throats (GSL-6, GSL-7).
* "Inflating" is distinct from ordinary expansion: the growing super-Nariai T-region accelerates and has anti-trapped spheres far outside its
  Hubble radius, while the Λ = 0 flip decelerates and recollapses within $c/6|\Theta_0|$.
* **The full metric is exact piecewise:**
  * near-extremal Reissner–Nordström exterior;
  * outgoing Bonnor–Vaidya–RN flip region;
  * RN–dS T-region $ds^2=-dr^2/(-f_C)+(-f_C)dx^2+r^2d\Omega^2$ for the growing daughter.
* Causal horizons:
  * the parent's event horizon has non-decreasing area (no flux crosses it);
  * the daughter observer's horizon stays non-compact while the daughter is in $D(\Sigma_0)$.
* *Open:* the junction dynamics at the bud's edges (P-R6f).

**Route S** (which premise of Wall's theorem fails?).
* For transient connected daughters, **none fails**: the theorem does not speak to them.
* For eternal daughters the premise that fails must be global hyperbolicity (E2) or semiclassicality (E3). A modified entropy (E4) does not
  help.

**Route Q.**
* Its semiclassical limit is Section 4, where QFC is exact and the duration bound holds.
* Its genuinely quantum regime $T<T_{\rm brk}$ is escape E3. There the semiclassical GSL and QNEC analysis loses validity, and a solvable
  quantum model is required.
  * SYK and JT baby-universe constructions (e.g. Sasieta–Swingle–Vilar López) are the natural laboratory for it.
  * In those constructions a "daughter" is a state of a closed universe, not an inflating region of a classical spacetime.
* *Open:* P-R6h.

**Route D — [Assessment].**
* Classical NEC violation from non-minimal couplings (Horndeski and beyond) is not QNEC-bounded, so GSL-4 does not apply. A genesis seed can
  open a throat classically.
* GSL-3 permits its daughters because genesis ends, i.e. they are transient. The appropriate entropy is the theory's dynamical entropy; its
  second law is established at linear order by Wall's method.
* Established facts:
  * Stable static wormholes exist beyond Horndeski and not in Horndeski (Franciolini–Hui–Penco–Santoni–Trincherini, arXiv:1811.05481;
    Mironov–Rubakov–Volkova, arXiv:1812.07022, 2212.05969).
  * Genesis with GR asymptotics exists (arXiv:1905.06249).
* Unresolved:
  * stability of a *time-dependent, inhomogeneous* throat;
  * strong coupling (Ageeva–Petrov–Rubakov 2020 for Horndeski genesis);
  * tree-level unitarity, which requires new physics below the cutoff (Cai et al., arXiv:2207.11772);
  * the GSL when $c_s\neq1$ near a horizon (see below).
* A pre-existing extremal throat would supply the neck that Proposition G7 demands, at budget $\propto f$.

**Direct GSL violation — [Literature + Assessment, kept separate].**
* Theories with several propagation speeds give different black-hole horizons to different fields, and admit energy extraction that lowers
  the generalised entropy:
  * Dubovsky–Sibiryakov (2006), via the Hawking effect;
  * Eling–Foster–Jacobson–Wall (hep-th/0702124), via a classical Penrose-type process.
* Benkel et al. (arXiv:1803.01624) show that the Penrose-splitting extraction fails in every theory where gravity is attractive.
* This is the only concrete GSL-violating class known, and it is contested. Beyond-Horndeski genesis with $c_s\ne1$ around a horizon could
  fall into it. That is a liability for Route D, not a resource.

## 7. Failed approaches of this round (detail in FAILED_APPROACHES.md)
* **F21.** "Route G must violate the GSL" (Round 6): withdrawn. Minimal-coupling classical NEC violation conflicts with black-hole
  thermodynamics as a theory, but genesis does not need to violate the GSL.
* **F22.** Sudden vacuum energy in a throat at rest: not valid initial data (KS constraint). Replaced by the T-region criterion.
* **F23.** A cheap QNEC flip of a static throat as a daughter seed: never grows ($C=g(a)<C_N$). It needs energy injection or a change of vacuum.
* **F24.** The electroweak corona as a natural false vacuum: a near miss ($\Lambda Q^2\le0.226<1/4$).
* **F25.** Extremality as a cure for the Cauchy-horizon route: the window closes as $\kappa_-\to0$ (GSL-9).

## 8. Open problems (OPEN_PROBLEMS.md)
* **P-R6f.** 1+1 Einstein–Maxwell–scalar evolution of a bud in an extremal throat: the edges, hiding (GSL-10) and persistence of the neck.
* **P-R6g.** QNEC and its strong 2D form for Landau-level fermions in a curved magnetic throat at large $q$ (MMP regime), with massive-mode
  corrections.
* **P-R6h.** The sub-$T_{\rm brk}$ regime: a Schwarzian/SYK treatment of anti-trapping and of the "bud" as a quantum state (escape E3).
* **P-R6i.** Theorem GSL-3 without spherical symmetry: the general version needs a mean-convexity comparison for the horizon's cross-section
  on normal data.
* **P-R6a (updated).** Beyond-Horndeski: dynamical-entropy second law around time-dependent scalar backgrounds with $c_s\ne1$.

## 9. Assessment

* **Problem B, GSL form.**
  * The GSL is respected by two daughter classes: transient connected daughters and daughters beyond a Cauchy horizon.
  * It forbids exactly the eternal, laboratory-determined daughter.
  * A good answer to "how to get around the GSL" is therefore: **do not fight it; make the daughter transient and keep it connected.** An
    extremal magnetic throat does both.
  * There the QNEC cost of anti-trapping falls by $\sqrt f$, and becomes payable with $c\sim N_fq$ light species when
    $T<(c/12)T_{\rm brk}/\sinh(\rho/Q)$.
  * Deep buds are automatically hidden.
* **What is not solved.**
  * The energy density: $\Lambda_{\rm loc}Q^2>3/16$ inside the throat.
  * Creating that false vacuum (G3/G5).
  * The existence of magnetic monopoles.
* The obstacle has moved from the second law to the vacuum-energy requirement $\rho\gtrsim m_P^6/M^2$. That requirement is the same for every
  route found so far.
