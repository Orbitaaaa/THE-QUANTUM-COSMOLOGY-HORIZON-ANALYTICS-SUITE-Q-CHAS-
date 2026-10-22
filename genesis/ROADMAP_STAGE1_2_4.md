# Roadmap Stages 1, 2 and 4: the mathematics, the quantum sector, and the search for another mechanism

**Programme (from the Roadmap to a Real Machine, v1.0):** close the open problems of Stage 1, settle the quantum-sector and
boundary-condition questions of Stage 2, and run the Stage 4 alternative-mechanism search in parallel. Stages 3 and 5–10
(machine, experiment, verification) are treated **theoretically only**, as derived requirements (Part F). **No experimental
or technological proposal is made or implied.**

Status labels as before: **[Theorem]** proved under stated hypotheses; **[Proposition-N]** established numerically with
stated coverage; **[Conjecture]**; **[Interpretation]**.

Code: `CODE/babyuniverse/bounds.py`, `stokes_exact.py`, `daughter.py`, `alternatives.py`; script
`CODE/scripts/roadmap_analysis.py`; data `DATA/roadmap_stage1_2_4.json`; figure `FIGURES/roadmap_stage1_2_4.png`;
tests `CODE/tests/test_bounds.py` (11), `test_stokes_exact.py` (4), `test_daughter.py` (4), `test_alternatives.py` (6); full suite 131.

---

## Summary of results

1. **A universal law for the laboratory exponent [Theorem, thin barrier].** For walls light compared with the vacuum
   energy ($\kappa\ll\chi$, the physically relevant case) the FGG exponent of every laboratory transition is
   $$B_{\rm lab}=B_{\rm CdL}\,\Phi\!\big((R_{\rm CdL}/R_M)^2\big)\,\big(1+O(\kappa^2/\chi^2)\big),\qquad
   \Phi(E)=\frac{16}{\pi E^2}\int z^2\sqrt{E-\frac{(z^3-1)^2}{z^4}}\,dz\ \ge\ \max\Big(1,\ \frac{8}{3E}+\frac23\Big),$$
   with $R_M=(3M/4\pi\rho)^{1/3}$ the radius of a false-vacuum ball of mass $M$. For bubbles larger than the CdL bubble,
   $B_{\rm lab}\simeq4\pi^2\sigma^2R_M^2/\rho$, a formula from which **Newton's constant drops out**. Checked to $6.5\times10^{-6}$; the lower bound on $\Phi$ is proved.
2. **A parameter-free frontier [Corollary].** Every thin-barrier laboratory transition obeys
   $B_{\rm lab}\ge\tfrac43\sqrt{\pi/2}\,\sqrt{B_{\rm CdL}\,(M/m_P)\,(R_M/\ell_P)}$. With metastability $B_{\rm CdL}\ge10$, a black hole of at
   least one Planck mass and a bubble of at least $100\,\ell_P$ (the least one can ask of a semiclassical description), $P\le10^{-23}$.
   This is reached only at $E_{\rm vac}\approx3\times10^{17}$ GeV.
3. **Where the FGG ambiguity lives [Theorem].** A laboratory transition is double-flip (B_an = B_sq, no ambiguity) **iff**
   $\kappa<\chi$ and $M<M_{\rm sd}=(\chi^2-\kappa^2)/(2G\chi^3)$. The boundary is the surface where the bounce turning point sits on
   the false-vacuum horizon, $R_2=1/\chi$. For light walls the ambiguous sliver is a fraction $\simeq\kappa^2/2\chi^2$ of the mass range.
4. **T5 is now mostly a theorem.** Rigorous bounds: $B_{\rm an}\ge B_{\rm sq}\ge$ the sine bound $\ge$ an elementary closed form that is
   asymptotically sharp for light walls. An explicit criterion for $\Delta t^W_E<0$ gives $B(M)\ge B_{\rm CdL}$ in the double-flip class.
5. **P-Q2 closed [Proposition-N, exact numerics].** Exact integration of the momentum-representation squared constraint at $M=0$,
   $\kappa>\chi$, gives the Stokes multiplier $\sqrt2$ in the limit $\hbar\to0$ for three orderings and four parameter sets. The
   corrections are $h(A\ln h+B)$. Controls with $\nu=0$ and $\nu=-1$ give 0 and $\sqrt{2\pi}(2r_q/a^2h)^{1/2}$ as predicted.
6. **P-Q9 reformulated [Theorem + Corollary].** In the full theory the orientation is configuration data (boost-invariant sign
   of $R'/\Lambda$). Wheeler–DeWitt squaring removes the time orientation, not the orientation: this corrects the Phase III
   heuristic (F14). The two families of the squared constraint are two phase conventions for the *same* daughter. They carry
   different conserved Casoratian charges. Which one is "probability" is the entire content of the B_an vs B_sq question.
7. **P-Q10 solved: no boundary condition is needed at the daughter edge [Theorem].** The finite-time edge of Theorem V-1 is a
   coordinate singularity of the static-frame momentum. $\tilde P=P\mp(R/2G)\ln|f_{\rm in}|$ is continuous, and the constraint and the true
   Hamiltonian are analytic across the interior horizon (F16 corrects V-1(ii)).
8. **The parent's clock folds back [Theorem].** Along every daughter, $dt_\infty/d\tau=A_{\rm out}/f_{\rm out}<0$, because the daughter wall lies
   in region III, where the Killing field is past-directed. Unitary evolution in the boundary clock therefore takes laboratory
   data to the *time-reversed* daughter. In that clock the physically expanding daughter is a negative-flux (Klein-type)
   channel. This is the minisuperspace face of Freivogel et al. (no unitary boundary process creates an inflating region) and of Wall's
   GSL theorem (no independent baby universes). It is also the Lorentzian image of the Euclidean fold of Theorem III-4.
9. **Stage 4: no mechanism found that evades the suppression within GR + quantum fields obeying QEIs.**
   * A seed black hole lowers $B$ by at most about 30%, and its data already contain a de Sitter equator (F15).
   * QEI-bounded NEC violation needs about $\tfrac83S_{\rm dS}$ species, which pushes the species cutoff below the daughter's Hubble rate.
   * Closed or de Sitter parents need horizon-scale bubbles; upward tunnelling from our vacuum costs about $e^{-3\times10^{122}}$.
   * The unsuppressed routes are cosmological (bubbles nucleated during inflation, GVZ), hypothetical (limiting curvature), or
     run beyond Cauchy horizons (charged and monopole exteriors).

---

## Part A — Analytic suppression bounds (Stage 4 task 30; Stage 1 P-T5)

Laboratory transitions: Schwarzschild exterior ($\eta=0$), regime A. Recall $C=(f_{\rm in}+f_{\rm out}-\kappa^2R^2)/2\sqrt{f_{\rm in}f_{\rm out}}$,
$\theta=\arccos C\in(0,\pi)$ (analytic branch), $\sin\theta=\kappa R\sqrt{V/f_{\rm in}f_{\rm out}}$ (proof of Theorem III-4),
$B_{\rm an}=\frac2G\int_{R_1}^{R_2}R\theta\,dR$ and $B_{\rm sq}=\frac2G\int R\arccos|C|\,dR$. Dimensionless variables: $R=sz$,
$s=(2GM/L^2)^{1/3}=R_0/\sqrt E$ with $R_0=2\kappa/L^2$, $V=1-g(z)/E$, and

$$g(z)=z^2-\frac{2\gamma}z+\frac1{z^4}=\frac{(z^3-1)^2}{z^4}+\frac{2(1-\gamma)}z,\qquad 1-\gamma=\frac{2\kappa^2}{L^2},\qquad E_{AB}=2-2\gamma .$$

### A1. Lemma S1 (sine bound) — [Theorem]
$$B_{\rm an}\ \ge\ B_{\rm sq}\ \ge\ B_{\rm mem}\equiv\frac{2\kappa}G\int_{R_1}^{R_2}R^2\sqrt{\frac V{f_{\rm in}f_{\rm out}}}\,dR\ \ge\ \frac{2\kappa}G\int_{R_1}^{R_2}R^2\sqrt V\,dR .$$
*Proof.* $\theta\ge\sin\theta$ on $[0,\pi]$. Also $\arccos|C|\ge\sqrt{1-C^2}=\sin\theta$. Finally $f\le1$ for $\eta=0$. $\square$
The last expression is the WKB exponent of a flat-space relativistic membrane of mass $4\pi\sigma R^2$ in the potential $V$.
Gravity can only increase the exponent relative to it.

### A2. Theorem S2 (elementary closed-form bound) — [Theorem]
$$B_{\rm an}\ \ge\ B_{\rm sq}\ \ge\ \frac{2\pi\kappa^2}{G}\,\frac{R_1^2\,f_{\rm out}(R_1)}{\big(L^2+4GM/R_1^3\big)\sqrt{f_{\rm in}(R_1)f_{\rm out}(R_2)}} .$$
*Proof.* $A_{\rm out}=GM/\kappa R^2-L^2R/2\kappa$ is strictly decreasing, with $|A_{\rm out}'|\le(L^2+4GM/R_1^3)/2\kappa$ on $[R_1,R_2]$.
$A_{\rm out}(R_1)=+\sqrt{f_{\rm out}(R_1)}$ (regime A) and $A_{\rm out}(R_2)=-\sqrt{f_{\rm out}(R_2)}\le-\sqrt{f_{\rm out}(R_1)}$ (Prop. G1; $f_{\rm out}$ increases),
so there is $R'\in(R_1,R_2]$ with $A_{\rm out}(R')=-\sqrt{f_1}$, $f_1=f_{\rm out}(R_1)$. On $[R_1,R']$:
* $R^2\ge R_1^2$;
* $V=f_{\rm out}-A_{\rm out}^2\ge f_1-A_{\rm out}^2$;
* $f_{\rm in}f_{\rm out}\le f_{\rm in}(R_1)f_{\rm out}(R_2)$.

Change variables to $w=A_{\rm out}$: $\int\sqrt{f_1-w^2}\,dR\ge(\pi f_1/2)/\max|A'_{\rm out}|$. Apply S1. $\square$

On the 54-model grid $B_{\rm sq}/B_{\rm S2}\ge1.95$. S2 is weak for heavy walls but **asymptotically sharp for light walls**:
$B_{\rm S2}/B=0.98$ at $\kappa/\chi=0.003$.

### A3. Theorem S3 (thin-barrier law) — [Theorem at fixed $E$, including the bound on $\Phi$; Proposition-N uniformly up to $M_{\rm cr}$]
Let $\epsilon=\kappa/\chi\to0$ at fixed $E$. Then
$$\frac{B_{\rm an}(M)}{B(0)}=\Phi(E)\,\big(1+O(\epsilon^2)\big),\qquad\Phi(E)=\frac{16}{\pi E^2}\int_{\{g_1<E\}}z^2\sqrt{E-g_1(z)}\,dz,\quad g_1=\frac{(z^3-1)^2}{z^4},$$
and $B(0)=B_{\rm CdL}=\pi\kappa R_0^3/8G\,(1+O(\epsilon^2))$ (T1, $q<0$).

*Proof.* $\kappa s=2\kappa^2/L^2\sqrt E=O(\epsilon^2)$. So $\sin\theta=\kappa sz\sqrt{1-g/E}/\sqrt{f_{\rm in}f_{\rm out}}=O(\epsilon^2)$ and
$\theta=\sin\theta\,(1+O(\epsilon^4))$. Next, $f_{\rm in}=1-\chi^2R_0^2z^2/E$ and $f_{\rm out}=1-L^2s^2/z$ are $1-O(\epsilon^2)$ on the barrier, and
$g=g_1+O(\epsilon^2)$. The turning-point perturbation changes $\int\sqrt{E-g}$ at relative order $\epsilon^2$, because $\partial_E\sqrt{E-g}$ is integrable.
Hence $B_{\rm an}=(2\kappa s^3/G\sqrt E)\int z^2\sqrt{E-g_1}\,(1+O(\epsilon^2))$ with $\kappa s^3/\sqrt E=\kappa R_0^3/E^2$. $\square$

**Uniformity.** As $E\to0$ at fixed $\epsilon$ the barrier narrows to $|R-R_f|<\tfrac23(\kappa/L^2)\sqrt F$, with $F=1-L^2R_f^2$. There
$\int\sqrt V\,dR=\pi\kappa F/3L^2$, while $1/\sqrt{f_{\rm in}f_{\rm out}}=1/F$. **The factors $F$ cancel**, so
$B\to2\pi\kappa^2R_f^2/3GL^2=\tfrac{8}{3E}B(0)$ up to $M_{\rm cr}$, where the naive expansion parameter $2(1-\gamma)/E$ is not small. Numerically the
relative error of $B=B(0)\Phi(E)$ is at most $6.5\times10^{-6}$ for $\kappa/\chi=0.003$ over $M/M_{\rm cr}\in[10^{-4},0.9]$. It scales as $\epsilon^2$.

**Properties of $\Phi$** (substitution $s=(z^3-1)/z^2$, monotone, $g_1=s^2$, $z^2dz=h(s)ds$, $h=z^5/(z^3+2)=\tfrac13(1+\tfrac43s+s^2)+O(s^3)$):
* $\Phi(E)=\dfrac{8}{3E}+\dfrac23+O(E)$ as $E\to0$;
* $\Phi(E)=1+\dfrac{16}{\pi}E^{-3/2}+o(E^{-3/2})$ as $E\to\infty$. The $-2/z$ term of $g$ contributes $\int z\,dz/\sqrt{E-z^2}=\sqrt E$. Hence
  $dB/dM|_0=2R_0=-\Delta t^W_E$, as the first law T2 requires.
* **$\Phi(E)\ge\max\big(1,\ \frac{8}{3E}+\frac23\big)$ for all $E$ [Theorem].** *Proof.* Write $h=\tfrac13\,d(z^3)/ds$ and let $x=z(s)>1$, $y=z(-s)<1$. Then
  $G(s)=z(s)^3-z(-s)^3=s(x^2+y^2)$. Integration by parts with the even weight $\sqrt{E-s^2}$ gives
  $\int(h-\tfrac13-\tfrac{s^2}3)\sqrt{E-s^2}\,ds=\tfrac13\int_0^{\sqrt E}K(s)\,s\,ds/\sqrt{E-s^2}$ with $K=G-2s-\tfrac23s^3$, so it suffices that
  $Q\equiv x^2+y^2\ge2+\tfrac23s^2$.
  From $x-1/x^2=s=1/y^2-y$ one finds, with $P=xy$:
  * $Q^2=P^4Q+2P^5$;
  * $Q\ge2P$, which forces $P\ge1$;
  * $s^2=(Q-2P)(Q+P)^2/Q^2$.
  The claim becomes $f=Q^3-6Q^2+6P^2Q+4P^3=c_1Q+c_0\ge0$, with $c_1=P^8+2P^5-6P^4+6P^2>0$ on $P\ge1$. Where $c_0<0$, squaring reduces it to
  $2c_1^2P^5-c_0^2-c_0c_1P^4=8P^6(P-1)^2(2P^7+7P^6+12P^5+20P^4+22P^3+6P^2-4P-2)\ge0$. The cofactor's only positive root is $0.428$. For $\Phi>1$ the
  same argument works with $K_2=G-s^3=s(x^2+y^2-s^2)>0$, since $0<s<x$. $\square$ (Symbolic certificate: `bounds.phi_bound_certificate`.)

**Physical form.** In the thin-barrier limit $R_f^3=2GM/L^2\to3M/4\pi\rho$, the radius $R_M$ of a ball of false vacuum of mass $M$. Also
$E=(R_0/R_f)^2=(R_{\rm CdL}/R_M)^2$ with $R_{\rm CdL}=3\sigma/\rho$. So
$$\boxed{\;B_{\rm lab}=B_{\rm CdL}\,\Phi\big(R_{\rm CdL}^2/R_M^2\big)\ \ge\ \frac{8}{3}\frac{R_M^2}{R_{\rm CdL}^2}B_{\rm CdL}=\frac{4\pi^2\sigma^2R_M^2}{\rho}\;}$$
**Newton's constant drops out.** The price of the transition is the flat-space membrane exponent; gravity only supplies the final
state, the $\varepsilon_{\rm out}=-1$ branch that has no counterpart in flat space. Equivalently
$B_{\rm lab}=4\pi\sqrt{2/27}\,\sqrt{B_{\rm CdL}}\,(3M/4\pi)^{2/3}\rho^{-1/6}$ in Planck units. This reproduces:
* the Phase VI optimum ($E_{\rm vac}=9\times10^{17}$ GeV, $\sigma^{1/3}=0.8\rho^{1/4}$, $M=10m_P$): 105 by the law, 107 exact;
* the T6 entry $4.4\times10^6$ at $0.5M_{\rm cr}$ for GUT parameters.

Why lab bubbles are harder to convert than CdL bubbles: a bubble of mass $M$ has radius $R_M\gg R_{\rm CdL}$, and the exponent grows like the area ratio.

### A4. Corollary S7 (the semiclassical frontier) — [Theorem, given S3]
Eliminate $\sigma$ and $\rho$ between $B_{\rm lab}\ge4\pi^2\sigma^2R_M^2/\rho$, $B_{\rm CdL}=27\pi^2\sigma^4/2\rho^3$ and $M=\tfrac{4\pi}3\rho R_M^3$:
$$\boxed{\;B_{\rm lab}\ \ge\ \tfrac43\sqrt{\pi/2}\ \sqrt{B_{\rm CdL}\cdot\frac{M}{m_P}\cdot\frac{R_M}{\ell_P}}\;}\qquad(\tfrac43\sqrt{\pi/2}=1.671).$$
No free parameters remain. Requiring (i) a false vacuum metastable at least at the level where the CdL formula means something,
(ii) a black hole of at least $M_{\min}$, and (iii) a bubble large enough to be semiclassical, the bound is attained at
$\rho_*=3M_{\min}/4\pi R_{\min}^3$:

| $B_{\rm CdL}\ge$ | $M\ge$ | $R_M\ge$ | $B_{\min}$ | $P_{\max}$ | $E_{\rm vac}$ at optimum | $\kappa/\chi$ |
|---|---|---|---|---|---|---|
| 10 | $m_P$ | $100\,\ell_P$ | 52.8 | $10^{-22.9}$ | $2.7\times10^{17}$ GeV | 0.050 |
| 10 | $m_P$ | $1000\,\ell_P$ | 167 | $10^{-72.6}$ | $4.8\times10^{16}$ GeV | 0.009 |
| 10 | $10m_P$ | $100\,\ell_P$ | 167 | $10^{-72.6}$ | $4.8\times10^{17}$ GeV | 0.089 |
| 100 | $m_P$ | $100\,\ell_P$ | 167 | $10^{-72.6}$ | $2.7\times10^{17}$ GeV | 0.089 |
| 100 | $10m_P$ | $1000\,\ell_P$ | 1671 | $10^{-726}$ | $8.5\times10^{16}$ GeV | 0.028 |

All optima have $\kappa/\chi\le0.16$, inside the validity of S3. **Check against the exact thin-wall exponent** (`observability.laboratory_case` at the optimum):

| case | frontier bound | exact $B_{\rm an}=B_{\rm sq}$ | gravitational $B_{\rm CdL}$ | $R_1$ | WKB $\epsilon$ |
|---|---|---|---|---|---|
| $(10, m_P, 100\ell_P)$ | 52.8 | 59.5 | 9.95 | $81\ell_P$ | $1.3\times10^{-5}$ |
| $(10, m_P, 1000\ell_P)$ | 167.1 | 173.8 | 10.0 | $883\ell_P$ | $1.2\times10^{-8}$ |
| $(10, 10m_P, 100\ell_P)$ | 167.1 | 171.9 | 9.84 | $89\ell_P$ | $1.4\times10^{-4}$ |

The excess over the bound is almost exactly $\tfrac23B_{\rm CdL}$, the next term of $\Phi$. All three are double-flip, so there is no FGG ambiguity.

**What sets the frontier (validity remarks).**
* *Quantum backreaction* is negligible there. The conformal-anomaly energy density in the false vacuum is $N H^4/960\pi^2$; relative to $\rho$ it is $NH^2/360\pi\approx2\times10^{-9}N$ at the
  first optimum.
* *Semiclassical geometry* is the binding constraint. That is why $R_{\min}$ is a parameter; the bound scales as $\sqrt{R_{\min}}$.
* *The thin-wall EFT* must also hold. A scalar wall has thickness $\delta\sim v^2/\sigma$ (vev $v$, mass $m\sim\sigma/v^2$), and $\delta\ll R_{\rm CdL}=3\sigma/\rho$ requires
  $v\ll\sqrt3\,\sigma/\sqrt\rho=m_P(\kappa/\chi)/\sqrt{2\pi}$. At the first optimum this is $v\ll0.02\,m_P$, with scalar mass $m\gtrsim0.014\,m_P$. Light walls
  ($\kappa\ll\chi$) are thin only for small vevs and near-Planckian scalar masses at these scales.
* So the controlled region ends close to where Corollary S7 puts it. Pushing $R_{\min}$ to a few $\ell_P$ would allow $P\sim e^{-10}$, but only in a regime where none of
  these approximations hold. **There is no semiclassically controlled laboratory transition with
$P\gtrsim10^{-23}$.** Higher probabilities require $R_M\lesssim20\,\ell_P$ or sub-Planckian black holes, where the theory used here is not trustworthy.

### A5. Theorem S4 (flip class; where the FGG ambiguity lives) — [Theorem]
For laboratory transitions ($\eta=0$, regime A):
$$\text{double flip}\iff\kappa<\chi\ \text{and}\ M<M_{\rm sd}=\frac{\chi^2-\kappa^2}{2G\chi^3}\iff E>E_{\rm sd}(\gamma)=\frac{1-\gamma^2}{\gamma^{2/3}} .$$
*Proof.* $C(R_2)=\varepsilon_{\rm in}\varepsilon_{\rm out}(R_2)=-\varepsilon_{\rm in}(R_2)$ (Prop. G1). $A_{\rm in}=(2GM/R-(\chi^2-\kappa^2)R^2)/2\kappa R$ changes sign once, at
$R_{\rm fin}=(2GM/(\chi^2-\kappa^2))^{1/3}$, iff $\chi>\kappa$, and $A_{\rm in}<0$ beyond it.
* $A_{\rm in}(R_1)>A_{\rm out}(R_1)>0$ in regime A, so $R_{\rm fin}>R_1$.
* $V(R_{\rm fin})=f_{\rm in}(R_{\rm fin})$ (since $A_{\rm in}=0$), which is positive iff $R_{\rm fin}<1/\chi$, i.e. iff $M<M_{\rm sd}$.
* $V>0$ only on $(R_1,R_2)$, so $M<M_{\rm sd}$ gives $R_{\rm fin}\in(R_1,R_2)$: double flip.
* $M>M_{\rm sd}$ gives $R_{\rm fin}\notin[R_1,R_2]$, so $R_{\rm fin}>R_2$: single flip.

At $M=M_{\rm sd}$, $R_2=R_{\rm fin}=1/\chi$. $\square$

**Geometric meaning.** The transition is unambiguous exactly when the daughter is born containing more than a de Sitter hemisphere of false
vacuum. Checked against the numerical flip class on all 54 grid models (`flip_class_theorem_agreement: true`). For light walls,
$(M_{\rm cr}-M_{\rm sd})/M_{\rm cr}\simeq\kappa^2/2\chi^2$. **The B_an/B_sq question (P-Q9) is moot for every transition in the physically relevant
regime except a sliver of relative width $\kappa^2/2\chi^2$ next to the critical mass.**

### A6. Theorem S5 (sign of the Euclidean wall time) — [Theorem; threshold location Proposition-N]
With $b=E_{AB}$, $W=E-b/z=Ef_{\rm out}$ and $s=(z^3-1)/z^2$, the exterior Euclidean Killing time along the wall is a positive multiple of
$$J=\int_{z_1}^{z_2}\frac{(1-z^3)\,dz}{z^2f_{\rm out}\sqrt{E-g}}=-E\int_{s_1}^{s_2}\frac{s\,q(z)}{W\sqrt{W-s^2}}\,ds,\qquad q=\frac{z^3}{z^3+2}.$$
Integrate by parts with $s/\sqrt{W-s^2}=-\partial_s\sqrt{W-s^2}+W'/2\sqrt{W-s^2}$; the boundary terms vanish at the turning points. This gives
$$-J/E=\int\frac{q'(W-s^2)W+qW'(s^2-W/2)}{W^2\sqrt{W-s^2}}\,ds,$$
with $q'=6z^5/(z^3+2)^3$ and $W'=bz/(z^3+2)$. So **$\Delta t^W_E<0$ whenever $D=6z(W-s^2)W-b(z^3+2)(W/2-s^2)\ge0$ on the barrier**. A
sufficient elementary condition is
$$E/E_{AB}\ \ge\ \max\big(\omega(z_1),\omega(z_2)\big),\qquad\omega(z)=\frac{z^2}{12}+\frac7{6z}.$$
It follows from $D\ge u(6zW-b(z^3+2))+bW(z^3+2)/2$ with $u=W-s^2\le W$.

Thresholds (bisection, verified above on a grid up to $10^6E_{AB}$):

| $\gamma$ | −0.99 | −0.5 | 0 | 0.2 | 0.4 | 0.6 | 0.8 | 0.99 |
|---|---|---|---|---|---|---|---|---|
| $\omega$-threshold $E/E_{AB}$ | 1.605 | 1.556 | 1.498 | 1.471 | 1.440 | 1.404 | 1.358 | 1.274 |
| $E_{\rm sd}/E_{AB}$ | — | — | — | 1.754 | 1.289 | 1.125 | 1.044 | 1.002 |

The pointwise $D$-threshold is lower, $1.25$–$1.33\,E_{AB}$, and also lies below $E_{\rm sd}$ for $0<\gamma\le0.4$. Across the full grid $J<0$
everywhere in regime A.

**Corollary (double-flip class: laboratory creation is never easier than vacuum decay).** In the double-flip class, $\mathcal A/4G=\pi/G\chi^2$
is constant in $M$, and the first law gives $dB/dM=-\Delta t^W_E$. Wherever the criterion holds along $[0,M]$, therefore
$$B(M)\ \ge\ B(0)=B_{\rm CdL}.$$
This covers the whole double-flip class for $0<\gamma\le0.4$ (pointwise $D$), and all masses below $E\ge1.61E_{AB}$ for every $\gamma$ ($\omega$ criterion).
The thin-barrier law S3 covers $\gamma\to1$ uniformly.

### A7. Status of T5
| statement | before | now |
|---|---|---|
| $B\ge B_{\rm CdL}$, light walls ($\kappa\ll\chi$) | Prop-N | **Theorem** (S3), with the explicit law $\Phi$ |
| $B\ge B_{\rm CdL}$, double flip, $E\ge E_{\rm thr}(\gamma)$ | Prop-N | **Theorem** (S5 + first law), threshold Prop-N |
| $B\ge0.564\,B(0)$, all regime A | Prop-N | Prop-N. Infimum at the A/B boundary as $\gamma\to-1$ (0.5640). Rigorous but weak: S1, S2 |

---

## Part B — P-Q2: the exact Stokes multiplier (Stage 1 task 5)

**Setup.** $M=0$, $\eta=0$, $\kappa>\chi$ (single-flip class at $M=0$). With $\psi(u)=\int e^{-u\theta/h}\phi(\theta)\,d\theta$, $h=\hbar G$ ($P_u=i\theta/G$,
$\hat u\to h\,\partial_\theta$), the squared constraint $g(\hat u)\cosh^2(G\hat P_u)-N(\hat u)^2$ becomes, for the ordering family "weight $\lambda$ of $g$ to the left",
the exact ODE
$$a^2h^2\phi''-2h(a-\chi^2\cos^2\theta)\phi'+\big(\sin^2\theta-2\lambda\chi^2h\sin2\theta\big)\phi=0,\qquad a=\chi^2+\kappa^2 .$$
* The characteristic roots $u_{1,2}=[(a-\chi^2c^2)\pm c\,r(\theta)]/a^2$, $r=\sqrt{\kappa^4-\chi^4\sin^2\theta}$, cross at $\theta_q=\pi/2$ ($C=0$). This is the coalescence
  point, with the crossing slope $r_q=\sqrt{\kappa^4-\chi^4}$.
* The transport equation $2\delta_iA_i'+\tau_iA_i=0$ ($\delta_i=\pm c\,r$, $\tau_i=a^2u_i'-2\lambda\chi^2\sin2\theta$) gives the local exponents
  $A_{1,2}\sim|\theta-\theta_q|^{-1/2}$ **for every $\lambda$**, because $\sin2\theta_q=0$. Hence $\nu=-\tfrac12$ (cf. Prop. III-8a).
* The Laplace transform maps $\nu\to-1-\nu$ (u-space recessive ↔ θ-space dominant), so $\nu=-\tfrac12$ is self-dual. The θ-space multiplier
  therefore tests the u-space one.

**Exact computation.** Integrate the ODE (DOP853, rtol $10^{-13}$) from $\theta_q-1$, starting on the recessive WKB solution $\phi_1$, to $\theta_q+1$.
Project onto the dominant solution with a Wronskian against first-order WKB data. The WKB solutions are normalised by $S_i(\theta_q)=0$ and
$\tilde A_i(\theta_q)=1$, for which Weber theory gives $s=\frac{\sqrt{2\pi}}{\Gamma(-\nu)}(2r_q/a^2h)^{-\nu-1/2}\to\sqrt2$.

| $h$ | 0.01 | 0.005 | 0.0025 | 0.00125 | 0.000625 |
|---|---|---|---|---|---|
| $s/\sqrt2-1$, $g$ left ($\lambda=1$) | $+6.0\times10^{-3}$ | $+5.0\times10^{-4}$ | $-1.0\times10^{-3}$ | $-1.1\times10^{-3}$ | $-8.8\times10^{-4}$ |
| Weyl ($\lambda=\frac12$) | $+9.2\times10^{-3}$ | $+4.3\times10^{-3}$ | $+2.0\times10^{-3}$ | $+9.2\times10^{-4}$ | $+4.2\times10^{-4}$ |
| $g$ right ($\lambda=0$) | $-2.2\times10^{-2}$ | $-1.4\times10^{-2}$ | $-8.1\times10^{-3}$ | $-4.7\times10^{-3}$ | $-2.7\times10^{-3}$ |

* The errors fit $h(A\ln h+B)$ with $A=0.73$, $0.09$ and $0.73$ for $\lambda=1,\tfrac12,0$; all tend to 0.
* Window independence: half-windows 0.7, 1.0, 1.3 agree to $2\times10^{-6}$.
* Other parameter sets: $(\chi,\kappa)=(1,1.2),(1,2.5),(0.5,1)$ give $s/\sqrt2=1.0093$, $0.99968$, $0.99974$ at $h=0.0025$.
* Controls with an added constant subprincipal term: $\nu=0$ gives $s=-2.7\times10^{-4}$ (no switching). $\nu=-1$ gives $22.03$ against the predicted $21.90$.

**Conclusion.** The universal multiplier $\sqrt2$ of Prop. III-8 is confirmed by an exact solution of the full (non-local-model) momentum-space
equation. It includes all cubic and higher terms of the symbol, over $\theta\in(0.27,2.87)$, i.e. most of the Euclidean range. The first correction
is $O(h\ln h)$ and depends on the ordering. **P-Q2 is closed.**

---

## Part C — The quantum sector: what the B_an / B_sq question really is (P-Q9, P-Q8; Stage 2 tasks 13–14)

### C1. Proposition Q9.1 (orientation is configuration data) — [Theorem]
On any spacelike slice through the wall in the static domain, the derivative of the area radius along the slice and the slice's
$\theta\theta$-momentum are the boosted wall-frame quantities:
$$\Big(\frac{R'}\Lambda,\ \frac{GP_\Lambda}R\Big)=\big(A\cosh\beta+\dot R\sinh\beta,\ \dot R\cosh\beta+A\sinh\beta\big),\qquad\Big(\frac{R'}\Lambda\Big)^2-\Big(\frac{GP_\Lambda}R\Big)^2=f>0 .$$
The vector is spacelike, so **$\mathrm{sign}(R'/\Lambda)=\mathrm{sign}(A)=\varepsilon$ on every slice**. The orientation is a property of the spatial geometry:
whether $r$ increases or decreases as one leaves the wall. The sign of $P_\Lambda$ is not invariant; it is the time-orientation datum.
*Proof.* Decompose $e=\cosh\beta\,n+\sinh\beta\,U$ and $u=\cosh\beta\,U+\sinh\beta\,n$. $\square$ (`daughter.orientation_is_configuration_data`)

**Consequence (correction F14).** The Wheeler–DeWitt constraint is quadratic in $P_\Lambda,P_R$, so it does not distinguish the *time orientation*. It
does distinguish $\varepsilon$, which enters through $R'$. Phase III motivated the sector-blind quantization by "the sign of $K_{\theta\theta}$ is exactly the
sign that squaring removes". **That is wrong in the full theory.** Squaring the bulk constraint does not lead to the squared reduced constraint.

### C2. Proposition Q9.2 (two conventions, one daughter) — [Theorem]
In the Lorentzian daughter region ($C<0$), the real-angle states ($\cosh GP_u=-C$, $P_u$ real) and the complex-angle states ($\cosh GP_u=C$,
$GP_u=i\pi+$real) describe the **same classical states**. They are related by the non-unitary intertwiner $e^{\pi u/\hbar G}$ (Theorem III-2c).
The squared constraint is invariant under the intertwiner, so its solution space contains each daughter state twice, with relative weight $e^{\pm\pi u/\hbar G}$.
$B_{\rm sq}$ is the exponent for arriving in the real copy, $B_{\rm an}$ for arriving in the complex copy. "Superselection of the orientation" is therefore a misnomer.
The physical question is **which copy's normalisation is probability**.

### C3. Theorem Q9.3 (Casoratian charges) — [Theorem]
For $\psi(u+ih)+\psi(u-ih)=W\psi$ with real-analytic $W$, the Casoratian $K[\psi,\psi^\sharp]$ is $ih$-periodic (Lemma III-6a). It is therefore a single
analytic function $\Phi(q)$ of $q=e^{2\pi u/h}$ on the static domain. A wave with $n$ flipped normals in the complex convention,
$\psi=a\,e^{-n\pi u/h}e^{iS/\hbar}$ with real $S$, has
$$K=(-1)^n\,2|a|^2\sinh(GS')\,e^{-2\pi nu/h}:$$
it contributes **only to the Laurent coefficient $q^{-n}$**. Exact check for constant $W$: `daughter.charge_of_wave`, periodicity and scaling to
$10^{-14}$. Hence:
* laboratory waves ($n=0$) and the real copy of the daughter carry the charge $q^0$, the Lebesgue current;
* the complex copy of the single-flip daughter carries $q^{-1}$;
* a mixture of the two determinations $\pm i\pi$ of the Sorkin angle carries $q^0$ only through cross terms, with weights fixed modulo
  quasi-constants.

No measure that is the same analytic function on both sectors turns "laboratory flux = daughter flux" into an identity in the complex
convention. Fixing the transition probability is **equivalent to choosing the relative weight of the $q^0$ and $q^{-1}$ charges**, which is the
anchor $u_*$ of Lemma III-6. For the squared constraint (step $2ih$) the Casoratian is $2ih$-periodic, the laboratory and the real-copy daughter
share the $Q^0$ charge ($Q=e^{\pi u/h}$), and the unitarity relation $|I|^2=|R|^2+|T_{\rm real}|^2$ is canonical.

### C4. Corollary Q9.4 (two principles, two answers) — [Theorem, conditional on the principle]
* **(i) Probability is one sector-independent current, the one that is probability in the laboratory.** Then probability reaches the daughter only
  through the real copy. The sector-blind theory is the one in which this is canonical:
  $$P_{\rm daughter}=2e^{-B_{\rm sq}}(1+O(\hbar))\quad\text{(single flip)},\qquad e^{-B_{\rm an}}\quad\text{(double flip)} .$$
* **(ii) Probability is the covariant Euclidean rate**, i.e. the pseudomanifold action with joint terms (Prop. III-3). Then $e^{-B_{\rm an}}$, which needs the
  entropy-weighted daughter measure $e^{(A-A_*)/4G}$ anchored at $u_*=u_2$.

Principle (i) is the standard probability interpretation of quantum mechanics extended across sectors. Principle (ii) is the standard
gravitational one, in which horizon entropy weights are built into the measure. Phase III adopted (ii) as primary. **On the evidence of C3, (i) is
the only one of the two with an anchor-free unitarity relation in the reduced theory.** The question is therefore no longer "is the orientation
superselected?" but "does the inner product of quantum gravity carry Bekenstein–Hawking weights on flipped sectors?". By Theorem S4 it matters
only in a sliver of the physical parameter space. Phase VI's bounds already used $\min(B_{\rm an},B_{\rm sq})$, so they hold under either principle.

### C5. P-Q8 (where the complex angles come from) — [Proposition, partial]
In the full theory the bulk WKB phase obeys $\delta S/\delta\Lambda=P_\Lambda=\pm(R/G)\sqrt{(R'/\Lambda)^2-F}$. Configurations with $R'_+>0$ and
$R'_+<0$ are separated in superspace by the classically forbidden zone $|R'_+|<\sqrt F\Lambda$, where $P_\Lambda$ is imaginary. The Wheeler–DeWitt
operator is regular at $R'=0$, so the phase continues analytically through the zone. With $R'/\Lambda=\sqrt F\cos\Theta$ and
$GP_\Lambda/R=i\sqrt F\sin\Theta$, the angle runs continuously from 0 to $\pi$, which is the atan2 determination.
**The complex-angle convention (C2) for the amplitude follows from analyticity of the WDW phase in superspace. The measure (C3) does not follow**;
it is the charge-weight question of C3.
*Open (P-Q9′):* evaluate the Klein–Gordon current of the full spherically symmetric WDW equation through the bridge. It is sector-independent
by construction, and its WKB evaluation requires the imaginary part of the Kuchař phase at the bridge, presumably $\pm A_{\rm BH}/8G$.
This is now the precise form of the last ambiguity.

---

## Part D — The daughter edge, the daughter clock, unitarity (P-Q10; Stage 2 tasks 15–21)

### D1. Theorem Q10.1 (the daughter edge is removable) — [Theorem]
Along the daughter the static-frame momentum diverges at the interior de Sitter horizon, $P\simeq\pm(R/2G)\ln|f_{\rm in}|$, because
$\dot R/A_{\rm in}\to\pm1$. Define
$$\tilde P=P-s_e\frac R{2G}\ln|f_{\rm in}|,\qquad s_e=\mathrm{sign}\,A_{\rm in}(1/\chi),$$
which is a canonical transformation generated by a function of $R$ alone. Then the real-angle constraint becomes
$$\sqrt{f_{\rm out}}\,\big(f_{\rm in}e^{G\tilde P_u}+e^{-G\tilde P_u}\big)=\varepsilon_{\rm in}\varepsilon_{\rm out}\,(f_{\rm in}+f_{\rm out}-\kappa^2R^2)\qquad(s_e=+1;\ \text{mirror form for }s_e=-1),$$
which is **analytic in $f_{\rm in}$ across 0**. The true Hamiltonian $H(u,\tilde P_u)=\frac R{2G}(1-s^2)$ solves a quadratic whose discriminant at $f_{\rm in}=0$ is
$c_1^2+4\kappa^2R^2>0$, so it too is analytic there. $dt_\infty/d\tau=A_{\rm out}/f_{\rm out}$ is finite through the crossing. As $R\to\infty$, $t_\infty\to-\infty$, so the reduced
flow is complete in both channels.

*Numerics* (`daughter.edge_report`; single flip $\chi=\kappa=1$, $M=0.1545$; double flip $\gamma=0.6$, $E=2.4$):
* $\tilde P$ is continuous across $R=(1\mp10^{-6})/\chi$ to $5\times10^{-5}$ (single) and $4\times10^{-7}$ (double), while $P$ diverges;
* the residual of the regularised constraint and $|H(\tilde P)-M|$ are below $6\times10^{-12}$ on both sides.

**Correction F16.** Theorem V-1(ii) ("the daughter channel ends at finite $t_\infty$; a quantum boundary condition is required") is wrong. The edge was a
property of the momentum coordinate, not of the reduced phase space. No boundary condition and no self-adjoint extension are needed there.

### D2. Theorem Q10.2 (the asymptotic clock folds back) — [Theorem]
On the parent branch ($R\le R_1<R_f$) $A_{\rm out}>0$, so $dt_\infty/d\tau>0$. On every daughter ($R\ge R_2>R_f$, $\varepsilon_{\rm out}=-1$, Prop. G1) $A_{\rm out}<0$, so
$$\frac{dt_\infty}{d\tau}=\frac{A_{\rm out}}{f_{\rm out}}<0 .$$
In Kruskal coordinates $\partial_tT_K=X/4GM$: positive in region I, negative in region III, where the daughter wall lives. The global Killing time is past-directed
there. The physical history (laboratory past → transition → expanding daughter) is therefore **not monotone in $t_\infty$**. It turns around at the
transition. This is the Lorentzian face of the Euclidean fold of Theorem III-4(b), where $dt_\infty/d\tau=0$ at the flip radius.

Checked on single- and double-flip models (`daughter.clock_orientation`): $dt/d\tau=+4.5$ and $+2.2$ on the parent, $-1.2$ to $-1.3$ on the daughter.

### D3. Corollary Q10.3 (Klein-type structure of creation in the parent's clock) — [Theorem in the reduced model; Interpretation flagged]
With $H=M$ generating $t_\infty$, the group velocity is $dR/dt_\infty=\dot Rf_{\rm out}/A_{\rm out}$.
* **Stationary picture.** The physically expanding daughter has $dR/dt_\infty<0$: it is an *incoming* wave, with negative $t_\infty$-flux. Flux conservation for a
  stationary state with a purely expanding daughter reads $|R|^2=|I|^2+|D|^2$ (reflected exceeds incident, as in the Klein paradox). Ordinary
  transmission $|I|^2=|R|^2+|T|^2$ describes laboratory → **time-reversed** daughter.
* **Wave-packet picture.** Unitary $t_\infty$-evolution of a laboratory packet populates the daughter channel's $t_\infty\to+\infty$ end. That end is the
  *contracting* leg of the bounce branch, which lies in the Kruskal past of the transition. A daughter that expands towards the parent's future is
  not in the image of unitary boundary-clock evolution of laboratory states.
* **Magnitudes.** The barrier problem is time-reversal invariant, so $|D|^2/|I|^2=|T|^2/|I|^2=e^{-B}(1+O(\hbar))$. The exponent is unaffected. What changes is
  the interpretation: in the parent's Hilbert space, creation of the future-directed daughter is a Bogoliubov-type process with $|\beta|^2\simeq e^{-B}$, not a
  unitary transition between positive-norm states.

**Consistency (Stage 2 task 21).** This is a minisuperspace realisation of Freivogel–Hubeny–Maloney–Myers–Rangamani–Shenker (2006): no unitary boundary
evolution of a pure state produces an inflating region. It is compatible with Wall's GSL theorem, which forbids baby universes that become independent
of the mother: in the boundary clock the daughter is never an independent positive-norm subsystem. [Interpretation] In the unitary picture the daughter's
expansion arrow points to the parent's past, a two-arrow configuration reminiscent of Carroll–Chen (2004) spontaneous inflation.

### D4. The daughter clock (Stage 2 task 15)
Any clock that is future-directed in both sectors must reverse orientation relative to $t_\infty$ across the transition. Examples are $-t_\infty$ on the
daughter sector, the wall's proper time, or the interior cosmic time. This replaces the (incorrect) edge argument of Theorem V-1 as the reason the
daughter "needs its own clock". With such a clock the transition probability is $e^{-B}$ per laboratory bubble, as a crossed (time-reversed) amplitude.

### D5. Probability per attempt (P-T3, reduced-theory part) — [Proposition]
The bounded branch is not periodic in $t_\infty$: it runs from the past horizon ($t=-\infty$) through $R_1$ to the future horizon ($t=+\infty$). **Each bubble
gets exactly one attempt.** There is no attempt frequency, and the 1-d uniform (Kemble) connection gives $P=1/(1+e^{B})\simeq e^{-B}$ (times 2 for single flip
under C4(i)). The field-theoretic one-loop determinant (fluctuations of the wall and the bulk) remains open (P-T3′).

### D6. Transition operator and global geometry (Stage 2 tasks 16–17)
In the boundary clock the scattering channels are:
* in: parent from the past horizon; daughter from $R=\infty$ at $t=-\infty$, which is the expanding leg's future end;
* out: parent to the future horizon; daughter to $R=\infty$ at $t=+\infty$, which is the contracting leg's past end.

The physical amplitude for laboratory → expanding daughter is the crossed element $\langle\Theta\,{\rm daughter_{in}}|S|{\rm lab_{in}}\rangle^*$, with $\Theta$ time reversal. The
global geometry (Prop. V-3) keeps spatial topology $\mathbb R^3$ and creates an Einstein–Rosen throat. Pinch-off still requires a singularity.

### D7. Backreaction, unitarity, information (Stage 2 tasks 18–19)
* Hawking evaporation competes only after the transition (Phase V Sec. 6).
* The which-branch entropy is $h(e^{-B})$.
* By D3, the daughter carries negative boundary-clock norm. In a unitary parent description its creation is accompanied by an equal
  enhancement of the parent's reflected norm ($|R|^2-|D|^2=|I|^2$). [Interpretation] No positive-norm probability is transferred to an independent
  subsystem at this order. The information question again reduces to the fate of Hawking partners, as in Phase V.

---

## Part E — The alternative-mechanism search (Stage 4 tasks 31–38)

| task | mechanism | assumption relaxed | quantitative cost | verdict |
|---|---|---|---|---|
| 30 | FGG tunnelling, light walls | classical GR | $B=B_{\rm CdL}\Phi\ge B_{\rm CdL}$; frontier S7 | suppressed, $P\le10^{-23}$ in controlled regime |
| 30 | FGG tunnelling, heavy walls ($\kappa>\chi$) | classical GR | $B\ge0.564B(0)$; $B(0)\sim1/G\chi^2$ | suppressed |
| 32 | seed black hole inside the bubble | regular data | $B$ reduced by $\le$ ~30% ($m\to M$ moves the data out of regime A) | not creation: data already contain a de Sitter equator (F15) |
| 32 | pre-existing throat / white hole ($r_u>0$ data) | past completeness | Cor. IV-2: needs a past singularity | not producible from laboratory data |
| 34 | NEC violation by quantum fields, QEI-bounded | NEC | anti-trapped radius $\le\sqrt{4\pi C_{FR}N}\,\ell_P=0.35\sqrt N\ell_P$; $N\ge\frac83S_{\rm dS}$ | self-defeating: species cutoff $m_P/\sqrt N=0.35\,H$ |
| 36 | closed / de Sitter parent | non-compact Cauchy surface | bubble beyond the parent horizon; upward tunnelling $\sim e^{-S_{\rm dS}}=e^{-3.3\times10^{122}}$ | cosmological only |
| 36 | bubbles nucleated during inflation (GVZ 2016) | laboratory data | $M>M_{\rm cr}\simeq1/2GL$; nucleation rate only | viable in the early universe; PBH signature (P-O1/O2) |
| 31 | monopoles, charged exteriors (Sakai et al. 2006) | global hyperbolicity | tunnelling, wall collision, or continuation through a Cauchy horizon (mass inflation); core inflation needs $\eta_{\rm sym}\sim0.3m_P$ | no predictable classical laboratory route |
| 33 | aspherical collapse, thick walls | symmetry, thin wall | Penrose/Farhi–Guth are symmetry-free; thick-wall CdL → thin-wall | no classical escape; spherical saddle expected dominant (open) |
| 35 | limiting curvature (Frolov–Markov–Mukhanov 1990) | GR at Planck curvature | any black hole with $M>0.93m_P$ reaches $K=\ell_P^{-4}$ inside the horizon ($r_*=(48M^2)^{1/6}\ell_P$) | hypothesis-dependent; daughter at Planck curvature, invisible outside |
| 37–38 | new mechanism | — | none found that beats $B\gtrsim B_{\rm CdL}$ within GR + QEI-bounded fields | negative result |

**E1. The seed (F15).** For $f_{\rm in}=1-2Gm/r-\chi^2r^2$ the junction formulas are unchanged and $V_{\rm seed}=V_{M-m}-2Gm/r$.
* At $\chi=\kappa=1$, $M=0.3M_{\rm cr}$: $B/B=1,0.977,0.929,0.875,0.816,0.783$ for $m/M=0,\dots,0.8$.
* At $m/M\ge0.9$ the initial wall is no longer in regime A ($\varepsilon_{\rm out}(R_1)=-1$).

The decisive point is geometric. With $m>0$ the interior side of the wall continues through the seed's bridge into a static SdS patch whose
time-symmetric slice reaches the cosmological bifurcation sphere $r_c$. That is a maximal sphere, i.e. a de Sitter equator, so the initial data already
contain the seed of an inflating region. A seed formed by collapse would have no such patch, and the static analysis would not apply. That case
belongs to the dynamical problem of Phase IV, where Prop. IV-1 forbids classical creation.

**E2. QEI species bound.** By Prop. IV-1, $\alpha^2\partial_u(r_u/\alpha^2)=-4\pi Gr\,T_{uu}$ along an ingoing ray from flat data ($r_u/\alpha^2=-\tfrac12$). Anti-trapping
needs $\int rT_{uu}\,du\le-1/8\pi G$. Sampled over $\tau\sim\Delta u\sim2r$, a Ford–Roman inequality $\langle\rho\rangle\ge-3N/32\pi^2\tau^4$ caps the left side at
$-C_{FR}N/2r^2$. Hence $r\le\sqrt{4\pi C_{FR}N}\,\ell_P$. An inflating daughter needs an anti-trapped de Sitter sphere $r>1/\chi$, so
$N\ge1/(4\pi C_{FR}\chi^2)=\tfrac83S_{\rm dS}$ for every vacuum scale. The species scale $m_P/\sqrt N$ then lies at $0.35H$, below the daughter's own Hubble rate.
The semiclassical description of the daughter fails exactly where the NEC violation becomes strong enough. [Proposition, scaling with O(1) constants.]
The Fewster–Eveson bound is tighter by 9/64, so the estimate is conservative.

**E4. The cosmological route quantified (P-O2, scaling level) [Proposition].** Take bubbles nucleated during inflation at rate $\lambda H_i^4$
(Garriga–Vilenkin–Zhang).
* Each bubble grows to the comoving horizon at nucleation, so $dn/d\ln r=\lambda/r^3$. A fraction $\frac{4\pi}3\lambda$ of space is filled per e-fold.
* Bubbles re-enter in the radiation era with $M\propto r^2$, and their energy fraction grows $\propto a\propto M^{1/2}$ until equality. Hence
  $$f(M)\equiv\frac1{\rho_{\rm DM}}\frac{d\rho}{d\ln M}\simeq\frac{4\pi}3\,\lambda\,\Big(\frac{M_{\rm eq}}M\Big)^{1/2},\qquad M_{\rm eq}\approx3\times10^{17}M_\odot ,$$
  with $O(1)$ formation efficiencies omitted.
* All of the dark matter in asteroid-mass ($10^{20}$ g) daughter-hosting black holes needs $\lambda\approx10^{-16}$, i.e. an inflationary tunnelling exponent
  $B_{\rm inf}\approx37$. The daughter interpretation also needs $M>M_{\rm cr}$, i.e. $E_{\rm vac}\approx2\times10^6$ GeV (Phase VI).
* Contrast: the laboratory frontier is $B\ge53$ per bubble, at $E_{\rm vac}\approx3\times10^{17}$ GeV and $R\approx100\,\ell_P$. The cosmological route is not
  exponentially harder, but it cannot be run on demand (`alternatives.gvz_fraction`).

**E3. The de Sitter parent: a non-flip window [Proposition-N].** In a Schwarzschild–de Sitter parent ($\eta>0$), Prop. G1 ($\varepsilon_{\rm out}(R_2)=-1$) can fail.
A laboratory-type transition ($\varepsilon_{\rm out}(R_1)=+1$, wall in the static region) can then land on a daughter born on the *near* side of the
bridge, $\varepsilon_{\rm out}(R_2)=+1$. There the exterior orientation does not flip, so the embedding obstruction NG3 does not apply and a manifold
instanton is not excluded. Bisection over $M$ shows the window opens only for

| $\kappa/\chi$ | 0.1 | 0.3 | 0.5 | 0.7 |
|---|---|---|---|---|
| $\eta_*/\chi$ | 0.863 | 0.799 | 0.764 | 0.813 |

i.e. for parents whose vacuum energy is at least about 58–75% of the false vacuum's. Inside the window $B_{\rm an}=0.07$–$0.24/G\chi^2$ (at $\kappa=0.3\chi$,
$\eta=0.95\chi$), still of order the de Sitter entropy scale $(m_P/H)^2$. For our universe $\eta/\chi\lesssim10^{-50}$, so the window is closed. It belongs to the
cosmological route, an inflating parent (`alternatives.sds_nonflip_threshold`).

---

## Part F — Theoretical requirements: the five numbers (Stages 3, 5–10, theory only)

For the most favourable controlled case (Corollary S7: $B_{\rm CdL}\ge10$, $M\ge m_P$, $R_M\ge100\,\ell_P$):

| number | value | comparison |
|---|---|---|
| $E_{\rm req}$ (false-vacuum scale) | $\approx3\times10^{17}$ GeV, with wall scale $\approx2\times10^{17}$ GeV | 13–14 orders of magnitude above any collider ($\sim10^4$ GeV) |
| $L_{\rm req}$ (bubble radius) | $\approx100\,\ell_P\approx1.6\times10^{-33}$ m | $10^{-15}$ of a proton radius |
| $\tau_{\rm req}$ (formation time) | $\lesssim R/c\approx5\times10^{-42}$ s | |
| $P_{\rm create}$ | $\le e^{-53}\approx10^{-23}$ per bubble (one attempt per bubble, D5) | falls as $\sqrt{R_M}$: $10^{-73}$ at $R_M=1000\,\ell_P$ |
| $S_{\rm verify}$ | 0 | the daughter is invisible to every exterior observable (Theorem V-2) and lies behind the horizon of a Planck-mass black hole that evaporates in $\sim10^{-40}$ s |

**Requirement for a target probability** (thin barrier, $M\ge m_P$): $E_{\rm vac}\ge m_P\,X^{-3/2}$ with
$X=\ln(1/P)/[4\pi\sqrt{2B_{\rm CdL}/27}\,(3/4\pi)^{2/3}]$. For example $P=10^{-100}$ requires $E_{\rm vac}\ge3\times10^{16}$ GeV ($B_{\rm CdL}\ge10$).
Rows with $\kappa/\chi>0.3$ or $R_M<10\,\ell_P$ are outside the theory's validity and are flagged in the data file.

**What would have to change** for a laboratory route (the Breakthroughs of the Roadmap, restated as conditions):
* (A) a quantum-gravity inner product with entropy-weighted flipped sectors *and* a UV completion that permits $R_M\sim\ell_P$;
* (B) NEC-violating matter beyond quantum energy inequalities, e.g. a stable phantom-like sector;
* (C) a parent with compact or horizon-bounded Cauchy surfaces on laboratory scales.

None of these is supported by known physics.

---

## Part G — Stage 1 remaining items

* **NG4 (P-N1).** Reformulated: $|\theta_D|=\frac{E^2}{2b}\int_{s_1}^{s_2}\frac{s\,q(z)}{W\sqrt{W-s^2}}ds$ (all $s>0$ in regime B). NG4 is equivalent to monotonicity
  of this integral in $E$ on $(g_{\min},E_{AB})$, given the critical-end theorem. This is still Proposition-N (Phase II: monotone on a 39×10 grid).
* **SdS (P-G1), done.** An explicit conformal atlas with eight Kruskal-glued blocks (`sds_atlas.py`, `FIGURES/sds_atlas.png`). Every wall worldline is
  timelike and future-directed and continuous across all horizons (84/84 branches). The finite null coordinate at each crossing matches the Kruskal
  prediction. Also new: the SdS non-flip window (E3), where laboratory-type transitions land on near-side daughters without an exterior flip; it
  opens only for $\eta\gtrsim0.76$–$0.86\chi$. Euclidean SdS with two horizons was treated in Phase II (SdS Smarr identity).
  Conical defects in the Euclidean SdS sections are not revisited here.
* **Literature (P-L1)**, from accessible abstracts only:
  * FMP (1990), PRD 42, 4042: Dirac + WKB; agreement with the Euclidean method requires degenerate metrics. This matches NG3 + Prop. III-3.
  * Ansoldi–Tanaka (2014/15), JETP 120, 460: the discrepancies come from the static-chart time and are removed by a smooth slicing. This is consistent
    with T0 (the slicing jump is a horizon monodromy) and with D1–D2 here.
  * Sakai et al. (2006): monopole route.
  * Garriga–Vilenkin–Zhang (2016): the cosmological route.
  * Full texts are still unavailable here.
* **Manuscript (task 12).** `PAPER/genesis_foundations.tex` now carries an addendum listing every correction:
  * "ambiguity removed" → located (III); unambiguous iff double flip (S4);
  * picogram → nanogram, sub-Planckian (V);
  * the daughter edge (F16);
  * T5 status (A7);
  * P-Q2 (B).

---

## Part H — Failed approaches and corrections of this round

* **F14** (Phase III heuristic): "squaring removes the sign of $K_{\theta\theta}$". In the full theory squaring removes the time orientation; the orientation
  is configuration data (Q9.1).
* **F15** (Stage 4): a seed black hole as a catalyst. The exponent drops only modestly, and the seeded initial data already contain a de Sitter equator.
* **F16** (Phase V): the finite-time daughter edge of Theorem V-1(ii) was a coordinate singularity of the momentum (Q10.1).

## Open problems after this round

* **P-Q9′.** The Klein–Gordon current of the full spherically symmetric WDW equation through the bridge (C5).
* **P-T3′.** The field-theoretic one-loop prefactor.
* **P-N1.** NG4 monotonicity.
* **P-T5′.** Close the strip $E_{AB}<E<1.6E_{AB}$ analytically. (The $\Phi$ bound is now proved.)
* **P-L1.** Full-text comparison.
* **P-O1/O2.** GVZ mass function.
