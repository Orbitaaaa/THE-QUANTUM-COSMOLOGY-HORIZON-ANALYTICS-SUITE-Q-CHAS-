# Phase III — The quantum transition

**Exit condition (from the programme):** *either a defensible new quantum transition construction or a precise
theorem identifying why a proposed construction fails.* Status: **both, in the reduced (minisuperspace) theory.**

* We prove precisely **why the FGG "slicing" (principal-branch) construction fails** (Theorem III-5). It is not the WKB
  limit of any solution of its own quantum constraint. It is also negative, i.e. "probability" > 1, in 33 of 36 models.
* We prove that **no single formulation removes the FGG ambiguity by itself; every formulation hides it in one and the same datum**
  (Theorems III-2, III-4, Lemma III-6). That datum is the convention for the boost angle of a flipped normal, which is
  equivalent to choosing which sheet of the two-valued true Hamiltonian carries the daughter. It is also equivalent to
  how the Hilbert spaces of the two sectors are identified.
* We give a **defensible construction**: complex (Sorkin) angles, one analytic quantum constraint, Lebesgue measure and
  the real-barrier contour class. Its assumptions are explicit (C1–C4) and its exponent is $B_{\rm an}$. Canonical,
  covariant-Euclidean and Lorentzian complex-saddle formulations agree on it.
* We construct the **sector-blind alternative** and settle it at leading order (Sec. 7). Quantizing the squared constraint, in
  any ordering, gives a universal Stokes multiplier $\sqrt2$ at the point where $C=0$ (Weber index $\nu=-\tfrac12$). Its exponent is
  $B_{\rm sq}=\frac2G\int R\arccos|C|\,dR$.
* **Classification.** The barrier contains a point with $C=0$ if and only if the transition is in the single-flip class
  (Lemma III-9, by concavity). So **in the double-flip class every consistent quantization gives $B_{\rm an}$.** The FGG ambiguity is
  fully resolved there. In the single-flip class exactly two consistent answers remain, $B_{\rm an}$ (sector-resolved) and
  $B_{\rm sq}<B_{\rm an}$ (sector-blind). They are separated by one physical question: **is the exterior orientation of the wall (the
  sign of its exterior extrinsic curvature) a superselected quantum number?**

The first draft of this report (earlier in this session) said the principal prescription was a "negative-tension
amplitude" and that the ambiguity was "eliminated, not hidden". **Both statements were wrong** and are corrected here.
FAILED_APPROACHES.md F10 records the error.

Code: `CODE/babyuniverse/contours.py`; script `CODE/scripts/phase3_analysis.py`; data
`DATA/phase3_quantum_transition.json`; figure `FIGURES/phase3_monodromy.png`; tests `CODE/tests/test_contours.py` (19).

---

## 1. The system, the states and what must be specified

**Reduced phase space** (QUANTIZATION.md): wall radius $R$, or the area variable $u=R^2/2$ (so $A=4\pi R^2=8\pi u$). Its
momentum is $P$, or $P_u=P/R$. The boundary clock is $t_\infty$ and the Hamiltonian is $H=M$. The **static domain** is the
set of states with $f_{\rm in},f_{\rm out}>0$ at the wall. The barrier always lies in it ($f\ge V>0$).

**Initial state $|{\rm in}\rangle$.** A packet of energy $M$ on the expanding leg of the regime-A bounded branch, sector
$(\varepsilon_{\rm in},\varepsilon_{\rm out})=(+,+)$, approaching the turning point $R_1$. This is regular laboratory data.

**Final states.** $|{\rm col}\rangle$: reflected at $R_1$ and collapsing to a black hole. $|{\rm dau}\rangle$: the outgoing bounce
branch for $R>R_2$, the daughter universe. Two classes occur (data: 24 and 12 of 36 grid models):

* **single flip:** $\varepsilon_{\rm out}$ changes sign on the barrier and $\varepsilon_{\rm in}$ does not, so the daughter is in sector $(+,-)$;
* **double flip:** both change sign, so the daughter is in sector $(-,-)$.

**What a construction must specify.** (i) a quantum constraint or Hamiltonian valid across the barrier; (ii) a solution
with the laboratory asymptotics; (iii) a probability, i.e. an inner product or conserved current, in *both* sectors. We
show that (i)–(iii) cannot all be read off from the classical reduced theory. The undetermined piece is exactly the FGG
ambiguity.

## 2. Theorem III-1 (horizon monodromy)

> The on-shell momentum $P(R)$ at fixed $M$ is analytic on a Riemann surface. Its finite singular points are the turning
> points (square root), the horizons $f_{\rm in}=0$, $f_{\rm out}=0$ (logarithmic) and $R=0$. Flip radii ($A=0$) are regular. A
> contour from $R_1$ to $R_2$ with winding numbers $n_h$ around the horizon points $r_h$ has
> $$B[\gamma]=B_{\rm an}+\sum_h n_h\,\sigma_h\,\frac{A(R_2)-A_h}{4G},\qquad\sigma_h=\pm1\ \text{(model dependent)},$$
> and an unchanged real part. The endpoint momentum is shifted by $i\pi(\sum_h\sigma_hn_h)R/G$.
>
> *Proof.* On shell $A^2-\dot R^2=f$, so $x=\dot R/A=\pm1$ exactly at $f=0$. There $\log\frac{1+x}{1-x}$ jumps by $\pm2\pi i$, so
> $P\to P\pm i\pi R/G$. A loop based at $R_c$ contributes $\pm\frac{i\pi}{2G}(R_c^2-r_h^2)$ and the shifted remainder
> contributes $\pm\frac{i\pi}{2G}(R_2^2-R_c^2)$. $\square$

**Zero-net periods.** When $\sum_h\sigma_hn_h=0$ the endpoint momentum is unchanged and
$\Delta B=\pm(A_{h_2}-A_{h_1})/4G$. For example, with the black-hole horizon and the interior de Sitter horizon,
$\Delta B=\pm(A_{\rm dS}-A_{\rm BH})/4G$.

**Verification.** Single windings $n=-1,1,2$ around the BH, dS and SdS horizons in four models agree with the prediction to
$10^{-4}$; ${\rm Re}$ is unchanged to $4\times10^{-15}$. Zero-net pairs agree to $10^{-6}$ (endpoint momentum unchanged to $10^{-15}$) in four
models including SdS. The cancelling orientation is $(1,-1)$ in single-flip models and $(1,1)$ in double-flip models.

## 3. Theorem III-2 (two conventions for the canonical momentum)

Write the momentum through the boost angle $\zeta$ between the wall normal $n$ and the static frame $e_r=\sqrt f\,\partial_r$ on each
side, $P=\frac RG(\zeta_{\rm out}-\zeta_{\rm in})$. Since $n\cdot e_r=A/\sqrt f$ and the tangential component is $\dot R/\sqrt f$:

* **real-angle convention:** $\zeta=\operatorname{artanh}(\dot R/A)$, real on the Lorentzian static domain. This is the angle modulo $i\pi$.
* **complex-angle (Sorkin) convention:** $\zeta_S=\log\frac{A+\dot R}{\sqrt f}=\operatorname{artanh}(\dot R/A)+i\pi[A<0]$. A normal pointing
  into the opposite wedge (flipped, $\varepsilon=-1$) has boost angle $i\pi$ plus real. This is the convention of Sorkin's Lorentzian
  trigonometry, used for corner and joint terms of the gravitational action (Sorkin, arXiv:1908.10022). Neiman showed that
  such imaginary parts of the Lorentzian action reproduce Bekenstein–Hawking entropy (arXiv:1301.7041).

> **Theorem III-2.**
> (a) $\zeta_S$ is analytic on the whole static domain, including the barrier, where $\zeta_S=i\Theta$ with $\Theta={\rm atan2}(\sqrt V,A)$.
> $A+\dot R$ vanishes only at $f=0$. **The analytic branch of Stages I–II is the complex-angle momentum**, and it satisfies
> the **sign-free analytic constraint**
> $$2\sqrt{f_{\rm in}f_{\rm out}}\,\cosh(GP_u)=f_{\rm in}+f_{\rm out}-2\kappa^2u\qquad(\equiv W_+(u)\sqrt{f_{\rm in}f_{\rm out}}).$$
> (b) **The principal branch is the real-angle momentum.** It satisfies the **sign-dressed constraint**
> $2\sqrt{f_{\rm in}f_{\rm out}}\cosh(GP_u)=\varepsilon_{\rm in}\varepsilon_{\rm out}(f_{\rm in}+f_{\rm out}-2\kappa^2u)$, whose right-hand side jumps at every flip radius.
> (c) The two are related by the complex canonical transformation generated by $F=\frac{i\pi}{2G}(n_{\rm out}-n_{\rm in})R^2$,
> with $n=[A<0]$, which is constant within a sector. Quantum mechanically this is multiplication by the anti-periodic factor
> $e^{-\pi(n_{\rm out}-n_{\rm in})u/\hbar G}$, which intertwines the $W_+$ and $\varepsilon\varepsilon W_+$ equations. It is **not unitary**.
> (d) At a turning point the complex-angle momentum equals $i\pi R(n_{\rm out}-n_{\rm in})/G$ and the real-angle momentum is $0$.
> Each convention connects continuously to its own Lorentzian momenta at *both* turning points: the analytic branch in
> the complex convention, and the principal branch (with jumps) in the real convention.
>
> *Proof.* (a) With $\cosh\zeta_S=A/\sqrt f$ and $\sinh\zeta_S=\dot R/\sqrt f$ (exact, including signs),
> $\cosh(\zeta_{\rm out}-\zeta_{\rm in})=(A_{\rm in}A_{\rm out}-\dot R^2)/\sqrt{f_{\rm in}f_{\rm out}}$, and $A_{\rm in}A_{\rm out}+V=(f_{\rm in}+f_{\rm out}-\kappa^2R^2)/2$.
> (b) ${\rm artanh}$ differs from $\zeta_S$ by $i\pi n$, and $\cosh(z-i\pi m)=(-1)^m\cosh z$ with $(-1)^{n_{\rm out}-n_{\rm in}}=\varepsilon_{\rm in}\varepsilon_{\rm out}$.
> (c) $P'=P+\partial_RF$, and $e^{iF/\hbar}=e^{-\pi(n_{\rm out}-n_{\rm in})u/\hbar G}$, which changes sign under $u\to u+i\hbar G$. (d) follows from
> (a), (b) and $\zeta_S=\log\varepsilon$ at $\dot R=0$. $\square$

Verification: analytic branch vs sign-free constraint, and principal branch vs sign-dressed constraint, both to $\le8\times10^{-16}$.
The cross-violations are $1.6$–$2.0$ (three models).

## 4. Theorem III-4 (the branch point of the true Hamiltonian is the flip is the fold)

Solve the sign-free constraint for the mass conjugate to $t_\infty$. With $F=\sqrt{f_{\rm in}}$, $c=\cosh GP_u$ and $s=\sqrt{f_{\rm out}}$:
$$s_\pm=Fc\pm\sqrt{F^2(c^2-1)+\kappa^2R^2},\qquad H_\pm(u,P_u)=\frac R{2G}\big(1-\eta^2R^2-s_\pm^2\big).$$

> **Theorem III-4.**
> (a) Along the complex classical trajectory on the barrier ($GP_u=i\theta$, analytic branch), the radicand equals
> $$\kappa^2R^2-f_{\rm in}\sin^2\theta=\kappa^2R^2\,\frac{A_{\rm out}^2}{f_{\rm out}},$$
> with a **double zero at the exterior flip radius**. The trajectory is tangent to the branch locus there. Along it the square root has
> the analytic determination $-\kappa RA_{\rm out}/\sqrt{f_{\rm out}}$, so $s=Fc-\kappa RA_{\rm out}/\sqrt{f_{\rm out}}$, and the sheet label is $-\varepsilon_{\rm out}$.
> (b) At that point $dt_\infty/d\tau=A_{\rm out}/f_{\rm out}=0$ and $\partial H/\partial P_u=\infty$. **The Euclidean wall turns around in Killing
> time**: this is the fold that makes FGG's exterior a pseudomanifold. The canonical branch point, the Euclidean fold and
> the slicing ambiguity are the same point.
> (c) On the real (Lorentzian) phase space the radicand is $\ge\kappa^2R^2>0$, so $H_+$ and $H_-$ are two distinct real-analytic
> Hamiltonians. Which one carries a given physical state depends on the convention of Theorem III-2:
>
> | state | real-angle convention | complex-angle convention |
> |---|---|---|
> | laboratory $(+,+)$ near $R_1$ | $H_-$ | $H_-$ |
> | single-flip daughter $(+,-)$ near $R_2$ | $H_-$ | $H_+$ |
> | double-flip daughter $(-,-)$ near $R_2$ | $H_+$ | $H_+$ |
>
> (d) Interior flips are not branch points of $H$. They are branch points of $\Lambda_{\rm in}(u,P_u)$, the function conjugate to the
> interior clock. This is consistent with the Stage One finding (R9) that interior flips do not obstruct embedding.
>
> *Proof.* (a) $F\sin\theta=F\sin(\Theta_{\rm out}-\Theta_{\rm in})=\sqrt V(A_{\rm in}-A_{\rm out})/\sqrt{f_{\rm out}}=\kappa R\sqrt{V/f_{\rm out}}$ and $f_{\rm out}-V=A_{\rm out}^2$.
> Also $s-F\cos\theta=(f_{\rm out}-A_{\rm in}A_{\rm out}-V)/\sqrt{f_{\rm out}}=-\kappa RA_{\rm out}/\sqrt{f_{\rm out}}$. (b) $f_{\rm out}\dot t=A_{\rm out}$, and Lorentzian/Euclidean alike
> $\partial H/\partial P_u\propto1/\sqrt{\rm radicand}$. (c) For a real $P_u$, $c\ge1$. The table follows from (a) at the turning points ($c=\pm1$), with
> $s=|A_{\rm out}|$. (d) The mass function does not involve the interior clock. $\square$

Verification: radicand at $R_f$ $\le7\times10^{-16}$ (relative), no sign change along the barrier, $|H_\pm(R_f)-M|\le4\times10^{-9}$. The
table holds to $10^{-15}$ in four models, including SdS.

**Corollary (sector decoupling).** Take any quantization that is diagonal in sectors: $\mathcal H=\bigoplus\mathcal H_\pm$,
$\hat H=\bigoplus\hat H_\pm$, each self-adjoint with a real symbol. This is the structure of a true-Hamiltonian (FMP-type) quantization in
the complex convention. Such a quantization gives **zero** laboratory → daughter amplitude in the single- and double-flip
classes. A non-zero amplitude needs a coupling of the sheets, which classically exists only at the tangency point (a) in the
*complexified* phase space.

In the real-angle convention the single-flip transition stays on $H_-$. But the complex level set $\{H_-=M\}$ under the barrier
then consists of two arcs ending at the distinct branch points $(u_f,\theta_f)$ and $(u_f,\pi-\theta_f)$. Connecting them on one
sheet needs a jump, and that jump is the principal prescription.

## 5. Theorem III-5 (why the FGG principal construction fails)

> **Theorem III-5.** Let $I$ be a real interval in the barrier containing an exterior flip radius $u_f$. Let $\psi$ be analytic on the strip
> $\{{\rm Re}\,u\in I,\ |{\rm Im}\,u|<\hbar G+\delta\}$; for $L^2$ functions this is what membership in the domain of $e^{\pm i\hbar GP_u}$
> means (Paley–Wiener). If
> $$\psi(u+i\hbar G)+\psi(u-i\hbar G)=\varepsilon_{\rm in}\varepsilon_{\rm out}(u)\,W_+(u)\,\psi(u)\quad(u\in I),$$
> then $\psi\equiv0$.
>
> *Proof.* At $R_f$, $W_+=2C=2V/\sqrt{f_{\rm in}f_{\rm out}}>0$, because $A_{\rm out}=0$. Put $g=W_+\psi$, which is analytic on $I$. The left side $L$ is analytic
> on $I$. $L=g$ on one side of $u_f$, so $L\equiv g$ on $I$. $L=-g$ on the other side, so $g=0$ there, hence $g\equiv0$, and
> $\psi\equiv0$ near $u_f$ (since $W_+\ne0$), hence everywhere. $\square$

So the principal-branch WKB function is a **splice of two different solutions at an ordinary point**. $R_f$ is not a turning
point, and $W_+$ is regular there, so no Stokes phenomenon can occur at $R_f$. In addition, $B_{\rm princ}<0$ in **33 of 36** grid
models (both flip classes). Read as a transition probability, $e^{-B_{\rm princ}}>1$. The construction can be rescued only by
reinterpreting $e^{-B}$ as an unnormalised rate that includes a final-state degeneracy $e^{+S}$ (detailed balance). That reading
already assumes the horizon-entropy weights that the construction was meant to produce.

## 6. Lemma III-6 (no conserved current; the measure is the hidden datum)

> **Lemma III-6.** Write $h=\hbar G$ and $(T\psi)(u)=\psi(u+ih)+\psi(u-ih)$.
> (a) The Casoratian $\mathcal C[\psi_1,\psi_2](u)=\psi_1(u)\psi_2(u+ih)-\psi_1(u+ih)\psi_2(u)$ of two solutions of $T\psi=W\psi$ is $ih$-periodic.
> It is *not* constant along real $u$, so the shift-form constraint has **no conserved probability current** by itself.
> (b) $T$ is symmetric in $L^2(\mathbb R,\mu\,du)$, for $\mu$ real-analytic and positive and on a dense set of entire functions decaying
> in the strip, **iff $\mu$ is $ih$-periodic**: $\mu=\Phi(e^{2\pi u/h})=\Phi(e^{A/4\hbar G})$.
> (c) If in addition $|\mu(u)|\le Ce^{c|{\rm Re}\,u|}$ with $c<2\pi/h$ on a period strip, then $\mu$ is constant. The same holds for any entire
> $ih$-periodic multiplier ("quasi-constant") of a solution.
>
> *Proof.* (a) $\mathcal C(u)-\mathcal C(u-ih)=\psi_2(u)[\psi_1(u+ih)+\psi_1(u-ih)]-\psi_1(u)[\psi_2(u+ih)+\psi_2(u-ih)]=0$. (b) Shifting contours with
> $\varphi^\sharp(z)=\overline{\varphi(\bar z)}$:
> $\langle\varphi,T\psi\rangle_\mu-\langle T\varphi,\psi\rangle_\mu=\int\varphi^\sharp(v-ih)\psi(v)[\mu(v-ih)-\mu(v)]dv+\int\varphi^\sharp(v+ih)\psi(v)[\mu(v+ih)-\mu(v)]dv$.
> (c) $\Phi$ is analytic on $\mathbb C^*$ with $|\Phi(q)|\le C\max(|q|,|q|^{-1})^{ch/2\pi}$, exponent $<1$, so all Laurent coefficients but $a_0$
> vanish (Cauchy estimates). $\square$

**Consequence.** Lebesgue measure is the unique symmetric measure that grows more slowly than $e^{A/4\hbar G}$, the exponential
of the wall's Bekenstein–Hawking entropy. Even with it, the probability of the daughter needs a normalisation in the daughter
sector. In the complex convention the daughter's WKB wave function carries $e^{-\pi u/\hbar G}$ relative to the real-angle WKB state,
whose flux is conserved. Relating the two needs an intertwiner $e^{\pm\pi(u-u_*)/\hbar G}$, which is anti-periodic and defined only up to the
constant $u_*$. Anchoring at the daughter's turning point ($u_*=u_2$) gives $B_{\rm an}$. Anchoring at the flip ($u_*=u_f$)
reproduces the principal value $B_{\rm an}-(A_2-A_f)/4G$. **This anchor is the FGG ambiguity in its canonical form.**

## 7. The sector-blind quantization: universal Stokes multiplier and the classification

A **sector-blind** quantization quantizes the squared constraint
$$\sigma_0(u,P_u)=g(u)\big[\cosh^2(GP_u)-C(u)^2\big],\qquad g=4f_{\rm in}f_{\rm out},$$
which contains both families $\cos\theta=\pm C$ and only the shifts $\pm2i\hbar G$. It is invariant under the intertwiner
$e^{\pi u/\hbar G}$ of Theorem III-2(c), so it does not distinguish orientation sectors. Dirac (Wheeler–DeWitt) quantization of the full
bulk-plus-shell system plausibly leads here: the bulk Hamiltonian constraint is quadratic in the gravitational momenta, and the sign of
the extrinsic curvature $K_{\theta\theta}=A/R$ is exactly the sign that squaring removes. This motivation is **heuristic**.

> **Roadmap note (F14).** The heuristic above is wrong in the full theory: the orientation is boost-invariant configuration data (sign of $R'/\Lambda$), and
> Wheeler–DeWitt squaring removes the time orientation instead (ROADMAP Prop. Q9.1). The two families of the squared constraint are two phase conventions
> for the same daughter, distinguished by their conserved Casoratian charge (Theorem Q9.3). The Stokes multiplier $\sqrt2$ below is confirmed by exact
> integration (ROADMAP Part B).

Its symbol has a **double point** wherever $C=0$, at $(u_q,\,P_q=i\pi/2G)$. There the branches $\theta$ and $\pi-\theta$ coalesce.

> **Lemma III-9 (where coalescence happens).** $N(R)\equiv f_{\rm in}+f_{\rm out}-\kappa^2R^2=2-(\chi^2+\eta^2+\kappa^2)R^2-2GM/R$ is strictly concave on
> $R>0$, and $C=N/2\sqrt{f_{\rm in}f_{\rm out}}$. On the barrier $C(R_1)=+1$ and $C(R_2)=\varepsilon_{\rm in}\varepsilon_{\rm out}(R_2)$. Hence:
> **double flip** ⇒ $C>0$ on the whole barrier (no coalescence point; $\theta<\pi/2$ throughout);
> **single flip** ⇒ exactly one coalescence point $R_q$, where $C$ changes sign with slope $c=C'(u_q)<0$.
>
> *Proof.* $N''=-2(\chi^2+\eta^2+\kappa^2)-4GM/R^3<0$. The set $\{N>0\}$ is an interval containing $R_1$; it contains $R_2$ iff $C(R_2)>0$. $\square$

Numerically: in all 12 double-flip grid models $\theta_{\max}\le0.30\pi$, and in all 24 single-flip models there is exactly one crossing.

> **Proposition III-8 (Weber index at the coalescence point; leading order).** Let $L$ be a quantization with Weyl symbol
> $\sigma_0+\hbar\sigma_1+O(\hbar^2)$. Conjugate by $e^{-\pi u/2\hbar G}$, which shifts $P_u$ by $P_q$, and put $x=u-u_q$, $p=P_u-P_q$. Then near the
> double point $\sigma_0=-g_q\,(G^2p^2+c^2x^2)+O(3)$, and $L\psi=0$ reduces to the Weber equation
> $-(\hbar G)^2\psi''+c^2x^2\psi=E\psi$ with $E=\hbar\sigma_1(u_q,P_q)/g_q$, i.e.
> $$\nu=-\tfrac12+\frac{\sigma_1(u_q,P_q)}{2G|c|\,g_q}.$$
> Locally the analytic branch $\theta$ is the Gaussian $e^{-|c|x^2/2\hbar G}$, recessive on the left. The solution carrying the laboratory
> state is therefore $D_\nu(-z)$, $z=x\sqrt{2|c|/\hbar G}$. To the right it contains the branch $\pi-\theta$ ($e^{+|c|x^2/2\hbar G}$) with
> multiplier
> $$s=\frac{\sqrt{2\pi}}{\Gamma(-\nu)}\,e^{i\pi\nu},$$
> which vanishes iff $\nu\in\{0,1,2,\dots\}$.
> (a) **Every ordering of the squared constraint** has $\nu=-\tfrac12$ and $|s|=\sqrt2$. "Ordering" means any operator obtained by ordering
> the factors $a(u)$ and $\cosh^2(G\hat P_u)$ of $\sigma_0$, including division by $g$ and symmetrisations. Each such reordering, and the
> standard-to-Weyl correction, contributes to $\sigma_1$ a multiple of $a'(u)\,\partial_P\cosh^2(GP)=G\,a'(u)\sinh(2GP)$, and
> $\sinh(2GP_q)=\sinh(i\pi)=0$.
> (b) The **sector-resolving** product $Q_-Q_+$, with $Q_\pm=T/2\mp C$, has $\sigma_1(u_q,P_q)=-cG\,g_q$ (from
> $[T,u]=-2i\hbar G\sinh GP_u$ and $\sinh GP_q=i$). Since $c<0$, $\nu=0$: **no switching**, consistent with $\ker Q_+\subset\ker Q_-Q_+$ (identity
> theorem). The other product $Q_+Q_-$ has $\nu=-1$ and $|s|=\sqrt{2\pi}$.
>
> *Status.* (a) and (b) are exact statements about the $O(\hbar)$ symbol. The reduction to the Weber equation is the standard normal form
> at a non-degenerate double point. Cubic and higher terms change the multiplier only at relative order $\hbar^{1/2}$ and higher. The
> connection coefficient $\sqrt{2\pi}/\Gamma(-\nu)$ is checked against mpmath's parabolic-cylinder function to $4\times10^{-4}$
> ($z=30$, finite-$z$ correction) for $\nu=-1,-\tfrac12,-\tfrac14,0,1$ (`contours.weber_multiplier`).

**Corollary (sector-blind exponent).** Beyond $R_q$ the branch $\pi-\theta$ decays more slowly than $\theta$. With $|s|=\sqrt2\ne0$ it dominates,
and the transmitted amplitude follows $\min(\theta,\pi-\theta)=\arccos|C|$:
$$B_{\rm sq}=\frac2G\int_{R_1}^{R_2}R\,\arccos|C|\,dR\ \le\ B_{\rm an}=\frac2G\int_{R_1}^{R_2}R\,\arccos C\,dR,$$
with probability prefactor enhanced by $|s|^2=2$. Equality holds iff there is no coalescence point, i.e. exactly in the double-flip class
(Lemma III-9). In the single-flip class $B_{\rm sq}=B_{\rm sw}$, where $B_{\rm sw}$ follows $\theta$ to $R_q$ and $\pi-\theta$ beyond. Its endpoint momentum is the
real-angle one. Numerically, $|B_{\rm sq}-B_{\rm sw}|/B_{\rm sw}\le5\times10^{-15}$ (single flip), $|B_{\rm sq}-B_{\rm an}|/B_{\rm an}\le4\times10^{-15}$ (double flip), and
$0<B_{\rm sq}<B_{\rm an}$ in all 24 single-flip models.

| γ | E/E_AB | class | $B_{\rm an}$ | $B_{\rm princ}$ | $B_{\rm sq}$ | $A_2/4G$ |
|---|---|---|---|---|---|---|
| −0.90 | 1.5 | single | 4.6467 | −2.9917 | 2.8602 | 9.6281 |
| −0.45 | 1.5 | single | 2.4034 | −1.8803 | 1.6891 | 5.8021 |
| 0.00 | 1.5 | single | 0.9357 | −1.0324 | 0.8032 | 3.0153 |
| 0.00 | 8.0 | single | 0.8128 | −2.1317 | 0.8083 | 3.1408 |
| 0.45 | 1.5 | double | 0.1828 | −0.2221 | 0.1828 | 1.1749 |
| 0.68 | 8.0 | double | 0.0206 | +0.0015 | 0.0206 | 0.3659 |

(units $G=\kappa=1$; full grid in the data file.)

## 8. Proposition III-3 (covariant action = complex-angle convention)

For regime A let $\Delta t_\infty=\Delta t^W_E+8\pi GM$, one half-turn of the exterior at the tip per half-bounce. Then, exactly,
$$
B_{\rm an}=-\tfrac{\Lambda_{\rm in}}{8\pi G}{\rm Vol}^{\rm in}_E-\tfrac\sigma2A^W_E+\tfrac M2\Delta t_\infty-M\Delta t_\infty+\tfrac{A_{BH}}{4G}
+\tfrac{n(R_2)A(R_2)}{4G}+\tfrac{A_\chi}{4G}\big|_{\varepsilon_{\rm in}(R_2)=-1},
$$
verified to $1.2\times10^{-9}$. The terms are bulk, wall, GHY at infinity, background, bolt, and fold or de Sitter bolt. Each fold or bolt
contributes $\pi A/8\pi G$ per half-instanton, a Hayward joint term with angle $\pi$. That is precisely the Euclidean image of the
$i\pi$ boost angle of a flipped normal. **So the covariant Euclidean action of the pseudomanifold with joint terms and the
canonical theory in the complex convention are the same prescription.** Their agreement is a consistency check, not
independent evidence for the convention.

The Lorentzian complex saddle is also computed on the analytic branch. Its imaginary asymptotic time is
${\rm Im}\,T^*=\tfrac12(-\Delta t^W_E+\partial_M\mathcal A/4G)$, verified to $5\times10^{-6}$, and each BH winding shifts it by $\beta/2$.

## 9. Verdict: does any resolution eliminate the FGG ambiguity, or merely hide it?

**No formulation eliminates it by its own rules. All formulations hide it in the same place, and that place is now explicit.**
The classification below then shows that this datum matters only in the single-flip class.

| formulation | where the FGG datum sits |
|---|---|
| canonical, true Hamiltonian | which sheet $H_\pm$ carries the daughter; how sheets couple at the tangency point (III-4) |
| canonical, shift constraint | which analytic constraint ($W_+$ vs sign-dressed); the intertwiner anchor $u_*$ (III-2, III-6) |
| covariant Euclidean | the action assigned to folds of the pseudomanifold (joint angle $\pi$ or $0$) (III-3) |
| Lorentzian | the boost angle of a flipped normal ($i\pi$ or $0$), i.e. Sorkin's convention (III-2) |
| sector-blind vs sector-resolved | whether the constraint distinguishes orientation sectors; not the ordering within the squared form (Prop. III-8) |

**Failure theorem.** The principal (FGG slicing) prescription solves no analytic quantum constraint (III-5). It is a splice at a
regular point and gives $e^{-B}>1$ in most of parameter space.

**Classification theorem (leading order, given C1, C3, C4).** Two consistent quantizations are constructed.
* **Double-flip class:** both give $B_{\rm an}$. Both conventions for flipped normals agree at the endpoints there. The FGG ambiguity is
  **eliminated** in this class: the only dissenting prescription is the principal splice, excluded by III-5.
* **Single-flip class:** sector-resolved gives $B_{\rm an}$, sector-blind gives $B_{\rm sq}<B_{\rm an}$ (Prop. III-8, Lemma III-9). The residual ambiguity is the
  single physical question whether the exterior orientation of the wall is superselected.

A third candidate, the single-sheet pseudo-differential quantization of $H_-$ (P-Q7), is not analysed.

**Defensible construction (sector-resolved, adopted as primary).**
* (C1) The reduced theory, with constraints solved classically and the boundary clock $t_\infty$, is the leading-order quantum theory.
* (C2) Boost angles are Sorkin's complex angles, equivalently Euclidean angles under continuation. The quantum constraint is then the single
  analytic equation $T\psi=W_+\psi$ on the static domain.
* (C3) The measure is Lebesgue, the unique symmetric measure growing more slowly than $e^{A/4\hbar G}$ (III-6c).
* (C4) The contour is in the real-barrier homotopy class: zero horizon winding and zero-net periods excluded. This is the Stokes segment between the two
  simple turning points.

Then the WKB connection from $|{\rm in}\rangle$ is unique and
$$P_{\rm daughter}=e^{-B_{\rm an}/\hbar}\,(1+O(\hbar))\quad\text{(daughter normalised at its moment of time symmetry, }u_*=u_2).$$
$B_{\rm an}$ is positive in every regime (Lemma P). It agrees with the covariant pseudomanifold action with joint terms and with the Lorentzian
complex saddle (Sec. 8). As $M\to0$ it tends to the closed form T1, which is the CdL exponent when the false vacuum is the larger cap (T1/T5).

**The sector-blind alternative (equally consistent).** Replace (C2) by the squared constraint. Then
$P_{\rm daughter}\simeq2\,e^{-B_{\rm sq}/\hbar}$ (single flip) or $e^{-B_{\rm an}/\hbar}$ (double flip).

**What decides between them.** Whether the full theory superselects the exterior orientation (Sorkin angles and the covariant joint terms
favour resolving it; Wheeler–DeWitt quantization of a constraint quadratic in the momenta favours blindness). Equivalently, whether
the inner product of the full theory anchors the daughter as in Sec. 6. Neither (C2) nor the squared alternative is derived here from a
gravitational path integral.

**Roadmap update.** (1) Theorem S4: the double-flip class is exactly $\kappa<\chi$, $M<M_{\rm sd}=(\chi^2-\kappa^2)/2G\chi^3$, which covers all but a sliver
$\simeq\kappa^2/2\chi^2$ of the physically relevant (light-wall) range. (2) Corollary Q9.4: if probability is one sector-independent current (the laboratory's),
the canonical answer is $2e^{-B_{\rm sq}}$. The primary choice $B_{\rm an}$ corresponds to an entropy-weighted daughter measure.

**For observational bounds (Phase VI)** we use the smaller exponent, $B_{\rm sq}=\min(B_{\rm an},B_{\rm sq})$. That is the optimistic choice, so any
inaccessibility conclusion holds under either quantization.

## 10. Open problems sharpened by Phase III

* **P-Q2 (reduced to a global check).** Prop. III-8 settles the multiplier at leading local order. Remaining: confirm it with an exact
  solution of the $M=0$, $\kappa>\chi$ squared constraint (a polynomial in $u$, hence an ODE in the momentum representation:
  $L^4\hbar^2\phi''+2i\hbar(L^2-\chi^2\mathcal C)\phi'+(\mathcal C-1-2i\hbar\chi^2\mathcal C')\phi=0$, $\mathcal C=\cosh^2Gk$), and compute the prefactors.
* **P-Q9.** Decide superselection of the exterior orientation in a full (field-theoretic or Wheeler–DeWitt) treatment of the wall. This is
  the one physical question left between $B_{\rm an}$ and $B_{\rm sq}$ in the single-flip class.
* **P-Q7.** A pseudo-differential quantization of the single sheet $H_-$ in the real-angle convention. Its tunnelling amplitude may be
  governed by the branch locus rather than by a saddle, as with eigenfunction decay of relativistic Schrödinger operators
  (Carmona–Masters–Simon 1990, J. Funct. Anal. 91, 117; used here only as an analogy).
* **P-Q3.** A Picard–Lefschetz analysis of the fixed-energy path integral: which thimbles (windings, zero-net periods) contribute.
* **P-Q8.** Derive (C2)–(C3) from a covariant inner product, e.g. a group-averaging inner product for the shift constraint, or the
  crossed-product construction of Phase V.
