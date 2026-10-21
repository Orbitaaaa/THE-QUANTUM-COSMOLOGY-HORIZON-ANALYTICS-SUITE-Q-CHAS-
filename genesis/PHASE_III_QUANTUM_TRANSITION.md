# Phase III — The quantum transition

**Exit condition (from the programme):** *either a defensible new quantum transition construction or a precise
theorem identifying why a proposed construction fails.* Status: **both, in the reduced (minisuperspace) theory.**

* We prove precisely **why the FGG "slicing" (principal-branch) construction fails** (Theorem III-5). It is not the WKB
  limit of any solution of its own quantum constraint. It is also negative, i.e. "probability" > 1, in 33 of 36 models.
* We prove that **no formulation removes the FGG ambiguity; every formulation hides it in one and the same datum**
  (Theorems III-2, III-4, Lemma III-6). That datum is the convention for the boost angle of a flipped normal, which is
  equivalent to choosing which sheet of the two-valued true Hamiltonian carries the daughter. It is also equivalent to
  how the Hilbert spaces of the two sectors are identified.
* We give a **defensible construction**: complex (Sorkin) angles, one analytic quantum constraint, Lebesgue measure and
  the real-barrier contour class. Its assumptions are explicit (C1–C4) and its exponent is $B_{\rm an}$. Canonical,
  covariant-Euclidean and Lorentzian complex-saddle formulations agree on it.
* We identify the one competitor not excluded: sector switching in the squared constraint, with exponent
  $0<B_{\rm sw}<B_{\rm an}$ (Sec. 7). It is decided by operator ordering, which is the sharpened open problem P-Q2.

The first draft of this report (earlier in this session) said the principal prescription was a "negative-tension
amplitude" and that the ambiguity was "eliminated, not hidden". **Both statements were wrong** and are corrected here.
FAILED_APPROACHES.md F10 records the error.

Code: `CODE/babyuniverse/contours.py`; script `CODE/scripts/phase3_analysis.py`; data
`DATA/phase3_quantum_transition.json`; figure `FIGURES/phase3_monodromy.png`; tests `CODE/tests/test_contours.py` (13).

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

## 7. Sector switching in the squared constraint (the competitor that survives)

A sector-blind quantization uses the squared constraint $\cosh^2(GP_u)=C^2$, with shifts $\pm2ih$. It has both families
$\cos\theta=\pm C$. Its symbol degenerates where $C=0$ ($\theta=\pi/2$), and there the branches $\theta$ and $\pi-\theta$ coalesce. In the
single-flip class the analytic $\theta$ runs from $0$ to $\pi$ and crosses $\pi/2$ at a point $R_q$, which is a genuine turning point
of the squared operator. A continuous connection then exists that follows $\theta$ up to $R_q$ and $\pi-\theta$ beyond it. It matches
the real-angle momentum at both ends, with
$$B_{\rm sw}=\frac2G\Big[\int_{R_1}^{R_q}R\,\theta\,dR+\int_{R_q}^{R_2}R\,(\pi-\theta)\,dR\Big],\qquad 0<B_{\rm sw}<B_{\rm an}\ \text{(all 24 single-flip models)}.$$
Beyond $R_q$ the branch $\pi-\theta$ decays more slowly than $\theta$. Any non-zero Stokes multiplier at $R_q$ therefore makes $B_{\rm sw}$
the dominant exponent.

> **Lemma III-7 (factorisation).** With $Q_\pm=T\mp W_+$, the ordering $Q_-Q_+$ contains $\ker Q_+$, and analytic continuation of an exact
> $\ker Q_+$ solution stays in $\ker Q_+$ (identity theorem). So in that ordering the laboratory solution can be taken in $\ker Q_+$ and does
> **not** switch: the exponent is $B_{\rm an}$. In other orderings ($T^2-W_+^2$ plus symmetrisations) the multiplier at $R_q$ comes from
> $O(\hbar)$ commutator terms and is not computed here.

| γ | E/E_AB | class | $B_{\rm an}$ | $B_{\rm princ}$ | $B_{\rm sw}$ | $A_2/4G$ |
|---|---|---|---|---|---|---|
| −0.90 | 1.5 | single | 4.6467 | −2.9917 | 2.8602 | 9.6281 |
| −0.45 | 1.5 | single | 2.4034 | −1.8803 | 1.6891 | 5.8021 |
| 0.00 | 1.5 | single | 0.9357 | −1.0324 | 0.8032 | 3.0153 |
| 0.00 | 8.0 | single | 0.8128 | −2.1317 | 0.8083 | 3.1408 |
| 0.45 | 1.5 | double | 0.1828 | −0.2221 | — | 1.1749 |
| 0.68 | 8.0 | double | 0.0206 | +0.0015 | — | 0.3659 |

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

**It is hidden, never eliminated. All formulations hide it in the same place, and that place is now explicit:**

| formulation | where the FGG datum sits |
|---|---|
| canonical, true Hamiltonian | which sheet $H_\pm$ carries the daughter; how sheets couple at the tangency point (III-4) |
| canonical, shift constraint | which analytic constraint ($W_+$ vs sign-dressed); the intertwiner anchor $u_*$ (III-2, III-6) |
| covariant Euclidean | the action assigned to folds of the pseudomanifold (joint angle $\pi$ or $0$) (III-3) |
| Lorentzian | the boost angle of a flipped normal ($i\pi$ or $0$), i.e. Sorkin's convention (III-2) |
| sector-blind | operator ordering of the squared constraint (Sec. 7) |

**Failure theorem.** The principal (FGG slicing) prescription solves no analytic quantum constraint (III-5). It is a splice at a
regular point and gives $e^{-B}>1$ in most of parameter space.

**Defensible construction (adopted).**
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

**What would overturn it.** A sector-blind quantization whose ordering gives a non-zero Stokes multiplier at $R_q$: then
$B_{\rm sw}<B_{\rm an}$ dominates (Sec. 7). Or a full-gravity inner product that anchors the daughter elsewhere (Sec. 6). Neither (C2) nor
(C3) is derived here from a gravitational path integral.

## 10. Open problems sharpened by Phase III

* **P-Q2 (sharpened).** At $M=0$ with $\kappa>\chi$ the squared constraint is polynomial in $u$, so the momentum representation turns it into an ODE
  ($L^4\hbar^2\phi''+2i\hbar(L^2-\chi^2\mathcal C)\phi'+(\mathcal C-1-2i\hbar\chi^2\mathcal C')\phi=0$, $\mathcal C=\cosh^2Gk$). The analytic $\theta$ runs from $0$ to $\pi$ there,
  so the coalescence point $R_q$ is present. **Compute the Stokes multiplier at $R_q$ for each ordering.** This decides between
  $B_{\rm an}$ and $B_{\rm sw}$.
* **P-Q7.** A pseudo-differential quantization of the single sheet $H_-$ in the real-angle convention. Its tunnelling amplitude may be
  governed by the branch locus rather than by a saddle, as with eigenfunction decay of relativistic Schrödinger operators
  (Carmona–Masters–Simon 1990, J. Funct. Anal. 91, 117; used here only as an analogy).
* **P-Q3.** A Picard–Lefschetz analysis of the fixed-energy path integral: which thimbles (windings, zero-net periods) contribute.
* **P-Q8.** Derive (C2)–(C3) from a covariant inner product, e.g. a group-averaging inner product for the shift constraint, or the
  crossed-product construction of Phase V.
