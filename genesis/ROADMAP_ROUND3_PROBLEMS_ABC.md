# Round 3: the three fundamental problems, the daughter's clock, and the Stage 1 closure

**Programme (updated master roadmap, 22 October 2026).** Answer the three fundamental problems, which are:

* **A.** Which quantum transition does the theory permit?
* **B.** Is there a more accessible formation mechanism?
* **C.** How would one verify that a daughter universe formed?

Also continue the Stage 1 publication work (NG4, the T5 strip, an independent check of S3, the SdS construction, the √2 domain). And address the critique
of the previous round: *"Reversing exterior Killing-time orientation is not by itself proof of backwards physical time."*
Stages 7–10 (experiment, engineering) are treated **theoretically only**. **No experimental or technological proposal is made or implied.**

Status labels as before:

* **[Theorem]**: proved under stated hypotheses;
* **[Proposition-N]**: established numerically, with stated coverage;
* **[Conjecture]**;
* **[Interpretation]**;
* **[Principle]**: an assumption stated as such, not derived.

**Files.**

* Code: `CODE/babyuniverse/selection.py` (Problem A), `signatures.py` (B, C), `crosscheck.py` (independent S3).
* New functions in `bounds.py`: `boundary_limit`, `boundary_certificate`, `strip_lp_certificate`, `theta_D_over_pi`, `ng4_ab_limit`, and a
  cancellation-free `delta_tE_dimless`.
* New functions elsewhere: `daughter.global_time_orientation`, `sds_atlas.euclidean_section`.
* Scripts: `CODE/scripts/round3_analysis.py` and `round3_stage1.py`.
* Data: `DATA/round3_abc.json`, `DATA/round3_stage1.json`.
* Figures: `FIGURES/round3_abc.png`, `FIGURES/round3_stage1.png`.
* Tests: `test_selection.py` (4), `test_signatures.py` (3), and new tests in `test_bounds.py` (+3), `test_daughter.py` (+1),
  `test_alternatives.py` (+1). The full suite has 143 tests.

---

## Summary

1. **Problem A, answered within the reduced theory [Theorem, conditional on three stated principles].** Assume:
   * P1, the reduced thin-wall theory;
   * P2, physical states are *geometric* and analytic in the shift strip, with each configuration counted once;
   * P3, probability is the constraint's own conserved charge, with no entropy weights.

   Then the transition amplitude is **uniquely determined**. A *charge selection rule* (Theorem A1) follows from the uniqueness of
   Laurent coefficients of the $ih$-periodic Casoratian:
   * **Double flip** (the whole physically relevant light-wall regime): $P=e^{-B_{\rm an}}(1+O(\hbar))$, with an ordinary, anchor-free flux balance.
     Here $B_{\rm an}=B_{\rm sq}$ (S4).
   * **Single flip:** the conserved current through the barrier vanishes identically. A laboratory wave is totally reflected. The daughter occurs only
     as a zero-flux standing wave.

   The last ambiguity of the previous rounds (the anchor $u_*$, B_an vs B_sq, the factor 2) disappears under P2. P3 is **equivalent** to the absence
   of a current weight $\propto e^{\pm A/4\hbar G}$ (Lemma A6). So the only remaining alternative is a Bekenstein–Hawking-weighted current, and it would
   matter only in the single-flip sliver. An exact single-flip toy, a Bessel
   difference equation (Theorem A2), confirms every step. The decay across the coalescence point equals the analytic-angle exponent $\pi a/h$ *exactly*,
   for every $h$.
2. **The daughter's clock, corrected [Theorem + correction F17].** The user's critique is right. On the expanding double-flip daughter the clock-free
   conserved charge is **positive**, with the same sign as the incident laboratory wave. Both Kruskal null coordinates increase along the daughter,
   on both sides of the wall (Proposition A5). Only the exterior Killing time $t_\infty$ runs backwards, because the daughter wall lies in region III,
   where the Killing field is past-directed. The previous round's "time-reversed daughter / Klein-type channel / Bogoliubov $|\beta|^2$" reading
   (Corollary Q10.3) is **withdrawn as physics**. It is bookkeeping in a clock that is not a global time function. The daughter's internal evolution is
   defined in the wall's proper time, in interior de Sitter time, or in any global time function. Its relation to the parent's description is through
   the clock-free charge, and through $dt_\infty/d\tau=A_{\rm out}/f_{\rm out}$, which folds once.
3. **Problem B: no accessible route [Proposition].**
   * Dynamical assistance (driving the wall) requires frequencies $\omega\gtrsim1/\Delta\tau_E$, which is $0.3$–$0.6$ of the false-vacuum energy scale
     itself ($1.6\times10^{17}$ GeV at the semiclassical frontier).
   * Theorem A1 closes the single-flip channel, the only one in which $B_{\rm sq}<B_{\rm an}$ could have helped.
   * Bubbles nucleated during inflation remain the only unsuppressed route. It is a *cosmological* route and cannot be initiated in a laboratory.
4. **Problem C: what can and cannot be verified [Theorem + Proposition].**
   * *State level:* after the event, no exterior observable distinguishes a black hole with a daughter from one without (Theorem V-2, Birkhoff).
   * *Event level:* the branches differ in their light curves (Proposition C1). On the collapse branch a radiating wall gives a late exponential tail
     $\propto e^{-u/4GM}$, carrying $\Gamma\tau_{\rm collapse}$ photons. On the daughter branch emission stops at the turning point, with a last redshift
     factor $\sqrt{f(R_1)}\approx0.74$. A late-photon count distinguishes the branches **if** $\Gamma\tau_{\rm collapse}\gg1$. It cannot show that the
     daughter inflates.
5. **T5 strip: the A/B end proved [Theorem]; the strip itself Proposition-N, with stated coverage.**
   * $\lim_{E\to E_{AB}^+}\Delta t_E^W\propto J_0(\gamma)+\pi<0$ for every $\gamma\in(-1,1)$. The proof reduces to the factorisation
     $(z^2+z+1)^2-z^2h'(z)=(z-1)(z+1)(z^2-z+3)(z^2+z+1)/z^2>0$.
   * Integration-by-parts certificates (sampling-based, re-verified on $2\times10^4$ points) hold from $1.02E_{AB}$ up.
   * Direct quadrature shows $J<0$ at all 2100 strip points down to $E/E_{AB}-1=10^{-7}$.
6. **NG4: both ends are now theorems; the interior is Proposition-N.**
   * At the A/B end, $|\theta_D|/\pi\to(|J_0|+\pi)/2\pi>1$, by the same certificate.
   * At the critical end, the Phase II closed form $F(t)>1$ holds.
   * $|\theta_D|$ is strictly increasing in $E$ at all 4800 grid points, with the cancellation-free quadrature.
7. **S3 checked independently [Proposition-N].** A from-scratch 25–35-digit implementation gives
   $B_{\rm an}=B(0)\Phi(E)\,[1+c(E)\kappa^2/\chi^2+O(\kappa^4/\chi^4)]$ with $-0.73<c(E)<0$ on all of regime A:
   * $c\to-0.722221$ as $E\to0$ (equal to $-13/18$ to six digits, [Conjecture]);
   * $c\simeq-5E^{-3/2}$ at large $E$.

   The correction is negative, so Corollary S7 carries the factor $1-0.73\kappa^2/\chi^2$. At the frontier optimum this is $-0.18\%$, and
   $P\le10^{-23}$ is unchanged.
8. **SdS (task 5) closed [Theorem].** No Euclidean period is smooth at both horizons ($\kappa_b>\kappa_c$). With both conical defects included, the
   action is $I(\beta)=-(A_b+A_c)/4G$ for **every** $\beta$. The $\beta$-dependence cancels through the SdS Smarr identity
   $\kappa_bA_b+\kappa_cA_c=4\pi\eta^2(r_c^3-r_b^3)$. Checked to $10^{-15}$, with the de Sitter limit $-3\pi/G\Lambda$.
9. **The √2 multiplier has a precise and limited domain.** It belongs to the squared (sector-blind) constraint, in the single-flip class only. P2
   excludes it, and the exact toy shows that the analytic constraint has no Stokes switching at the coalescence point.

---

## Part A — Problem A: which quantum transition does the reduced theory permit? (Stage 2 tasks 8–18)

### A0. What is assumed, and why

The three previous rounds reduced the B_an/B_sq question to one undetermined datum: *the relative weight of the Laurent sectors $q^0$ and $q^{-1}$ of the
Casoratian charge* (Q9.3), equivalently the anchor $u_*$ of Lemma III-6. Here the datum is fixed by three principles. **They are principles, not
derivations.** Deriving them from the full Wheeler–DeWitt theory remains open (A8).

* **P1 (reduced theory).** Spherically symmetric thin wall, interior de Sitter, exterior Schwarzschild. Area variable $u=R^2/2$, $h=\hbar G$. The constraint
  acts by complex shifts: $\psi(u+ih)+\psi(u-ih)=W_+(u)\psi(u)$ with $W_+=2C(u)$ and $C=(f_{\rm in}+f_{\rm out}-\kappa^2R^2)/(2\sqrt{f_{\rm in}f_{\rm out}})$.
  This is the *analytic* (complex-angle, Sorkin) constraint: in a flipped region $GP_u=i\pi n+p$.
* **P2 (geometric states, each configuration once).**
  * A physical state is analytic in the strip $|{\rm Im}\,u|\le2h$ over the static domain.
  * In every classically allowed region it is *geometric*: a superposition of WKB waves $a\,e^{-n\pi u/h}e^{\pm i\int p\,du/h}$ whose amplitude varies slowly
    across a full period strip. Here $n=n_{\rm out}-n_{\rm in}$ is that region's net flip number.
  * Each classical configuration is represented once.
  * Motivation: C5 (the complex angle is the analytic continuation of the WDW phase through the forbidden zone) and Q9.2 (the squared constraint
    contains every daughter twice, in two non-unitarily related copies).
* **P3 (probability is the constraint's own charge).**
  * Probability flux is the conserved Casoratian charge $K[\psi,\psi^\sharp](u)=\psi(u)\psi^\sharp(u+ih)-\psi(u+ih)\psi^\sharp(u)$, with
    $\psi^\sharp(u)=\overline{\psi(\bar u)}$.
  * No further measure is inserted, in particular no Bekenstein–Hawking weight $e^{A/4G}$ on flipped sectors.
  * Motivation: Lemma III-6(c). Lebesgue measure is the unique symmetric measure growing more slowly than $e^{A/4\hbar G}$.

### A1. Lemma A0 (what "geometric" excludes) — [Theorem at WKB level, to all orders in $h$]
In a region with net flip number $n$, the WKB solutions of the analytic constraint are $a\,e^{-n\pi u/h}e^{\pm iS/\hbar}$ with real $S$
($\cosh GS'=|C|$). A solution of the **squared** constraint that is regular there can be written as $\psi_++\psi_-$:

* $\psi_+$ solves the analytic constraint ($\cosh GP=C$);
* $\psi_-$ solves the sign-dressed one ($\cosh GP=-C$), whose WKB waves carry $e^{-(n\pm1)\pi u/h}$.

So the squared constraint's extra solutions are exactly the **non-geometric** waves, with the Sorkin angle mismatched by $i\pi$. On the laboratory side
this is the family $T=-W$, with $GP=i\pi\pm p$. P2 removes them. *Proof.* The symbol of the squared operator factorises as
$[\cosh GP-C][\cosh GP+C]$, and the two factors have disjoint WKB phases. (The operators themselves do not factorise, because $W$ does not commute with the
shifts; hence the WKB qualifier.) $\square$

### A2. Theorem A1 (charge selection rule) — [Theorem, under P1–P3]
*Let $\psi$ satisfy P2 on a static domain containing a laboratory region $\mathcal L$ (net flip number 0) and a daughter region $\mathcal D$ (net flip
number $n_D\in\{0,1\}$). Then:*

* *$n_D=0$ (double flip):* $K\equiv\varphi_0$ *is constant on the whole domain, and*
  $|I|^2-|R|^2=(k_D/k_L)\,|D|^2$ *with $k=\sinh GS'$ on each side;*
* *$n_D=1$ (single flip):* $K\equiv0$. *The laboratory wave is totally reflected ($|R|=|I|$). The daughter component is a standing wave with zero charge.*

*Proof.*
1. By Lemma III-6(a), $K(u+ih)=K(u)$. Analyticity of $\psi$ in $|{\rm Im}\,u|\le2h$ makes $K$ analytic for $0\le{\rm Im}\,u\le h$ (it involves
   $\psi(u+ih)$ and $\psi^\sharp(u+ih)=\overline{\psi(\bar u-ih)}$). Periodicity continues it to the whole vertical band over the static domain
   $u_{\min}<{\rm Re}\,u<u_{\max}$. Hence $K=\Phi(q)$, $q=e^{2\pi u/h}$, with $\Phi$ analytic on the annulus $e^{2\pi u_{\min}/h}<|q|<e^{2\pi u_{\max}/h}$.
   Its Laurent series $\Phi=\sum_kc_kq^k$ converges there and is unique.
2. Each $c_k=\frac1{2\pi i}\oint_{|q|=\rho}\Phi\,q^{-k-1}dq$ may be computed on *any* circle. A circle $|q|=\rho$ is the image of a horizontal segment
   ${\rm Re}\,u=u_0$, $0\le{\rm Im}\,u<h$.
3. **The cross terms vanish identically.** In each region pick one geometric solution $\psi_+$ and set $\psi_-=\psi_+^\sharp$, which is again a solution
   because $W_+$ is real-analytic. By P2, $\psi=I\psi_++R\psi_-$ on $\mathcal L$ with *constant* $I$ and $R$. Since $\psi^\sharp=\bar I\psi_-+\bar R\psi_+$
   and $K[f,f]\equiv0$,
   $$K=(|I|^2-|R|^2)\,\kappa_L,\qquad\kappa_L=K[\psi_+,\psi_+^\sharp]\quad\text{exactly on }\mathcal L,$$
   and likewise $K=(|D_+|^2-|D_-|^2)\,\kappa_D$ on $\mathcal D$. For a single geometric wave (Q9.3),
   $\kappa=(-1)^n2\sinh(GS')\,e^{-2\pi nu/h}\,[1+\delta]$, with $|\delta|\le\delta_0$ on full period segments. P2 gives $\delta_0=O(h)$; we need only $\delta_0<\frac12$.
4. Write $A_L=(|I|^2-|R|^2)\,2\sinh GS'_L$ and $A_D=(|D_+|^2-|D_-|^2)\,2\sinh GS'_D$. Taking Fourier coefficients on a laboratory circle $\rho_L$ gives
   $|c_0-A_L|\le\delta_0|A_L|$ and $|c_k|\rho_L^k\le\delta_0|A_L|$ for $k\ne0$. On a daughter circle $\rho_D$,
   $|c_{-n_D}\rho_D^{-n_D}\mp A_D\rho_D^{-n_D}|\le\delta_0|A_D|\rho_D^{-n_D}$ and $|c_k|\rho_D^k\le\delta_0|A_D|\rho_D^{-n_D}$ for $k\ne-n_D$.
5. **Case $n_D=1$.** From $c_0$: $(1-\delta_0)|A_L|\le|c_0|\le\delta_0|A_D|\rho_D^{-1}$. From $c_{-1}$: $(1-\delta_0)|A_D|\rho_D^{-1}\le|c_{-1}|\rho_D^{-1}\le\delta_0|A_L|\,\rho_L/\rho_D$.
   Together these give $|A_L|\,[(1-\delta_0)^2-\delta_0^2\rho_L/\rho_D]\le0$. The laboratory lies at smaller $u$, so $\rho_L/\rho_D=e^{-2\pi(u_D-u_L)/h}<1$, and with
   $\delta_0<\frac12$ the bracket is positive. Hence $A_L=0$, then $A_D=0$, so $K=A_L\kappa_L/(2\sinh GS'_L)\equiv0$ on $\mathcal L$, and everywhere by analyticity.
6. **Case $n_D=0$.** Both regions carry $q^0$. Uniqueness of $c_0$ gives $A_L\langle1+\delta_L\rangle=A_D\langle1+\delta_D\rangle$, where $\langle\cdot\rangle$ is the
   period average. This is the flux balance with charges $2\sinh GS'(1+O(h))$. $\square$

**Remarks.**
* The proof never chooses an anchor $u_*$ or an intertwiner. P2 forbids multiplying the daughter component by the anti-periodic factor
  $e^{\pm\pi(u-u_*)/h}$, because that would make the state non-analytic or non-geometric somewhere in the strip.
* The rule is a *superselection of the net flip number in the flux*, not of the orientation. Orientation remains configuration data (Q9.1).
* The squared constraint escapes the rule only through its non-geometric solutions (Lemma A0). That is why the previous rounds found a
  canonical unitarity relation $|I|^2=|R|^2+|T_{\rm real}|^2$ there (C3).

### A2′. Lemma A6 (what P3 actually assumes) — [Theorem]
*Let $Q[\psi,\varphi](u)$ be bilinear in two solutions and conserved, $Q(u+ih)=Q(u)$, for every pair of solutions.*

* *(i) If $Q=a\,\psi(u)\varphi(u+ih)+b\,\psi(u+ih)\varphi(u)$ (nearest-shift, off-diagonal: the form of every discretised Wronskian), then $b=-a$ and
  $a(u+ih)=a(u)$. That is, $Q=\Phi(q)\,K[\psi,\varphi]$.*
* *(ii) In general $Q=\sum_{ij}\Phi_{ij}(q)K[\psi,f_i]K[\varphi,f_j]$ for a basis $f_{1,2}$ of solutions. $Q$ is independent of the choice of reference
  solutions (invariant under periodic $SL_2$ changes of basis) iff $Q=\Phi(q)K[\psi,\varphi]$.*
* *(iii) By Lemma III-6(c), $\Phi$ is constant unless it grows at least like $q^{\pm1}$. And $q=e^{2\pi u/h}=e^{A/4\hbar G}$ is exactly the exponentiated
  Bekenstein–Hawking entropy of the wall's area.*

*Proof.*
* (i) Eliminate $\psi(u+2ih)=W\psi(u+ih)-\psi(u)$ and equate coefficients of the four independent data products. This gives $b=-a$, $a(u+ih)=a(u)$
  (`selection.current_uniqueness_certificate`).
* (ii) Bilinear forms on the two-dimensional solution space over periodic functions are spanned by $K[\cdot,f_i]K[\cdot,f_j]$. The only $SL_2$-invariant one
  is the antisymmetric (symplectic) form, which is $K[f_1,f_2]^{-1}\cdot\big(K[\psi,f_1]K[\varphi,f_2]-K[\psi,f_2]K[\varphi,f_1]\big)=K[\psi,\varphi]$. $\square$

**Consequence.** Within the reduced theory, P3 is **equivalent** to one statement: the probability current carries no weight growing like
$e^{\pm A/4\hbar G}$. The *only* alternative to Corollary A3 is a current weighted by the exponentiated horizon entropy of the wall. With such a weight,
$\Phi=q$ maps the single-flip daughter's $q^{-1}$ sector into $q^0$ and reopens that channel with $e^{-B_{\rm an}}$ (Phase III principle (ii)).
P-Q9′ is thus a yes/no question: does the reduction of the full Wheeler–DeWitt Klein–Gordon current through the bridge produce a factor $e^{A/4\hbar G}$?
Either answer leaves the double-flip result, and hence the frontier, unchanged.

### A3. Theorem A2 (exact single-flip toy) — [Theorem]
Take $C(u)=-\lambda u$, which gives laboratory $C>1$ for $u<-1/\lambda$, barrier $|C|<1$, and daughter $C<-1$ for $u>1/\lambda$:
$$\psi(u+ih)+\psi(u-ih)=-2\lambda u\,\psi(u).$$
The Bessel recurrence $Z_{\nu+1}+Z_{\nu-1}=(2\nu/x)Z_\nu$ with $\nu=-iu/h$ and $x=i/(\lambda h)$ solves it exactly. Using the order asymptotics
$J_\nu(x)\sim(x/2)^\nu/\Gamma(1+\nu)$:

| solution | daughter ($u\to+\infty$) | laboratory ($u\to-\infty$) | Laurent charge |
|---|---|---|---|
| $J_\nu$, $Y_\nu$ | grow $\sim e^{+\pi u/h}$ (non-geometric) | geometric | $q^{+1}$ in $\mathcal D$ |
| $e^{-i\pi\nu}J_{-\nu}$ | decays $\sim e^{-\pi u/h}$ (geometric, $n=1$) | grows $\sim e^{2\pi\lvert u\rvert/h}$ (non-geometric) | $q^{-1}$ in $\mathcal D$, $q^{-2}$ in $\mathcal L$ |
| $H^{(1)}_\nu$ | bounded, $\sim e^{-\pi u/h}\times$ standing wave | geometric | $K\equiv0$ |

**$H^{(1)}_\nu$ is the unique solution, up to a constant, that is geometric in both regions, and its Casoratian vanishes identically.** This is
Theorem A1 in an exactly solvable case. Numerically $|K|<10^{-55}$ at all $u$, against $|K|=0.032$ for $J_\nu$ in the laboratory region
(`selection.toy_report`).

**The exponent is the analytic angle, exactly.** The identity $H^{(1)}_{-\nu}(x)=e^{i\pi\nu}H^{(1)}_\nu(x)$ gives, for real $u$,
$$|H^{(1)}(-u)|=e^{\pi u/h}\,|H^{(1)}(u)| .$$
So the decay across the segment $[-a,a]$ is $e^{-\pi a/(\lambda h)}$ for **every** $h$. And $\int_{-a}^a\arccos(-\lambda u)\,du/h=\pi a/(\lambda h)$ is
precisely the analytic-angle ($B_{\rm an}$) exponent. The sector-blind exponent $\int\arccos|\lambda u|\,du/h$ is smaller. There is no Stokes switching
at the coalescence point $C=0$.

| $h$ | $a$ | measured $-\ln\lvert\psi(a)/\psi(-a)\rvert$ | $\pi a/h$ (analytic angle; $\frac12B_{\rm an}$ on the segment) | $\int\arccos\lvert C\rvert/h$ (sector-blind; $\frac12B_{\rm sq}$) |
|---|---|---|---|---|
| 0.1 | 0.3 | 9.42478 | 9.42478 | 8.518 |
| 0.1 | 0.5 | 15.70796 | 15.70796 | 13.151 |
| 0.1 | 0.8 | 25.13274 | 25.13274 | 18.296 |
| 0.05 | 0.8 | 50.26548 | 50.26548 | 36.592 |

Figure: `FIGURES/round3_abc.png`, panels 1–2.

### A4. Corollary A3 (the physical amplitude) — [Theorem, under P1–P3]
* **Double flip** ($\kappa<\chi$, $M<M_{\rm sd}$; S4): $P_{\rm daughter}=|D|^2k_D/(|I|^2k_L)=e^{-B_{\rm an}}(1+O(\hbar))$, with $B_{\rm an}=B_{\rm sq}$. No principle
  beyond P3 is needed. This is the whole light-wall regime except a sliver of relative width $\simeq\kappa^2/2\chi^2$ in mass (S4). It includes every
  frontier optimum of Corollary S7.
* **Single flip** ($M_{\rm sd}<M<M_{\rm cr}$, or $\kappa>\chi$): $P_{\rm daughter}=0$ as a flux. The time-symmetric state ($K=0$) has a daughter component of
  relative magnitude $e^{-B_{\rm an}/2}$, but no probability *flows* into it.
* **Consequences for earlier rounds:**
  * C4(i) "$2e^{-B_{\rm sq}}$" is superseded. It assumed the squared constraint, which P2 excludes.
  * The Phase VI bounds, which used $\min(B_{\rm an},B_{\rm sq})$, remain valid and are now conservative.
  * The frontier $P\le10^{-23}$ stands (double flip).

### A5. Proposition A4 (charge signs of physical daughters) — [Theorem + numerics]
For a WKB wave with net flip number $n$ the charge is $(-1)^n2|a|^2\sinh(GS')\,q^{-n}$, with $S'$ the real-angle momentum
(`daughter.lorentzian_momentum`). For $n=0$ its sign is the sign of $S'$. For $n=1$ it lies in the $q^{-1}$ sector and cannot be compared with the laboratory charge. `selection.charge_signs` gives:

| model | $P_{\rm lab}$ | $P_{\rm daughter}$ (expanding) | $n_D$ | $dR/dt_\infty$ on the daughter |
|---|---|---|---|---|
| single flip, $\chi=\kappa=1$, $M=0.1545$ | +0.147 | −0.897 | 1 | −0.243 |
| double flip, $\gamma=0.6$, $E=2.4$ | +0.025 | **+0.241** | 0 | −0.305 |
| double flip, $\gamma=0.3$, $E=5$ | +0.025 | **+0.460** | 0 | −0.292 |
| double flip, light wall $\kappa=0.05\chi$ | +0.056 | **+0.076** | 0 | −0.524 |

**The expanding double-flip daughter carries positive charge, the same sign as the incident laboratory wave, although $dR/dt_\infty<0$.** In the
constraint's own conserved current it is an ordinary transmitted wave. The negative $t_\infty$-velocity is a property of the clock.

### A6. Proposition A5 (global time orientation) and the daughter's internal evolution — [Theorem]
**The critique, stated precisely.** Theorem Q10.2 ($dt_\infty/d\tau=A_{\rm out}/f_{\rm out}<0$ on every daughter) is a correct statement about the *exterior Killing
time*. Corollary Q10.3 read it as "unitary boundary evolution produces the time-reversed daughter" and "the expanding daughter is a Klein-type channel".
That inference needs $t_\infty$ to be a global time function. It is not: in Kruskal coordinates $\partial_tT_K=X/4GM$ changes sign across the bifurcation
sphere.

**Proposition A5.**
*At any static point of the wall ($f_{\rm in},f_{\rm out}>0$) on any branch, in either direction of motion, both Kruskal null coordinates of the exterior,
and both of the interior de Sitter patch, increase with the wall's proper time.*

*Proof.*
1. On each side, $du/d\tau=(A-\dot R)/f$, $dv/d\tau=(A+\dot R)/f$ and $dt/d\tau=A/f$.
2. The wall lies in the static patch whose Killing field has orientation ${\rm sign}(A)$. For the exterior this is region I for $A_{\rm out}>0$ and region
   III for $A_{\rm out}<0$ (Prop. G1). For the interior it is the patch containing the centre for $A_{\rm in}>0$, otherwise the antipodal patch.
3. There the Kruskal coordinates are $U=-{\rm sign}(A)e^{-ku}$, $V={\rm sign}(A)e^{kv}$ (black-hole type), or $U={\rm sign}(A)e^{ku}$,
   $V=-{\rm sign}(A)e^{-kv}$ (cosmological type). So $dU/d\tau=k|U|\,{\rm sign}(A)(A-\dot R)/f$, and similarly for $V$.
4. The constraint $\dot R^2=A^2-f$ with $f>0$ gives $|\dot R|<|A|$, hence ${\rm sign}(A\mp\dot R)={\rm sign}(A)$, and both derivatives are positive. $\square$

Checked on single- and double-flip models (`daughter.global_time_orientation`, `test_global_time_increases_along_parent_and_daughter`). The global SdS
atlas of the previous round gives the same result numerically across all horizons (84/84 worldlines future-directed).

**The daughter's internal evolution, and its relation to the parent (the task set by the critique).**
* *Internal clock.* Use the wall's proper time $\tau$, the interior de Sitter time, or any global time function such as Kruskal $T$ or the conformal
  atlas time. All are future-directed along the whole history: laboratory → barrier → expanding daughter. Along the expanding daughter $R$ grows
  monotonically, so $R$ (or the area $8\pi u$) is also a good relational clock there. On the laboratory branch it is not, since $R$ turns around at $R_1$.
* *Relation to the parent's description.* The parent observer's $t_\infty$ is related to $\tau$ by $dt_\infty/d\tau=A_{\rm out}/f_{\rm out}$. This is positive on the parent,
  negative on the daughter, and folds once at the transition. In a stationary $t_\infty$ description, a daughter that expands in its own (global) future
  appears as $t_\infty$-incoming. That is a statement about the clock. The transition probability is defined by the clock-free charge $K$ (P3), which is
  positive for the expanding double-flip daughter (A5).
* *What survives from the previous round.* Q10.1 (the edge is removable), Q10.2 (as a statement about $t_\infty$), and D4 (a clock that is future-directed in both
  sectors cannot be $t_\infty$).
* *What is withdrawn.* The physical reading of Q10.3, together with D7's "negative boundary-clock norm". This is **F17**.
* The consistency remarks citing Freivogel et al. and Wall remain literature facts. They are no longer claimed as consequences of this model.

### A7. Stage 2 tasks 8–18: status after this round

| task | statement | status | where |
|---|---|---|---|
| 8 | full WDW current through the bridge | **open** (P-Q9′); the reduced charge is fixed by P3 | C5, A8 |
| 9 | physical inner product | reduced theory: Lebesgue $\times$ Casoratian charge, unique under P2–P3 | A0, Lemma III-6 |
| 10 | orientation sectors | orientation = configuration data; net flip number superselected in the *flux* | Q9.1, Theorem A1 |
| 11 | select $B_{\rm an}$ vs $B_{\rm sq}$ | $B_{\rm an}$ (P2); exact toy confirms | Theorems A1, A2 |
| 12 | Dirac quantization | open (P-Q1); Lemma A0 identifies what squaring adds | A1 |
| 13 | single-sheet Hamiltonian | the analytic constraint is the single-sheet object; $t_\infty$-Hamiltonian is clock-dependent (F17) | A6 |
| 14 | Picard–Lefschetz | open (P-Q3); the toy's exact solution supplies a test case | A3 |
| 15 | one-loop | open (P-T3′) | — |
| 16 | daughter clock and observables | proper time / interior dS time / global $T$; $t_\infty$ folds once | A6 |
| 17 | quantum-information accounting | which-branch entropy $h(e^{-B})$; single flip: no flux, no transfer | D7 (revised), A4 |
| 18 | GSL / singularity checks | Proposition A5 is compatible with Wall's GSL argument (no independent positive-flux baby universe without a past singularity, Prop. IV-1); no new violation | A6 |

### A8. What is not established
* **P2 and P3 are principles.** A derivation of the Casoratian charge as the reduction of the full Klein–Gordon (WDW) current through the bridge
  (P-Q9′) would turn Corollary A3 into a theorem of quantum gravity in spherical symmetry. If instead the full theory carries Bekenstein–Hawking weights
  on flipped sectors, the single-flip channel reopens with $e^{-B_{\rm an}}$ (Phase III principle (ii)). It would still be suppressed by S1–S3.
* The selection rule uses analyticity in a strip of width $2h$ over the *whole* static domain. Near the exterior horizon, where $f_{\rm out}\to0$, the WKB form
  degenerates. The laboratory circle must be taken away from it, which costs nothing because any circle in $\mathcal L$ will do.
* The toy is linear in $u$. The real $C(u)$ has turning points at finite $u$ and horizons at the ends of the static domain. The proof does not use the
  toy's special form; the toy only exhibits the mechanism exactly.

---

## Part B — Problem B: is there a more accessible formation mechanism?

### B1. What Theorem A1 removes
The single-flip channel was the only place where a smaller exponent ($B_{\rm sq}<B_{\rm an}$) or an $O(1)$ enhancement ($|\sqrt2|^2$) could appear. Under P1–P3 it
carries no flux. The double-flip channel has $B=B_{\rm an}=B_{\rm sq}\ge B_{\rm CdL}$ (S3, S5) and $\ge$ the frontier S7.

### B2. Proposition B1 (dynamical assistance needs vacuum-scale frequencies) — [Proposition]
A time-dependent drive (a field coupled to the wall, or a modulated false-vacuum energy) assists tunnelling only if it has Fourier components with
$\omega\gtrsim\omega_c=1/\Delta\tau_E$. This is the standard Keldysh/Büttiker–Landauer criterion: slower drives are adiabatic and change $B$ only through the
instantaneous parameters. Here $\Delta\tau_E=\int_{R_1}^{R_2}dR/\sqrt V$ is the proper Euclidean crossing time (`signatures.euclidean_crossing_time`). In the
thin-barrier regime $\Delta\tau_E\simeq(\pi/3)R_{\rm CdL}$.

| case | $\Delta\tau_E$ | thin estimate | $\omega_c$ | $\omega_c/E_{\rm vac}$ |
|---|---|---|---|---|
| semiclassical frontier ($B_{\rm CdL}=10$, $M=m_P$, $R=100\ell_P$) | $78.1\,\ell_P$ | $74.4\,\ell_P$ | $1.6\times10^{17}$ GeV | 0.58 |
| GUT ($10^{16}$ GeV, $M=10m_P$) | $3851\,\ell_P$ | $3833\,\ell_P$ | $3.2\times10^{15}$ GeV | 0.32 |

**Assistance requires driving the wall at a large fraction of the false-vacuum energy scale itself.** A drive at that scale is itself a source of the
energy density whose absence made the problem hard. The exponent reduction available from an assisting field of frequency $\omega\ll\omega_c$ is
exponentially small in $\omega_c/\omega$.

### B3. Ledger update (adds to Part E of the previous round)

| mechanism | assumption relaxed | cost | verdict |
|---|---|---|---|
| single-flip channel, $B_{\rm sq}$ | P2 | closed by Theorem A1 | no flux |
| dynamically assisted tunnelling | adiabaticity | $\omega\gtrsim0.3$–$0.6\,E_{\rm vac}$ | not accessible |
| bubbles nucleated during inflation (GVZ) | laboratory initial data | $f(M)\simeq\frac{4\pi}3\lambda(M_{\rm eq}/M)^{1/2}$; $B_{\rm inf}\approx37$ for all dark matter at $10^{20}$ g | viable only cosmologically |

**On the user's remark.** A mechanism that works during inflation is not a laboratory mechanism, and the analysis agrees. In the GVZ route the
*parent* is an inflating de Sitter space, and the bubble's mass exceeds $M_{\rm cr}$ by cosmological expansion, not by preparation. Neither can be
reproduced in an asymptotically flat parent: Corollary IV-2 and the dS-parent cost $e^{-3\times10^{122}}$ (Part E).

### B4. Verdict on Problem B
Within GR plus quantum fields that obey QEIs, with the reduced quantum theory of Part A, **no mechanism more accessible than $P\le10^{-23}$ per
semiclassically controlled laboratory bubble has been found**. The remaining unsuppressed routes are cosmological (GVZ), hypothetical (limiting
curvature, new NEC-violating sectors), or depend on Planckian physics ($R\lesssim20\,\ell_P$).

---

## Part C — Problem C: how would one verify that a daughter formed?

### C1. State-level no-go — [Theorem] (Theorem V-2, Birkhoff)
After the transition the exterior is Schwarzschild of mass $M$ on both branches. Every exterior observable (the abelian algebra generated by $M$ and
fields outside the horizon) has the same expectation value in the two final states, up to the Hawking-partner question. **No measurement made after the
event on the exterior can certify the daughter.**

### C2. Proposition C1 (event-level signature: light-curve truncation) — [Proposition]
The branches differ in their *history*. Let the wall emit radially outward at a constant proper rate $\Gamma$ and rest-frame frequency $\nu_e$. With
$u=t-r_*$ and $du/d\tau=(A_{\rm out}-\dot R)/f_{\rm out}$:

* **Collapse branch:** received rate $\Gamma\,d\tau/du$ and frequency $\nu_ed\tau/du$. Both decay as $e^{-u/4GM}$ near the horizon (fit: e-fold 0.1716 against
  $4GM=0.1721$). The number of photons arriving after the turning-point signal is $N_{\rm late}=\Gamma\tau_{\rm collapse}$.
* **Daughter branch:** after the transition the wall is in region III, and no future-directed null ray from it reaches region I. Emission stops at
  $u_0$, the arrival time of the turning-point signal, with a finite last redshift factor $\sqrt{f_{\rm out}(R_1)}$. Hence $N_{\rm late}=0$.

| model | $\Gamma\tau_{\rm collapse}$ per unit $\Gamma/\chi$ | last redshift factor | tail e-fold / $4GM$ |
|---|---|---|---|
| $\gamma=0.6$, $E=2.4$ | 0.096 | 0.737 | 0.997 |
| $\gamma=0.3$, $E=5$ | 0.115 | 0.756 | 1.004 |

`signatures.light_curve`, `late_photon_fraction`. Figure: `FIGURES/round3_abc.png`, panel 3. The light curve is continuous at the turning point
(`test_light_curve_is_continuous_at_turning_point`).

**What it establishes.**
* It is an exterior, event-level discriminator *between the two branches of one prepared bubble*, with signal-to-noise set by
  $\Gamma\tau_{\rm collapse}=\Gamma\times O(R_1)$. It requires a wall that radiates at a known rate, and $\Gamma\tau_{\rm collapse}\gg1$.
* It does **not** establish that the region behind the horizon inflates. A wall that stops radiating for any other reason, or a collapse whose
  tail is absorbed, mimics it.
* In the controlled regime ($R\sim100\,\ell_P$, $\Gamma\lesssim1/\ell_P$), $\Gamma\tau_{\rm collapse}\lesssim10^2$, while the event rate is $\le10^{-23}$.

### C3. Success criterion (theoretical, Stage 6 task 46)
A claim "a daughter universe formed" would need, jointly:

1. a prepared bubble of known $M$, $\sigma$ and $\rho$, inside the double-flip class;
2. a radiating wall with $\Gamma\tau_{\rm collapse}\gg1$, and an observed truncation, i.e. zero late photons where $\Gamma\tau_{\rm collapse}$ are expected;
3. an unchanged exterior mass (Birkhoff), and no burst beyond the transition's own field-theoretic emission (C4);
4. a frequency of truncation events per prepared bubble consistent with $e^{-B_{\rm an}}$, which requires $\gtrsim e^{B}$ trials.

Item 4 alone makes the criterion unattainable at $B\ge53$. **Stating this is not a proposal to attempt it.**

### C4. Candidate transient (open)
The transition changes the wall's state discontinuously in the parent's description: the momentum jumps between the branches. Quantum fields coupled to the
wall then see a sudden change, and their Bogoliubov coefficients would give a burst at the transition, with a spectrum set by $1/\Delta\tau_E$. This is
**not computed here** (it needs the field-theoretic one-loop analysis, P-T3′). It is listed because it is the only candidate signal that would be
present on the daughter branch and absent on the collapse branch, rather than the reverse.

### C5. Cosmological populations (Stage 6 tasks 40–45)
In the GVZ route the observable is the mass function of the hosting black holes, $f(M)\propto M^{-1/2}$ below $M_{\rm eq}$, with daughters only for
$M>M_{\rm cr}$. A population signal tests the inflationary bubble-nucleation hypothesis, not the presence of a daughter in any particular hole. By C1,
individual daughters remain unverifiable.

---

## Part D — Stage 1 closure (tasks 1–7)

### D1. Task 2: the T5 strip $E_{AB}<E<1.6E_{AB}$
**Proposition S5-b (A/B end) — [Theorem, with the pinch limit derived by matched asymptotics and checked to $10^{-7}$].** For every $\gamma\in(-1,1)$,
$$\lim_{E\to E_{AB}^+}J(\gamma,E)=J_0(\gamma)+\pi<0,\qquad J_0=-\int_1^{z_2}\frac{(z^2+z+1)\,dz}{z\sqrt{E_{AB}-g}}=-\int_0^{\pi/2}F(z(\varphi))\,d\varphi .$$
Here $h(z)=(z^3-1)(z^2+z+1)/z^3$ is defined by $E_{AB}-g=(z-1)(E_{AB}-h)/z$, the substitution is $h(z(\varphi))=E_{AB}\sin^2\varphi$, and
$F=2(z^2+z+1)^2/(z^2h'(z))$.

* **The $+\pi$.** At $E=E_{AB}$ the outer horizon $W=Ef_{\rm out}=0$ pinches the turning point $s=0$ at $z=1$. For $E=E_{AB}+\epsilon$ the integrand
  $-E\,sq/(W\sqrt{W-s^2})$ has a spike of width $O(\epsilon^2)$ at $z_1$. With $W\simeq\epsilon+bx$, $s\simeq3x$ and $q\simeq\frac13$, it contributes
  $3\cdot\frac{\epsilon}{b}\cdot\frac{\pi}{c}=\pi$ with $c=3\epsilon/b$.
* **The sign.** $F(1)=2$ and
  $$(z^2+z+1)^2-z^2h'(z)=\frac{(z-1)(z+1)(z^2-z+3)(z^2+z+1)}{z^2}>0\quad(z>1),$$
  so $F>2$ on $z>1$, $|J_0|>\pi$, and the limit is negative. $\square$ (`bounds.boundary_certificate`)
* Values: $\gamma=-0.99,-0.5,0,0.5,0.9,0.99$ give $-1.461,-1.134,-0.760,-0.370,-0.071,-0.0070$.

**Certificates in the interior of the strip [Proposition-N].** Theorem S5 is integration by parts with one multiplier, $\varphi=k=q/W$. Allowing
$\varphi=\sum_j(a_jk+b_jq+c_jq\sqrt W)x^j$ (degree 5) and maximising the minimum of the modified integrand by linear programming
(`bounds.strip_lp_certificate`; 800 points, re-verified on $2\times10^4$):

| $\gamma$ | $1.01E_{AB}$ | $1.02E_{AB}$ | $1.05E_{AB}$ | $1.3E_{AB}$ | $1.6E_{AB}$ |
|---|---|---|---|---|---|
| −0.9 | −0.002 | **0.017** | 0.054 | 0.247 | 0.257 |
| 0.0 | −0.005 | **0.017** | 0.050 | 0.207 | 0.223 |
| 0.9 | **0.005** | 0.013 | 0.026 | 0.079 | 0.100 |

Positive margins certify $J<0$ on the samples. No pointwise multiplier works arbitrarily close to $E_{AB}$: the positive spike (+π) must be beaten by
the negative bulk, which is a non-local comparison. That is exactly what Proposition S5-b handles. **Direct quadrature** (cancellation-free, below)
gives $J<0$ at all 2100 points of $\gamma\in[-0.98,0.98]$ × $E/E_{AB}-1\in[10^{-7},0.62]$, with maximum $-0.0140$ at $\gamma=0.98$, equal to the boundary limit.

**Domain statement for T5 (double-flip class, $\Delta t^W_E<0\Rightarrow B\ge B_{\rm CdL}$):**
* *Theorem* for $E\ge E_\omega(\gamma)$ ($1.27$–$1.61\,E_{AB}$), wherever the pointwise $D\ge0$ holds, and in the limit $E\to E_{AB}^+$;
* *Theorem* for light walls at fixed $E$ (S3);
* *Proposition-N* on $E_{AB}<E<E_\omega$: sampled certificates from $1.02E_{AB}$, and direct quadrature on the full strip.

**Numerical correction (F18).** The previous double-precision quadrature lost accuracy at the pinch. $E-g$ and $f_{\rm out}$ were computed as
differences of $O(1)$ numbers where both are $O(\epsilon^2)$. `delta_tE_dimless` now integrates in $y=z-z_1$ with
$E-g=-y[2z_1+y+2\gamma/(z_1z)-(z_1+z)(z_1^2+z^2)/(z_1^4z^4)]$ and $W=s(z_1)^2+by/(z_1z)$, breaking the range geometrically on the spike width. A 30-digit
mpmath evaluation agrees to $10^{-10}$. A preliminary LP run in this round had reported certificates at $1.0001E_{AB}$ with margin 1.0 (the cap). That was
an undersampling artefact and is withdrawn.

### D2. Task 1: NG4
In regime B, $|\theta_D|/\pi=-EJ/(2\pi E_{AB})$ (`bounds.theta_D_over_pi`).
* **Critical end [Theorem, Phase II]:** $\to F(t)>1$.
* **A/B end [Theorem, new]:** from below, the pinch contributes $-\pi$ to $J$ (the spike lies at $s>0$ and adds to the negative bulk). So
  $|\theta_D|/\pi\to(|J_0|+\pi)/2\pi$, which is $>1$ exactly because $|J_0|>\pi$, the certificate of D1 (`bounds.ng4_ab_limit`).

| $\gamma$ | −0.9 | −0.5 | 0 | 0.5 | 0.9 | 0.99 |
|---|---|---|---|---|---|---|
| critical end $F(t)$ | 1.1864 | 1.1554 | 1.1088 | 1.0557 | 1.0111 | 1.00111 |
| A/B end $(\lvert J_0\rvert+\pi)/2\pi$ | 1.2235 | 1.1805 | 1.1209 | 1.0588 | 1.0113 | 1.00111 |

* **Interior [Proposition-N]:** strictly increasing in $E$ at all 4800 points (60 values of $\gamma$ × 80 Chebyshev-spaced $E$). The minimum, at the critical
  end, is 1.0022. The endpoint values match the two limits to $1.5\times10^{-5}$.

NG4 is therefore reduced to a monotonicity statement with both boundary values proved. The previous double-precision scan briefly showed spurious
non-monotone steps near $E_{AB}$, which were the F18 artefact. A 30-digit check at $\gamma=0.9$ confirms monotone approach to the limit 1.01126177.

### D3. Task 3: independent check of S3 and its exact range
`crosscheck.py` rebuilds the junction data, the barrier, $B_{\rm an}$ (atan2 angle), $B_{\rm sq}$ and $\Phi$ from scratch in mpmath (25–35 digits). It imports
nothing from `model.py`, `tunnelling.py` or `bounds.py`. The angle identity $\cos^2\theta+\sin^2\theta=1$ holds symbolically.

| $E$ | $c(E)=(B_{\rm an}/B_0\Phi-1)/(\kappa/\chi)^2$ at $\kappa/\chi=0.1$ | 0.03 | 0.01 |
|---|---|---|---|
| $1.5E_{AB}$ | — | −0.7210 | −0.7221 |
| 0.3 | −0.6641 | −0.6692 | −0.6697 |
| 1 | −0.5565 | −0.5608 | −0.5612 |
| 3 | −0.3479 | −0.3507 | −0.3510 |
| 30 | −0.0276 | −0.0278 | −0.0278 |

**Result [Proposition-N].** The law holds uniformly on all of regime A ($E>E_{AB}$), with
$$B_{\rm an}=B(0)\,\Phi(E)\,\big[1+c(E)\,\kappa^2/\chi^2+O(\kappa^4/\chi^4)\big],\qquad-0.73<c(E)<0 .$$
* $c(0^+)=-0.7222205$ at $\kappa/\chi=10^{-3}$. **[Conjecture]** $c(0^+)=-13/18$. Note that $B(0)=\frac{\pi\kappa R_0^3}{8G}(1+\kappa^2/\chi^2+\dots)$, so this is not
  simply the correction to $B(0)$.
* $c(E)\simeq-5E^{-3/2}$ for $E\gtrsim30$.
* All sampled points are double flip ($B_{\rm an}=B_{\rm sq}$ to 12 digits).
* The single-flip sliver $E_{AB}<E<E_{\rm sd}\simeq E_{AB}(1+\kappa^2/3\chi^2)$ is not resolved at these $\kappa$ and carries no flux anyway (Theorem A1).

**Consequence for Corollary S7.** The correction is negative, so the rigorous form is
$B_{\rm lab}\ge1.671\sqrt{B_{\rm CdL}MR_M}\,(1-0.73\kappa^2/\chi^2)$. At the frontier optimum ($\kappa/\chi=0.050$) the bound moves from 52.84 to 52.75.
The exact exponent there is 59.5, so $P\le10^{-23}$ is unaffected.

### D4. Task 5: the SdS construction — complete
* **Both horizons and the causal structure:** the eight-block Kruskal atlas of the previous round (`sds_atlas.py`). Every branch is timelike,
  future-directed and continuous (84/84).
* **Euclidean periods and conical defects [Theorem].** $\beta_b=2\pi/\kappa_b<\beta_c=2\pi/\kappa_c$ for every sub-Nariai SdS, so every Euclidean section carries
  at least one conical singularity. With deficits $\delta_h=2\pi-\kappa_h\beta$, $\int\sqrt gR=2\delta_hA_h$ on each defect, and $R=4\Lambda$ in the bulk:
  $$I(\beta)=-\frac{\Lambda\beta}{8\pi G}\frac{4\pi}{3}(r_c^3-r_b^3)-\frac1{16\pi G}\sum_h2\delta_hA_h=-\frac{A_b+A_c}{4G}+\frac{\beta}{8\pi G}\Big[\kappa_bA_b+\kappa_cA_c-4\pi\eta^2(r_c^3-r_b^3)\Big],$$
  and the bracket vanishes identically (SdS Smarr). **$I=-(A_b+A_c)/4G$ for every $\beta$**, smooth or not.
  * Checks: three $(M,\eta)$ and four $\beta$ each, with $\beta$-independence to $10^{-15}$. $R$ is recomputed from the metric. The de Sitter limit is
    $-\pi/(G\eta^2)=-3\pi/G\Lambda$ (`sds_atlas.euclidean_section`, `test_sds_euclidean_action_is_period_independent`).
  * Consequence: the SdS tunnelling exponents of Phase II (T3′), built from horizon-area terms, do not depend on a choice of period. The conical-defect
    prescription and the Hamiltonian one agree.

### D5. Task 6: the domain and significance of the √2 Stokes result
* **Domain.** The multiplier $\sqrt2$ (Phase III Prop. III-8, confirmed by exact integration in Part B of the previous round) belongs to the **squared**
  constraint. It enters only where $C$ crosses 0 inside the barrier, which is the **single-flip class** ($\theta$ runs from 0 to $\pi$). In the
  double-flip class $C\ge0$ on the barrier, and there is no coalescence point.
* **Under P2** the squared constraint is excluded (Lemma A0). The exact toy (A3) shows that the analytic constraint has no switching at $C=0$: the decay
  is $e^{-\pi a/h}$ exactly, i.e. $\nu=0$ in the language of the previous round's controls.
* **Global significance.** In the squared theory $\sqrt2$ produced the factor 2 in $2e^{-B_{\rm sq}}$, an $O(1)$ prefactor in a sliver of relative mass width
  $\kappa^2/2\chi^2$. In the analytic theory the same sliver carries no flux. In neither case does it affect a frontier or an exponent.

### D6. Task 4: literature
Full texts remain unavailable from this environment. The environment's network policy rejects arXiv, INSPIRE, Semantic Scholar, OSTI, CORE and the
publishers. The comparison of the previous round (abstract level: FGG, FMP 1990, Ansoldi–Tanaka 2014/15, Sakai et al. 2006, GVZ 2016) stands. No literature
result is claimed beyond it. The line-by-line comparison stays open (P-L1).

### D7. Task 7: manuscript
`PAPER/genesis_foundations.tex` gains a second addendum: Theorem A1 and Corollary A3, the clock correction F17, Proposition S5-b, NG4 at both ends, the S3
correction coefficient, and the Euclidean SdS action.

---

## Part E — Stages 3 and 5–10 (theory only): what changed
* **Stage 3 (tasks 19–25).** Unchanged requirements. The thick-wall first law and Smarr relation (P-T2) and aspherical stability (P-D2) remain open.
  Proposition A5 gives the time orientation needed to pose the backreaction problem on the daughter in a global time.
* **Stage 5 (five numbers).** Unchanged: $E_{\rm req}\approx3\times10^{17}$ GeV, $L_{\rm req}\approx100\,\ell_P$, $\tau_{\rm req}\lesssim5\times10^{-42}$ s, $P\le10^{-23}$.
  $S_{\rm verify}$ is refined from 0 at state level to the event-level truncation of C2, at signal level $\Gamma\tau_{\rm collapse}$. These remain
  conditional predictions of the investigated model.
* **Stages 7–10.** Nothing to add. Every requirement is set by Stage 5 and by the success criterion C3, and none is within reach of any conceivable technology.
  **No experimental or technological proposal is made or implied.**

---

## Part F — Corrections and failed approaches of this round
* **F17 (previous round, Corollary Q10.3 and D7).** "Unitary boundary-clock evolution yields the time-reversed daughter; the physically expanding daughter
  is a Klein-type, negative-norm channel; creation is a Bogoliubov process with $|\beta|^2\simeq e^{-B}$." *Failure:* this inferred physical time
  reversal from the orientation of the exterior Killing time, which is not a global time function. The clock-free conserved charge of the expanding
  double-flip daughter is positive (A5), and every global time function increases along it (Proposition A5). The $t_\infty$-flux bookkeeping remains
  correct as bookkeeping. **Type:** interpretational error, identified by the user's critique.
* **F18 (this round, numerics).** Double-precision quadrature at the A/B pinch gave spurious non-monotone NG4 steps, $10^{-3}$-level errors at
  $\epsilon\lesssim10^{-5}$, and a spurious LP certificate at $1.0001E_{AB}$. Fixed by the cancellation-free integrand and a 30-digit cross-check.
* **Superseded:** C4(i) "$P=2e^{-B_{\rm sq}}$ (single flip)" (by Theorem A1 under P2); the summary item 8 and D3 of the previous round (by F17).

## Open problems after this round
* **P-Q9′.** Derive P3, the Casoratian charge as the reduction of the full WDW Klein–Gordon current through the bridge. This is the last step to a theorem
  of spherically symmetric quantum gravity.
* **P-N1′.** NG4: monotonicity of $EJ(E)$ on regime B (both boundary values are proved).
* **P-T5″.** The strip $E_{AB}<E<E_\omega$: a rigorous (interval-arithmetic or analytic) version of the sampled certificates, plus a quantitative remainder for
  the boundary limit.
* **P-S3.** Derive $c(E)$ (conjecture $c(0^+)=-13/18$; tail $-5E^{-3/2}$).
* **P-C4.** The field-theoretic burst at the transition (one-loop, P-T3′).
* **P-L1.** Full-text literature comparison (requires network access to the literature hosts).
* **P-O1/O2.** GVZ formation efficiencies and the supercritical fraction.
