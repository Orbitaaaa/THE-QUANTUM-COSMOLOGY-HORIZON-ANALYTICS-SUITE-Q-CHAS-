# The Farhi–Guth–Guven transition: reconstruction and new results

Status labels as in NO_GO_THEOREMS.md. Code: `CODE/babyuniverse/tunnelling.py`; data:
`DATA/B_of_M.json`, `DATA/M0_closed_form_grid.json`, `DATA/smarr_identity_grid.json`,
`DATA/first_law_checks.json`, `DATA/laboratory_regime_bounds.json`, `DATA/physical_scales.json`;
figures: `FIGURES/tunnelling_exponent.png`, `FIGURES/euclidean_wall_exterior.png`,
`FIGURES/pseudomanifold_multiplicity.png`, `FIGURES/laboratory_regime_bounds.png`.

---

## 1. The problem

* **Initial state:** the bounded branch of regime A at its turning point $R_1$, with
  $\varepsilon_{\rm in}=\varepsilon_{\rm out}=+1$. The time-symmetric slice is regular laboratory data: no extremal spheres,
  asymptotically flat, wall in region I (GLOBAL_GEOMETRY.md, NO_GO_THEOREMS.md Sec. 2).
  Classically the bubble recollapses into a black hole.
* **Final state:** the bounce branch at $R_2$ with $\varepsilon_{\rm out}=-1$ (always, Proposition G1). This is
  bag-of-gold data with an Einstein–Rosen bridge, and it evolves into an inflating, causally
  disconnected daughter (criteria D1–D3).
* **Forbidden region:** $R_1<R<R_2$, where $V>0$. Same exterior mass $M$; energy is conserved.
* **Question:** is there a consistent amplitude for $R_1\to R_2$, and what is it?

The classical obstruction (Farhi–Guth) forbids a Lorentzian interpolation. NG2 shows that the final
state's classical past is singular, so a transition must replace that past.

## 2. Canonical (reduced phase space) formulation

Solve the constraints classically, keep $(R,P)$, and take the asymptotic Killing time $t$ as the clock,
with Hamiltonian $H=M$ (the ADM mass).

**Momentum.** Hamilton's equation $\partial H/\partial P=dR/dt=\dot R f_{\rm out}/A_{\rm out}$ is equivalent to
$\partial P/\partial M|_R=A_{\rm out}/(\dot Rf_{\rm out})$. It is solved by
$$
P(R,\dot R)=\frac RG\Big[\operatorname{artanh}\frac{\dot R}{A_{\rm out}}-\operatorname{artanh}\frac{\dot R}{A_{\rm in}}\Big]+F(R),
$$
with $F$ an arbitrary real function, i.e. a canonical transformation that does not affect any exponent.
Verified symbolically with exact rational arithmetic (`symbolic.reduced_hamiltonian`). The inverse
Hamiltonian $H(R,P)$ was checked numerically in all three static sign sectors that occur,
$(+,+)$, $(+,-)$, $(-,-)$, with $dR/dt=\partial H/\partial P$ to $6\times10^{-10}$
(`quantum.hamilton_equation_check`, `DATA/quantum_checks.json`).

**Coupling-constant Hamilton–Jacobi relations** (derived by implicit differentiation at fixed $R,M$; these
are new in this form as far as we know):
$$
\frac{\partial P}{\partial M}\Big|_R=\frac{dt_{\rm out}}{dR},\qquad
\frac{\partial P}{\partial\Lambda_{\rm in}}\Big|_{R,M}=-\frac{1}{8\pi G}\,\frac{4\pi R^3}3\,\frac{dt_{\rm in}}{dR},\qquad
\frac{\partial P}{\partial\sigma}\Big|_{R,M}=-4\pi R^2\,\frac{d\tau}{dR}.
$$
(Derivation: with $\delta A_{\rm in}=\delta A_{\rm out}$ from the junction and $\partial_u\operatorname{artanh}(u/A)=1/A$,
$\partial_f\operatorname{artanh}(u/A)=-u/(2Af)$, one finds $\delta P=(R/G)\,\delta f_{\rm in}A_{\rm in}/(2uf_{\rm in})$ for
$\delta\chi^2$ and $\delta P=-(R^2/Gu)\,\delta\kappa$.)

**Lorentzian slicing problem.** Inside horizons ($f<0$), $|\dot R|>|A|$ and the artanh's become complex. The
asymptotic-time clock is singular wherever the wall crosses a horizon. This is the canonical face of
the FGG "slicing" difficulty. In the barrier, by contrast, $f\ge V>0$ on both sides (since $f-V=A^2\ge0$),
so the clock is regular throughout the tunnelling region.

## 3. Euclidean momentum and the branch problem

Under the barrier $\dot R=i\sqrt V$ and $P=iP_E$ with
$$
P_E(R)=\frac RG\big[\Theta_{\rm out}(R)-\Theta_{\rm in}(R)\big],\qquad (A,\sqrt V)=\sqrt f\,(\cos\Theta,\sin\Theta),
$$
and $\cos(\Theta_{\rm out}-\Theta_{\rm in})=(f_{\rm in}+f_{\rm out}-\kappa^2R^2)/2\sqrt{f_{\rm in}f_{\rm out}}$.
The angle is defined only up to its branch. There are two candidate prescriptions:

* **analytic:** $\Theta={\rm atan2}(\sqrt V,A)\in[0,\pi]$;
* **principal:** $\Theta=\arctan(\sqrt V/A)\in(-\pi/2,\pi/2)$, which jumps by $\pi$ where $A$ changes sign.

**[Proposition] T0 (selection of the branch).** Three independent arguments select the analytic prescription.

1. *Analyticity.* As a function of complex $R$, $P$ has square-root branch points at the turning points
   and logarithmic branch points where $\dot R/A=\pm1$, i.e. where $f=0$ (the horizons; on shell $A^2-\dot R^2=f$).
   The real segment $(R_1,R_2)$ contains none of them, since $f>0$ there. The flip radius $R_f$, where $A=0$, is a regular
   point: $\arctan y$ is analytic through $y=\infty$, and its branch points are $y=\pm i\Leftrightarrow f=0$. The WKB phase
   is the analytic continuation of the classical Hamilton–Jacobi function along the barrier, and that
   continuation is the atan2 branch. The principal prescription equals the analytic one plus a horizon
   monodromy inserted at $R_f$. Around $f=0$, $\operatorname{artanh}$ shifts by $i\pi$, so $P\to P\pm i\pi R/G$, and the Killing
   time, which has a simple pole of residue $\propto1/f'(r_h)$, shifts by $\pm i\beta/2$ ($\beta=2\pi/\kappa_s$). **The
   FGG slicing jump, a slice rotated by $\pi$ around the Euclidean horizon, is exactly this horizon monodromy.**
2. *Positivity.* For tunnelling out of a state localised on the bounded branch, the WKB solution must decay
   under the barrier, so $B>0$. The principal prescription gives $B<0$ in 41 of 45 grid models; the analytic one
   gives $B>0$ in all.
3. *$M\to0$ limit.* Only the analytic prescription tends to the action of a genuine manifold instanton (T1).

**Phase III refinement of T0** (PHASE_III_QUANTUM_TRANSITION.md).
* The analytic branch is the momentum written with Sorkin's complex boost angles: a flipped normal has angle $i\pi$ plus real.
  The principal branch is the real-angle (artanh) momentum.
* Analytic continuation of the *laboratory* momentum through the barrier lands on the complex-angle momentum of the daughter
  sector. So argument 1 shows that analyticity plus laboratory initial data imply the complex convention on the daughter side.
* What argument 1 does not fix is the normalisation of the daughter sector relative to the laboratory sector (Lemma III-6).
  The FGG ambiguity survives there, as the anchor $u_*$ of an anti-periodic intertwiner $e^{\pm\pi(u-u_*)/\hbar G}$.
* The principal prescription admits no analytic solution of its own quantum constraint across the flip (Theorem III-5).
* The true Hamiltonian has a square-root branch point exactly at the flip radius, where the trajectory touches the branch
  locus tangentially. That point is also the fold of the Euclidean wall in Killing time (Theorem III-4).

**Definition.** $B\equiv2\int_{R_1}^{R_2}P_E\,dR$ (full bounce), and the transition probability per attempt
is $\Gamma\propto e^{-B/\hbar}$. The prefactor is not computed here (OPEN_PROBLEMS.md P-T3).

## 4. Results

### T1. Closed form at $M=0$ and a Morse-theoretic identity — [Theorem]

At $M=0$, $V=1-R^2/R_0^2$ with $R_0=2\kappa/L^2$. Let $p=-A_{\rm out}(R_0)$ and $q=A_{\rm in}(R_0)$ be the signed
cosines of the polar angles of the two 4-sphere caps of the compact $O(4)$ instanton (a flat 4-ball if
$\eta=0$), with $p+q=\kappa R_0>0$. Then
$$
S_E[\text{instanton}]=-\frac{\pi R_0^2}{2G}\Big[\frac1{1+p}+\frac1{1+q}\Big],\qquad
B=\frac{\pi R_0^2}{2G}\,\frac{p+q+2p_+q_+}{(1+|p|)(1+|q|)},
$$
$p_+=\max(p,0)$, $q_+=\max(q,0)$, and equivalently
$$
\boxed{\;B=S_E[\text{instanton}]+\frac1{4G}\Big(\sum_{\rm maximal\ spheres}A-\sum_{\rm minimal\ spheres}A\Big)\;}
$$
over the extremal 2-spheres of the instanton's time-symmetric equatorial slice. These are the de
Sitter equator of a cap larger than a hemisphere, or the wall when both caps are smaller than
hemispheres.

*Proof.* $B=(2R_0^2/G)[J(-p)-J(q)]$ with $J(s)=\int_0^1x\,{\rm atan2}(\sqrt{1-x^2},sx)\,dx$, evaluated by parts:
$J(s)=\pi/4(1+s)$ for $s\ge0$ and $\pi/2-\pi/4(1-s)$ for $s<0$ (checked against quadrature to $10^{-14}$).
$S_E$ follows from the on-shell action $-\frac{\Lambda}{8\pi G}{\rm Vol}-\frac\sigma2A_W$ with cap volumes
$\frac{2\pi^2}{H^4}(\frac23-c+\frac{c^3}3)$ and $\kappa R_0=p+q$. Comparing case by case gives the Morse form; the case
$p,q<0$ is impossible. $\square$

Checked on 60 parameter sets including de Sitter exteriors: numerical and closed form agree to
$7\times10^{-15}$ (`DATA/M0_closed_form_grid.json`). Consequences:

* If the false vacuum is the larger cap ($q<0$; typical for $\kappa<\chi$), then
  $B(M\to0)=S_E+\pi/G\chi^2=B_{\rm CdL}$, the **Coleman–De Luccia decay exponent of the false vacuum itself**.
  The thick-wall CdL bounces of Sec. 6 converge to this value.
* If both caps are small ($\kappa\gtrsim\chi$), then $B(M\to0)=S_E+A_W/4G<B_{\rm CdL}$.
* $B(M\to0)>0$ always.

### T2. First law for the FGG exponent — [Theorem]

Let $\tilde B=B-\mathcal A/4G$, where $\mathcal A$ collects the endpoint areas: $\pm4\pi R_i^2$ for an endpoint at which
the wall is a maximal/minimal sphere of its slice, plus $4\pi/\chi^2$ if the final slice contains the de Sitter
equator (`tunnelling.endpoint_area_term`). Within a regime,
$$
\boxed{\;d\tilde B=\frac{{\rm Vol}_E^{\rm in}}{8\pi G}\,d\Lambda_{\rm in}+A_E^W\,d\sigma-\Delta t_E^{W}\,dM\;}
$$
where ${\rm Vol}^{\rm in}_E$ is the Euclidean 4-volume of the interior region swept by the full bounce,
$A^W_E$ the Euclidean area of the wall worldvolume, and $\Delta t^W_E$ the exterior Euclidean Killing time
elapsed along the wall ($dt_E/d\tau_E=A_{\rm out}/f_{\rm out}$; it is negative in all models studied).

*Proof.* Differentiate $B=2\int P_E\,dR$ with respect to each parameter. Endpoint terms are
$2P_E(R_i)\,dR_i$ with $P_E(R_i)=\pi R_i n_i/G$, $n_i=n_{\rm out}-n_{\rm in}\in\{0,\pm1\}$; these are exactly
$d(\mathcal A/4G)$. The integrands are the Euclidean continuations of the three Hamilton–Jacobi
relations of Sec. 2. For the interior volume, the Stokes form $(r^3-r_h^3)/3\,dt_E$ (regular at the de Sitter
Euclidean horizon) plus the axis contribution reproduces the $4\pi/\chi^2$ endpoint term. $\square$

Finite-difference check (`DATA/first_law_checks.json`): all three relations hold to 8 significant
digits in four models covering $\gamma<0$, $\gamma=0$, $\gamma>0$ with and without interior flips.

### T3. Smarr identity: a covariant meaning for the canonical exponent — [Theorem]

$$
\boxed{\;B=-\frac{\Lambda_{\rm in}}{8\pi G}{\rm Vol}^{\rm in}_E-\frac\sigma2A^W_E-\frac M2\,\Delta t^W_E+\frac{\mathcal A}{4G}\;}
$$

*Proof.* At fixed $G$, $\tilde B(\lambda^{-2}\Lambda,\lambda^{-1}\sigma,\lambda M)=\lambda^2\tilde B(\Lambda,\sigma,M)$
(lengths scale as $\lambda$). Euler's relation $2\tilde B=-2\Lambda\partial_\Lambda\tilde B-\sigma\partial_\sigma\tilde B+M\partial_M\tilde B$
combined with T2 gives the result. $\square$

Numerical check over a 45-model grid covering regimes A and B, $\gamma\in[-0.85,0.85]$, $M/M_{\rm cr}$ from $\sim10^{-2}$
to $0.99$: maximum relative residual $2.5\times10^{-8}$ (`DATA/smarr_identity_grid.json`). This includes models whose
exterior wraps the Euclidean horizon dozens of times. The first two terms are the on-shell Euclidean action
of the interior region and of the wall (bulk $\mathcal R=4\Lambda$; distributional curvature
$\frac1{8\pi G}[K]=-\frac32\sigma$). The exterior enters **only** through $-\tfrac M2\Delta t^W_E$.

**Why this matters for FGG.** FGG had to invent a rule for the action of an object that is not a
manifold. T3 shows that the canonical exponent, fixed by T0, is reproduced by covariant data that
never require the exterior region: the interior Euclidean region (which *is* a manifold region in all
models studied), the wall worldvolume, and the exterior Killing time along the wall. **The
pseudomanifold enters the answer only through one boundary term.** Any proposed Euclidean prescription
must reproduce $-\tfrac M2\Delta t^W_E+\mathcal A/4G$ from its exterior to agree with the canonical theory.

**Integral representation.** Integrating T2 in $M$ from the manifold instanton at $M=0$:
$$
B(M)=S_E[\text{compact instanton}]+\frac{\mathcal A(M)}{4G}-\int_0^M\Delta t^W_E(M')\,dM' ,
$$
a pseudomanifold-free definition of the FGG exponent in terms of a genuine instanton and the
Euclidean wall trajectory.

### T4. The Euclidean geometry: obstruction and winding

* **NG3** (NO_GO_THEOREMS.md): for every regime-A transition the exterior Euclidean region is not embedded.
  This is the pseudomanifold, now a theorem with explicit hypotheses.
* **Winding.** The exterior angle is $\theta=\Delta t_E/4GM=2\pi\,\Delta t_E/\beta$: the wall's Euclidean history lasts
  $|\Delta t^W_E|/\beta$ thermal periods of the black hole. For small $M$, $\beta=8\pi GM\to0$ while $\Delta t^W_E$ stays finite,
  so the winding diverges ($\theta_D=-73$ at $E=20\,g(1)$; $-851$ at $E=100\,g(1)$, `FIGURES/euclidean_wall_exterior.png`).
  By T3 its contribution to $B$, $-\tfrac M2\Delta t^W_E$, **vanishes** as $M\to0$. The winding is a property of
  the exterior's multiple cover, not of the exponent. A fuller scrutiny is in VERIFICATION.md, entry V-E5.
* The interior region is embedded in all models studied, including those with an $\varepsilon_{\rm in}$ flip.

### T5. Laboratory transitions are not easier than vacuum decay — [Proposition-N]

Over regime A (the transitions that start from laboratory data), the ratio $B(M)/B(M\to0)$ ranges over
(`DATA/laboratory_regime_bounds.json`, `FIGURES/laboratory_regime_bounds.png`):

| $\gamma$ | $-0.95$ | $-0.65$ | $-0.35$ | $-0.05$ | $0.25$ | $0.55$ | $0.85$ |
|---|---|---|---|---|---|---|---|
| $\min_A B/B(0)$ | 0.564 | 0.572 | 0.613 | 0.766 | 1.000 | 1.001 | 1.003 |
| $\max_A B/B(0)$ | 1.000 | 1.000 | 1.002 | 1.128 | 1.68 | 2.89 | 8.10 |

For $\gamma\gtrsim0.2$ ($\chi\gtrsim1.2\kappa$, the false vacuum larger than a hemisphere at $R_2$),
$B\ge B(0)=B_{\rm CdL}$: **creating a daughter universe from laboratory data is never exponentially easier than
the CdL decay of the same false vacuum**. For all $\gamma$, $B\ge0.56\,B(0)$. The minimum over regime A sits at the
A/B boundary, where the initial wall touches the horizon ($R_1=2GM$, Proposition G1(iii)).

### T6. Physical scales

Using $B(\chi,\kappa,M,G)=b(\kappa/\chi,GM\chi)/G\chi^2$ (`DATA/physical_scales.json`; Planck units, $\rho=E_{\rm vac}^4$, $\sigma=E_{\rm wall}^3$):

| $E_{\rm vac}$ | $E_{\rm wall}$ | $\kappa/\chi$ | $M_{\rm cr}$ | $B(M\to0)$ | $B(10^{-3}M_{\rm cr})$ | $B(0.5M_{\rm cr})$ | $B_{\rm CdL}$ | $S_{\rm dS}$ |
|---|---|---|---|---|---|---|---|---|
| $10^{16}$ GeV | $10^{16}$ GeV | $3.6\times10^{-3}$ | 5.6 g | 133 | $7.0\times10^4$ | $4.4\times10^6$ | 133 | $8.3\times10^{11}$ |
| $10^{16}$ GeV | $3\times10^{16}$ GeV | 0.096 | 5.6 g | $7.0\times10^7$ | $1.1\times10^8$ | $3.2\times10^9$ | $7.0\times10^7$ | $8.3\times10^{11}$ |
| $10^{16}$ GeV | $10^{17}$ GeV | 3.6 | 1.8 g | $1.1\times10^{11}$ | $1.1\times10^{11}$ | $9.0\times10^{10}$ | $7.1\times10^{11}$ | $8.3\times10^{11}$ |
| $10^{15}$ GeV | $10^{16}$ GeV | 0.36 | 0.54 kg | $1.1\times10^{14}$ | $1.1\times10^{14}$ | $4.4\times10^{14}$ | $1.1\times10^{14}$ | $8.3\times10^{15}$ |
| $10^{18}$ GeV | $10^{18}$ GeV | 0.36 | 0.54 mg | 105 | 107 | 438 | 105 | $8.3\times10^{3}$ |

For small $\kappa/\chi$ ($\gamma\to1$) the exponent rises steeply with $M$. In the first row $B/B(0)=1.01$ at $M=10^{-9}M_{\rm cr}$,
6 at $10^{-6}$, $5\times10^2$ at $10^{-3}$ and $3\times10^4$ at $0.5$. The value $B(M\to0)=B_{\rm CdL}$ is reached only for picogram
masses; for milligram-to-gram laboratory masses $B\sim10^5$–$10^7$. At $M=0.5M_{\rm cr}$ the initial wall radius is
$R_1\approx6.6\times10^{-30}$ m, about $1.6$ Schwarzschild radii of the bubble. The last row is near the Planck scale, where
the semiclassical treatment is not trustworthy.

**Phase V correction.** "Picogram masses" should read nanogram: $10^{-9}M_{\rm cr}=5.6\times10^{-12}$ kg. More importantly, these masses are
**below the Planck mass** ($2.2\times10^{-8}$ kg). The tunnelling barrier is still semiclassical ($R_{1,2}\sim10^2$–$10^3\,\ell_P$, WKB parameter about $10^{-11}$).
But the black hole the parent would form has $2GM\ll\ell_P$ and is not semiclassical. For semiclassical black holes ($M\gtrsim10\,m_{\rm Pl}$) the GUT
exponent is $B\gtrsim7\times10^4$ (PHASE_V_QUANTUM_SYSTEM.md Sec. 6).

The exponent depends steeply on $\sigma^4/\rho^3$ (flat-space CdL scaling) and is **not** generically of order
$(M_{\rm Pl}/H)^2$. That scaling holds only for $\kappa\sim\chi$. The thin-wall approximation needs the wall
thickness $\ll R_1$; the GUT-scale numbers are indicative only.

## 5. Relation to the literature

* **FGG (1990)**: canonical attempt flawed by slicing; Euclidean interpolation a "pseudomanifold" with a
  conjectured action rule. Here the slicing jump is identified as a horizon monodromy (T0), the pseudomanifold
  is proved unavoidable (NG3), and its action is fixed by the Smarr identity (T3) without any rule for the exterior.
* **FMP (1990)**: Hamiltonian (Dirac + WKB) treatment; consistency with the Euclidean method requires "metrics
  with degeneracies". Our reduced canonical exponent is of the FMP type. We could not consult the full text in this
  environment (arXiv is blocked by the network policy); the identification of our analytic branch with FMP's choice
  rests on secondary sources stating FGG–FMP agreement.
* **de Alwis–Muia–Pasquarella–Quevedo (2019/2020)** and **Céspedes–de Alwis–Muia–Quevedo (2023)**: Hamiltonian
  transitions between vacua; Minkowski→dS allowed; transitions understood as $M\to0$ limits of black-hole
  transitions; detailed balance with the Hartle–Hawking sign. T1 is consistent with the $M\to0$ picture. A
  term-by-term comparison of their formulae with T1 is pending (P-L1).
* **Chen–Hu–Yeom (2016)**: "something-to-something" vs "nothing-to-something" readings of thin-shell instantons;
  consistent rates. T1/T3 give an explicit covariant decomposition of the something-to-something rate.
* **Ulvestad–Albrecht (2012)**: thick-wall FGG numerics, mass lower bounds.
* **Wang–Wu–Li–Wang (2023/2025)**: wormhole-assisted FGG instantons.

Novelty assessment (to be confirmed against full texts, P-L1): T0's monodromy identification, T1's Morse form,
T2–T3 (first law and Smarr identity for the FGG exponent), NG3 as a theorem, and T5 are, to our knowledge,
not stated in this form in the literature we could access.

## 6. What is still ambiguous

1. **Choice of reduced theory.** The exponent is that of the quantum theory obtained by solving the constraints
   classically with the asymptotic clock. Quantizing first (Dirac/WDW) could differ beyond WKB.
2. **Prefactor and interpretation.** $\Gamma\propto e^{-B/\hbar}$ per bounce of the bounded branch; neither the prefactor
   nor the meaning of "probability to create a causally disconnected region" (unitarity, Freivogel et al.) is settled.
3. **Euclidean side.** T3 fixes what the exterior must contribute but does not derive it from a path integral over
   degenerate or complex geometries. That derivation is the Phase III target.
4. **Thin wall.** T1–T3 are thin-wall statements; thick walls converge in the cases tested, but a thick-wall Smarr
   identity is open.
