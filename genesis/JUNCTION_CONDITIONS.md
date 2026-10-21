# Junction conditions for a spherical thin wall

Status: **derived by hand and verified symbolically** (`CODE/babyuniverse/symbolic.py`,
functions `extrinsic_curvature_static`, `israel_equations`, `squaring_and_signs`,
`dimensionless_form`, `mass_formula`, `wall_equation_of_motion`; log in
`DATA/derivation_log.txt`; tests `CODE/tests/test_symbolic.py`, `test_model.py`).

Conventions (fixed for the whole programme):

* signature $(-,+,+,+)$, $c=1$; $G$ and $\hbar$ explicit where they matter;
* $K_{ab}=e_a^{\mu}e_b^{\nu}\nabla_\mu n_\nu$, with the unit normal $n$ pointing **from the interior (–) to the exterior (+)**;
* $[X]=X_+-X_-=X_{\rm out}-X_{\rm in}$;
* surface stress tensor $S_{ab}$ defined by $\delta S_{\rm wall}=-\tfrac12\int_W\sqrt{-h}\,S_{ab}\,\delta h^{ab}$.

These are the conventions of Israel (1966) and of Poisson, *A Relativist's Toolkit*, Sec. 3.7.

---

## 1. Geometry on each side

On each side the metric is static and spherically symmetric,
$$
ds^2_\pm=-f_\pm(r)\,dt_\pm^2+\frac{dr^2}{f_\pm(r)}+r^2d\Omega^2 ,
$$
with, in the baseline model,
$$
f_{\rm in}(r)=1-\chi^2r^2,\qquad f_{\rm out}(r)=1-\frac{2GM}{r}-\eta^2r^2 ,
$$
$\chi^2=\Lambda_{\rm in}/3=8\pi G\rho_{\rm in}/3$, $\eta^2=\Lambda_{\rm out}/3=8\pi G\rho_{\rm out}/3$.
By Birkhoff's theorem (with cosmological constant) these are the most general
spherically symmetric vacuum solutions on each side, so the only freedom is
$(\chi,\eta,M)$ and the wall trajectory.

The wall worldvolume $W$ is $r=R(\tau)$, $t_\pm=T_\pm(\tau)$, with $\tau$ its proper time;
its induced metric is the same from both sides (first junction condition),
$$
ds^2_W=-d\tau^2+R(\tau)^2d\Omega^2 ,
$$
which is guaranteed because both sides use the area radius $r$ and proper time:
$-f_\pm\dot T_\pm^2+\dot R^2/f_\pm=-1$.

## 2. Extrinsic curvature

Tangent $u=\dot T\partial_t+\dot R\partial_r$, unit normal $n_\mu=\varepsilon'(-\dot R,\dot T,0,0)$
($n\cdot n=1$ follows from the normalisation of $u$). Define the **signed root**
$$
A\equiv n^r=\varepsilon\sqrt{f+\dot R^2},\qquad \varepsilon={\rm sgn}\,n^r .
$$

*Angular components.* $K_{\theta\theta}=-n_\mu\Gamma^\mu_{\theta\theta}=-n_r\Gamma^r_{\theta\theta}=rfn_r=r\,n^r$:
$$
K_{\theta\theta}=R\,A,\qquad K^\theta{}_\theta=K^\phi{}_\phi=\frac{A}{R}.
$$

*Time component* (a coordinate-free derivation). Let $\xi=\partial_t$ be the Killing field and
$e=-\xi\cdot u=f\dot T$. Since $\nabla_{(\mu}\xi_{\nu)}=0$, $\dot e=-\xi_\mu a^\mu$ with
$a=u\cdot\nabla u$. For radial motion $a=\alpha\, n$ with $\alpha=n\cdot a$, and
$\xi\cdot n=-\varepsilon'\dot R$, $e=\varepsilon' A$. Hence $\dot A=\alpha\dot R$. Differentiating
$A^2=f+\dot R^2$ gives $\dot A=\dot R(\ddot R+f'/2)/A$, so
$$
K_{\tau\tau}=-n\cdot a=-\frac{\ddot R+\tfrac12f'(R)}{A},\qquad
K^\tau{}_\tau=\frac{\ddot R+\tfrac12 f'(R)}{A}.
$$

Both results are reproduced from the Christoffel symbols for a *generic* $f$ by
`symbolic.extrinsic_curvature_static()`.

## 3. Israel conditions

For a domain wall $S_{ab}=-\sigma h_{ab}$ (surface energy density $\sigma$, surface pressure $-\sigma$).
Then $S=-3\sigma$, and the trace-reversed Israel condition
$[K_{ab}]=-8\pi G\,(S_{ab}-\tfrac12h_{ab}S)$ gives
$$
[K_{ab}]=-4\pi G\sigma\,h_{ab}\equiv-\kappa\, h_{ab},\qquad \kappa\equiv4\pi G\sigma .
$$

**$\theta\theta$ component (first-order equation of motion):**
$$
\boxed{\;\varepsilon_{\rm in}\sqrt{f_{\rm in}+\dot R^2}-\varepsilon_{\rm out}\sqrt{f_{\rm out}+\dot R^2}=\kappa R\;}
\qquad\text{i.e.}\qquad A_{\rm in}-A_{\rm out}=\kappa R .
$$
This is the equation quoted in the programme brief, now *derived*, with the sign of $\kappa R$
fixed by the conventions above (check: a positive-energy shell decreases $K_{\theta\theta}$
outward, as it must; for dust $S_{ab}=\mu u_au_b$ the same formula gives $[K_{\theta\theta}]=-4\pi G\mu R^2<0$).

**$\tau\tau$ component (second-order equation):**
$$
\frac{\ddot R+\tfrac12f_{\rm out}'}{A_{\rm out}}-\frac{\ddot R+\tfrac12f_{\rm in}'}{A_{\rm in}}=-\kappa .
$$
It is not independent: $\frac{d}{d\tau}(A_{\rm out}-A_{\rm in}+\kappa R)=\dot R\times(\tau\tau\ \text{equation})$
identically (verified symbolically). This is the conservation law $D_aS^{ab}=-[T^{b}{}_n]$
with constant $\sigma$: vacuum energy on both sides gives no flux $T_{\tau n}=0$.

**Normal equation of motion.** The remaining component of the wall's equations of motion,
$S^{ab}\{K_{ab}\}=[T_{nn}]$ with $\{X\}=\tfrac12(X_++X_-)$, reads for a domain wall in vacuum energies
$$
\sigma\,\{K\}=\rho_{\rm out}-\rho_{\rm in},\qquad K=K^\tau{}_\tau+2K^\theta{}_\theta ,
$$
and is implied by the $\theta\theta$ equation together with $\ddot R=-V'(R)/2$
(`symbolic.wall_equation_of_motion`, exact arithmetic).

## 4. Exact resolution of the signed roots (no squaring ambiguity)

Write $X=A_{\rm in}$, $Y=A_{\rm out}$. The $\theta\theta$ equation is $X-Y=\kappa R$, and
*identically* $X^2-Y^2=(f_{\rm in}+\dot R^2)-(f_{\rm out}+\dot R^2)=f_{\rm in}-f_{\rm out}$.
Therefore $X+Y=(f_{\rm in}-f_{\rm out})/\kappa R$ and
$$
\boxed{\;A_{\rm in}(R)=\frac{f_{\rm in}-f_{\rm out}+\kappa^2R^2}{2\kappa R},\qquad
A_{\rm out}(R)=\frac{f_{\rm in}-f_{\rm out}-\kappa^2R^2}{2\kappa R}\;}
$$
are rational functions of $R$ alone on every solution, and
$$
\dot R^2+V(R)=0,\qquad V(R)=f_{\rm out}-A_{\rm out}^2=f_{\rm in}-A_{\rm in}^2 .
$$

> **Lemma J1 (equivalence).** Let $\sigma>0$. A curve $R(\tau)$ with signs
> $(\varepsilon_{\rm in},\varepsilon_{\rm out})$ solves the unsquared $\theta\theta$ equation **iff**
> $\dot R^2+V(R)=0$ and $\varepsilon_{\rm in}={\rm sgn}\,A_{\rm in}(R)$, $\varepsilon_{\rm out}={\rm sgn}\,A_{\rm out}(R)$
> wherever these are non-zero.
>
> *Proof.* ($\Rightarrow$) shown above. ($\Leftarrow$) If $\dot R^2=A_{\rm out}(R)^2-f_{\rm out}$ then
> $\sqrt{f_{\rm out}+\dot R^2}=|A_{\rm out}(R)|$; with $\varepsilon_{\rm out}={\rm sgn}A_{\rm out}$ this is
> $A_{\rm out}(R)$, and likewise inside; $A_{\rm in}-A_{\rm out}=\kappa R$ holds identically. $\square$

Consequences, which answer the brief's warning about squaring:

1. Squaring introduces **no spurious solution** provided the signs are *read off* from $A_{\rm in},A_{\rm out}$
   rather than *imposed*. The test `test_unsquared_junction_holds_on_shell` checks the unsquared
   equation to $10^{-9}$ on every allowed point of five models and that the opposite sign fails.
2. The squared potential is invariant under $\sigma\to-\sigma$ (verified); the sign-reversed
   assignment $(\varepsilon_{\rm in},\varepsilon_{\rm out})\to(-\varepsilon_{\rm in},-\varepsilon_{\rm out})$ is the *negative-tension*
   solution. These are the only "spurious" solutions of the squared equation.
3. Branches that would be *discarded* by imposing, e.g., $\varepsilon_{\rm out}=+1$ (a bubble "in our universe")
   are not spurious: they are the global geometries with the wall beyond the Einstein–Rosen bridge,
   which is exactly where the inflating solutions live (GLOBAL_GEOMETRY.md).

> **Lemma J2 (where orientation can change).** Along a Lorentzian solution $\varepsilon_{\rm out}$ can change
> sign only at $R_f$ with $A_{\rm out}(R_f)=0$, and there $f_{\rm out}(R_f)=-\dot R^2\le0$: the flip happens inside
> a dynamical (trapped or anti-trapped) region, or at a turning point lying exactly on a horizon.
> Along a Euclidean solution ($\dot R\to iR'$, $R'^2=V>0$) the flip happens where $f_{\rm out}=V>0$,
> i.e. in a static region. The same holds for $\varepsilon_{\rm in}$.

For the baseline functions
$$
A_{\rm out}=0\iff R_f^3=\frac{2GM}{\chi^2+\kappa^2-\eta^2},\qquad
A_{\rm in}=0\iff R_{f,\rm in}^3=\frac{2GM}{\chi^2-\eta^2-\kappa^2},
$$
each existing only if the denominator is positive ($M>0$).

**Meaning of the signs.** On a spatial slice through the wall:

| sign | geometry on that side |
|---|---|
| $\varepsilon_{\rm in}=+1$ | area radius increases from the regular centre to the wall (interior less than a hemisphere of the de Sitter 3-sphere when $\chi>0$) |
| $\varepsilon_{\rm in}=-1$ | the interior contains the de Sitter equator $r=1/\chi$ (more than a hemisphere) |
| $\varepsilon_{\rm out}=+1$ | area radius increases away from the wall: in static regions the wall is in region I of Kruskal |
| $\varepsilon_{\rm out}=-1$ | area radius decreases away from the wall: wall in region III; the exterior contains the Einstein–Rosen bridge |

## 5. The effective potential

Substituting the baseline functions,
$$
V(R)=1-\frac{L^4}{4\kappa^2}R^2+\frac{GM(\chi^2-\eta^2-\kappa^2)}{\kappa^2R}-\frac{G^2M^2}{\kappa^2R^4},
\qquad L^4\equiv(\chi^2-\eta^2+\kappa^2)^2+4\kappa^2\eta^2 .
$$
With $R=(2GM/L^2)^{1/3}z$ the potential collapses to a one-parameter shape (valid for all $\chi,\eta,\kappa$, $M>0$):
$$
V=1-\frac{g(z)}{E},\qquad g(z)=z^2-\frac{2\gamma}{z}+\frac1{z^4},\qquad
\gamma=\frac{\chi^2-\eta^2-\kappa^2}{L^2},\qquad E=\frac{4\kappa^2}{L^{8/3}(2GM)^{2/3}},
$$
and the signs become
$$
\varepsilon_{\rm out}={\rm sgn}(1-\beta_{\rm out}z^3),\quad \beta_{\rm out}=\frac{\chi^2+\kappa^2-\eta^2}{L^2};\qquad
\varepsilon_{\rm in}={\rm sgn}(1-\gamma z^3).
$$
Identity: $L^4-(\chi^2-\eta^2-\kappa^2)^2=4\kappa^2\chi^2\ge0$, so $|\gamma|\le1$.
$g$ has a single minimum at $z_*^3=\tfrac12(-\gamma+\sqrt{\gamma^2+8})\in[1,2]$, which gives the
critical mass at which the barrier disappears,
$$
M_{\rm cr}=\frac{4\kappa^3}{G\,g(z_*)^{3/2}L^4}.
$$
For $\eta=0$: $\beta_{\rm out}=1$, $\gamma=(\chi^2-\kappa^2)/(\chi^2+\kappa^2)$,
$g(1)=2-2\gamma$, $f_{\rm out}=1-2(1-\gamma)/(Ez)$, $f_{\rm in}=1-(1-\gamma^2)z^2/E$. The
classification that follows from these formulas is in GLOBAL_GEOMETRY.md.

For $M=0$: $V=1-R^2/R_0^2$ with $R_0=2\kappa/L^2$ (the thin-wall Coleman–De Luccia/Parke radius).

## 6. The mass formula

Solving the junction condition for $M$:
$$
M=\frac{4\pi}{3}(\rho_{\rm in}-\rho_{\rm out})R^3+\varepsilon_{\rm in}\,4\pi\sigma R^2\sqrt{1-\chi^2R^2+\dot R^2}-8\pi^2G\sigma^2R^3 ,
$$
volume energy + (redshifted, boosted) wall energy − wall self-binding. For $G\to0$ and $\varepsilon_{\rm in}=+1$
it is the flat-space bubble energy $\tfrac{4\pi}3\rho R^3+4\pi\sigma R^2\gamma_{\rm Lorentz}$. The sign
$\varepsilon_{\rm in}=-1$ makes the wall term negative; this is how large false-vacuum regions can have
small exterior mass. Verified by substitution (`symbolic.mass_formula`) and numerically
(`test_mass_function_inverts_junction`).

## 7. Euclidean continuation

Under $\tau=-i\tau_E$, $t=-it_E$: $\dot R=iR'$, $R'^2=V(R)$, and
$f\,dt_E/d\tau_E=A$ with the **same** rational $A(R)$. Writing $(A,\sqrt V)=\sqrt f(\cos\Theta,\sin\Theta)$,
$$
\cos(\Theta_{\rm out}-\Theta_{\rm in})=\frac{f_{\rm in}+f_{\rm out}-\kappa^2R^2}{2\sqrt{f_{\rm in}f_{\rm out}}}
$$
(`symbolic.euclidean_angle_identity`). This angle relation is the classical symbol of the
quantum constraint of QUANTIZATION.md, and the angle difference is the Euclidean canonical
momentum of TUNNELLING.md.
