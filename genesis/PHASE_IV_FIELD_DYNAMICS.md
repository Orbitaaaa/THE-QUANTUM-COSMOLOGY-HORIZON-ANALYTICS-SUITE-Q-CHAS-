# Phase IV — Replacing the idealised wall with actual fields

**Exit condition (from the programme):** *a self-consistent field-theoretic model in which the relevant parent and
daughter configurations can be characterised.* Status: **met, in spherical symmetry.**

* A physically specified scalar potential (quartic with a tilt, false vacuum $U=\rho_f$ and true vacuum $U=0$) is coupled to
  gravity. The full Einstein–scalar equations are solved in double-null coordinates, through horizon formation and up to the
  singularity.
* **Parent:** laboratory false-vacuum bubbles collapse to black holes. As the wall thins, the field-theory evolution converges
  to the Israel thin-wall dynamics: $R_{\max}\to R_1$, $r_{AH}\to2GM$, $A_{\rm in}-A_{\rm out}\to\kappa R$. No anti-trapped region ever forms,
  in 10/10 laboratory configurations, and Proposition IV-1 proves this in general for spherical NEC matter.
* **Daughters:** both inflating thin-wall branches exist as solutions of the field equations and converge to the thin-wall
  trajectories.
  * The **bounce branch** has $\varepsilon_{\rm out}=-1$ throughout and lies beyond the Einstein–Rosen bridge. It is reached from data that are
    regular at the centre but contain a "reversed normal" exterior ($r_v<0<r_u$).
  * The **white-hole branch** (regime C) is reached from anti-trapped data.
  * Both develop an anti-trapped (inflating) region. Corollary IV-2: in spherical symmetry such a region is always in the future
    of initial data with $r_u>0$, i.e. data that already lie beyond a bridge, in a white hole or beyond a cosmological horizon.
* **Thin → thick:** thin-wall predictions hold to about 1% for $w/R\lesssim0.05$. At $w/R=0.4$ the thin-wall model fails ($r_{AH}=0.79\cdot2GM$) but
  the qualitative outcome, a black hole and no anti-trapped region, persists.
* **Stability:** the outcomes are robust to radial perturbations of the data. A probe-wall (Jacobi) analysis gives an aspherical
  amplification of at most about 2 before the collapsing wall reaches $2GM$.

Code: `CODE/babyuniverse/doublenull.py`, `CODE/babyuniverse/stability.py`; script `CODE/scripts/phase4_dynamics.py`
(about 3 min); data `DATA/phase4_dynamics.json`; figures `FIGURES/phase4_bubbles.png`, `FIGURES/phase4_causal_regions.png`;
tests `CODE/tests/test_doublenull.py` (8).

---

## 1. Model

Matter: $S_m=\int\sqrt{-g}\,[-\tfrac12(\nabla\phi)^2-U(\phi)]$ with the potential of CLASSICAL_MODEL.md Sec. 7,
$$U(\phi)=\tfrac\lambda4(\phi^2-v^2)^2+\rho_f\big(1-s(\phi/v)\big),\qquad s(x)=\tfrac12+\tfrac34x-\tfrac14x^3 .$$
The false vacuum is $\phi=-v$ ($U=\rho_f$) and the true vacuum $\phi=+v$ ($U=0$, Minkowski/Schwarzschild). The kink width is $w=\sqrt2/(\sqrt\lambda v)$ and the
tension $\sigma_0=\tfrac43v^2/w$. Units are $G=1$, $\chi=\sqrt{8\pi G\rho_f/3}=1$, and $\kappa=4\pi G\sigma_0$ is a parameter (mostly $\kappa=1$).

Metric $ds^2=-\alpha^2du\,dv+r^2d\Omega^2$. Field equations (verified symbolically, `doublenull.symbolic_check`: each is a
non-vanishing multiple of a component of $G_{\mu\nu}-8\pi GT_{\mu\nu}$ or of the wave equation):
$$
\begin{aligned}
r_{,uv}&=-\frac{r_{,u}r_{,v}}r-\frac{\alpha^2}{4r}+2\pi G\alpha^2rU, &
(\ln\alpha)_{,uv}&=\frac{r_{,u}r_{,v}}{r^2}+\frac{\alpha^2}{4r^2}-4\pi G\phi_{,u}\phi_{,v},\\
(r\phi)_{,uv}&=\phi\,r_{,uv}-\tfrac{\alpha^2}4rU'(\phi), &
r_{,uu}-2\frac{\alpha_{,u}}\alpha r_{,u}&=-4\pi Gr\phi_{,u}^2,\quad r_{,vv}-2\frac{\alpha_{,v}}\alpha r_{,v}=-4\pi Gr\phi_{,v}^2 .
\end{aligned}
$$
Misner–Sharp mass $m=\frac r{2G}(1+4r_{,u}r_{,v}/\alpha^2)$. Causal character of a sphere: normal $r_u<0<r_v$, "reversed normal"
$r_v<0<r_u$ (region III, beyond the bridge), trapped $r_u,r_v<0$, anti-trapped $r_u,r_v>0$.

## 2. Numerical method (NUMERICAL_METHODS.md, new section)

* Second-order diamond (characteristic) scheme with a three-pass predictor–corrector at the cell centre.
* The scalar is evolved as $\psi=r\phi$. The $\phi$-form, with the term $-(r_u\phi_v+r_v\phi_u)/r$, was **unstable at the centre**: a massless pulse
  blew up to $|\phi|\sim10^4$ (FAILED_APPROACHES F11). The $\psi$-form is exact in flat space.
* Centre: $r=\psi=0$ and $\alpha=2r_v$ (regularity). Within two points of the centre, $r_ur_v+\alpha^2/4=(\alpha^2/4)(2Gm/r)$ is replaced by its
  regular limit $m=\tfrac{4\pi}3\rho r^3$, which is exact in de Sitter. Without this the $1/r^2$ term magnified errors by about 8 per ray (F11).
* Two characteristic problems:
  * (a) a regular centre with data on the light cone of the centre, in gauge $r=v/2$ or $\alpha=1$;
  * (b) data on two null rays with no centre, for anti-trapped data.
* Excision:
  * $r<2h$ in a trapped region (the singularity);
  * $r>r_\infty$ in an anti-trapped region (beyond future null infinity of the inflating region);
  * $r$ increasing along an already trapped outgoing ray. This is forbidden by Raychaudhuri under the NEC; it occurred only through
    numerical breakdown next to the singularity at under-resolved settings (F12).

**Validation.**

| test | result |
|---|---|
| flat-space massless pulse ($G\to0$), exact solution $r\phi=F(v)-F(u)$ | error $1.9\times10^{-14}$ |
| exact de Sitter in the code's gauge (closed form, `de_sitter_exact`) | $\max\lvert\delta r\rvert$: $5.9,2.3,0.79,0.26\times10^{-4}$ for $h=0.02\dots0.0025$, about second order; $\ln\alpha$ second order away from the centre, first order at it |
| self-gravitating pulse ($A=0.15$), self-convergence away from the centre and trapped region | ratios $r$: 3.94, 3.98; $\phi$: 4.35, 4.49; $\ln\alpha$: 3.6, 3.1 (second order) |
| critical collapse | $A=0.03, 0.04$ disperse; $A=0.05$ forms a black hole ($r_{AH}=0.0630$ vs $0.0627$ at $h=0.002,0.001$); $A=0.15$: $r_{AH}=0.4591=1.98\,M_{\rm initial}$ |
| symbolic derivation | all five equations confirmed by sympy (test) |

## 3. Parent: the laboratory bubble

Data on the light cone of a regular centre: $\phi=v\tanh((r-0.4)/w)$, i.e. false vacuum inside $r\approx0.4$, $\kappa=1$.

| $w$ | $M$ | $R_{\max}$ | thin-wall $R_1(M)$ | dev. | $r_{AH}/2GM$ | Israel $(A_{\rm in}-A_{\rm out})/\kappa R$ | max dev. of $R(\tau)$ | anti-trapped |
|---|---|---|---|---|---|---|---|---|
| 0.04 | 0.1560 | 0.4055 | 0.4139 | −2.0% | 0.974 | 1.056 | — | 0 |
| 0.02 | 0.1546 | 0.4084 | 0.4118 | −0.84% | 0.987 | 0.998 | 0.92% | 0 |
| 0.01 | 0.1539 | 0.4087 | 0.4108 | −0.51% | 0.994 | 0.983 | 0.55% | 0 |
| 0.005 | 0.1536 | 0.4089 | 0.4103 | −0.35% | 0.996 | 0.992 | 0.37% | 0 |

The wall rises slightly to $R_{\max}$, then collapses. The apparent horizon appears at the wall as it crosses $2GM$. The false-vacuum
interior is crushed into the singularity with the wall, and the singularity reaches the centre. All of this is the thin-wall bounded
branch of regime A ($E=3.47>2-2\gamma$). The deviations shrink roughly as $w^{0.6\text{–}0.7}$, slower than linearly. This matches the slow
convergence of the tilt contribution to the tension already found in Stage One (F8).

### Proposition IV-1 (Farhi–Guth, forward in time, in field theory)

> Let the matter satisfy the NEC; for a scalar field $T_{uu}=\phi_u^2\ge0$. In spherical symmetry the $E_{uu}$ constraint reads
> $\alpha^2\,\partial_u(r_u/\alpha^2)=-4\pi G\,r\,\phi_u^2\le0$. So $r_u/\alpha^2$ is non-increasing along every ingoing null ray. If the initial outgoing cone
> has $r_u<0$ everywhere, then $r_u<0$ in its whole future development, and **no anti-trapped sphere ever forms**.
>
> For laboratory data, meaning a regular centre, $r_v>0$ and $r>2Gm(r)$ everywhere on the cone (no trapped, marginally trapped or
> anti-trapped sphere), one has $r_u=-(1-2Gm/r)\alpha^2/4r_v<0$ on the cone. Hence an inflating (anti-trapped) false-vacuum region can never
> develop from laboratory data, for any potential and any wall thickness. $\square$

**Numerical confirmation.** The scan over $\kappa\in\{0.5,1,2\}$ × initial wall radius $r_w\in[0.15,0.85]$ (18 runs) found:
* **10 laboratory configurations:** 0 anti-trapped points, and **10/10 form black holes**. Even the smallest bubble ($M=0.013$, $2GM/r_w=0.17$)
  collapses to a black hole rather than dispersing.
* **8 non-laboratory configurations:** $2GM>r_w$ on the cone, so a white-hole exterior is present from the start. 7 of 8 contain
  anti-trapped points from the outset. The marginal case $2GM/r_w=1.03$ has none.

Robustness: replacing tanh by erf, or adding small pulses inside or outside, changes nothing qualitatively (`stability.parent_robustness`).

### Corollary IV-2 (where daughters can come from)

> Under the same hypotheses, every anti-trapped point lies to the future, along its ingoing null ray, of initial data with $r_u>0$.
> In spherical symmetry with NEC matter there are three possibilities: a "reversed normal" exterior beyond an Einstein–Rosen bridge
> (region III), a white-hole region, or a false vacuum already beyond its own cosmological horizon. This is the field-theoretic
> counterpart of the Stage One statements (GLOBAL_GEOMETRY.md, NO_GO_THEOREMS.md Sec. 4). The only way past it is quantum: Phase III.

## 4. Daughters: the inflating branches as solutions of the field equations

### 4.1 Bounce branch ($\varepsilon_{\rm out}=-1$), from data with a regular centre

In the gauge $\alpha=1$ on the cone, $r$ follows from $r_{vv}=-4\pi Gr\phi_v^2$ and **drops** across the wall. On each side of the wall
$r_v=\lambda(\dot R+A_\pm)$, with $\lambda=w/w_v$ fixed by the null width $w_v$ of the profile, and the thin-wall jump
$\Delta r_v=-\lambda\kappa R$ is produced automatically by the constraint. Choosing $w_v=2w\sqrt{1-\chi^2R_2^2}$ puts the wall at the thin-wall bounce
radius $R_2$ with $\dot R=0$. Outside, $r_v<0<r_u$: the exterior is region III, and the outgoing cone falls into the black hole.

| $w$ | $M$ | thin-wall $R_2(M)$ | $R$: start → end | $\varepsilon_{\rm in},\varepsilon_{\rm out}$ | Israel ratio (median; p10–p90) | $\dot R$ vs thin wall ($\sigma_0$; $\sigma_{\rm thin}$) | anti-trapped |
|---|---|---|---|---|---|---|---|
| 0.02 | 0.1622 | 0.9860 | 0.983 → 3.17 | +, − throughout | 1.027; 1.004–1.034 | +1.8%; −1.2% | yes |
| 0.01 | 0.1583 | 0.9867 | 0.985 → 2.17 | +, − throughout | 1.013; 1.004–1.017 | +1.0%; −0.5% | yes |
| 0.005 | 0.1563 | 0.9870 | 0.986 → 1.91 | +, − throughout | 1.006; 1.004–1.008 | +0.7%; −0.14% | yes |

Agreement with the thin-wall bounce branch improves as $w\to0$ (figure `phase4_bubbles.png`). The causal map (`phase4_causal_regions.png`)
shows the geometry claimed by the thin-wall analysis:
* an anti-trapped (inflating) false-vacuum region forms behind the wall;
* the wall lies beyond the Einstein–Rosen bridge (region III);
* the outer exterior is a black hole.

Robustness: initial $\dot R_0=\pm0.05$ still bounces and inflates with $\varepsilon_{\rm out}<0$ throughout.

### 4.2 White-hole branch (regime C), from anti-trapped data on two null rays

The ingoing ray carries exact anti-trapped de Sitter ($r_0=1.2/\chi$). The outgoing ray crosses the wall into a white-hole exterior. The
result: $M=5.08\gg M_{\rm cr}=0.385$ (regime C). The wall inflates from $R=1.70$ to $R=16.5$, with $\dot R$ matching thin wall to
**0.4%** (median), Israel ratio 0.98, and $\varepsilon_{\rm in}=\varepsilon_{\rm out}=+1$ initially (region IV exterior).

### 4.3 Characterisation of the parent and daughter configurations

| | parent (laboratory) | daughter, bounce branch | daughter, white-hole branch |
|---|---|---|---|
| initial data | regular centre, no trapped or anti-trapped sphere | regular centre; reversed-normal exterior (beyond a bridge) | anti-trapped false vacuum; white-hole exterior |
| wall | $\varepsilon_{\rm out}=+1$, rises to $R_1$, collapses | $\varepsilon_{\rm out}=-1$ throughout, bounces at $R_2$, inflates | $\varepsilon_{\rm out}=+1$ at start, inflates |
| outcome | black hole $M$; false vacuum crushed | black hole $M$ outside; inflating region behind the horizon | inflating region; exterior black hole/white hole of mass $M$ |
| anti-trapped region | never (Prop. IV-1) | forms, fed by $r_u>0$ beyond the bridge (Cor. IV-2) | present from the start |
| thin-wall counterpart | bounded branch, regime A | bounce branch, regime A | monotonic branch, regime C |

## 5. Thin to thick walls

Laboratory bubble with $r_w=0.4$ and increasing $w$:

| $w/R$ | 0.0125 | 0.025 | 0.05 | 0.1 | 0.2 | 0.4 |
|---|---|---|---|---|---|---|
| $R_{\max}/R_1^{\rm thin}$ | 0.996 | 0.994 | 0.990 | 0.977 | 0.956 (no rise) | 0.931 (no rise) |
| $r_{AH}/2GM$ | 0.995 | 0.991 | 0.983 | 0.973 | 0.951 | 0.790 |
| black hole / anti-trapped | yes / 0 | yes / 0 | yes / 0 | yes / 0 | yes / 0 | yes / 0 |

For $w/R\lesssim0.05$ the thin-wall model is accurate to about 1–2%. By $w/R\sim0.2$ the wall no longer rises at all (the initial data are already
past the maximum). At $w/R=0.4$ the matter is spread over the bubble and the first apparent horizon forms well inside $2GM$. The
qualitative statements (collapse, black hole, no anti-trapped region) do not depend on the thin-wall idealisation, as
Proposition IV-1 requires.

## 6. Stability

* **Radial, non-linear:** robust to the profile shape and to added radiation (parent); robust to the initial wall velocity (daughter). See above.
* **Aspherical, probe approximation** (`stability.py`). Normal deformations $\xi_\ell$ of a constant-mean-curvature wall obey the Jacobi equation
  $\ddot\xi+2(\dot R/R)\dot\xi+[\ell(\ell+1)/R^2-K_{ab}K^{ab}-{\rm Ric}(n,n)]\xi=0$. The $\ell=1$ translation mode $\xi=\gamma$ and the $\ell=0$ mode $\xi=\dot R$ are reproduced to
  $5\times10^{-6}$. Along the collapse from $R_{\max}$, the amplification of the shape distortion $\delta_\ell=\xi_\ell/R$ is:
  * to $R=2GM$: $\ell=2$: 2.1; $\ell=3$: 1.3; $\ell=4$: 0.54; $\ell=6$: 1.2; $\ell=10$: 0.83;
  * to $R=GM$: $\ell=2$: 15.7; $\ell=3$: 7.1; $\ell=4$: 3.6; $\ell=6$: 8.2; $\ell=10$: 5.5.

  So small asphericities stay small until the wall is inside its gravitational radius. The spherical treatment of horizon
  formation is self-consistent. Growth inside the horizon is irrelevant to the exterior.
* **Limitation.** The probe approximation neglects the wall's own gravity. A self-gravitating aspherical thin-shell or field analysis, and
  the aspherical stability of the inflating daughter, remain open (P-D2).

## 7. What Phase IV did not do

* No aspherical non-linear evolution (no 2+1 or 3+1 code).
* No evolution of the white-hole past or the past singularity (Corollary IV-2 only locates the required data).
* No quantum fields on these backgrounds (Hawking radiation, backreaction on the transition): Phase V.
* One potential family (quartic with tilt), with $\chi=1$ and $\kappa\in\{0.5,1,2\}$.
