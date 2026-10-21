# Global geometry: branch classification and causal structure

Status: **analytic classification + two independent numerical constructions that agree**
(chart-free itinerary `CODE/babyuniverse/branches.py`; explicit Kruskal / global de Sitter
integration `CODE/babyuniverse/geometry.py`; agreement asserted in `tests/test_branches.py`;
atlas in `DATA/branch_atlas_cases.json`, figures `FIGURES/penrose_*.png`,
`FIGURES/phase_diagram_schwarzschild.png`).

The brief warns that an effective-potential curve is not a spacetime. Accordingly no diagram here is
drawn from $V(R)$. Each one is built from an integrated trajectory in global null coordinates on each
side. A second, coordinate-free algorithm must reproduce the same sequence of causal regions.

---

## 1. Branches of $\dot R^2+V(R)=0$

For $M>0$, $V\to-\infty$ as $R\to0$ (term $-G^2M^2/\kappa^2R^4$) and as $R\to\infty$
(term $-L^4R^2/4\kappa^2$), and $V$ has one maximum. Hence there are 0, 1 (double) or 2 turning points.

| case | branches |
|---|---|
| $M<M_{\rm cr}$ (two roots $R_1<R_2$) | **bounded**: $0\to R_1\to0$; **bounce**: $\infty\to R_2\to\infty$ |
| $M=M_{\rm cr}$ | unstable static wall at $R_*=z_*R_{\rm unit}$, plus solutions asymptotic to it |
| $M>M_{\rm cr}$ | **monotonic**: $0\to\infty$ and its time reverse |
| $M=0$ | single root $R_0=2\kappa/L^2$: bounce $\infty\to R_0\to\infty$ (Lorentzian CdL wall) |

The forbidden interval $(R_1,R_2)$ is the tunnelling barrier.

## 2. Orientation at the turning points (Schwarzschild exterior, $\eta=0$)

In $z$ variables $\varepsilon_{\rm out}={\rm sgn}(1-z^3)$, and the minimum of $g$ is at $z_*^3=\tfrac12(-\gamma+\sqrt{\gamma^2+8})\ge1$.

> **Proposition G1.** (i) $z_2>z_*\ge1$, so **every bounce branch has $\varepsilon_{\rm out}(R_2)=-1$**:
> the wall of every inflating solution bounces in region III, beyond the Einstein–Rosen bridge.
> (ii) $\varepsilon_{\rm out}(R_1)=+1$ iff $E>g(1)=2-2\gamma$ (regime **A**), and $-1$ iff $g_{\min}<E<g(1)$ (regime **B**).
> (iii) At the A/B boundary $E=g(1)$ the bounded wall's maximal sphere sits exactly on the
> bifurcation sphere, $R_1=2GM$.
> (iv) In regime A the $\varepsilon_{\rm out}$-flip radius $R_f$ ($z=1$) lies inside the barrier, $R_1<R_f<R_2$.
>
> *Proof.* $g$ decreases on $(0,z_*)$ and increases on $(z_*,\infty)$; $z_2>z_*\ge1$ gives (i).
> $z_1<1\iff g(z_1)=E>g(1)$ gives (ii) and (iv). At a turning point $f_{\rm out}=A_{\rm out}^2$, and
> $f_{\rm out}(z=1)=1-2(1-\gamma)/E$ vanishes iff $E=2-2\gamma$, which gives (iii). $\square$

Regimes in the $(\gamma,\log E)$ plane: `FIGURES/phase_diagram_schwarzschild.png`, left. Regime B
is a thin sliver between $g_{\min}$ and $g(1)$; in mass, $M_{AB}/M_{\rm cr}\approx0.85$–$0.99$
(`DATA/laboratory_regime_bounds.json`). The right panel shows the interior orientation at $R_2$:
for $\gamma>0$ ($\chi>\kappa$, the physically common case of a tension small compared with the vacuum
energy scale) and moderate $E$, $\varepsilon_{\rm in}$ flips inside the barrier and the inflating
bubble at $R_2$ contains the de Sitter equator.

## 3. Chart-free causal itinerary

> **Proposition G2 (itinerary).** On either side, the causal region occupied by the wall at proper
> time $\tau$ is determined by the triple $({\rm sgn}f(R),\ \varepsilon(R),\ {\rm sgn}\dot R)$:
> * $f>0$: a static region; $\varepsilon=+1$ is the static region in which the normal points to increasing
>   $r$ (region I for Schwarzschild, "N" patch for the de Sitter interior), $\varepsilon=-1$ the other one
>   (III, "S");
> * $f<0$: $r$ is a time function; future-directed motion has ${\rm sgn}\dot R=-1$ in the black-hole (II)
>   or contracting de Sitter (P) regions and $+1$ in the white-hole (IV) or expanding (F) regions.
>
> Between events the wall stays in one region, and by Lemma J2 $\varepsilon$ changes only inside dynamical
> regions. The sequence of regions therefore follows from the ordered list of horizon radii and
> flip radii crossed along each monotonic leg.
>
> *Proof sketch.* In a static region $\varepsilon\sqrt{f+\dot R^2}=f\dot t=-\xi\cdot u$ is the Killing energy
> relative to the Killing field $\xi$ that is future-directed in the "right" static region. A
> future-directed $u$ has positive Killing energy where $\xi$ is future-directed and negative where it is
> past-directed, so $\varepsilon$ labels the static region. In a dynamical region
> $g^{ab}\nabla_ar\nabla_br=f<0$, so $\nabla r$ is timelike and its time orientation is ${\rm sgn}\dot R$ on any
> future-directed curve. $\square$

This algorithm (`branches._itinerary`) needs no coordinates and works for any $f$ with any number of
horizons (Schwarzschild–de Sitter included). It is cross-checked by integrating the wall in

* **Kruskal coordinates** for Schwarzschild, $UV=(1-r/2GM)e^{r/2GM}$, using the two algebraically
  equivalent forms $\dot U=-U/(4GMs_+)=s_-/(2kV)$, $\dot V=V/(4GMs_-)=-s_+/(2kU)$ with
  $s_\pm=A\pm\dot R$ ($s_+s_-=f$) and $k=(4G^2M^2/r)e^{-r/2GM}$, switching to the better-conditioned form
  so that horizon crossings are integrated without loss of precision;
* **global conformal coordinates** for de Sitter, $r=\sin\psi/(H\cos T)$,
  $\dot p=H\cos p/s_-=Hs_+\cos^2T/\cos q$, $\dot q=H\cos q/s_+=Hs_-\cos^2T/\cos p$, with $p,q=T\pm\psi$.

The integrated trajectories are tested independently: finite differences of the integrated null
coordinates reproduce $u\cdot u=-1$ and the *unsquared* junction value $n^r=A(R)$ to $<10^{-3}$
median error (`tests/test_geometry.py`). **In all 15 atlas branches the chart-free itinerary and the
explicit integration give identical region sequences.**

Penrose compactification: $\tilde U=\arctan U$, $\tilde V=\arctan V$. With this choice the singularities
$UV=1$ are exactly the straight lines $\tilde T=\pm\pi/4$, because $\arctan U+\arctan(1/U)=\pi/2$. de Sitter is
drawn as the square $(\psi,T)\in[0,\pi]\times(-\pi/2,\pi/2)$. Interior charts show the false vacuum to
the **left** of the wall, exterior charts the true vacuum to the **right**.

## 4. The branch atlas (baseline, $\eta=0$)

Exterior itinerary / interior itinerary (from `DATA/branch_atlas_cases.json`; the two test models
are $\gamma=+0.6$ ($\chi=2\kappa$) and $\gamma=-0.5$ ($\chi=\kappa/\sqrt3$)):

| regime | branch | exterior | interior | inflates? | seen from region I? | past |
|---|---|---|---|---|---|---|
| A | bounded | IV → I → II | N | no | **yes** | white-hole singularity |
| A | bounce | III | P → S → F ($\gamma=0.6$), P → N → F ($\gamma=-0.5$) | **yes** | no | eternal; exterior contains the WH singularity |
| B | bounded | IV → III → II | N | no | no | white-hole singularity |
| B | bounce | III | P → N → F | **yes** | no | as above |
| C | monotonic (expanding) | IV → III | N → F | **yes** | no | white-hole singularity |
| C | monotonic (collapsing) | III → II | P → N | no | no | — |

Figures: `FIGURES/penrose_gamma{+0.6,-0.5}_{A,B,C}_{branch}.png`. Each shows $V(R)$ (with horizons
and flip radii), the interior de Sitter diagram and the exterior Kruskal diagram, with the wall and the
retained side shaded.

**What the atlas says.**

1. The only branch visible from the parent's asymptotic region I is the **bounded** branch of regime A.
   It is the classical evolution of a "laboratory" false-vacuum bubble: it expands to $R_1$ and
   collapses into a black hole. Its maximal development emerges from a white-hole singularity, but
   the time-symmetric slice at $R_1$ is regular laboratory data (CLASSICAL_MODEL.md Sec. 5). Its
   future is regular until the black-hole singularity.
2. **Every inflating branch** ends with $\varepsilon_{\rm out}=-1$ in region III. Its interior reaches the expanding
   region F (anti-trapped spheres), and its exterior contains the white-hole singularity: either the
   wall emerges from it (C) or the whole of region IV lies on the exterior side (bounces). This is
   the Blau–Guendelman–Guth picture, reproduced here independently.
3. No branch is an inflating bubble that remains an ordinary, visible region of the parent. "Inflating"
   and "visible from region I" are mutually exclusive in this model.

## 5. What counts as a daughter universe here

We use an explicit criterion rather than "the interior grows". A branch produces a **daughter universe**
if

* (D1) the false-vacuum region contains complete expanding de Sitter regions (the interior reaches F
  and $R\to\infty$), so its spatial volume grows without bound;
* (D2) no future-directed causal curve from the inflating region reaches the parent's future null
  infinity $\mathscr I^+_{\rm I}$. The wall is in region III; future-directed curves from the interior must
  cross the wall into III, and from III they can only enter II and end on the singularity;
* (D3) the parent remains asymptotically flat with ADM mass $M$: to an outside observer the event is
  the formation of a black hole of mass $M$.

All inflating branches satisfy D1–D3 (flag `daughter_causally_disconnected` in the atlas).

**Topology.** Every Cauchy slice of these spacetimes is topologically $\mathbb R^3$, including the bag-of-gold
slice at $R_2$: ball ∪ ($S^2\times[0,\infty)$) through the bridge. There is **no change of spatial topology**
on regular Cauchy slices. The "separation" of the daughter is causal: after the bridge pinches off at the
region-II singularity, the only way to split late-time slices into two components is through the
singularity. Any topology-change interpretation is therefore singular in this classical model
(cf. Geroch, Tipler, Borde; NO_GO_THEOREMS.md Sec. 5).

**Curvature vs coordinate singularities.** $r=0$ in II/IV is a curvature singularity
($R_{abcd}R^{abcd}=48G^2M^2/r^6$). All horizons are coordinate singularities of the static charts,
resolved by the Kruskal / global charts used here. When the interior bubble shrinks to $R\to0$ at the
de Sitter pole, the false-vacuum region simply disappears into the exterior singularity.

## 6. $M=0$ and the de Sitter parent

* **$\eta=0$, $M=0$:** $\varepsilon_{\rm out}=-1$ for all $R$, so the "exterior" is a flat **ball**. The solution is a closed
  universe (de Sitter cap + flat ball, the Lorentzian continuation of the compact $O(4)$ instanton),
  not attached to any asymptotic region. It is the $M\to0^+$ limit of the bounce branch: region III
  becomes the flat ball and the bridge to region I shrinks to zero size.
* **$\eta>0$, $M=0$ (de Sitter parent):** the exterior itinerary is dS-past → dS-static(−) → dS-future.
  The spacetime is **globally regular**: no horizon singularities, complete. It contains an inflating
  false-vacuum region beyond the cosmological horizon of the observer at the antipode
  (`FIGURES/penrose_M0_dS_exterior_bounce-M0.png`). This is the Lorentzian continuation of the CdL
  (upward / Lee–Weinberg) instanton. The Farhi–Guth theorem does not apply because the Cauchy
  surfaces are compact ($S^3$). This is a classical, singularity-free "baby universe" inside a
  *closed* parent. It is not created from regular data: the wall contracts from $\mathscr I^-$.

## 7. Schwarzschild–de Sitter exterior

The potential, signs, turning points and chart-free itinerary are implemented for general
$(\chi,\eta,\kappa,M)$ (labels BH/WH, SdS-static(±), cosmo-future/past). Explicit conformal diagrams for
SdS need a patched Kruskal atlas around both horizons, which is not yet implemented. The atlas here is
restricted to $\eta=0$ and $M=0$ (OPEN_PROBLEMS.md P-G1).
