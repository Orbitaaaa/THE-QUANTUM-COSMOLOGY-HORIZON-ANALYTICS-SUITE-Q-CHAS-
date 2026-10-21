# The classical model: action, field equations, admissible data

This file specifies the physical system of Stage One completely: fields, action, boundary
terms, variational principle, field equations, admissible initial data, and the
energy conditions they satisfy. Derived consequences are in JUNCTION_CONDITIONS.md
(wall equations) and GLOBAL_GEOMETRY.md (solutions).

## 1. The system

A spherically symmetric spacetime $(\mathcal M,g)$ divided by a timelike hypersurface $W$
(the wall) into

* an **interior** $\mathcal M_-$ in a false vacuum: energy density $\rho_{\rm in}=3\chi^2/8\pi G>0$
  (cosmological constant $\Lambda_{\rm in}=3\chi^2$), containing a regular centre;
* an **exterior** $\mathcal M_+$ in a lower vacuum $\rho_{\rm out}=3\eta^2/8\pi G\ge0$
  ($\eta=0$: asymptotically flat "laboratory" parent; $\eta>0$: de Sitter parent);
* a **domain wall** $W$ of constant tension $\sigma>0$ ($\kappa=4\pi G\sigma$).

Units $c=1$; signature $(-,+,+,+)$; $G$, $\hbar$ kept where they matter. The baseline is
$D=4$. Parameters: $(\chi,\eta,\kappa,M)$, where $M$ is the exterior Schwarzschild(-de Sitter)
mass parameter, which Birkhoff's theorem makes the only exterior datum.

## 2. Action

$$
\begin{aligned}
S={}&\frac1{16\pi G}\int_{\mathcal M_-}\!\!\sqrt{-g}\,(\mathcal R-2\Lambda_{\rm in})
+\frac1{16\pi G}\int_{\mathcal M_+}\!\!\sqrt{-g}\,(\mathcal R-2\Lambda_{\rm out})
-\sigma\int_W\!\sqrt{-h}\\
&+\frac1{8\pi G}\int_{W}\!\sqrt{-h}\,\big(K^{(-)}-K^{(+)}\big)
+\frac1{8\pi G}\int_{\partial\mathcal M_\infty}\!\!\sqrt{|h|}\,(K-K_0)
+\frac1{8\pi G}\sum_{\rm joints}\int\!\sqrt{\gamma}\;\Theta .
\end{aligned}
$$

* $K^{(\pm)}$ is the trace of the extrinsic curvature of $W$ computed from either side with the
  common normal $n$ (pointing from $\mathcal M_-$ to $\mathcal M_+$). The two Gibbons–Hawking–York
  terms on $W$, one for each region with its own outward normal, combine into $-\frac1{8\pi G}\int_W[K]$.
  They are **not optional**: they equal the distributional part of the Einstein–Hilbert term,
  $\mathcal R\supset-2[K]\,\delta(\ell)$, when the metric is $C^0$ across $W$ with a jump in its normal
  derivative. Hence $S$ is the Einstein–Hilbert action of the glued $C^0$ metric.
* $\partial\mathcal M_\infty$ is an outer timelike boundary (large sphere) with Dirichlet data;
  $K_0$ is the extrinsic curvature of the same boundary isometrically embedded in flat space,
  which makes the action finite for asymptotically flat exteriors. For $\eta>0$ the exterior is closed
  or bounded by a cosmological horizon and this term is replaced as described in TUNNELLING.md.
* Joint (Hayward) terms appear where spacelike initial/final slices meet $W$ or
  $\partial\mathcal M_\infty$ at boost angle $\Theta$. They are needed for the Hamiltonian
  (fixed-slice) variational problem and for the Euclidean actions with corners. In the
  time-symmetric configurations used for tunnelling the slices meet the wall orthogonally.
* Vacuum energy is represented by cosmological constants. Equivalently, it can be represented by
  scalar fields at the minima of $U(\phi)$ (Sec. 6), with $\Lambda=8\pi G\,U(\phi_{\min})$.

## 3. Variational principle and field equations

Varying with the induced metric fixed on $\partial\mathcal M_\infty$ and on the initial and final slices:

1. **Bulk.** $G_{ab}+\Lambda_\pm g_{ab}=0$ in $\mathcal M_\pm$. Birkhoff's theorem gives
   $f_{\rm in}=1-\chi^2r^2$ (a centre-regular de Sitter patch; a mass term would make the centre singular)
   and $f_{\rm out}=1-2GM/r-\eta^2r^2$.
2. **Metric variations on $W$.** On each side
   $\delta(S_{\rm EH}+S_{\rm GHY})\supset\frac1{16\pi G}\int\sqrt{|h|}\,(K_{ab}-Kh_{ab})\,\delta h^{ab}$,
   with the side's outward normal. The two sides combine with opposite normals into
   $-\frac1{16\pi G}\int_W\sqrt{-h}\,([K_{ab}]-h_{ab}[K])\,\delta h^{ab}$, and the wall gives
   $-\tfrac12\int_W\sqrt{-h}\,S_{ab}\delta h^{ab}$ with $S_{ab}=-\sigma h_{ab}$. Stationarity gives the
   **Israel condition** $[K_{ab}]-h_{ab}[K]=-8\pi G\,S_{ab}$. Continuity of $h_{ab}$ is the first junction
   condition, imposed by construction.
3. **Embedding variations of $W$** give the normal equation of motion $S^{ab}\{K_{ab}\}=[T_{nn}]$, i.e.
   $\sigma\{K\}=\rho_{\rm out}-\rho_{\rm in}$. It is implied by 2. together with the bulk equations
   (verified symbolically, JUNCTION_CONDITIONS.md Sec. 3).

All wall equations reduce to
$$
\dot R^2+V(R)=0,\qquad \varepsilon_{\rm in}={\rm sgn}A_{\rm in}(R),\quad\varepsilon_{\rm out}={\rm sgn}A_{\rm out}(R),
$$
with $V$, $A_{\rm in}$, $A_{\rm out}$ given in JUNCTION_CONDITIONS.md.

## 4. Conserved quantities and the Hamiltonian boundary term

The ADM mass of the exterior is $M$ (for $\eta=0$). On solutions of the constraints, the
reduced action takes the form $\int(P\,dR-M\,dt_\infty)$ plus horizon/joint terms
(QUANTIZATION.md). $M$ is the Hamiltonian conjugate to the asymptotic Killing time
$t_\infty$. The mass formula (JUNCTION_CONDITIONS.md Sec. 6) expresses $M(R,\dot R)$.

## 5. Admissible initial data and the notion of "laboratory" data

A Cauchy slice $\Sigma$ with data $(h_{ij},K_{ij})$ satisfying the constraints. We single out:

* **time-symmetric data** ($K_{ij}=0$): slices through a turning point $\dot R=0$. In spherical symmetry
  the constraint reads $2rr''=1-r'^2-8\pi G r^2\rho$ along proper radial distance, with
  $r'=\varepsilon\sqrt f$ on either side of the wall.
* **regular laboratory data**: asymptotically flat, complete, no closed trapped, anti-trapped or
  marginally trapped surfaces. For time-symmetric data this means **no minimal or maximal 2-spheres**,
  i.e. $r'>0$ everywhere. In the thin-wall model this is exactly
  $\varepsilon_{\rm in}=\varepsilon_{\rm out}=+1$ at the turning point (`causal.time_symmetric_slice`).
  These are the initial data an experimenter could in principle prepare, and the ones relevant to
  "creating a universe in the laboratory".
* **bag-of-gold data**: $\varepsilon_{\rm out}=-1$. The slice contains a minimal sphere (the
  Einstein–Rosen bridge at $r=2GM$), which is marginally trapped in both time directions.

## 6. Field-theoretic (thick-wall) completion

To test the idealisation of item 3 of Sec. 1, the wall is replaced by a scalar field:
$$
S_\phi=\int\sqrt{-g}\Big[-\tfrac12(\nabla\phi)^2-U(\phi)\Big],\qquad
U(\phi)=\frac\lambda4(\phi^2-v^2)^2+\rho_f\Big(1-s(\phi/v)\Big),\quad s(x)=\tfrac12+\tfrac34x-\tfrac14x^3 .
$$
$s'(\pm1)=0$ keeps the extrema at $\phi=-v$ (false vacuum, $U=\rho_f$) and $\phi=+v$
(true vacuum, $U=0$) for every $\rho_f$. The flat kink has width $w=\sqrt2/(\sqrt\lambda v)$ and
tension $\sigma_0=\tfrac{2\sqrt2}{3}\sqrt\lambda v^3=\tfrac43v^2/w$. The thin-wall limit is
$w\to0$ at fixed $\sigma_0$ (so $v^2=3\sigma_0 w/4$, $\lambda=2/(w^2v^2)$). Implemented in
`CODE/babyuniverse/thickwall.py`:

* time-symmetric initial data from the Hamiltonian constraint (any profile is admissible data);
* the $O(4)$ Coleman–De Luccia bounce with gravity,
  $\phi''+3(a'/a)\phi'=U'$, $a''=-\tfrac{8\pi G}3a(\phi'^2+U)$, constraint $a'^2=1+\tfrac{8\pi G}3a^2(\tfrac12\phi'^2-U)$.

Full Lorentzian evolution of the coupled system is Phase IV (OPEN_PROBLEMS.md).

## 7. Energy conditions satisfied by the model

| source | stress tensor | NEC $T_{ab}k^ak^b\ge0$ | WEC | SEC |
|---|---|---|---|---|
| vacuum energy | $T_{ab}=-\rho g_{ab}$ | $=0$ (saturated) | yes ($\rho\ge0$) | violated ($\rho>0$) |
| domain wall | $T_{ab}=-\sigma h_{ab}\delta(\ell)$ | $\sigma(n\cdot k)^2\delta\ge0$ | yes | violated |
| canonical scalar | $\nabla\phi\nabla\phi-g(\tfrac12(\nabla\phi)^2+U)$ | $(k\cdot\nabla\phi)^2\ge0$ | if $U\ge0$ | can fail |

For null $k$ the wall's $h(k,k)=-(n\cdot k)^2$, so the NEC holds for every $\sigma\ge0$ (checked
numerically, `causal.nec_wall`). **Every source in the model satisfies the NEC** and violates the SEC.
The SEC violation is what allows inflation and wall repulsion. The NEC is what the Penrose and
Farhi–Guth theorems require (NO_GO_THEOREMS.md).
