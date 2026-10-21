# Obstructions and no-go results

Status labels: **[Theorem]** proved under the stated hypotheses (proof given or cited);
**[Proposition-N]** established numerically with stated coverage, proof open; **[Conjecture]**.

---

## 1. The Farhi–Guth obstruction, stated precisely

**Penrose (1965), time-reversed.** Let $(\mathcal M,g)$ be a spacetime such that

* (P1) $R_{ab}k^ak^b\ge0$ for all null $k$ (the null convergence condition; with the Einstein equations this is the NEC $T_{ab}k^ak^b\ge0$);
* (P2) $\mathcal M$ contains a **non-compact** Cauchy hypersurface;
* (P3) $\mathcal M$ contains a closed (compact, boundaryless) **anti-trapped** 2-surface $S$: both future-directed null expansions are strictly positive.

Then $\mathcal M$ is **past** null geodesically incomplete. The proof shows that some past-directed null
generator of $\partial J^-(S)$ has finite affine length.

**[Theorem] Farhi–Guth (1987).** Assume (P1), (P2), classical GR, and that a false-vacuum region
*inflates*, i.e. contains de Sitter spheres of radius $r>1/\chi$ in its expanding part. These spheres
are anti-trapped, so (P3) holds and the spacetime is past null incomplete. If the spacetime is past
complete (for example, a regular laboratory history whose past is ordinary flat space), no inflating
region can form. The incomplete generator lies on $\partial J^-(S)$, so the singularity is in the causal
past of the inflating region: it is an *initial* (white-hole-like) singularity.

Precisely what this does **not** say:

* It does not use spherical symmetry, a thin wall, or any property of the matter beyond the NEC.
  Thick walls with canonical scalar fields satisfy the NEC, so **thickening the wall cannot evade it**.
* It says nothing about quantum transitions: a tunnelling amplitude is not a Lorentzian spacetime
  satisfying (P1)–(P3).
* It does not apply to closed parents: if the Cauchy surfaces are compact, (P2) fails. GLOBAL_GEOMETRY.md
  Sec. 6 exhibits an explicit globally regular inflating bubble in a de Sitter parent.
* It does not forbid incompleteness to the past of the laboratory slice $\Sigma_0$. If one demands only
  that $J^+(\Sigma_0)\cap J^-(S)$ be regular, the theorem pushes the singularity into $J^-(\Sigma_0)$: the data on
  $\Sigma_0$ must themselves have a singular past. For data an experimenter prepares, the relevant
  condition is past completeness of the whole history, and then the obstruction is absolute.

**Check in the model.** `causal.farhi_guth_audit` over 42 branches on a $(\gamma,E)$ grid covering
regimes A, B and C (`DATA/farhi_guth_audit.json`): all 21 inflating branches reach anti-trapped
spheres in the interior, and all 21 have a past-singular exterior. The NEC holds for the wall
($\min T(k,k)=0$). The model never contradicts the theorem.

## 2. Time-symmetric version: the bag-of-gold data

**[Theorem] NG2.** Let the time-symmetric data at the bounce turning point $R_2$ be given, so that
$\varepsilon_{\rm out}(R_2)=-1$ (Proposition G1). The slice contains a minimal 2-sphere, the Einstein–Rosen bridge at
$r=2GM$, which is marginally trapped in both time directions. Its maximal globally hyperbolic
development contains region II (closed trapped spheres) and region IV (closed anti-trapped spheres).
With (P1) and (P2), it is therefore both future and past null incomplete.

*Proof.* In the thin-wall model the exterior of the development is explicitly Kruskal–Schwarzschild
to the exterior side of the wall, which contains II and IV whenever $\varepsilon_{\rm out}=-1$ at the time-symmetric
slice (GLOBAL_GEOMETRY.md). Then apply Penrose's theorem and its time reverse. More generally, a
minimal surface in NEC-satisfying time-symmetric data has $\theta_\pm=0$, and the Raychaudhuri equation
makes nearby displaced spheres strictly trapped/anti-trapped whenever shear or $R_{ab}k^ak^b$ is non-zero
along the generators. $\square$

Consequence: any history passing through the final state of a "universe-creation" transition has a
singular classical past. A quantum transition can produce such a state only by **replacing the past
half of its classical development**. This is the exact sense in which tunnelling is the only remaining
classical-GR route in an asymptotically flat parent.

By contrast, the initial state at $R_1$ in regime A is regular laboratory data. It has no extremal
spheres (`tests/test_causal.py::test_fgg_endpoint_slices`); see `FIGURES/fgg_slices_embedding.png`.

## 3. A semiclassical no-go: the exterior embedding obstruction

Setting. Euclidean section of a static exterior $f\,dt_E^2+dr^2/f+r^2d\Omega^2$ with a single
non-degenerate horizon $r_h$ ($f(r_h)=0$, $f'>0$ on $(r_h,\infty)$, $f\to1$). With the period fixed by
regularity, the $(r,t_E)$ section is a smooth "cigar" $\mathcal C\cong\mathbb R^2$. Use polar coordinates
$(\rho,\theta)$ with $\rho$ any increasing function of $r$ vanishing at $r_h$, and $\theta=\kappa_s t_E$.
The Euclidean wall solves $R'^2=V$, $f\,dt_E/d\tau_E=A_{\rm out}(R)$.

Hypotheses:

* (H1) the half-bounce runs from a turning point $R_1$ to a turning point $R_2>R_1$ with $R'>0$ in between;
* (H2) $A_{\rm out}$ has exactly one simple zero $R_f\in(R_1,R_2)$, with $A_{\rm out}(R_1)>0>A_{\rm out}(R_2)$.
  Here $\theta$ increases on $(R_1,R_f)$ and decreases on $(R_f,R_2)$;
* (H3) the initial and final exterior slices are the totally geodesic (time-symmetric) slices:
  the ray $\theta=0$, $\rho\ge\rho_1$ (since $\varepsilon_{\rm out}(R_1)=+1$), and the broken line through the tip
  $\{\theta=\theta_D,\ 0\le\rho\le\rho_2\}\cup\{\theta=\theta_D+\pi,\ \rho\ge0\}$ (since $\varepsilon_{\rm out}(R_2)=-1$);
* (H4) the exterior region is an oriented surface $\tilde E$ with a local isometry $\iota:\tilde E\to\mathcal C$
  (Birkhoff: the Euclidean vacuum exterior is locally Euclidean Schwarzschild), whose boundary maps
  onto initial slice ∪ wall ∪ final slice and whose end is the asymptotic end of $\mathcal C$.

**[Theorem] NG3 (exterior embedding obstruction).** Under (H1)–(H4), $\iota$ is not injective. The
Euclidean exterior region interpolating between regular laboratory data and an inflating bounce is
not a domain in Euclidean Schwarzschild.

*Proof.* Orientation (from $n=$ tangent rotated clockwise, JUNCTION_CONDITIONS.md Sec. 7) places
$\iota(\tilde E)$ locally:

* on the counterclockwise side of the initial ray;
* outside (larger $\rho$) the $\varepsilon=+1$ arc $\gamma_+$;
* inside the $\varepsilon=-1$ arc $\gamma_-$;
* on the side $\theta>\theta_D$ of the final segment, and on the side $\theta<\theta_D+\pi$ of the final ray.

Because $R$ increases monotonically, every point of $\gamma_-$ has larger $\rho$ than every point of $\gamma_+$.
Suppose $\iota$ is injective, so the multiplicity $m$ is 1 on the image and 0 off it, and jumps by exactly
one across each boundary piece.

1. *No wrapping.* If $\gamma_-$ spans an angle $\ge2\pi$, then on a generic radial line it has two
   consecutive crossings $\rho_a<\rho_b$. Just outside $\rho_a$, $m=0$; just inside $\rho_b$, $m=1$; nothing
   lies between them. Contradiction. The same argument applies to $\gamma_+$ with the sides reversed.
   So both arcs span less than $2\pi$, with $\theta_{\max}>0$ and $\theta_{\max}-\theta_D<2\pi$.
2. *At infinity.* $m(\infty,\theta)$ changes only across the initial ray ($\theta=0$) and the final ray
   ($\theta\equiv\theta_D+\pi$), each by one with the sides above. Hence $m(\infty,\cdot)=1$ exactly on the
   counterclockwise arc $(0,a)$, $a\in(0,2\pi)$, $a\equiv\theta_D+\pi$, and $0$ on $(a,2\pi)$.
3. *Outside $\gamma_-$.* For a generic $\theta^*\in(\theta_D,\theta_{\max})$ the radial line meets $\gamma_-$ once.
   Just outside it $m=0$, and nothing lies further out, so $m(\infty,\theta^*)=0$.
4. *Cases.*
   * If $\theta_D<0$: angles $\theta^*\in(0,\epsilon)$ satisfy 3., but lie in $(0,a)$. Contradiction.
   * If $0\le\theta_D<\pi$: then $a=\theta_D+\pi$ and $(\theta_D,\min(\theta_{\max},a))\subset(0,a)$ is non-empty. Contradiction.
   * If $\pi\le\theta_D<2\pi$: then $a=\theta_D-\pi$. Take $\theta^*\in(\max(a,\theta_{\max}-2\pi),\theta_D)$, which is
     non-empty and lies in the angular range of $\gamma_+$ but not of $\gamma_-$. Just outside $\gamma_+$, $m=1$, and
     the outward radial path crosses nothing, so $m(\infty,\theta^*)=1$. But $\theta^*\in(a,2\pi)$. Contradiction.
   * $\theta_D\ge2\pi$ would make $\gamma_+$ wrap, which is excluded by 1.

$\square$

For the baseline model, (H1)–(H4) hold for **every regime-A transition**, i.e. every transition from
regular laboratory data to the inflating bounce (Proposition G1(iv); $A_{\rm out}\propto2GM-(\chi^2+\kappa^2)R^3$ is
monotone). The argument uses no property of the interior, so it extends to any interior geometry and
any asymptotically flat exterior with one non-degenerate horizon that satisfies (H1)–(H2).

**Independent numerical confirmation.** The multiplicity is computed by degree theory: the winding
number of the boundary loop around grid points (`tunnelling.multiplicity_map`). In every regime-A model
of the 45-model grid it takes values outside $\{0,1\}$ (`DATA/smarr_identity_grid.json`,
`tests/test_tunnelling.py::test_embedding_obstruction_regime_A`). Typical values are 2 for moderate $M$
(`FIGURES/pseudomanifold_multiplicity.png`). For small $M$ the wall winds many times around the horizon:
$\theta_D=-73$ at $M/M_{\rm cr}\approx10^{-2}$, and multiplicities range from $-3$ to $2$, so
orientation-reversed sheets appear.

This makes precise the **pseudomanifold** of Farhi–Guth–Guven. It is not a failure to find the right
instanton: within the thin-wall dS|Schwarzschild model no manifold instanton exists for the laboratory
transition.

**[Proposition-N] NG4 (regime B).** For transitions from the hidden bounded branch ($\varepsilon_{\rm out}=-1$
throughout, no flip), the slices pass through the tip and the half-bounce region is embedded only if
$|\theta_D|\le\pi$. Across a scan of 100 models covering $\gamma\in[-0.95,0.95]$ and the full width of
regime B, $\min|\theta_D|=3.159>\pi$, and multiplicity maps give 2. Regime-B transitions are therefore
also non-manifold, through a different mechanism (horizon-angle overshoot). A proof is open (P-N1).

**Counter-example to the naive interior analogue.** An $\varepsilon_{\rm in}$ flip inside the barrier does **not** obstruct
an embedding of the interior. In the de Sitter section the wall moves *towards* the Euclidean horizon
as $R$ grows, so step 3 of the proof fails. The computed interior multiplicities are $\{0,1\}$ whenever
the interior does not wrap (`test_interior_flip_can_be_embedded`). The obstruction is tied to a sign flip
**while the wall recedes from a Euclidean horizon**.

## 4. The unifying statement

The three obstructions share one origin, the **creation of a closed marginally trapped surface** (the
bridge) in an asymptotically flat exterior:

| level | obstruction | origin |
|---|---|---|
| classical evolution | Farhi–Guth: past singularity | Penrose focusing of the anti-trapped sphere |
| initial data | NG2: the inflating state has a minimal sphere, so its development is singular both ways | the bridge |
| semiclassical (Euclidean) | NG3: no manifold instanton (pseudomanifold) | $\varepsilon_{\rm out}$ flip ⇔ creation of the bridge |

**[Conjecture] NG5.** For any spherically symmetric NEC-satisfying matter (thin or thick wall), no
smooth Riemannian manifold with the boundary data of regular laboratory data and of minimal-surface
data in an asymptotically flat exterior is a stationary point of the Euclidean action. That is, the
FGG transition never has a real Euclidean saddle that is a manifold, and complex or degenerate
saddles are unavoidable. Evidence: NG3 (thin wall), Fu–Marolf (2019) for AdS exteriors
("smooth Euclidean saddles never produce inflating bubbles"), and Freivogel et al. (2006).

## 5. Topology change and causal structure

* **Geroch (1967):** a compact, time-orientable Lorentzian interpolation between spacelike slices of
  different topology contains closed timelike curves.
* **Tipler (1977); Borde (1994):** with energy conditions, classical topology change requires singularities.
* **Gibbons–Hawking (1992):** selection rules (spin structure, kink number) for Lorentzian cobordisms.

In the present model the Cauchy slices never change topology (GLOBAL_GEOMETRY.md Sec. 5). The daughter's
"separation" is causal plus singular. A genuine topology-changing description would need non-Hausdorff or
degenerate-metric structures. The FGG pseudomanifold is a degenerate structure of exactly this kind,
on the Euclidean side.

## 6. Quantum-gravitational obstructions in the literature

* **Wall (2013),** *The generalized second law implies a quantum singularity theorem*: with the GSL assumed
  on all causal horizons, quantum (generalized-entropy) trapped surfaces imply incompleteness, and it is
  impossible to form baby universes that eventually become independent of the mother universe. Caveat
  (Wall's erratum): the GSL is assumed on causal horizons only. Relation to tunnelling: Wall's
  theorem concerns semiclassical spacetimes with quantum matter, not transition amplitudes between
  classical geometries. Its applicability to the FGG amplitude is an open question (OPEN_PROBLEMS.md P-Q4).
* **Freivogel, Hubeny, Maloney, Myers, Rangamani, Shenker (2006):** in AdS/CFT a pure state cannot evolve
  into the mixed state that describes the exterior of an inflating bubble, so no process (tunnelling
  included) creates an inflating region from a pure state.
* **Fu & Marolf (2019):** smooth Euclidean saddles in AdS never produce inflating bubbles; good saddles do
  produce bags of gold.

## 7. Map of assumptions and escape routes

| assumption | used by | escape route | status / cost |
|---|---|---|---|
| classical GR | FG, NG2 | quantum tunnelling (FGG, FMP) | Euclidean side not a manifold (NG3); canonical side well defined (TUNNELLING.md) |
| NEC | FG, NG2 | NEC-violating matter; quantum stress tensors | QNEC/ANEC, Fewster–Galloway / Fewster–Kontou semiclassical singularity theorems, Wall's GSL theorem constrain it; no controlled model here |
| non-compact Cauchy surface | FG | closed parent (de Sitter) | realised: regular inflating bubble in a dS parent (GLOBAL_GEOMETRY.md Sec. 6), but not *created* from regular data |
| global hyperbolicity | FG | Cauchy horizons (charged/monopole exteriors), topology change | Sakai et al. (2006): monopoles need a collision or tunnelling to inflate; topology change costs CTCs or singularities |
| asymptotic flatness | NG3 | different exterior asymptotics | AdS: Fu–Marolf obstruction; dS: open (P-G1) |
| thin wall | NG3 (as stated) | thick walls | FG unaffected; NG3 thick-wall analogue open (P-T2); Ulvestad–Albrecht (2012) studied FGG with thick walls numerically |
| spherical symmetry | atlas, NG3 | non-spherical collapse | FG unaffected; Aguirre–Johnson (2005): aspherical instabilities of bound bubbles |
