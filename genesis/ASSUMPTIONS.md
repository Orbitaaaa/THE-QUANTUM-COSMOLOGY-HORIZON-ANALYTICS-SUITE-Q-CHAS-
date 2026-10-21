# Assumptions register

Every result in this programme is conditional on the assumptions below. Status:
**adopted** (a modelling choice, kept fixed in Stage One), **tested** (checked against a less idealised
model), **derived** (follows from earlier assumptions), **relaxed** (removed in a later result).

| id | assumption | status | where it matters | what relaxing it would require |
|---|---|---|---|---|
| A1 | Classical general relativity in 4 dimensions, Einstein–Hilbert action, signature $(-,+,+,+)$ | adopted | everything | higher dimensions, higher-curvature terms |
| A2 | Spherical symmetry | adopted | atlas, tunnelling, NG3 | non-spherical perturbations (Aguirre–Johnson instability); FG theorem does not need it |
| A3 | Vacuum regions with constant energy densities $\rho_{\rm in}>\rho_{\rm out}\ge0$ | adopted | $f_{\rm in},f_{\rm out}$ | rolling fields inside the bubble (inflaton dynamics) |
| A4 | Thin wall with constant tension, $S_{ab}=-\sigma h_{ab}$, $\sigma>0$ | adopted; **tested** | junction conditions, tunnelling | thick walls: time-symmetric data and CdL bounces converge to the thin wall as $O(w)$ (thickwall.py); full Lorentzian evolution is Phase IV |
| A5 | Interior is a centre-regular de Sitter patch (no interior mass) | derived (regularity + Birkhoff) | $f_{\rm in}$ | interior black holes / charged cores |
| A6 | Exterior: Schwarzschild ($\eta=0$) or SdS/dS; for $\eta=0$ asymptotically flat | adopted | atlas, NG3, T1–T5 | AdS (Fu–Marolf), charged (RN: Cauchy horizons) |
| A7 | NEC holds for all sources | derived for the model (CLASSICAL_MODEL.md Sec. 7) | Farhi–Guth, NG2 | quantum stress tensors; exotic matter |
| A8 | Globally hyperbolic, time-orientable spacetimes with non-compact Cauchy surfaces (laboratory parent) | adopted | Farhi–Guth | closed parents (relaxed in GLOBAL_GEOMETRY.md Sec. 6), topology change |
| A9 | "Laboratory" initial data = time-symmetric, asymptotically flat, no extremal spheres | adopted (definition) | choice of initial state for tunnelling | moving initial data; non-time-symmetric preparation |
| A10 | Semiclassical (WKB) approximation, leading exponential order | adopted | $B$ | prefactors, loop corrections, validity $B\gg\hbar$ |
| A11 | Reduced phase space: constraints solved classically, asymptotic Killing time as clock, $H=M$ | adopted | T0–T3 | Dirac quantization (Route B, QUANTIZATION.md) |
| A12 | WKB phase defined by analytic continuation along the real barrier interval | **derived** from analyticity, positivity and the $M\to0$ limit (T0) | value of $B$ | the alternative branch is excluded by T0 |
| A13 | The initial quantum state is localised on the bounded branch (decaying under the barrier) | adopted | sign of the exponent | superpositions of branches |
| A14 | Euclidean continuation through moments of time symmetry, exterior period fixed by regularity at the horizon | adopted | NG3, T3 pieces | complex (non-Euclidean) saddles, Lorentzian contours (Phase III) |
| A15 | Planck-unit estimates assume $\rho=E_{\rm vac}^4$, $\sigma=E_{\rm wall}^3$ | adopted | T6 table | specific particle-physics models |
| A16 | The crossed-product manuscript is **not** assumed | adopted | — | its claims must be verified separately |

Assumptions that turned out unnecessary or false during Stage One:

* *"Squaring the junction condition creates spurious branches"*: false once the signs are read off
  (Lemma J1). Only negative-tension solutions are spurious.
* *"An interior orientation flip obstructs an embedded interior instanton"*: false. Computed counter-example;
  the obstruction is specific to the exterior (NO_GO_THEOREMS.md Sec. 3).
* *"The FGG exponent scales as $(M_{\rm Pl}/H)^2$"*: only for $\kappa\sim\chi$ (TUNNELLING.md T6).
