# New approaches, ranked

Ranking criteria (brief Sec. 16): relevance, mathematical clarity, novelty potential, testability, foundational
significance, tractability, connection to later stages. Score 1–5 each; the total decides the order.

| rank | approach | rel | clar | nov | test | found | tract | conn | total |
|---|---|---|---|---|---|---|---|---|---|
| 1 | N3 Monodromy-fixed WKB / complex-contour formulation | 5 | 5 | 4 | 5 | 4 | 4 | 5 | 32 |
| 2 | N1 Branch-aware quantization of the squared constraint | 5 | 4 | 4 | 4 | 4 | 4 | 5 | 30 |
| 3 | N4 Covariant derivation of T3 from a degenerate-metric path integral | 5 | 3 | 5 | 4 | 5 | 3 | 4 | 29 |
| 4 | N5 Thick-wall Smarr identity and thick-wall NG3 | 4 | 4 | 3 | 5 | 3 | 4 | 4 | 27 |
| 5 | N2 Relational-clock / crossed-product Hilbert space for the exterior | 4 | 3 | 4 | 3 | 5 | 3 | 4 | 26 |
| 6 | N6 dS parent: creation from regular data in a closed universe | 4 | 4 | 3 | 4 | 3 | 4 | 4 | 26 |
| 7 | N7 Non-Hausdorff / cobordism description of the separation | 3 | 2 | 4 | 2 | 5 | 2 | 3 | 21 |
| 8 | N8 Singularity resolution in effective quantum cosmology applied to the white-hole past | 3 | 2 | 3 | 2 | 4 | 2 | 3 | 19 |

## N3. Monodromy-fixed WKB and complex contours (Phase III, first target)

The canonical exponent is the integral of an analytic function $P(R)$ whose only branch points are the turning
points and the horizons $f=0$ (T0). *Proposal:* define the amplitude by the contour class of the real barrier
interval. Prove that every Lorentzian / Picard–Lefschetz evaluation of the minisuperspace path integral, at
fixed energy $M$, deforms to this class unless it encircles a horizon branch point. Then show that encircling
one corresponds to inserting a half-period $i\beta/2$ of Euclidean time, i.e. a *thermal* contribution. This
would turn the FGG ambiguity into a precise statement: the alternative prescriptions are thermally weighted
contributions of a different initial state, not alternative definitions of the same amplitude.
*Test:* M = 0, where the exact solution (N1) is available.

## N1. Branch-aware quantization

The squared constraint $\mathcal C^2$ (QUANTIZATION.md Sec. 2) contains all four orientation sectors as sheets.
Quantize $\mathcal C^2$ as a second-order difference operator with imaginary step in the area variable.
Solutions then carry sector labels automatically, and the sectors communicate exactly where the classical theory
says they do (under the barrier). At $M=0$ the problem is an ODE in the momentum representation. *Test:*
asymptotics vs T1; the self-adjoint extension and its inner product; the sector content of the physical
states.

## N4. Covariant meaning of the Smarr identity

T3 says what the exterior pseudomanifold must contribute: $-\tfrac M2\Delta t^W_E+\mathcal A/4G$. *Proposal:* show that
this is the action of the multiply covered exterior counted with the degree of the covering map, with
horizon (bolt) terms $\propto(\text{covering angle}-2\pi)A_h/8\pi G$ and Hayward corner terms. If so, the FGG
"pseudomanifold action" becomes a theorem about branched coverings, and complex or degenerate saddles
acquire a definite weight. *Test:* the Smarr grid (any proposed rule must reproduce it to $10^{-8}$).

## N5. Thick walls

Extend the first law and the Smarr identity to scalar-field walls. Euler homogeneity still holds; the first
law needs $\partial S/\partial\lambda$ for the potential couplings. Test whether the thick-wall Euclidean exterior is
still non-embedded: Ulvestad–Albrecht's thick-wall FGG solutions are the natural data.

## N2. Relational clocks and crossed products

The exterior algebra with boundary Hamiltonian $M$ and the boundary clock $t_\infty$ is a crossed-product
setting. *Questions:* is the reduced Hilbert space of Route A the physical Hilbert space of a crossed product?
How does the daughter sector appear (as a different superselection sector, or as a mixed state)? This connects
to Freivogel et al. and to Marolf–Maxfield. **The manuscript "Quantum Gravity from the Crossed Product" is not
assumed**; its Theorems T and L must be checked in this specific setting.

## N6. de Sitter parent

In a closed de Sitter parent the Farhi–Guth theorem does not apply, and a regular inflating bubble exists
(GLOBAL_GEOMETRY.md Sec. 6). *Question:* can such a configuration be reached from regular data with a localised
false-vacuum seed in SdS ($M>0$, $\eta>0$)? Classically, quantum mechanically (Pasmatsiou 2019), and how does
B compare with Lee–Weinberg up-tunnelling?

## N7 / N8. Longer-term

Non-Hausdorff or cobordism descriptions of the causal separation; quantum-gravity resolution of the white-hole
past (limiting curvature, effective bounces). Low tractability now; to be revisited after Phase III.
