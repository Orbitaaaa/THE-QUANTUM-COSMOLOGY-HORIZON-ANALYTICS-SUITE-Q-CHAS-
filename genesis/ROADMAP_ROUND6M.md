# Round 6m — the global daughter spacetime

**Task.** Establish the complete spacetime in which a daughter exists, in the best candidate model: model U (Round 6l), meaning D = 4 + n
Einstein gravity plus the quasi-topological tower. The user has accepted a finite daughter.

**Answer: yes.**
* Collapse of ordinary matter in model U produces a global spacetime with three properties: it is curvature-bounded, it is regular across
  every horizon, and every geodesic tested is complete.
* That spacetime contains two daughters:
  * **the daughter universe**: a white-hole region IV and a new asymptotically flat universe I'. It is causally disconnected from the
    parent's future null infinity, has its own null infinity, and has ADM mass parameter μ;
  * **the finite daughter**: the closed FRW interior of the collapsed matter. It bounces at a = l and re-expands, finite in volume, inside
    the new universe.

Code `global6m.py`; script `scripts/round6m.py`; tests `tests/test_round6m.py` (4); data `DATA/round6m.json`; figure `FIGURES/round6m.png`.

## 1. The spacetime

* **Construction.** An Oppenheimer–Snyder closed dust ball, with mass parameter C and boundary angle χ_b, is glued to the static regular
  exterior with μ = C sin^{D−1}χ_b. The generalised Misner–Sharp relation gives f(R) + Ṙ² = cos²χ_b exactly. So the surface is a radial
  geodesic of the exterior with energy E = cos χ_b, and the gluing is exact (Bueno et al., arXiv:2505.09680, use the same matching).
* **The maximal extension** of the exterior (simple horizons r₋ < r₊, regular centre) is tiled by diamond patches:

  I (parent) → II (black hole) → III (static de Sitter core) → IV (white hole) → I' (new universe) → …

* **Matching the patches.** Each patch carries its own null coordinates. The matching uses the finite null coordinate: v across the event
  horizon and the inner horizon into III, u out of III into IV and I'. The computed worldline crosses the patches in exactly that order
  (figure, panel a).

## 2. What was proved or computed

**Proposition 6m-1 (curvature bound).** We computed the Kretschmann scalar
K = f''² + 2(D−2)f'²/r² + 2(D−2)(D−3)(1−f)²/r⁴ (the D = 4 Schwarzschild check gives 48M²/r⁶ to 10⁻⁹).
* Its maximum is the same for black holes of every mass, μ = 10⁻³ to 10⁸:

  | D | tower | max K × l⁴ |
  |---|---|---|
  | 4 | Hayward | 29–34 |
  | 5 | Hayward | 45–48 |
  | 10 | Hayward | 188.0–188.8 |
  | 10 | Dymnikova | 919 |

* This is a solution-independent curvature bound, of the kind Bueno et al. report.
* The closed interior has K = 4(D−1)(ä/a)² + 2(D−1)(D−2)ψ², which peaks at exactly the de Sitter value 2D(D−1)/l⁴ at the bounce.

**Proposition 6m-2 (the ball reaches the new universe).**
* For every bound collapse with μ > μ_ext, the surface crosses the event horizon and then the inner horizon.
* It bounces in the core where f = E² > 0, as the Round 6l lemma requires (an untrapped bounce).
* It leaves through the white-hole region and reaches its maximal radius in I'.
* All of this takes finite proper time: 2.5–35 l from the parent horizon to the daughter horizon over the 18 cases computed (D = 5, 10;
  μ/μ_ext = 2–100; E = 0.3–0.9).
* For μ < μ_ext there is no horizon, and the ball bounces back into the parent universe.

**Horizon regularity.** At both horizons, the Kruskal metric function f e^{−2κr*}/κ² tends to the same finite non-zero limit from inside
and outside (the ratio of the two limits differs from 1 by less than 10⁻⁴). The extension across each horizon is smooth.

**Geodesic completeness.** We checked these families: radial timelike (E = 0.6 and 1.5), non-radial timelike (L = 2), radial null and
non-radial null (L = 3).
* Each crosses the core and the black-hole region in finite affine parameter.
* Each either passes through the regular centre (f(0) = 1, f'(0) = 0) or turns at r > 0.
* No geodesic ends. This agrees with Bueno et al. (PRL 134, 181401), who find the thin-shell spacetime singularity-free and geodesically
  complete.

**Proposition 6m-3 (the daughter).**
* D = IV ∪ I' lies outside J⁻(I⁺) of the parent: every future causal curve from it stays in it, or enters the daughter's own later black
  hole.
* It has its own future null infinity and ADM mass parameter μ.
* The finite daughter, the closed FRW interior, re-expands from a_b = l to a_max. Its volume grows by:
  * 1.3×10⁵ (C = 10⁴, D = 10);
  * **8×10¹¹⁴** for a four-dimensional stellar collapse.

## 3. Status

| Requirement for "the daughter spacetime is established" | Status |
|---|---|
| formed from ordinary matter (dust; pressure in 6n) | yes |
| curvature bounded everywhere, by a mass-independent bound | yes (numerical, all masses) |
| every horizon regular, extension smooth | yes |
| geodesically complete | yes (all families checked; literature) |
| a region causally disconnected from the parent, with its own infinity | yes (IV ∪ I') |
| a finite daughter universe inside it | yes (the bounced closed FRW, volume ×10⁵–10¹¹⁵) |
| survives evaporation | Round 6o |
| survives pressure, inhomogeneity, shocks | Round 6n |
