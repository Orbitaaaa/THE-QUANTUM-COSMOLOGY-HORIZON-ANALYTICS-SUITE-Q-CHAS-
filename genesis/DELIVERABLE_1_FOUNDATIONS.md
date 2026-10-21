# Project Baby Universe — Genesis Branch
## Deliverable 1: Mathematical Foundations Report (Stage One)

*Date:* 21 October 2026. *Scope:* tasks M1–M8 of the initial campaign. Detailed derivations are in the
companion files cited in each section; every number is reproducible with `CODE/scripts/run_all.py`
and checked by `CODE/tests` (67 tests).

---

## 1. Precise formulation of the problem

**System.** A spherically symmetric spacetime made of a false-vacuum region (de Sitter, $H_{\rm in}=\chi$), a domain
wall of tension $\sigma$ ($\kappa=4\pi G\sigma$), and an exterior that is Schwarzschild with mass $M$ (the laboratory parent),
or Schwarzschild–de Sitter/de Sitter ($H_{\rm out}=\eta$). Full specification, action and admissible data:
CLASSICAL_MODEL.md.

**Questions, made precise.**

* (Q-classical) Does a solution exist in which regular laboratory data evolve into a configuration containing a
  daughter universe, defined by criteria D1–D3 (inflation; causal disconnection from the parent's future null
  infinity; parent remains asymptotically flat)?
* (Q-quantum) For data that cannot do so classically, is there a consistent transition amplitude to the daughter
  branch? What is it, and is it formulation-independent?

"Regular laboratory data" means time-symmetric, asymptotically flat, complete data without extremal 2-spheres
(CLASSICAL_MODEL.md §5). In the thin-wall model this is exactly $\varepsilon_{\rm in}=\varepsilon_{\rm out}=+1$ at the turning point.

## 2. Literature review (summary)

Full annotated list, with verification status and a novelty table: LITERATURE.md. Main threads:

* **Classical dynamics:** Israel (1966); Blau–Guendelman–Guth (1987); Berezin–Kuzmin–Tkachev (1987);
  Aguirre–Johnson (2005).
* **Obstructions:** Penrose (1965); Farhi–Guth (1987); Vachaspati–Trodden (1999); topology-change theorems
  (Geroch, Tipler, Borde, Gibbons–Hawking); Wall's GSL theorem (2013).
* **Tunnelling:** Coleman–De Luccia (1980); Parke (1983); Lee–Weinberg (1987); Farhi–Guth–Guven (1990);
  Fischler–Morgan–Polchinski (1990); Chen–Hu–Yeom (2016); de Alwis et al. (2019/2020);
  Céspedes et al. (2023/2024); Pasmatsiou (2019); Ulvestad–Albrecht (2012); Wang et al. (2023/2025).
* **Holography and information:** Freivogel et al. (2006); Fu–Marolf (2019); Sahu–Simidzija–Van Raamsdonk (2023);
  Marolf–Maxfield (2020).

Status of the field in one paragraph. The classical obstruction is a theorem. The quantum route has a canonical
calculation that the literature considers consistent (FMP; de Alwis et al.). Its Euclidean counterpart requires an
object that is not a manifold (FGG), and holographic arguments (Freivogel et al.; Fu–Marolf) suggest that inflating
bubbles do not arise from smooth saddles or unitary processes in AdS. What was missing: a precise statement of
*why* the Euclidean object must be singular, and a formulation-independent definition of the exponent. Stage One
supplies both (Sections 8–10).

## 3. Baseline gravitational action

$$
S=\sum_\pm\frac1{16\pi G}\int_{\mathcal M_\pm}\sqrt{-g}(\mathcal R-2\Lambda_\pm)-\sigma\int_W\sqrt{-h}
-\frac1{8\pi G}\int_W\sqrt{-h}\,[K]+\frac1{8\pi G}\int_{\partial\mathcal M_\infty}\sqrt{|h|}(K-K_0)+\text{joints}.
$$
The wall GHY pair equals the distributional Einstein–Hilbert curvature $\mathcal R\supset-2[K]\delta$. Variations give the
Einstein equations, the Israel condition, and the wall's normal equation of motion (CLASSICAL_MODEL.md §2–3).
Energy conditions: every source satisfies the NEC and violates the SEC (§7).

## 4. Israel junction conditions (derived)

With $n$ from interior to exterior, $A\equiv n^r=\varepsilon\sqrt{f+\dot R^2}$:
$$
K_{\theta\theta}=RA,\qquad K^\tau{}_\tau=\frac{\ddot R+\tfrac12f'}{A},\qquad
A_{\rm in}-A_{\rm out}=\kappa R,\qquad \Big[\frac{\ddot R+\tfrac12f'}{A}\Big]=-\kappa .
$$
The ττ equation is the τ-derivative of the θθ one. Exact sign resolution:
$$
A_{\rm in}=\frac{f_{\rm in}-f_{\rm out}+\kappa^2R^2}{2\kappa R},\quad A_{\rm out}=\frac{f_{\rm in}-f_{\rm out}-\kappa^2R^2}{2\kappa R},
$$
and **Lemma J1**: no spurious solutions arise when the signs are read off. Signs flip only inside dynamical
regions (Lorentzian) or static regions (Euclidean), as stated in Lemma J2. All of this is verified symbolically
(JUNCTION_CONDITIONS.md).

## 5. Effective bubble-radius dynamics

$$
\dot R^2+V=0,\qquad V=1-\frac{L^4}{4\kappa^2}R^2+\frac{GM(\chi^2-\eta^2-\kappa^2)}{\kappa^2R}-\frac{G^2M^2}{\kappa^2R^4}
=1-\frac{g(z)}{E},\quad g=z^2-\frac{2\gamma}z+\frac1{z^4}.
$$
Barrier for $M<M_{\rm cr}=4\kappa^3/(Gg_{\min}^{3/2}L^4)$; bounded and bounce branches; monotonic above $M_{\rm cr}$.
At $M=0$ there is a single turning point $R_0=2\kappa/L^2$.

## 6. Classification of global branches

Regimes A, B, C in the $(\gamma,E)$ plane (`FIGURES/phase_diagram_schwarzschild.png`). Proposition G1:

* every inflating bounce has $\varepsilon_{\rm out}=-1$ (it lies beyond the Einstein–Rosen bridge);
* the bounded branch is visible from region I iff $E>2-2\gamma$ (regime A);
* at the A/B boundary $R_1=2GM$.

The chart-free itinerary (Proposition G2) and explicit Kruskal / global de Sitter integration agree on all
15 atlas branches (`FIGURES/penrose_*.png`). Inflating and visible-from-the-parent are mutually exclusive. Every
inflating branch is a daughter in the sense D1–D3, and every one has a singular past. A de Sitter parent admits a
globally regular inflating bubble (compact Cauchy surfaces). There is no spatial topology change on regular slices.
Details: GLOBAL_GEOMETRY.md.

## 7. Assumptions behind the singularity obstruction

Penrose (time-reversed): NEC + non-compact Cauchy surface + closed anti-trapped surface ⇒ past null incompleteness.
Farhi–Guth: inflation supplies the anti-trapped surface. The theorem uses neither spherical symmetry nor a thin
wall, so thick walls cannot evade it. It does not apply to closed parents or to quantum transitions. **NG2**: the
inflating state's time-symmetric slice contains a minimal sphere, so any classical history through it is singular
both to the past and to the future. A full table of assumptions and escape routes is in NO_GO_THEOREMS.md §7.

## 8. Mathematical structure of the Farhi–Guth–Guven proposal

Initial state: the regime-A bounded branch at $R_1$ (regular laboratory data). Final state: the bounce at $R_2$
(bag of gold + bridge). Forbidden interval $(R_1,R_2)$, with the same $M$. Canonical variables $(R,P)$ with the
asymptotic clock and $H=M$:
$$
P=\frac RG\Big[\operatorname{artanh}\frac{\dot R}{A_{\rm out}}-\operatorname{artanh}\frac{\dot R}{A_{\rm in}}\Big],\qquad
2\sqrt{f_{\rm in}f_{\rm out}}\cosh\frac{GP}R=\varepsilon_{\rm in}\varepsilon_{\rm out}(f_{\rm in}+f_{\rm out}-\kappa^2R^2).
$$
Under the barrier $P=iP_E$ with $P_E=\frac RG(\Theta_{\rm out}-\Theta_{\rm in})$, and the WKB exponent is $B=2\int_{R_1}^{R_2}P_E\,dR$
(TUNNELLING.md, QUANTIZATION.md).

## 9. Critical assessment of the canonical and Euclidean difficulties

* **Canonical ("slicing").** Constant-Killing-time slices are singular at horizons. In the Lorentzian region the
  momentum is complex inside horizons. In the barrier, $f>0$ on both sides, so the reduced clock is regular. The only
  remaining issue is the branch of $\Theta$. **T0** resolves it: the analytic continuation along the barrier is the WKB
  phase. The alternative branch is the analytic one plus a horizon monodromy ($P\to P\pm i\pi R/G$, $t\to t\pm i\beta/2$),
  and it gives $B<0$ in 41/45 models. Assessment: **resolved within the reduced theory**. Whether the full theory
  selects the same contour is Phase III's problem P\*.
* **Euclidean ("pseudomanifold").** **NG3**: for every laboratory transition the exterior Euclidean region is not a
  domain in Euclidean Schwarzschild. Degree-theory multiplicities are $\ge2$, and negative (orientation-reversed sheets)
  at small $M$ (`FIGURES/pseudomanifold_multiplicity.png`). Regime B is also non-manifold, numerically (NG4). The
  interior is not obstructed (counter-example). Assessment: **no manifold instanton exists in this model**. The
  FGG object is forced, not a failure of search.
* **The bridge between them.** **T3 (Smarr identity)** expresses the canonical exponent through the interior
  Euclidean volume, the wall area and the exterior Killing time along the wall. The pseudomanifold therefore
  enters only through $-\tfrac M2\Delta t^W_E$. The winding of the Euclidean wall around the horizon (dozens of turns
  at small $M$) is invariant, but it is a property of the multiple cover, not of $B$ (VERIFICATION.md V-E).

## 10. Ranked new research approaches

N3 (monodromy-fixed contours / Picard–Lefschetz) > N1 (branch-aware quantization of the squared constraint) >
N4 (covariant derivation of the Smarr identity as an action of branched coverings) > N5 (thick-wall Smarr and NG3) >
N2 (crossed-product / relational clock) ≈ N6 (de Sitter parent) > N7, N8. Scores and rationale: NEW_APPROACHES.md.

## 11. First independently checked calculations

| calculation | independent check | agreement |
|---|---|---|
| extrinsic curvature, Israel, potential, mass, normal EOM | sympy vs hand derivation | exact |
| global itineraries | chart-free vs Kruskal/dS integration | 15/15 branches |
| trajectories | finite differences of integrated coordinates vs unsquared junction | $<10^{-3}$ median |
| $M=0$ exponent | numerical quadrature vs closed form | $7\times10^{-15}$ (60 models) |
| $M=0$ exponent | thick-wall CdL bounces with gravity | 0.18 % at $w=0.005$, $O(w)$ |
| thin-wall junction | scalar time-symmetric data | $O(w)$, both orientations |
| first law | finite differences vs Euclidean pieces | 8 digits |
| Smarr identity | canonical $B$ vs covariant pieces | $2.5\times10^{-8}$ (45 models) |
| Hamilton's equations | numerical $\partial H/\partial P$ vs $dR/dt$ | $6\times10^{-10}$, three sectors |
| pseudomanifold | multiplicity maps, grid-converged | values $\notin\{0,1\}$ in every regime-A model |

## 12. Path toward an original theorem / consistent quantum construction

Already obtained: T1 (closed form + Morse identity), T2 (first law), T3 (Smarr identity), NG3 (embedding
obstruction), plus numerical propositions NG4 and T5 (RESULTS.md). The path to a *consistent quantum
construction* is problem P\* (RESEARCH_QUESTIONS.md):

1. solve the $M=0$ area-shift constraint exactly and recover $e^{-B/2\hbar}$ with $B$ from T1;
2. show that the fixed-energy Lorentzian path integral's contour reduces to the real barrier class, with
   horizon-encircling contours identified as thermal (half-period) contributions;
3. derive T3 from a path integral over branched coverings (degree-weighted exterior action).

Success would yield a theorem: *the FGG amplitude is uniquely defined at leading order, equals $e^{-B/2\hbar}$ with
$B$ given covariantly by T3, and the pseudomanifold is the Euclidean shadow of a horizon monodromy.* Failure in step 2
(non-unique Stokes structure) would itself be a sharp result: the ambiguity is physical and must be fixed by a boundary
condition at the horizon.

## 13. Verification plan

Already in place: VERIFICATION.md (every claim with type of check), a 67-test suite run in CI, and the
`run_all.py` reproduction. Next steps:

* **(V1)** proofs of NG4 and T5 or sharp counter-examples;
* **(V2)** term-by-term comparison with FGG/FMP and de Alwis et al. once full texts are available;
* **(V3)** high-precision (mpmath) re-evaluation of the Smarr identity at extreme parameters;
* **(V4)** a second implementation of $B$ written without reference to the canonical code. Partly in place: the
  Smarr check evaluates the interior volume, wall area and Killing time by separate quadratures. A fully separate
  Euclidean action code (including corner terms) remains to be written;
* **(V5)** SdS extension with explicit two-horizon charts.

## 14. Continuation plan

Following the user's Phase II–VI programme:

* **Phase II** (independent verification): V1–V5 above; scrutiny of the sign sectors (done, F2) and of the winding
  (VERIFICATION.md V-E, done at the Stage One level).
* **Phase III** (quantum transition): problem P\*; exact $M=0$ amplitude; contour analysis; covariant derivation of
  the Smarr identity.
* **Phase IV** (fields): spherically symmetric Einstein–scalar evolution (double-null), thick-wall Smarr identity,
  aspherical stability.
* **Phase V** (full quantum system): physical Hilbert space with sectors; crossed-product exterior algebra; information
  bookkeeping of the Morse/endpoint areas.
* **Phase VI** (observability): early-universe daughter production versus laboratory rates; primordial black holes.

### What has actually been established, and what has not

* Established: the classical picture (reproduced and sharpened); that the semiclassical laboratory transition has
  **no manifold instanton** in this model (theorem); a **unique, covariantly expressible** leading-order exponent
  within the reduced canonical theory (theorems T1–T3 with T0); and that the exponent is **not smaller than** the false
  vacuum's own decay exponent for γ ≳ 0.2 (numerical).
* Not established: that the full quantum gravity amplitude agrees with the reduced theory; prefactors; any statement
  beyond spherical symmetry and thin walls (except where tested); that a daughter universe could be created by any
  practical means. Nothing here suggests a laboratory route. The required initial state is a gram-scale lump of
  GUT-scale false vacuum confined to about $10^{-29}$ m (1–2 of its own Schwarzschild radii), and even then
  $B\sim10^6$ (TUNNELLING.md T6).
