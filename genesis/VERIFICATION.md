# Verification dossier

For every claim: how it was verified, by which independent route, and what remains unproven.
Verification types: **S** symbolic (sympy, exact), **A** analytic proof (by hand, in the text),
**N** numerical, **X** cross-check between two independent methods, **L** literature agreement.

## V-J. Junction conditions and potential

| id | claim | type | evidence |
|---|---|---|---|
| V-J1 | $K_{\theta\theta}=RA$, $K^\tau{}_\tau=(\ddot R+f'/2)/A$ for generic $f$ | S + A | `symbolic.extrinsic_curvature_static` (Christoffels); Killing-energy derivation (JUNCTION_CONDITIONS.md Sec. 2) |
| V-J2 | Israel: $A_{\rm in}-A_{\rm out}=\kappa R$; ττ = d/dτ(θθ)/Ṙ | S | `symbolic.israel_equations` |
| V-J3 | Lemma J1 (no spurious branches) | A + N | proof; `test_unsquared_junction_holds_on_shell` (5 models, residual $<10^{-9}$, wrong signs fail) |
| V-J4 | explicit $V$, the two squared forms equal, $V(\sigma)=V(-\sigma)$, $\lvert\gamma\rvert\le1$ | S + N | `symbolic.squaring_and_signs`; `test_potential_forms_agree` |
| V-J5 | dimensionless reduction $V=1-g(z)/E$ | S + N | `symbolic.dimensionless_form`; `test_dimensionless_reduction` |
| V-J6 | mass formula | S + N | `symbolic.mass_formula`; `test_mass_function_inverts_junction` |
| V-J7 | normal EOM $\sigma\{K\}=\rho_{\rm out}-\rho_{\rm in}$ implied | S | `symbolic.wall_equation_of_motion` |
| V-J8 | Ipser–Sikivie limit ($\chi=\eta=M=0$): $R=(2/\kappa)\cosh(\kappa\tau/2)$ | A | from $V=1-\kappa^2R^2/4$ |
| V-J9 | flat-space limit of $H$: $\sqrt{P^2+m^2}$ | S | `symbolic.flat_space_limit` |

## V-G. Global geometry

| id | claim | type | evidence |
|---|---|---|---|
| V-G1 | Proposition G1 | A + N | proof; `test_sign_flip_radii`; A/B boundary $R_1=2GM$ to $10^{-6}$ (`laboratory_regime_bounds.json`) |
| V-G2 | chart-free itinerary = explicit charts | X | 15/15 atlas branches identical (`branch_atlas_cases.json`, `test_bgg_itineraries_and_chart_agreement`) |
| V-G3 | integrated trajectories satisfy $u\cdot u=-1$ and $n^r=A$ (unsquared junction) | N (independent finite differences) | `test_geometry.py` |
| V-G4 | BGG global structure (inflating walls beyond the bridge) | L + X | agrees with Blau–Guendelman–Guth (1987) |
| V-G5 | every inflating branch is past-singular / reaches anti-trapped spheres | N + A | 21/21 (`farhi_guth_audit.json`); consistent with Farhi–Guth |

## V-T. Tunnelling

| id | claim | type | evidence |
|---|---|---|---|
| V-T1 | $P$ satisfies Hamilton's equation with $H=M$ | S (exact rationals) + N | `symbolic.reduced_hamiltonian`; `test_hamilton_equation_all_sectors` |
| V-T2 | T0 branch selection | A + N | analyticity argument; $B_{\rm principal}<0$ in 41/45 models (`smarr_identity_grid.json`) |
| V-T3 | T1 closed form and Morse identity | A + S + N + X | proof; mpmath $J(s)$; 60 models to $7\times10^{-15}$; **thick-wall CdL bounces converge to it** (0.18 % at $w=0.005$, $O(w)$) |
| V-T4 | T2 first law | A + N | proof; finite differences to 8 digits in 4 models |
| V-T5 | T3 Smarr identity | A + N | Euler + T2; 45 models, residual $\le2.5\times10^{-8}$ |
| V-T6 | NG3 exterior embedding obstruction | A + N | proof; degree-theory multiplicity outside $\{0,1\}$ in every regime-A model (grid-converged) |
| V-T7 | NG4 regime-B overshoot | N only | $\min\lvert\theta_D\rvert=3.159>\pi$ over 100 models; **no proof** |
| V-T8 | interior flip embeddable (counter-example) | N | `test_interior_flip_can_be_embedded` |
| V-T9 | T5 laboratory bound | N only | 39 values of $\gamma$, 50 masses each; **no proof** |
| V-T10 | continuity $B(M\to0^+)=B(0)$ | N | $10^{-3}$ relative at $M=10^{-6}M_{\rm cr}$ |

## V-W. Thick walls

| id | claim | type | evidence |
|---|---|---|---|
| V-W1 | scalar time-symmetric data → thin-wall junction and mass formula, $O(w)$ | N | 6 widths, both orientations (`thickwall_convergence.json`, `test_time_symmetric_data_approach_thin_wall`) |
| V-W2 | bag-of-gold data form an Einstein–Rosen bridge at $r=2GM_{\rm ADM}$ | N | $r_{\min}/2GM-1<10^{-3}$ |
| V-W3 | gravitational CdL bounce → thin-wall exponent | N | rel. error 3.1 %, 1.4 %, 0.58 %, 0.18 % at $w=0.04\ldots0.005$; constraint $<3\times10^{-10}$ |

## V-E. Euclidean winding: physical or artefact? (Phase II item)

| id | question | finding |
|---|---|---|
| V-E1 | is the angle coordinate-dependent? | No. $\theta=2\pi\Delta t_E/\beta$ uses the Killing time normalised at infinity and the period fixed by smoothness at the horizon; both are invariant. |
| V-E2 | is the winding a numerical artefact? | No. $\Delta t^W_E=\int A_{\rm out}/(f_{\rm out}\sqrt V)\,dR$ is a convergent integral of smooth functions on the barrier, computed with endpoint-exact quadrature. The large values come from $\beta=8\pi GM\to0$. |
| V-E3 | does it affect physics? | Its only contribution to the exponent is $-\tfrac M2\Delta t^W_E\to0$ as $M\to0$ (T3). It is a property of the exterior's multiple cover (the pseudomanifold), not of the amplitude. |
| V-E4 | interpretation | $\lvert\Delta t^W_E\rvert/\beta$ = number of thermal periods of the black hole spanned by the wall's Euclidean history. With a horizon monodromy (T0) each half-turn is a shift $i\beta/2$ of Lorentzian time. |
| V-E5 | conclusion | **meaningful but representation-specific**: real and invariant, irrelevant to $B$, and it disappears from any formulation that does not build the exterior as a region (canonical, Smarr). |

## V-Q. Quantum transition (Phase III)

| id | check | method | result |
|---|---|---|---|
| V-Q1 | horizon monodromy | complex-contour integration with branch tracking (`contours.continue_P`), 4 models × all horizons × $n=-1,1,2$ | ratio $=n$ to $1.1\times10^{-4}$, $\Delta{\rm Re}\le4\times10^{-15}$ |
| V-Q2 | zero-net periods | two-loop contours | $\Delta B=(A_{\rm dS}-A_{\rm BH})/4G$ to $10^{-6}$, endpoint $P$ unchanged to $10^{-15}$ |
| V-Q3 | conventions (III-2) | analytic vs sign-free, principal vs sign-dressed constraint | $\le8\times10^{-16}$ |
| V-Q4 | branch locus (III-4) | radicand at $R_f$; identity radicand $=\kappa^2R^2A_{\rm out}^2/f_{\rm out}$ along the barrier; sheets of $H$ at both turning points in both conventions | $\le10^{-12}$; table of III-4 holds to $10^{-15}$ |
| V-Q5 | Lorentzian saddle | finite-difference ${\rm Im}\,T^*$ vs Euclidean pieces | $5\times10^{-6}$ |
| V-Q6 | covariant action rule | Smarr pieces + bolt/fold terms | $\le1.2\times10^{-9}$ |
| V-Q7 | competitor exponents | $B_{\rm princ}$, $B_{\rm sw}$ on a 9×4 grid | $B_{\rm princ}<0$ in 33/36; $0<B_{\rm sw}<B_{\rm an}$ in 24/24 single-flip models |
| V-Q8 | coalescence lemma and sector-blind exponent | zeros of $N$ on the barrier; $B_{\rm sq}$ vs $B_{\rm sw}$, $B_{\rm an}$ | one zero in 24/24 single-flip, none in 12/12 double-flip; $B_{\rm sq}=B_{\rm sw}$ resp. $B_{\rm an}$ to $5\times10^{-15}$ |
| V-Q9 | Weber connection coefficient | mpmath `pcfd` at $z=30$ vs $\sqrt{2\pi}/\Gamma(-\nu)$ | $\le4\times10^{-4}$ (finite-$z$), exact zeros at $\nu=0,1$ |

Analytic proofs (III-1, III-2, III-4, III-5, III-6, III-7) are in PHASE_III_QUANTUM_TRANSITION.md. Tests: `tests/test_contours.py`.

## V-D. Field dynamics (Phase IV)

| id | check | result |
|---|---|---|
| V-D1 | symbolic derivation of the double-null equations (sympy) | 5/5 confirmed |
| V-D2 | flat-space pulse vs exact $r\phi=F(v)-F(u)$ | $1.9\times10^{-14}$ |
| V-D3 | exact de Sitter | $r$: about second order; $\ln\alpha$: second order away from the centre, first order at it |
| V-D4 | self-convergence (self-gravitating pulse) | ratios 3.1–4.5 |
| V-D5 | critical collapse bracket, black-hole mass from $r_{AH}$ | $A\in(0.04,0.05)$; $r_{AH}=1.98M$ |
| V-D6 | thin-wall limit of the parent ($R_{\max}$, $r_{AH}$, Israel, $R(\tau)$) | converges as $w\to0$ (Sec. 3 table) |
| V-D7 | thin-wall limit of the daughters ($\dot R(R)$, Israel, signs) | converges as $w\to0$ (Sec. 4) |
| V-D8 | Prop. IV-1 against numerics | 0 anti-trapped points in all laboratory runs (after fix F12) |
| V-D9 | probe Jacobi zero modes | $5\times10^{-6}$ |

## Unproven or unverified items (proof gaps)

1. NG4 (regime B non-embeddability): numerical only.
2. T5 (laboratory bound): numerical only, over the scanned range.
3. Thick-wall analogues of T2, T3, NG3: open (the Lorentzian thick-wall dynamics is done in Phase IV).
4. Agreement of our analytic-branch exponent with FMP/FGG numbers: unverified (full texts inaccessible).
5. SdS exterior: no explicit conformal diagrams; tunnelling formulas are coded for SdS but the Euclidean pieces
   (two horizons, conical defects) are not.
6. All statements are leading-order WKB; prefactors are not computed.
7. Phase III: assumptions (C2) complex angles and (C3) Lebesgue measure are motivated (covariance, analyticity from laboratory data,
   growth bound) but not derived from a gravitational path integral. The Stokes multiplier at the coalescence point $R_q$ of the
   squared constraint is derived at leading local order (Prop. III-8); an exact global confirmation is open (P-Q2). Superselection of the
   exterior orientation (P-Q9) is open.
