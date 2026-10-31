# Failed approaches and dead ends (kept, with reasons)

Each entry records what was tried, why it failed, what the failure taught us, and whether the obstruction is
fundamental, model-specific, approximation-dependent or methodological.

### F1. Principal-branch canonical momentum
*Tried:* $P_E=(R/G)[\arctan(\sqrt V/A_{\rm out})-\arctan(\sqrt V/A_{\rm in})]$, the formula one writes first.
*Failure:* $B<0$ in 41/45 grid models (exponential enhancement, impossible for decay from a localised state). It
jumps by $\pi R/G$ at the flip radius, and it does not tend to any manifold instanton at $M\to0$.
*Lesson:* the jump is a horizon monodromy inserted at a regular point (T0). **Type:** methodological (wrong
branch), now resolved.

### F2. Bracketed inversion of $H(R,P)$ across sign sectors
*Tried:* find $M$ with $P(R,\dot R(R,M);M)=P$ by bracketing over a global $M$ grid.
*Failure:* "no root in bracket" in the $(+,-)$ sector. First cause: the branch sign of $\dot R$ was inferred from
${\rm sgn}\,P$, which is wrong where $\varepsilon_{\rm out}=-1$ (there $P<0$ for $\dot R>0$). Second cause: the admissible $M$-interval
at fixed $R$, which requires the same sector, a static region and $|\dot R|<|A|$, can be narrow, so global brackets miss it.
*Fix:* damped Newton iteration in $M$ using the exact $\partial P/\partial M$, constrained to the sector.
*Lesson:* at fixed $R$, $H$ is defined on a sector-dependent domain. The sectors are genuinely different sheets.
**Type:** methodological, resolved (Phase II item).

### F3. Quantizing the effective potential
*Tried:* treat $\dot R^2+V=0$ as a particle in a potential, with mass $m(R)$.
*Failure:* WKB exponents disagree with the gravitational one by factors between 1.06 and 13, and depend on the
arbitrary $m(R)$.
*Lesson:* orbits do not determine the symplectic structure; only the reduced gravitational $P$ is admissible.
**Type:** fundamental for this approach.

### F4. "An interior orientation flip also obstructs an embedding"
*Tried:* extend NG3 to $\varepsilon_{\rm in}$ by symmetry.
*Failure:* computed counter-examples ($\gamma=0.3,0.6$): the interior region with an $\varepsilon_{\rm in}$ flip has multiplicity in
$\{0,1\}$.
*Lesson:* the obstruction needs the wall to **recede** from the Euclidean horizon at the flip. In the de Sitter section
it approaches it. **Type:** the conjecture was false; the theorem's hypotheses are now sharper.

### F5. Multiplicity maps with a badly sized window (two successive errors)
*Tried (1):* degree-theory multiplicity on a fixed window of radius 2.5 in the cigar chart.
*Failure (1):* two small-$M$ regime-A models were reported "embedded" because the wall extended beyond the window
($\rho$ up to 7).
*Tried (2):* a window sized from the largest radius of the boundary loop.
*Failure (2):* that radius is set by the radial rays, which run to the far boundary, so the window became huge and
the grid too coarse to resolve the multiplicity-2 region near the wall. Seven regime-A and all nine regime-B
models were then reported "embedded".
*Fix:* size the window from the wall alone (1.25 × its largest $\rho$); beyond the wall the multiplicity can change
only across the two rays. Result: all 36 regime-A and 9 regime-B grid models are non-embedded and all interiors
embedded, with multiplicities $-3\ldots2$. A grid-refinement test is now part of the suite.
*Lesson:* a numerical certificate must cover the whole object **and** resolve it. **Type:** methodological.

### F6. Endpoint quadrature of the Euclidean trajectory
*Tried:* $\int dR/\sqrt V$ with $V$ evaluated near the turning points.
*Failure:* round-off made $V$ slightly negative within $10^{-17}$ of $R_1$, the integrand blew up ($\theta\sim10^{133}$),
and a Smarr residual of $10^{284}$ appeared.
*Fix:* parametrise $R=R_a+(R_b-R_a)\sin^2(s/2)$, compute the distances to the endpoints without cancellation,
and use the exact endpoint limit $\sqrt{(R_b-R_a)/|V'|}$.
**Type:** numerical, resolved; recorded because the error was silent in most runs.

### F7. Closed-form CdL exponent by subtraction
*Tried:* $B_{\rm CdL}=S_E-S_E[S^4]$ evaluated by subtracting two numbers of size $\pi/G\chi^2$.
*Failure:* for $\kappa/\chi=3.6\times10^{-3}$ the subtraction lost all digits ($1.7\times10^4$ instead of 133).
*Fix:* algebraically cancellation-free forms (`tunnelling._stable_cap_cosines`).
**Type:** numerical, resolved.

### F8. A naive thin-wall tension for the tilted potential
*Tried:* $\sigma=\int\sqrt{2U}\,d\phi$ including the vacuum-energy tilt, as input to the thin-wall formula.
*Outcome:* it converges to $\sigma_0$, but slowly. Excess $+18\%,+10\%,+5.6\%,+3.1\%$ at
$w=0.04,0.02,0.01,0.005$, roughly $\propto w^{0.85}$. The CdL exponent built from it is off by 9 % at $w=0.005$,
against 0.18 % with $\sigma_0$.
*Lesson:* the thin-wall tension is that of the degenerate part; the tilt's contribution vanishes only as $w\to0$
and is better treated as part of the volume term. **Type:** approximation-dependent (not a failure of the
thin-wall limit itself).

### F9. Shooting the CdL bounce from the centre
*Tried:* shoot on $\phi(0)$ directly.
*Failure:* for thin walls the required $|\phi(0)-v|\sim e^{-mR_0}$ is far below double precision. Also, near the far
pole the Hubble term becomes anti-friction and amplifies residual errors.
*Fix:* start from the linearised Bessel solution at the radius where the deviation is $10^{-6}v$; cut the numerical
solution at closest approach to the false vacuum and complete it with the exact de Sitter cap.
**Type:** numerical, resolved.

### F10. "The principal prescription is a negative-tension amplitude" (Phase III, first draft)
*Claimed:* beyond the flip radius the principal-branch momentum satisfies $\cos(GP_E/R)=-C$. This is the "other sheet" of the
squared constraint, which the first draft identified via Lemma J1 with the negative-tension sector. On that basis the FGG
ambiguity was called "eliminated, not hidden".
*Why it is wrong:* the constraint $2\sqrt{f_{\rm in}f_{\rm out}}\cosh(GP_u)=\pm(f_{\rm in}+f_{\rm out}-2\kappa^2u)$ is even in $\kappa$. Also $P$ and
$\varepsilon_{\rm in}\varepsilon_{\rm out}$ are invariant under $\kappa\to-\kappa$ ($A_{\rm in}\to-A_{\rm out}$, $A_{\rm out}\to-A_{\rm in}$). So the sign of the right-hand side
does not encode the tension. In fact the principal momentum satisfies the constraint with the *actual local* $\varepsilon_{\rm in}\varepsilon_{\rm out}$: it is
the real-angle (artanh) momentum. The analytic momentum satisfies the sign-free constraint: it is the complex-angle (Sorkin)
momentum (Theorem III-2).
*Lesson:* the two FGG prescriptions are two conventions for the boost angle of a flipped normal. Neither is a different
physical sector. What excludes the principal prescription is Theorem III-5 (no analytic solution across the flip) together
with $B<0$. The ambiguity is relocated to the convention and the sector normalisation, not eliminated (PHASE_III Sec. 9).
*Related slip, caught by the numerics:* a zero-net double winding was first assumed to need opposite orientations around the
black-hole and de Sitter horizons. The cancelling orientation is model dependent ($\sigma_h$ in Theorem III-1).
**Type:** interpretational error in a draft, corrected before release.

### F11. Evolving phi and ln alpha with their raw equations at the centre (Phase IV)
*Tried:* the diamond scheme for $\phi_{uv}=-(r_u\phi_v+r_v\phi_u)/r-\tfrac{\alpha^2}4U'$ and $(\ln\alpha)_{uv}=(r_ur_v+\alpha^2/4)/r^2-4\pi G\phi_u\phi_v$, with a
reflection-symmetric extrapolation to the centre.
*Failure:* (i) the $1/r^2$ term magnifies the even-parity truncation error of $\alpha$ by about $2/n^2$ per step at $n$ points from the centre.
The run blew up within three rays ($\ln\alpha$ error 360 against exact de Sitter). (ii) Even with $\alpha$ regularised, a massless pulse blew up
to $|\phi|\sim10^4$ through the $1/r$ term of the wave equation.
*Fix:* evolve $\psi=r\phi$ ($\psi_{uv}=\phi r_{uv}-\tfrac{\alpha^2}4rU'$, $\psi=0$ at the centre, exact in flat space); impose regularity $\alpha=2r_v$ at the
centre; replace $r_ur_v+\alpha^2/4$ by its regular limit $(\alpha^2/4)(8\pi G/3)\rho r^2$ within two points of the centre.
**Type:** numerical, resolved.

### F12. Excision by small r only (Phase IV)
*Tried:* excising points only when $r<2h$ in a trapped region.
*Failure:* at under-resolved settings ($h=w/16$, $w=0.005$) one ray next to the singularity jumped to $r=559$, $\phi=-14v$ instead of $r\to0$.
That produced 171 spurious "anti-trapped" points in a laboratory bubble, contradicting Proposition IV-1.
*Fix:* also excise where $r$ would increase along an already trapped outgoing ray (forbidden by Raychaudhuri under the NEC). The
spurious points disappear. At $h\le w/24$ they never occurred.
*Lesson:* Proposition IV-1 doubles as a numerical diagnostic. **Type:** numerical, resolved.


### F13. Absolute quadrature tolerances for the exponent (Phase VI)
*Tried:* `quad(..., epsabs=1e-13)` and angle differences $\Theta_{\rm out}-\Theta_{\rm in}$ (or $\arccos C$) for $B$.
*Failure:* at low vacuum scales and Planckian masses the dimensionless exponent is about $10^{-15}$ and $\theta\sim10^{-8}$. Then $1-C\sim10^{-16}$, so the angles lose
every digit and the absolute tolerance swamps the integral. One scenario showed $B_{\rm sq}>B_{\rm an}$, which violates $\arccos|C|\le\arccos C$.
*Fix:* relative tolerances, and the cancellation-free angle $\theta={\rm atan2}(\kappa R\sqrt{V/f_{\rm in}f_{\rm out}},\,C)$ (`tunnelling.theta_stable`, from the proof of
Theorem III-4). Stage One and Phase II–III results are unchanged (all tests pass). Their exponents were about 5000× above the old tolerance.
**Type:** numerical, resolved.


### F14. "Squaring removes the sign of the extrinsic curvature" (Phase III heuristic, corrected in the roadmap round)
*Claim:* Phase III motivated the sector-blind (squared) quantization by saying that a Wheeler–DeWitt constraint, quadratic in the momenta,
removes exactly the sign of $K_{\theta\theta}=A/R$.
*Failure:* on any slice through the wall in the static domain, $(R'/\Lambda,\ GP_\Lambda/R)$ is the boosted $(A,\dot R)$ with invariant $f>0$. So
$\mathrm{sign}(R'/\Lambda)=\mathrm{sign}(A)$ is boost-invariant **configuration** data, while the sign of $P_\Lambda$ (time orientation) is what squaring removes
(ROADMAP Prop. Q9.1). The squared reduced constraint is therefore not the image of the squared bulk constraint. The real content of
$B_{\rm an}$ vs $B_{\rm sq}$ is which conserved Casoratian charge is probability (ROADMAP Theorem Q9.3).
**Type:** heuristic error, corrected; no Phase III theorem is affected.

### F15. A seed black hole as a catalyst for daughter creation (Stage 4)
*Tried:* a false-vacuum bubble around a black hole of mass $m$ (interior Schwarzschild–de Sitter), by analogy with black-hole catalysis of vacuum decay.
*Result:* $V_{\rm seed}=V_{M-m}-2Gm/R$. The exponent falls only to $0.78B$ at $m=0.8M$ ($\chi=\kappa=1$, $M=0.3M_{\rm cr}$), and for $m\ge0.9M$ the initial wall leaves
regime A. More importantly, with an eternal seed the interior side of the wall continues through the seed's bridge into a static SdS patch whose
time-symmetric slice reaches the cosmological bifurcation sphere. That is a de Sitter equator: **the initial data already contain the daughter.**
A collapse-formed seed has no such patch and falls under Prop. IV-1.
**Type:** dead end with an instructive reason (`alternatives.SeededModel`).

### F16. The finite-time "daughter edge" (Phase V Theorem V-1(ii), corrected in the roadmap round)
*Claim:* the daughter channel ends at the interior de Sitter horizon at finite $t_\infty$, where $P\to\mp\infty$, so a quantum boundary condition is required.
*Failure:* the divergence is $P\simeq\pm(R/2G)\ln|f_{\rm in}|$, a property of the static-frame momentum coordinate. After the canonical transformation
$\tilde P=P\mp(R/2G)\ln|f_{\rm in}|$, the constraint and the true Hamiltonian are analytic across $f_{\rm in}=0$. $\tilde P$ is continuous and the flow is complete
($t_\infty\to-\infty$ as $R\to\infty$). No boundary condition is needed (ROADMAP Theorem Q10.1). The daughter still needs its own clock, but for a different
reason: the asymptotic clock runs backwards along every daughter (Theorem Q10.2).
**Type:** interpretational error, corrected.

### F17. "Time-reversed daughter", "Klein-type channel" (roadmap round, Corollary Q10.3 and D7; corrected in round 3)
*Claim:* because $dt_\infty/d\tau=A_{\rm out}/f_{\rm out}<0$ on every daughter, unitary boundary-clock evolution produces the time-reversed daughter. The physically
expanding daughter would then be a negative-flux, Klein-type channel, created as in a Bogoliubov process with $|\beta|^2\simeq e^{-B}$.
*Failure* (raised by the user: "reversing exterior Killing-time orientation is not by itself proof of backwards physical time"): $t_\infty$ is not a global
time function. Its Killing field is past-directed in region III, where every daughter wall lies.
* Both Kruskal null coordinates increase along the whole history on both sides of the wall (ROADMAP_ROUND3 Proposition A5, `daughter.global_time_orientation`).
* The clock-free conserved charge of the expanding double-flip daughter is positive, the same sign as the incident laboratory wave (A5,
  `selection.charge_signs`).

The $t_\infty$-flux relations remain correct as bookkeeping.
**Type:** interpretational error (a clock artefact read as physics), corrected.

### F18. Double-precision quadrature at the A/B pinch (round 3)
*Symptom:* spurious non-monotone steps of $|\theta_D|$ near $E_{AB}$, $10^{-3}$-level errors in $J$ for $E/E_{AB}-1\lesssim10^{-5}$, and a spurious
integration-by-parts certificate at $1.0001E_{AB}$ (margin at the LP cap, from undersampling a spike of width $O(\epsilon^2)$).
*Cause:* $E-g$ and $f_{\rm out}$ were evaluated as differences of $O(1)$ numbers where both are $O(\epsilon^2)$.
*Fix:* integrate in $y=z-z_1$ with cancellation-free expressions and geometric breakpoints on the spike width. Cross-checked against 30-digit mpmath.
**Type:** numerical pitfall, fixed before any conclusion was drawn from it.

### F19. A ghost-sign scalar pulse as a universe seed (Round 6)
*Tried:* the Phase IV double-null solver with $G\to-G$ in the matter couplings (massless ghost scalar, zero potential). Gaussian pulses of amplitude up
to 2 (Planck units), widths 0.15–0.3, imploded from flat data.
*Result:* no anti-trapped point before the evolution ended near the centre ($\max r_u<0$).
*Lesson:* the budget G1 is necessary, not sufficient. NEC violation must be coherent and focused on the ingoing rays through one sphere. Ghost
matter destabilises the vacuum in any case, so the line was not pursued.
**Type:** inconclusive numerical experiment.

### F20. Breathing seeds in a coarse scan (Round 6)
*Tried:* wall-plus-dust surface energies of opposite signs and negative-tension walls (96 parameter sets), looking for a stably oscillating orbit in
our region that could later reach a daughter.
*Result:* none found. Orbits either reach $\sigma=0$ (the wall becomes null) or are separated from the daughter by a barrier.
*Consistency:* Theorem G2 forbids a classical bridge crossing for any wall matter. The Guendelman–Sakai breathing orbits may lie outside the scanned range.
**Type:** negative scan, not a theorem.

### F21. "Route G must violate the GSL" (Round 6)
*Claim:* a sector that creates a daughter classically must violate the GSL, because laboratory creation is "restarting inflation".
*Error:* NEC violation is not GSL violation (Gao–Jafferis–Wall). Properly defined, the GSL forbids only eternal (de Sitter-ending)
daughters determined by the data, and disconnected ones (Theorem GSL-3). A transient daughter is not "restarting inflation" in Wall's sense.
*What survives:* a minimally coupled *classical* NEC-violating sector conflicts with black-hole thermodynamics as a theory, because its stress
can be dropped into a black hole.
**Type:** overstated conclusion, corrected (user feedback).

### F22. Vacuum energy switched on in a throat at rest (Round 6b)
*Tried:* KS dynamics of a homogeneous throat segment with $r=Q$, $\dot r=0$ and $\Lambda$ switched on, using only the $xx$ equation.
*Error:* the radial null constraint $\ddot r=(\dot a/a)\dot r$ ($T_{uu}=0$ for vacuum energy and radial magnetic field) forbids acceleration
from marginal homogeneous data; such data are not valid. Replaced by the RN–dS T-region criterion $C>C_N$ (Proposition GSL-8).
**Type:** incomplete use of the constraints, caught before any conclusion.

### F23. A cheap QNEC flip of a static throat as a daughter seed (Round 6b)
*Result:* a static throat in a fixed vacuum has $C=g(a)<C_N$ (local minimum vs maximum of $g$), so the flipped region recollapses. Growth
needs energy injection or a change of vacuum (GSL-8).
**Type:** structural negative result.

### F24. The electroweak corona as a natural false vacuum in a throat (Round 6b)
*Result:* self-consistently $\Lambda Q^2\le0.075$ (full restoration at $eB=m_H^2$), $\le0.226$ even with the overestimate $eB=m_W^2$ and full
$\rho_{\rm EW}$: below $1/4$, and a corona that is part of the static solution never buds anyway.
**Type:** near miss.

### F25. Extremality as a cure for the Cauchy-horizon route (Round 6b)
*Result:* the classical window survives as $Q/M\to1$ but its width is $\le r_+-r_-=2\kappa_-r_-^2$: weakening the blueshift closes the window,
and even the extremal Cauchy horizon is extendible only non-uniquely (Gajic–Luk).
**Type:** trade-off, not a cure.

### F26. "0.6–0.9% of the hole mass" and "about 400× cheaper" (Round 6c checkpoint)
*Claim:* the Nariai-throat tipping pulse is 0.6–0.9% of $M_N$, about 400 times cheaper than flat-space genesis.
*Error:* the quoted number was $\int\lvert N\rvert du$, the energy measured by the geodesic observer at the corner (a local, frame-dependent
energy), not a mass. The gauge-invariant Misner–Sharp jump is $\sqrt{f_0}$ times it ($f_0=2\delta/3$): $5\times10^{-7}M_N$ at $\delta=10^{-8}$,
scaling as $\sim10^{-2}\sqrt\delta M_N$. In local energy the threshold is $0.55\%M_N$ and the flat-space flip is 270 times larger.
**Type:** mislabelled quantity, corrected (ROADMAP_ROUND6C Part I.5).

### F27. QNEC-saturating CFT pulses as the tipping source (Round 6c)
*Tried:* conformal-map states of $c$ fields, first a step in $y'$ (which leaves a permanent positive flux $(c/48\pi)Y^2$), then optimal rise/slow-fall
profiles with the positive energy spread over up to $150\gg\tau$, $c$ up to $500=1.7r_N^2$.
*Result:* no tip. The weighted deficit of any CFT state is $\le c/48\pi\tau$ (Proposition 6c-10), so tipping needs $c\ge0.89r_N^2$ (Corollary 6c-11),
and the Bunch–Davies back-reaction of the same fields raises the threshold.
**Type:** theorem-backed negative result.

### F28. Vacuum polarisation (one-loop anti-evaporation) as the source (Round 6c)
*Tried:* Bunch–Davies back-reaction of $c=30$–$1000$ conformal fields on the near-Nariai configuration, without any pulse, domains up to $400$.
*Result:* no growth; the polarisation is stabilising (it raises the classical threshold above $W=0.071$ at $c=300$).
**Type:** negative numerical result.

### F29. Positive-energy perturbations (Round 6c)
*Result:* outgoing positive energy ($W<0$) at three amplitudes gives the trapped (collapsing) branch, never a daughter: Farhi–Guth, and Proposition
6c-13 (the minimal anti-trapped radius cannot decrease under the NEC, in any parent).
**Type:** consistency check.

### F30. A hole formed in true vacuum with a thin false-vacuum layer (Round 6c)
*Tried:* true vacuum deeper than $u_{\rm inner}=6$–$60$ along the ingoing ray (the black-hole horizon is crossed at $u_{\rm hole}=34.6$).
*Result:* no growth, even at four times the pulse; with $u_{\rm inner}\ge100$ the daughter grows and is hidden. The false vacuum must continue inside
the horizon.
**Type:** requirement found (Proposition 6c-15).

### F31. Non-minimal coupling as the classical NEC violator (Round 6c)
*Result:* $\xi\phi^2R$ is conformally equivalent to GR with NEC-satisfying matter; Farhi–Guth holds in the Einstein frame, so a Jordan-frame
daughter with $\phi\to$ const is impossible unless $\Omega^2=1-8\pi G\xi\phi^2\to0$ (strong coupling).
**Type:** structural negative result.

### F32. Hiding with an exterior already trapped on the data (Round 6c)
*Problem:* the heavy wall (vphi = 0.005, ε = 2) focused the outgoing data ray (Δr_v = −0.017 against r_v = +4e-5), so the true vacuum beyond
the wall was trapped on the data surface, with r decreasing outward into the parent's singularity: a configuration beyond the Einstein–Rosen
bridge, with no asymptotic parent on the data.
*Fix:* light walls (ε = 4, vphi = 1.5e-3 at δ = 1e-4) keep the data normal and bridge-free; the daughter still grows and hides (Round 6d).
**Type:** defect in the Round 6c set-up, corrected.

### F33. Re-tuning the corner radius against Bunch–Davies back-reaction (Round 6d)
*Tried:* moving the corner to the semiclassical static point r_N(1 + k/2r_N²) and around it (η from −0.026 to +0.02) at c = 300.
*Result:* the threshold stays between W = 0.07 and 0.14. The back-reaction is not a shift of the Nariai point.
**Type:** negative.

### F34. Method of lines without u-substeps (Round 6d)
*Result:* explicit u-steps of the Volterra operator (norm ~m²L_v/4 ≈ 10³ for the scalar's mass m² = 21 on L_v = 200) are unstable; 64
substeps per output line fix it.
**Type:** numerical.

### F35. Regular-centre collapse in the r = v/2 cone gauge (Round 6d)
*Result:* the lapse jumps by e³ across a massive shell and the stiff scalar becomes unstable; the α = 1 gauge (profiles in areal radius,
the same physical data) fixes it, but the near-horizon band is then too compressed to follow.
**Type:** numerical.

### F36. Emitting the push inside the false-vacuum ball (Round 6d)
*Result:* the complete route closes (hidden daughter), but the push sits where f ≈ 0.75 and its invariant energy is ~1.5 M_N: the Nariai
discount needs the NEC source in the band.
**Type:** partial success.

### F37. Condensate data at rest in the static frame (Round 6e)
*Result:* near both horizons the static observers are boosted without limit, and the condensate's null fluxes blow up (W_v = 10³–10⁶).
The physical frame is the corner's geodesic frame.
**Type:** wrong set-up.

### F38. Splitting the tilt into "charge source only" or "vacuum drift only" (Round 6e)
*Result:* each half violates energy–momentum conservation. The relative constraint residual is 0.04–0.39 and does not converge; the
conserving pair gives 2e-6 → 2e-7. Only the combination is physical.
**Type:** inconsistent model.

### F39. Exits that do not work (Round 6e)
*Result:* three variants fail.
* A fixed hilltop false vacuum drains the band before it tips.
* An exactly stationary false vacuum never exits without a seed.
* The unclipped tilted quartic can roll into a second positive-energy minimum.

The working exit is driven by the condensate clock (Proposition-N 6e-7).
**Type:** dead end / artefact.

### F40. Coherent classical lift of the ball, without dissipation (Round 6e)
*Result:* an ingoing massive-field packet crosses the barrier only transiently.
* Late false-basin points cover ≤ 10% of the interior.
* With full gravity the packets collapse to black holes (2m/r up to 0.9).
* Capture into a metastable basin needs a sink for the kinetic energy.

**Type:** dead end.

### F41. Confined-radiation fireball (Round 6e)
*Result:* a normal closed ball turns around with ρ_rad ≥ ρ_f at a < l_H/√2 (Proposition 6e-9), radius ≤ 0.42 l_H. It cannot be the
route's empty ball of radius 0.6 l_H.
**Type:** no-go (homogeneous model).

### F42. Explicit stiff mass term in the diamond scheme (Round 6e)
*Result:* the explicit potential force is unstable wherever m h e^s/2 > 1 (s ≳ 2 inside the parent hole).
* A field 10⁻⁶ v off the false vacuum grows into a spurious true-vacuum intrusion. That removed the daughter in every seeded run, and it
  contaminated the first transient-daughter run (withdrawn).
* Fixed in `condensate.py` with an implicit (Newton) potential step. Seeded and unseeded daughters are now identical.
* Results with the field exactly at the minimum (Rounds 6c and 6d, and the unseeded Round 6e runs) are unaffected, because that state is an
  exact solution.
  * Residual effect: the implicit step is stable but does not damp grid-scale oscillations. Where the light wall is crushed against the
    forming parent singularity (s ≈ 4–6), ϕ shows a checkerboard (≈ +3v alternating with −v) that does not resolve as h → 0 (it would
    need h ≲ 10⁻³). Hiding is therefore judged by points where the field actually sits at the true vacuum (|ϕ/v − 1| < 0.5): there are
    none at h = 0.08 and 0.04. The raw ϕ > 0 count is kept in the data as `phi_positive_flags`.

**Type:** numerical (fixed).

### F43. Thermal creation of the false-vacuum ball (Round 6f; withdraws Proposition-N 6e-10)
*Result:* heating restores the false vacuum only above E_min = 1.74–4.7 ρ_f (stall equilibrium; R118). The largest cone-normal hot ball
is about the size the route needs (17.5–18.25 against 17.85). In continuous GR evolutions (plasma, field and gravity, from a regular
centre) the false-vacuum region shrinks or is trapped before the radiation can leave (R119).
* Proposition-N 6e-10 used the high-temperature free energy at m/T ≈ 15, where it is invalid. Its "free-streaming emptying" was not a GR
  evolution.

**Type:** no-go (model, O(1) margin, robust over ε = 1–5.2, Γ/m = 0.03–3, f_H = 0.95–0.999, x_c = 5–30).

### F44. Numerical problems in the formation solver (Round 6f)
* **Box-scheme checkerboard.** The temperature alternated with the floor, η reached ±12 and ϕ ran to −6v. Fixed with Preissmann
  weighting θ = 0.6 and T_floor = 0.02 T_c. The plasma sector then converges at first order.
* **An invalid homogeneous test** (constant data on null rays). Replaced by ODE data on the rays.
* **The PDE–ODE mismatch** came from clip = 1 against the unclipped ODE. With clip = 0 they agree at second order.
* **Centre-stage axis blow-up** (s → −15 by row 117). The axis clock was NaN during the correction step, which silently skipped it. Fixed
  by setting the clock first.
* **Residual near-axis temperature error** (2.4% within 1–2 cells). Accepted.

**Type:** numerical (fixed).

### F45. Corner-gauge grid inside an inflating daughter (Round 6f)
*Result:* after the daughter's false vacuum passes the spinodal, h²e^{2s}|U″|/16 reaches ~4×10³. The implicit step cannot follow the
tachyonic roll, so rows break off at T ≈ 250–270, before the exit.
* A homogeneous-FRW hand-off failed: surfaces of constant condensate clock are not homogeneous slices (charge spread ~100%).
* Exit and recollapse remain established only in the band geometry (R114).

**Type:** numerical limitation (open, P-R6f-b).

### F46. Misleading l-mode diagnostics (Round 6f)
*Result:* two mistakes, both corrected (R123).
* Maximising the amplification over the whole domain measures the crushing near r → 0 at the singularities (up to 10⁷, equally in the
  controls).
* Using the Misner–Sharp mass in the Regge–Wheeler term adds a spurious tachyonic −Λ in vacuum-energy regions.

**Type:** artefact (fixed).

### F47. Stability maps from white-noise rms amplification (Round 6f)
*Result:* stable fast modes grow algebraically (~ω₀^{1/4}) in dS₂, which contaminates the rms measure, and the per-mode integration was
far too slow at large μ. Replaced by the WKB growth exponent with k and l modes, cross-checked against the full ODE (R124).

**Type:** method (replaced).

### F48. The route with ε = 4 and a seeded exit (Round 6f)
*Result:* at λ_f = 10⁻⁶ the false vacuum reaches its spinodal (T = 44) before the throat forms (T ≈ 125). It exits near the centre and no
daughter forms. The exit must follow throat formation (R127). Unseeded runs "work" only because ϕ = −v exactly is an unstable
equilibrium.

**Type:** physics constraint.

### F49. Conformal FRW gauge for an inflating patch (Round 6g)
*Result:* with s = ln a(η) of the radiation background on the initial cone, the null coordinates do not follow the patch's own
exponential expansion. Inside a super-Hubble patch, r tripled within one cell and the explicit predictor diverged (r ~ 10²² at the next
cell).
* The affine cone gauge (s = 0 on the cone, comoving fluid at η_f = ln a) and h = 0.05 fix this for the first 10²–10⁴ of growth.
* The inflating interior still becomes stiff later (as in F45), which ends the computed domain.

**Type:** numerical (gauge; partly fixed).

### F50. A machine humanity could build (Round 6g)
*Result:* every candidate fails Theorem 6g-2, including invented ones:
* colliders, lasers, explosions;
* an ideal focus of one year of world energy in 100 TeV quanta;
* a parsec-long collider;
* a Dyson sphere.

The required power (~10⁵² W) or force (~10⁴³ N) is a fixed fraction of c⁵/G and c⁴/G, whatever the vacuum scale. Only physics
beyond 4D GR would change that (gravity strong at the TeV scale). It is excluded below 10.1 TeV, and would itself be a set of new
problems.

**Type:** no-go (theorem).

### F51. "No machine of any kind" (Round 6g, withdrawn in Round 6h after review)
*Result:* the claim does not follow from Proposition 6g-2, for three reasons.
* Rapid external delivery and static holding are two preparation strategies, not every possible history.
* The Planck force and luminosity are conjectures with counterexamples in GR.
* The computed stress (1.3×10⁴³ N) is below c⁴/4G, and an assembly time of 1.5 light-crossings puts the power below c⁵/4G.

The rigorous strategy-independent statement is Theorem 6h-1 (Penrose / Farhi–Guth). It applies to machines of NEC-respecting matter.

**Type:** overclaim (corrected).

### F52. Curvature activation as the trigger (Round 6i, route 1)
*Result:* a field the Gauss–Bonnet coupling destabilises in a vacuum-energy throat (GB > 0) rolls so as to raise F = λ²f, and the weighted
Gauss–Bonnet flux is then pro-collapse (Proposition 6i-1). Scalarised Nariai throats collapse or lose predictivity (degenerate kinetic matrix);
the reversed sign with a positive-energy minimum gave 0 daughters in 144 runs. Laboratories fall short of activation by 19 orders.

**Type:** no-go (proposition + numerical).

### F53. Unbounded-potential "daughters" in the Gauss–Bonnet scan (Round 6i, discarded in the same round)
*Result:* runs with V = 1 − μ²φ²/2 reached b = 3 r_N and were first counted as daughters. The potential had already turned negative
(t ≈ 4–8 r_N) and the solution ran away (φ ~ 10²⁷, constraint lost). With a positive-energy minimum no daughter occurs.

**Type:** numerical artefact (caught).

### F54. Monopole to daughter at the monopole's own mass (Round 6i, route 2)
*Result:* the thin-shell exponent is only ≈ 650/λ, but M ≪ Q: the exterior has no bridge, and the β_out < 0 branch is a universe
disconnected from ours, which Gauss's law forbids as a process in our universe. The connected channel needs a magnetic black hole
(M > Q ≈ 7 m_P) and has B ≥ 3×10³. Negative dust lowers B by at most 20%.

**Type:** no-go (assessment + numerical).

### F55. Breathing bubbles, gauged ghosts, resonance, collisions, membranes (Round 6i, routes 3–7)
*Result:*
* Breathing: the hidden case needs NEC-violating dust and changes side only through zero net surface energy at light speed.
* Gauged universal condensate: the window is empty by ≥ 5 orders (screening).
* Resonance: never removes the last barrier.
* Collision pockets: cannot reach the hidden branch in a normal parent.
* Membranes: Farhi–Guth walls; B ≈ B_CdL.

All are instances of Theorem 6i-11 or of "no free lunch".

**Type:** no-go (propositions).

### F56. The Round 6c partner-vacuum formula (found in Round 6j)
*Result:* `nariai.partner_vacuum` used ρ_f = 1/(24πM²). The Nariai condition 9M²Λ = 1 gives ρ_f = Λ/8π = 1/(72πM²), the relation
`accessibility.py` always used (ρ_f = 1/8πr_N², M_N = r_N/3). The partner scale drops by 3^{1/4}: 0.329 GeV (M_☉/M)^{1/2}, a GUT vacuum
pairs with 2.2 g (not 3.7 g), an electroweak one with 3.6 (not 6) Earth masses. No conclusion depended on the factor.

**Type:** error (corrected).

### F57. A field-dependent gravitational coupling as a shortcut (Round 6j, #3)
*Result:* the F^{3/4} reduction of the requirement is paid for by the core itself. Its canonical excursion Δχ = √(3/16π)|ln F| M_P has
gradient energy ≥ 2πLΔχ², so the core is horizon-class for F < 0.32, long before the needed F ≈ 10⁻⁵³. f(R) is the same theory.

**Type:** no-go (proposition).

### F58. Environment-prepared states (Round 6j, #6–#9)
*Result:*
* Accretion cannot reach the partner density: Eddington flows are 10¹⁹ short, merger debris 0.25–1.8%.
* The most compact bag star holds a sub-Hubble vacuum and collapses to a black hole.
* Evaporating-hole catalysis is Planck-scale only.
* Critical cosmological patches existed only in the early universe.
* Small supplied energies need precursors within ~10⁻⁸ of threshold that are stable and older than the universe: none exists.

**Type:** no-go (propositions + numerical).

### F59. A connected quantum channel without a throat (Round 6j, #10–#11)
*Result:* in Euclidean axion wormholes the cost of staying connected is ΔS = 3π³a₀² (8πG = 1), proportional to the throat area. The
throat-free limit is the disconnected no-boundary universe. Together with Q* (Round 6i): suppression disappears exactly when the
connection does.

**Type:** no-go (numerical, exact fit).

### F60. A four-dimensional ECSK daughter at collider scale; the Planet-Nine/electroweak coincidence (Round 6k)
*Result:* ECSK puts a bounce inside every hole above M_min ≈ 0.15 M_P²/m (4.3×10¹⁰ kg for nucleons). A heavier fermion lowers it only as 1/m
(top quark: ~2×10⁸ kg), so in four dimensions no collider-scale hole holds an ECSK bounce. A collider route needs extra dimensions (P1).
Separately, a 5–10 Earth-mass primordial hole is the Nariai partner of a 60–85 GeV vacuum, close to the electroweak V(0)^{1/4} = 104 GeV. But
the symmetric phase is not metastable at zero temperature, and restoring it thermally costs ~200× its vacuum energy (Round 6f). The
closeness is a coincidence, not a route.

**Type:** scaling bound + assessment.

### F61. A homogeneous Einstein-Cartan ball re-expanding inside a vacuum exterior (Round 6k picture, corrected in Round 6l)
*Result:* the spin pressure does work, so the homogeneous ball's Misner–Sharp mass falls from (sin³χ/2)A to (sin³χ/2)a_b by the bounce:
99 % of it for S = 10⁻⁶. That is incompatible with a vacuum exterior of fixed mass. Kinematically, any bounce sphere is untrapped
(Proposition 6l-3). So the ECSK bounce must be inhomogeneous: an untrapped core inside an inner apparent horizon, with the outer layers
still falling. The Round 6k statement that ECSK interiors re-expand into a new universe survives only in this form.

**Type:** consistency correction (proposition + numerical).

### F62. Inflating daughters, and extra daughters by adding matter, in the unified model (Round 6l)
*Result:* in model U, the one consistent action found that contains both P1 and P2:
* vacuum collapse leaves a static de Sitter core;
* matter cores accelerate only for N ≈ ((D−3)/(D−1)) ln(a_max/l) e-folds and then recollapse;
* lasting inflation needs a vacuum as dense as the hole.

Later infall meets the bounced matter inside the same core and, after 7–30 l of delay, saturates the curvature bound: it does not quietly
create a second daughter. The daughter that does exist is a new asymptotic region beyond a white hole, not an inflating cosmos.

**Type:** no-go within model U (propositions + numerical); open beyond thin shells (P-R6l-a).

### F63. Crossing null thin shells at the inner horizon of model U (Round 6l, Proposition-N 6l-6b; withdrawn in Round 6n)
*Result:* the junction conditions of pure-gravity quasi-topological theories admit no null thin shells in these regular black holes (Di
Filippo, Kubizňák & Srinivasan, PRD 114, 044061 (2026)). The Dray–'t Hooft–Redmount argument we used assumed such shells, so its
derivation does not hold in model U. Round 6n treats the same encounter, later matter meeting the bounced core at the inner horizon, with
continuous timelike matter instead: the curvature stays bounded, matter piles up at R ≈ l, and the daughter's white-hole region begins to
form.

**Type:** withdrawn derivation (literature).

### F64. Comparing M_D across conventions (Rounds 6k-6l; corrected in Round 8)
*Result:* our formulas use the Dimopoulos–Landsberg horizon radius, i.e. G_D = 1/M_D^{n+2}. Experimental bounds use the PDG convention,
G_D = (2π)^n/(8πM_D^{n+2}), so M_D(DL) = 0.89–0.38 M_D(PDG) for n = 2–6. Rounds 6k–6l quoted "M_D ≳ 5–11 TeV allowed" without converting.
Statements made in terms of M_ext (TeV) are unaffected. The allowed range of M_D(DL) is wider: M_D(DL) > 2.0–8.8 TeV.

**Type:** convention error (corrected).

### F65. Flat-efficiency LHC bound and crude FCC-hh reach (Rounds 8-9; corrected in Round 10)
*Result:* Round 8 bounded model U with a single inclusive ATLAS limit (σ × A × ε < 0.3 fb) and an assumed 50 % efficiency. That gave
M_ext ≳ 8–9 TeV (full trapping) or 4.7–5.5 TeV (head-on trapping). Round 9 used an order-of-magnitude Z(→νν) background with a
leading-log ISR spectrum.

*What Round 10 found:*
* The full chain (Pythia shower, detector, pile-up, exact selection) and a 13-bin likelihood give a weaker LHC exclusion,
  M_ext ≲ 3.5–4.2 TeV (head-on) and ≲ 6.7–7.2 TeV (full trapping).
* The realistic 100 TeV discovery reach is ≈ 32.5 TeV, not ≈ 40 TeV.
* The reason is physical: a heavy invisible remnant only recoils against initial-state radiation, so its MET (≈ 0.5–1 TeV) does not
  grow with M_ext. Near the top of the window the signal sits on the Z(→νν)+jet background at low MET.

The window stays open and testable. Its numbers move.

**Type:** method too crude (corrected with full simulation).

### F66. Double-null evolution in (r, σ) near a regular Cauchy horizon (Round 10)
*Result:* in the standard (r, σ, φ) variables, truncation errors near the inner horizon are blueshifted, exactly as physical
perturbations are. ψ is derived from r_u r_v e^{−2σ}, so it overshoots the tower's pole ψ = 1/α, and even the static hole drifts
for v ≳ 20 l. Evolving the mass function instead, (r, μ, φ), fails at apparent horizons because of the division by f.

*What works:* a hybrid. μ is carried with the first-law fluxes and ψ = h⁻¹(μ/r^{D−1}) feeds the metric equations. The static hole
is then exact to 10⁻⁷, and constraint violation converges.

**Type:** numerical formulation (replaced).

### F67. Profile-likelihood fits on raw counts with linear nuisances (Round 10, caught before reporting)
*Result:* with up to 5×10⁷ events per bin, the raw −ln L is ~10⁷. L-BFGS-B then stopped at non-minima: a "free" fit came out worse
than the constrained one, giving Z = 11.5 where the correct value is 4.0. Linear normalisation nuisances (1 + σθ) also let a signal be
switched off at a bounded cost.

*Fix:*
* use the saturated-model-subtracted likelihood;
* use log-normal normalisation nuisances;
* fit with Fisher scoring (damped Gauss–Newton) on the analytic Jacobian, gradients checked to 2×10⁻⁷;
* validate against the analytic Asimov formula (agreement to 10⁻¹⁰) and against pseudo-experiments.

**Type:** numerical method (fixed).
