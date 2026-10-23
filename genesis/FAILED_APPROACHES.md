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
