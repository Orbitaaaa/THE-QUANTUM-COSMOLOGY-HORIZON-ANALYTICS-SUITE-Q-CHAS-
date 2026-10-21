# Canonical structure and quantization

Status: classical reduced phase space **derived and verified**; quantum operator formulations
**formulated**, with their semiclassical (eikonal) limit verified; full quantum solutions are
**Phase III/V work**. Code: `CODE/babyuniverse/quantum.py`, `symbolic.reduced_hamiltonian`;
data `DATA/quantum_checks.json`.

---

## 1. Degrees of freedom

Spherical symmetry plus Birkhoff leaves, after the constraints are solved:

* the wall radius $R$ and its conjugate momentum $P$ (one physical degree of freedom);
* the exterior mass $M$, conjugate to the asymptotic Killing time $t_\infty$ (the boundary Hamiltonian);
* no local degrees of freedom in either vacuum region.

The reduced action is $S=\int(P\,dR-M\,dt_\infty)$ plus joint/horizon terms. Physical observables are
functions of $(R,P)$ on the constraint surface, with $t_\infty$ as relational clock: the time read by an
asymptotic observer (a boundary clock). The radius $R$ is a gauge-invariant **relational observable**:
the area radius of the wall sphere, evaluated when the boundary clock reads $t_\infty$.

## 2. The Hamiltonian is multivalued: orientation sectors

From $P=\frac RG[\operatorname{artanh}(\dot R/A_{\rm out})-\operatorname{artanh}(\dot R/A_{\rm in})]$, write $A=\varepsilon\sqrt f\cosh\tilde\alpha$,
$\dot R=\sqrt f\sinh\tilde\alpha$ in static regions. One finds $A_{\rm in}A_{\rm out}-\dot R^2=\varepsilon_{\rm in}\varepsilon_{\rm out}\sqrt{f_{\rm in}f_{\rm out}}\cosh(GP/R)$.
Since also $A_{\rm in}A_{\rm out}-\dot R^2=\tfrac12(f_{\rm in}+f_{\rm out}-\kappa^2R^2)$, the constraint in every static sector is
$$
\boxed{\;2\sqrt{f_{\rm in}f_{\rm out}}\,\cosh\frac{GP}R=\varepsilon_{\rm in}\varepsilon_{\rm out}\,\big(f_{\rm in}+f_{\rm out}-\kappa^2R^2\big)\;}
$$
Because $\cosh>0$, the sign of $f_{\rm in}+f_{\rm out}-\kappa^2R^2$ fixes the product $\varepsilon_{\rm in}\varepsilon_{\rm out}$. In the $(+,+)$ sector the
constraint can be solved for $M$:
$$
M=H(R,P)=\frac R{2G}\Big[1-\Big(\sqrt{f_{\rm in}}\cosh q-\sqrt{f_{\rm in}\sinh^2q+\kappa^2R^2}\Big)^2\Big],\qquad q=\frac{GP}R\quad(\eta=0).
$$
Its $G\to0$ limit is $\sqrt{P^2+m^2}$ with $m=4\pi\sigma R^2$, a relativistic membrane (`symbolic.flat_space_limit`).

The four sign sectors are four sheets of one constraint surface. The **squared constraint**
$$
\mathcal C^2\equiv4f_{\rm in}f_{\rm out}\cosh^2\frac{GP}R-\big(f_{\rm in}+f_{\rm out}-\kappa^2R^2\big)^2=0
$$
contains all of them. Under the barrier $\cosh\to\cos(GP_E/R)$, and the $\varepsilon_{\rm in}\varepsilon_{\rm out}=\pm1$ sectors join
analytically where $\Theta_{\rm out}-\Theta_{\rm in}$ crosses $\pi/2$. The Euclidean regime is where the sectors communicate.
This is the starting point of the "branch-aware" quantization of NEW_APPROACHES.md (N1).

**Verification.** $dR/dt_\infty=\partial H/\partial P$ holds to $6\times10^{-10}$ in sectors $(+,+)$, $(+,-)$, $(-,-)$ at 22 phase-space
points (`DATA/quantum_checks.json`). $H$ is inverted by a damped Newton iteration in $M$ using the exact
derivative $\partial P/\partial M=A_{\rm out}/\dot Rf_{\rm out}$. An earlier bracket search could not stay inside a sector on
narrow $M$-intervals; this is the "sign-sector failure" noted in Phase II, resolved here.

## 3. Area variables and imaginary shift operators

With $u=R^2/2$ and $P_u=P/R$, $(u,P_u)$ is canonical ($P\,dR=P_u\,du$) and the momentum enters only through
$\cosh(GP_u)$. Quantizing $P_u\to-i\hbar\,\partial_u$ gives
$$
\cosh(GP_u)\to\tfrac12(T_++T_-),\qquad (T_\pm\psi)(u)=\psi(u\mp i\hbar G).
$$
The Hamiltonian generates **translations of the wall area in the imaginary direction**: $\Delta A=8\pi\Delta u=\mp8\pi i\hbar G=\mp8\pi i\ell_P^2$.
The quantum constraint is a second-order difference equation with imaginary step. In the squared form,
with $\cosh^2=\tfrac14(T_+^2+2+T_-^2)$:
$$
f_{\rm in}(\hat u)f_{\rm out}(\hat u)\,\big[\psi(u-2i\hbar G)+2\psi(u)+\psi(u+2i\hbar G)\big]=\big(f_{\rm in}+f_{\rm out}-2\kappa^2u\big)^2\psi(u)
$$
(one ordering among several). Wavefunctions must be analytic in a strip $|{\rm Im}\,u|\le2\hbar G$, which is the
natural domain for questions of self-adjointness. The eikonal limit $\psi\sim e^{iS/\hbar}$ reproduces the
classical constraint, and under the barrier the Euclidean angle relation. The symbol vanishes on the
analytic-branch Euclidean orbit to $10^{-12}$ (`quantum.eikonal_check`).

**$M=0$: an exactly polynomial problem.** For $M=0$, $f_{\rm in}=1-2\chi^2u$ and $f_{\rm out}=1-2\eta^2u$ are linear in $u$,
and the squared constraint is quadratic in $u$:
$$
\cosh^2(GP_u)\,\big(1-2bu+cu^2\big)-(1-au)^2=0,\quad a=\chi^2+\eta^2+\kappa^2,\ b=\chi^2+\eta^2,\ c=4\chi^2\eta^2 .
$$
(The classical turning point $u_0=2\kappa^2/L^4=R_0^2/2$ follows from $a^2-c=L^4$.) In the momentum representation
$u\to i\hbar\,d/dk$ this becomes a **second-order linear ODE in $k$** with $\cosh^2(Gk)$ coefficients. That makes the
$M=0$ transition exactly solvable in principle, which is the Phase III benchmark for comparing the exact quantum
amplitude with $e^{-B/2\hbar}$ (T1).

## 4. Two inequivalent routes, and a dead end

**Route A: reduce, then quantize** (this file). Physical Hilbert space: functions of $u$ (or $R$) analytic
in the strip, with the asymptotic clock. The WKB exponent is $B$ of TUNNELLING.md.

**Route B: quantize, then reduce** (Dirac/Wheeler–DeWitt). Canonical variables: Kuchař-type
$(M,p_M)$ for each Schwarzschild region (Kuchař 1994), the wall $(R,P_R)$, and the Hamiltonian and
momentum constraints of spherically symmetric gravity, imposed as operators. The reduced constraint of
Route A is the semiclassical image of the Dirac constraints after the gauge $t=t_\infty$. The two routes
are expected to agree at leading WKB order and to differ in ordering, measure and prefactor. Making this
comparison concrete is OPEN_PROBLEMS.md P-Q1.

**Dead end: quantizing the effective potential.** $\dot R^2+V(R)=0$ arises from many Lagrangians with the same
orbits but different symplectic forms, e.g. $L=\tfrac12m(R)(\dot R^2-V)$. Their WKB exponents
$2\int m\sqrt V\,dR$ disagree with the gravitational one by factors between 1.06 and 13, with no fixed relation (`DATA/quantum_checks.json`):

| $\gamma$ | $E$ | $B_{\rm can}$ | $B_{\rm naive}$ ($m=4\pi\sigma R^2$) | $B_{\rm naive}$ ($m=1/G$) |
|---|---|---|---|---|
| −0.6 | 3.84 | 2.70 | 0.62 | 0.54 |
| 0.0 | 6.0 | 0.881 | 0.362 | 0.938 |
| 0.6 | 8.0 | 0.0361 | 0.0275 | 0.480 |

The effective potential determines the classical orbits, **not** the quantum theory. Only the gravitational
symplectic structure, $P$ above, is admissible (FAILED_APPROACHES.md F3).

## 5. Wheeler–DeWitt-type questions, answered or posed

| question | status |
|---|---|
| operator domain | functions analytic in $\lvert{\rm Im}\,u\rvert\le2\hbar G$ (Route A); to be made precise |
| ordering | open; the coefficient functions are not analytic at $u=0$ when $M>0$ ($f_{\rm out}\ni u^{-1/2}$) |
| self-adjointness | open; known results for shell difference operators (Berezin; Hájíček) give discrete mass spectra for bound shells |
| inner product | Route A: $L^2$ in $u$ with the boundary clock; Route B: group averaging / relational |
| probability | per attempt $e^{-B/\hbar}$ (fixed energy, decaying under the barrier); prefactor open |
| sectors / global geometry | the four orientation sectors are sheets of one constraint ($\mathcal C^2$). A superposition across sectors corresponds to different *global* geometries (wall in I vs III), not to different spatial topologies |

## 6. Boundary conditions

For the laboratory problem the boundary condition is **fixed by the initial state**: a packet localised on the
bounded branch decays into the barrier. No cosmological boundary proposal (no-boundary or tunnelling
wavefunction) needs to be chosen, unlike in creation from nothing. The final state is the outgoing (expanding)
WKB branch for $R>R_2$. At $M=0$ the configuration is a closed universe; there the choice between growing and
decaying branches is the Hartle–Hawking vs Vilenkin question, and T1 shows the decaying branch gives
$S_E+\text{Morse}$. The literature (Céspedes et al. 2023) finds that detailed balance favours the Hartle–Hawking
sign for nucleated spacetimes. Reconciling the two readings at $M\to0^+$ is P-Q3.

## 7. Connection to the crossed-product manuscript (to be examined, not assumed)

In Route A the asymptotic Killing time is the clock and $M$ its generator. This is the structure of a crossed
product of the exterior observable algebra by the boundary time translation. If the manuscript's Theorem T/L
(relational time from a Hamiltonian constraint) applies to the exterior algebra with $H=M$, it would supply the
physical Hilbert space of Route A. The daughter sector is invisible to that algebra after pinch-off. Whether the
construction survives the sector structure of $\mathcal C^2$ is NEW_APPROACHES.md N2. None of the manuscript's claims
are used in any result here.
