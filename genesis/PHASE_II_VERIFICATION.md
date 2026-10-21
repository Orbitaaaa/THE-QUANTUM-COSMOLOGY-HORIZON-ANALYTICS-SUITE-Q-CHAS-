# Phase II — Independent mathematical verification

**Exit condition (from the programme):** *a reproducible and independently checked mathematical foundation,
with unresolved assumptions documented.* Status: **met**, with the explicit residual gaps in Sec. 8.

Scripts: `CODE/scripts/phase2_verification.py` (independent 40-digit mpmath re-implementation, imports no
tunnelling code), `CODE/scripts/phase2_regimeB_scan.py`. Data: `DATA/phase2_high_precision.json`,
`DATA/phase2_regimeB_scan.json`. New tests in `CODE/tests/test_tunnelling.py`.

---

## 1. The sign-sector failure: resolved and classified

**Symptom (Stage One).** Numerical inversion of $H(R,P)$ failed in the $(\varepsilon_{\rm in},\varepsilon_{\rm out})=(+,-)$ sector.

**Diagnosis.** Two independent defects (FAILED_APPROACHES.md F2):

1. The orbit branch ${\rm sgn}\,\dot R$ was inferred from ${\rm sgn}\,P$. In the $\varepsilon_{\rm out}=-1$ sector
   $P=\frac RG[\operatorname{artanh}(\dot R/A_{\rm out})-\dots]$ has $A_{\rm out}<0$, so $P<0$ for $\dot R>0$ (the asymptotic
   clock runs backwards in region III). The inference was wrong.
2. At fixed $R$, $P(M)$ is defined only on the set of $M$ for which (i) $V(R;M)<0$, (ii) both $f>0$,
   (iii) the sector is unchanged, and (iv) $|\dot R|<|A_{\rm in}|,|A_{\rm out}|$. This set can be a narrow interval that global
   brackets miss.

**Fix.** Damped Newton iteration in $M$ using the exact Hamilton–Jacobi derivative $\partial P/\partial M=A_{\rm out}/(\dot Rf_{\rm out})$,
constrained to the sector and the branch. Verified in all three static sectors that occur, $(+,+),(+,-),(-,-)$:
$\partial H/\partial P=dR/dt_\infty$ to $6\times10^{-10}$ at 22 points.

**Classification (not a bug but physics).** Where (iv) fails, the wall is inside a horizon on at least one side.
There the artanh's are complex ($P$ acquires $\pm i\pi R/2G$) and the asymptotic Killing-time clock is singular.
This is the canonical face of the FGG slicing problem (TUNNELLING.md Sec. 2). The reduced Hamiltonian with the
asymptotic clock is a well-defined real function exactly on the static–static domain. The barrier always lies in
that domain ($f\ge V>0$).

## 2. Independent high-precision re-implementation

All key quantities were re-coded from the formulas in the documents with mpmath at 40 digits, sharing no code
with `babyuniverse.tunnelling`:

| quantity | method | result |
|---|---|---|
| $B$ | arccos formula (Sec. 3), tanh–sinh quadrature | — |
| Smarr identity, 9 cases incl. $M/M_{\rm cr}\sim10^{-6}$, $\gamma=\pm0.95$, SdS, $M=0$ | Stokes volume, wall area, Killing time, excised exterior volume | residual $\le3.4\times10^{-21}$ (8 cases); $2\times10^{-15}$ for $M/M_{\rm cr}\sim10^{-6}$ |
| $M=0$ closed form (5 cases incl. $\kappa/\chi=3.6\times10^{-3}$) | closed form vs quadrature | relative $\le1.5\times10^{-37}$ |
| regime-B critical angle (Sec. 4) | closed form vs quadrature at $E=g_{\min}(1+10^{-8})$ | agreement to $3\times10^{-9}$ (= the $10^{-8}$ offset) |

The double-precision production code agrees with these values to its own accuracy (Smarr $\le2.5\times10^{-8}$).

## 3. New lemma and an independent formula for the exponent

> **Lemma P (positivity).** On the analytic branch, $0<\Theta_{\rm out}-\Theta_{\rm in}<\pi$ throughout the barrier. Hence
> $P_E>0$ and $B>0$ for every transition in every regime.
>
> *Proof.* $\Theta\in(0,\pi)$ with $\cot\Theta=A/\sqrt V$, and $A_{\rm out}=A_{\rm in}-\kappa R<A_{\rm in}$. Since $\cot$ decreases on $(0,\pi)$,
> $\Theta_{\rm out}>\Theta_{\rm in}$. $\square$

> **Corollary (arccos formula).**
> $$B=\frac2G\int_{R_1}^{R_2}R\,\arccos\frac{f_{\rm in}+f_{\rm out}-\kappa^2R^2}{2\sqrt{f_{\rm in}f_{\rm out}}}\,dR,$$
> with the principal $\arccos\in[0,\pi]$. No angles and no branch logic are needed.

Agreement with the angle formula: $\le10^{-12}$ in double precision, and used as the 40-digit reference above. This
gives a second derivation of T0's positivity criterion: the analytic branch is positive *by theorem*, and the
principal branch is negative in 41/45 models.

**Exact form of the FGG ambiguity.** For every model tested,
$$B_{\rm analytic}-B_{\rm principal}=\frac{A(R_2)-A(R_f)}{4G}-\Big[\frac{A(R_2)-A(R_{f,\rm in})}{4G}\Big]_{\text{if }\varepsilon_{\rm in}\text{ flips}},$$
with $A=4\pi R^2$ and $R_f$, $R_{f,\rm in}$ the flip radii, to 10 digits. The two prescriptions differ by **horizon-type area terms
over $4G$**. This is the horizon monodromy of T0 made explicit (Phase III, PHASE_III_QUANTUM_TRANSITION.md).

## 4. NG4 (regime B): the critical limit proved

In the limit $M\to M_{\rm cr}$ the barrier shrinks to the maximum $z_*$ of $V$, and the half-bounce takes Euclidean time $\pi/\omega$
with $\omega^2=-V''/2$. Using $E=g(z_*)=3(z_*^6-1)/z_*^4$, $1-\gamma=(z_*^3+2)(z_*^3-1)/z_*^3$,
$f_{\rm out}(z_*)=(z_*^3-1)/3(z_*^3+1)$ and $g''(z_*)=6(z_*^6+2)/z_*^6$, everything collapses to

> **Theorem NG4-crit.** As $M\to M_{\rm cr}$ in regime B,
> $$\frac{|\theta_D|}{\pi}\;\longrightarrow\;F(t)=\frac{3\sqrt3\,(t+1)^2}{4\,(t+2)\sqrt{t^2+2}},\qquad t=z_*^3=\tfrac12\big(-\gamma+\sqrt{\gamma^2+8}\big)\in(1,2],$$
> and $F(t)>1$ for all $t>1$, i.e. for all $\gamma<1$. Hence $|\theta_D|>\pi$ near criticality and the half-instanton is
> not embedded.
>
> *Proof.* $F^2>1\iff27(t+1)^4-16(t+2)^2(t^2+2)>0$, and this polynomial equals $(t-1)(11t^3+55t^2+121t+101)$. $\square$

Near $\gamma=1$: $F=1+\delta/9+\delta^2/324+O(\delta^3)$, $\delta=1-\gamma$. So the margin vanishes only in the degenerate
limit $\kappa\to0$. The formula agrees with direct quadrature to $3\times10^{-9}$ (Sec. 2).

**Away from criticality** (numerical): the scan in `DATA/phase2_regimeB_scan.json` covers 39 values of
$\gamma\in[-0.95,0.95]$ × 10 values of $E$ spanning regime B. Results:
$\min|\theta_D|/\pi=1.0056$, attained at $\gamma=0.95$ at the critical end, where the closed form gives $F=1.0056$; and
**$|\theta_D|$ increases monotonically with $E$ for every $\gamma$**. Together with the theorem at the critical end,
$|\theta_D|>\pi$ throughout regime B. A proof of monotonicity in $E$ would make NG4 a full theorem.

## 5. Extension of the first law and Smarr identity to de Sitter parents

For $\eta>0$, the same implicit differentiation gives
$\partial P/\partial\Lambda_{\rm out}|_{R,M}=+\frac1{8\pi G}\frac{4\pi R^3}3\frac{dt_{\rm out}}{dR}$, the mirror image of the interior relation.
Euler homogeneity now includes $\Lambda_{\rm out}$:

> **Theorem T3′ (SdS exterior).**
> $$B=-\frac{\Lambda_{\rm in}}{8\pi G}{\rm Vol}^{\rm in}_E+\frac{\Lambda_{\rm out}}{8\pi G}W^{\rm out}_E-\frac\sigma2A^W_E-\frac M2\Delta t^W_E+\frac{\mathcal A}{4G},
> \qquad W^{\rm out}_E=4\pi\oint\frac{R^3}3\,dt_E^{\rm out},$$
> with $\partial\tilde B/\partial\Lambda_{\rm out}=-W^{\rm out}_E/8\pi G$.

$W^{\rm out}_E$ is the 4-volume of the region $r<R$ excised from the exterior static chart along the wall. The exterior vacuum
energy enters as the bulk action of the *missing* exterior relative to the background. Verified: Smarr to
$5\times10^{-10}$ (double precision) and $\le3\times10^{-21}$ (40 digits) for SdS models with $\eta/\chi$ up to 0.8. The
$\Lambda_{\rm out}$ first law holds to 7 digits.

## 6. The Euclidean winding: physical, coordinate effect, or artefact?

Answer: **physical and invariant, but representation-specific and irrelevant to the exponent.**

* *Invariance.* $\theta=2\pi\Delta t_E/\beta$ uses the Killing time normalised at infinity and the period fixed by smoothness
  at the horizon. It is not a coordinate effect.
* *Not numerical.* The winding is confirmed **independently of any time integration** by the first law:
  $\Delta t^W_E=-\partial\tilde B/\partial M$, computed by differentiating $B$ (built from $P_E$ alone), agrees with the direct
  integral to 8 digits (`DATA/first_law_checks.json`). The 40-digit code reproduces the integrals.
* *Origin.* $\beta=8\pi GM\to0$ as $M\to0$ while $\Delta t^W_E$ stays finite. The wall's Euclidean history spans many thermal
  periods of a tiny black hole.
* *Effect.* It contributes $-\tfrac M2\Delta t^W_E\to0$ to $B$ (Smarr). It is a property of the multiply covered exterior,
  which canonical and Smarr formulations never construct.
* *Interpretation.* Each half-turn around the Euclidean horizon is a horizon monodromy, a shift $i\beta/2$ of Lorentzian
  Killing time. Section 3 shows its effect on the exponent is an area term over $4G$, so the winding is the Euclidean
  shadow of horizon-entropy bookkeeping.

## 7. Consistency of conventions across the codebase

| convention | checked by |
|---|---|
| sign of κR in the Israel condition | symbolic derivation; dust-shell sanity check; thick-wall data (jump $r'_+-r'_-=-4\pi G\sigma R$) |
| ε as ${\rm sgn}\,n^r$, region labels | itinerary vs Kruskal/dS charts (15/15); finite-difference $n^r=A$ |
| Euclidean orientation (exterior to the right) | multiplicity maps consistent with NG3; boundary loops close with the correct winding |
| $B$ normalisation (full bounce) | $M=0$ closed form = CdL exponent; thick-wall CdL bounces |
| Hamiltonian clock direction | $\partial H/\partial P=dR/dt_\infty$ in region III ($dR/dt<0$ for $\dot R>0$) |

## 8. Unresolved assumptions and gaps (documented)

1. **NG4 away from criticality.** Numerical: monotone in $E$ on a 39×10 grid, $\min|\theta_D|/\pi=1.0056$. A proof of
   monotonicity in $E$ is open.
2. **T5** (laboratory bound): numerical only.
3. **Literature comparison** (P-L1): blocked by network policy.
4. **Thick-wall analogues** of T2/T3/NG3: Phase IV.
5. **SdS conformal diagrams**: not needed for the identities above. Open for the atlas (P-G1).
6. All results are leading-order WKB in the reduced theory (A10, A11).
