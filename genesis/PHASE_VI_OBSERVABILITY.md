# Phase VI — What nature might allow us to test

**Exit condition (from the programme):** *a set of physically meaningful predictions, or a rigorous explanation of why the
proposed transition is observationally inaccessible.* Status: **both, with the scope stated.**

* **Individually inaccessible (Proposition VI-1).** No observer outside the black hole can tell whether it hosts a daughter. This is rigorous in
  spherical symmetry: Birkhoff, plus the daughter region lying in the black-hole region.
* **Loophole search (Sec. 2b–2c).** We optimised over all vacuum and wall scales under the physical constraints. The best conceivable
  laboratory case gives $P\le10^{-46}$ per bubble, at the edge of validity near the Planck scale. The only unsuppressed route is cosmological.
* **Laboratory creation exponentially suppressed wherever the calculation is valid (Proposition VI-2).** The bound $B_{\rm lab}\ge0.56\,B(0)$ (T5) is small
  only where $B_{\rm CdL}\lesssim2$, i.e. where the false vacuum is not even metastable and no laboratory bubble can be prepared ("no free lunch").
* **The only potentially testable imprint is population-level.** In an inflating early universe the Farhi–Guth obstruction (Prop. IV-1) does not
  apply, and daughters can form classically. Their black holes have a characteristic mass
  $$M_*\simeq M_{\rm cr}=\frac1{2GL}\approx5.6\ {\rm g}\Big(\frac{10^{16}\ {\rm GeV}}{E_{\rm vac}}\Big)^2,$$
  which maps vacuum energy scales onto primordial-black-hole windows. A daughter-hosting PBH population is constrained exactly like an
  ordinary PBH population of the same masses, and is indistinguishable from it.

No experimental or technological proposal is made or implied.

Code: `CODE/babyuniverse/observability.py`; script `CODE/scripts/phase5_6_analysis.py`; data `DATA/phase6_observability.json`;
figure `FIGURES/phase6_observability.png`.

---

## 1. Proposition VI-1 (exterior indistinguishability)

> In spherical symmetry, whatever its interior, a vacuum exterior is Schwarzschild with mass $M$ (Birkhoff). In every daughter configuration
> found in Stage One and Phase IV, the inflating region is anti-trapped. By Corollary IV-2 it lies to the future of data beyond an
> Einstein–Rosen bridge, in a white hole, or beyond a cosmological horizon; in every case it lies outside the causal past of the exterior's
> future null infinity (it is inside the black-hole region). Hence no observation at $\mathscr I^+$, and no exterior measurement at any time, depends on
> whether a daughter exists. The two outcomes of Phase III, $|{\rm col}\rangle$ and $|{\rm dau}\rangle$, have identical exterior observables (Theorem V-2).

**Scope.** Beyond spherical symmetry, the daughter region remains inside the event horizon by the same causal argument: anti-trapped regions are
not in the domain of outer communication. But the exterior is then Kerr-like only asymptotically, and a general no-hair argument for the
formation transient is not attempted.

## 2. Proposition VI-2 (laboratory: no free lunch)

The laboratory exponent obeys $B_{\rm lab}\ge0.56\,B(0)$ for all γ (T5, numerical). $B(0)$ is the closed form T1, which equals the CdL decay exponent of
the false vacuum when the false vacuum is the larger cap. Evaluated over vacuum and wall scales from $10^3$ to $10^{18.5}$ GeV (figure, left):

* $B_{\rm lab}\lesssim1$ only in the band where $B(0)\approx B_{\rm CdL}\lesssim2$ (red dashed line). There the false vacuum decays essentially instantly: it is not
  metastable, and a laboratory bubble of it cannot be prepared.
* Elsewhere the bound grows to astronomically large values: $\log_{10}B$ up to about 60 for low vacuum and high wall scales.
* GUT reference ($\rho^{1/4}=\sigma^{1/3}=10^{16}$ GeV): $B(0)=133$ (flat CdL 133). For semiclassical black holes ($M>10\,m_{\rm Pl}$), $B\ge7\times10^4$. The total
  probability, even allowing repeated attempts until evaporation, is below $10^{-30490}$ (Phase V Sec. 6).

> **Proposition VI-2 (conditional on T5).** If a false vacuum is metastable enough to be prepared in a laboratory bubble, the probability that
> the bubble becomes a daughter universe is at most $e^{-0.56\,B_{\rm CdL}}$. Since metastability requires $B_{\rm CdL}\gg1$, laboratory daughter creation is
> exponentially suppressed in every regime where the semiclassical calculation applies. In the single-flip class, the sector-blind
> alternative $B_{\rm sq}<B_{\rm an}$ (Phase III) lowers the exponent. Phase II/III data show $B_{\rm sq}/B_{\rm an}\in[0.5,1]$, so the conclusion is unchanged. Re-proving T5 for
> $B_{\rm sq}$ is open.

## 2b. Loophole search: is there any way in the laboratory?

"The probability is small" is a statement about typical parameters. We therefore searched for the **most favourable** laboratory scenario
the calculation can describe (`observability.best_laboratory_case`). We minimised the exponent, taking the smaller of $B_{\rm an}$ and $B_{\rm sq}$, over
vacuum scale $\rho^{1/4}\in[10^{10},5\times10^{18}]$ GeV and wall/vacuum ratio $\sigma^{1/3}/\rho^{1/4}\in[0.5,1.6]$: 1620 scenarios. Three physical requirements apply:

1. **Metastability.** The false vacuum must not decay during the bubble's life: $B_{\rm CdL}>\ln(\rho\cdot\tfrac{4\pi}3R_1^4)+3$.
2. **Semiclassical black hole.** $M\ge10\,m_{\rm Pl}$, or the more conservative $100\,m_{\rm Pl}$. Smaller "black holes" are not geometric objects (Phase V Sec. 6).
3. **WKB validity.** $\hbar G|d\theta/du|<0.1$.

| black-hole mass | admissible scenarios | best case | $B_{\min}$ | $P_{\max}=e^{-B}$ |
|---|---|---|---|---|
| $10\,m_{\rm Pl}$ (0.2 mg) | 1002/1620 | $\rho^{1/4}=9\times10^{17}$ GeV, $\sigma^{1/3}=0.80\rho^{1/4}$ ($B_{\rm CdL}=8.7$), $R_1=39\,\ell_P$, $2GM=20\,\ell_P$ | 107 | $10^{-46.5}$ |
| $100\,m_{\rm Pl}$ | 925/1620 | $\rho^{1/4}=2.9\times10^{17}$ GeV, ratio 0.85 ($B_{\rm CdL}=18.8$), $R_1=410\,\ell_P$ | 1510 | $10^{-656}$ |

The structure of the optimum explains why there is no way round it:
* **Minimise $B_{\rm CdL}$.** The wall must be tuned slightly lighter than the vacuum (ratio about 0.8–0.9). That keeps the vacuum only
  marginally metastable, with $B_{\rm CdL}\approx9$–40. Going further makes the false vacuum decay before the bubble exists.
* **Use the highest scale.** At fixed ratio, $B_{\min}$ falls roughly as $E_{\rm vac}^{-0.9}$ (figure, right panel). The best case is therefore always at
  the edge of validity, near the Planck scale ($\varepsilon\approx2.6\times10^{-3}$ at the optimum). Beyond it the semiclassical calculation, and hence any
  claim either way, fails.
* **Use the smallest mass.** $B$ grows steeply with $M$, so the optimum is at the smallest mass that still makes a geometric black hole.

**Result.** Within the domain where the theory can be trusted, the best conceivable laboratory scenario has probability at most about
$10^{-46}$ per bubble. It needs Planck-scale energy densities confined to about 40 Planck lengths, and it yields a 0.2 mg black hole that evaporates
in about $10^{-40}$ s. The daughter it might contain would be undetectable (Prop. VI-1). Outside that domain (trans-Planckian, $\varepsilon\gtrsim0.1$, or
sub-Planckian masses) the question is outside semiclassical gravity. No conclusion, positive or negative, can be drawn there without a
theory of quantum gravity.

## 2c. Which assumption would have to fail?

The obstructions, collected:

| route | obstruction | assumption it rests on | can it fail? |
|---|---|---|---|
| classical, from laboratory data | Farhi–Guth (Penrose's theorem); Prop. IV-1 in spherical symmetry. The proof uses only $T_{uu}\ge0$, so it holds for **any NEC matter** (scalars, gauge fields, monopoles: "topological inflation" from laboratory data is excluded too) | the null energy condition; initial data without anti-trapped spheres | NEC is violated by quantum fields. But Wall's generalized-second-law singularity theorem (Wall 2013; cited, not re-checked) reaches the same conclusion assuming only the GSL |
| quantum tunnelling from laboratory data | $P\le e^{-B}$, best case $10^{-46}$ (Sec. 2b) | semiclassical gravity (C1–C4 of Phase III or the sector-blind alternative) | only beyond the Planck scale, where nothing is computable here |
| classical, from cosmological data | none: Cor. IV-2 allows it once anti-trapped regions exist | the data are already inflating (super-Hubble false vacuum) | this is the route nature can use (Sec. 3) |
| from data with an existing throat ("bag of gold") | Phase IV: inflation follows | such data already contain a region beyond an Einstein–Rosen bridge; they cannot be prepared from laboratory data (Prop. IV-1) | — |

**So the way a daughter universe can actually arise is not a laboratory process.** It is the classical, unsuppressed route in an inflating
early universe (Cor. IV-2), where anti-trapped regions are generic.

## 3. The early universe: where daughters are not suppressed

Proposition IV-1 assumes an asymptotically flat exterior and initial data without anti-trapped spheres. During inflation every super-Hubble
region is anti-trapped, so these assumptions fail and Corollary IV-2 allows **classical** daughter formation with no tunnelling suppression.
Examples are bubbles nucleated during inflation that later collapse into black holes with inflating interiors. This possibility was discussed
by Garriga, Vilenkin and Zhang (2016) and Deng, Garriga and Vilenkin (2017); cited from knowledge, the full texts were not re-checked here.

**Characteristic mass.** The thin-wall regime boundary $M_{\rm cr}$ (regime A ↔ C) tends to $1/(2GL)$ as $\kappa/\chi\to0$. This was verified numerically to $2\times10^{-6}$ at
$\kappa/\chi=0.0036$, and analytically from $g_{\min}\simeq4\kappa^2/L^2$. It is the mass whose Schwarzschild radius equals the Hubble radius of the false vacuum.
For $\sigma^{1/3}\approx\rho^{1/4}=E$:

| $E_{\rm vac}$ (GeV) | $M_*$ (g) | Hawking lifetime (s, one species) | PBH window |
|---|---|---|---|
| $10^{18}$ | $5\times10^{-4}$ | $10^{-35}$ | evaporated before BBN |
| $10^{16}$ | 5.6 | $10^{-23}$ | evaporated before BBN |
| $10^{13}$ | $5.6\times10^{6}$ | $10^{-5}$ | evaporated before BBN |
| $10^{11}$ | $5.6\times10^{10}$ | $10^{7}$ | BBN/CMB constrained |
| $10^{9}$ | $5.6\times10^{14}$ | $10^{19}$ | gamma-ray constrained |
| $10^{7}$ | $5.6\times10^{18}$ | $10^{31}$ | asteroid-mass window |
| $10^{5}$ | $5.6\times10^{22}$ | $10^{43}$ | microlensing constrained |
| $10^{3}$ | $5.6\times10^{26}$ | $10^{55}$ | microlensing constrained |

Inverting: vacuum scales $E_{\rm vac}\approx2\times10^5$–$7\times10^7$ GeV put $M_*$ in the asteroid-mass window, where PBHs could be all of the dark matter.
$7\times10^9$–$7\times10^{11}$ GeV give holes evaporating during BBN/CMB, which are tightly constrained. Above about $10^{12}$ GeV, including GUT scales, the holes evaporate
before BBN and leave no relic. Window boundaries are orders of magnitude taken from the PBH literature, and the one-species lifetime
underestimates the evaporation rate by a factor of order 10.

**What this is and is not.** $M_*$ is the characteristic scale of the regime structure. It is not a derived mass function: that needs the
bubble-nucleation statistics during inflation and the embedding of thin-wall dynamics in FRW, neither of which was done here. It is
"physically meaningful" in a limited sense. It ties an observable (PBH abundance at mass $M_*$) to a microphysical parameter (the vacuum
energy scale), under the assumption that daughter-forming bubbles dominate PBH formation at that scale.

## 4. Distinguishable consequences, systematically

| candidate signature | verdict |
|---|---|
| exterior of an individual black hole (metric, multipoles, Hawking spectrum) | identical to an ordinary black hole of the same mass (Prop. VI-1, Theorem V-2) |
| gravitational waves from the formation itself | none in spherical symmetry (Birkhoff). Aspherical transients come from the exterior formation dynamics and are not sourced by the causally disconnected inflating region. Not computed here |
| gravitational waves from the vacuum transition (bubble collisions) | standard first-order phase-transition background; not specific to daughter formation |
| PBH population (abundance vs mass) | the only handle: constrains bubble abundance at $M_*(E_{\rm vac})$, identically for daughter-hosting and ordinary PBHs |
| laboratory creation | excluded in practice; exponentially suppressed wherever the calculation applies (Prop. VI-2) |

## 5. Answer to the programme's question for Phase VI

*Can nature let us test the transition?* **Not for an individual event.** A daughter universe, created by the laboratory quantum route or
formed classically in the early universe, cannot be distinguished from an ordinary black hole by any exterior measurement (Prop. VI-1). The
laboratory route is exponentially suppressed wherever the calculation can be trusted (Prop. VI-2). The only conceivable empirical contact is
statistical. If vacuum bubbles of energy scale $E_{\rm vac}$ formed in the early universe, the holes they leave have masses around
$M_*\approx5.6\,{\rm g}\,(10^{16}\,{\rm GeV}/E_{\rm vac})^2$, and PBH constraints at that mass constrain those bubbles. Such constraints test the vacuum-bubble scenario, not the
presence of daughter universes inside the holes.
