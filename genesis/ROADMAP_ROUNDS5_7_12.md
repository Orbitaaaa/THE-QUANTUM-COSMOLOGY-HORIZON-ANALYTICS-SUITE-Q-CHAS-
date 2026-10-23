# Rounds 5, 7 and 8–12: field theory, observables, and what the later rounds would need

Companion to [ROADMAP_ROUND6_MECHANISM_ATLAS.md](ROADMAP_ROUND6_MECHANISM_ATLAS.md).

Round 5 asks for a field theory that *derives* the bubble's material properties; Round 7 asks for radiation from an explicit quantum
field. Rounds 8–12 (experiment, engineering, integration) are treated **as theory only**: what each would have to establish, derived from
Rounds 5–7. **No apparatus is designed and no experimental procedure is proposed.**

**Files.**

* Code: `CODE/babyuniverse/field_theory.py` (Round 5), `radiation.py` (Round 7).
* Script: `CODE/scripts/round5_7.py`.
* Data: `DATA/round5_7.json`.
* Figure: `FIGURES/round5_7.png`.
* Tests: `test_field_theory.py` (3), `test_radiation.py` (2).

---

## Round 5 — GENESIS FIELD THEORY

### R5.1 Theory F1: the false vacuum, fully specified
Take one real scalar, minimally coupled, with the tilted quartic potential
$$U(\phi)=\frac\lambda4(\phi^2-v^2)^2+\rho_f\big(1-s(\phi/v)\big),\qquad s(x)=\tfrac12+\tfrac34x-\tfrac14x^3,\qquad \rho_f=\epsilon\,\frac{\lambda v^4}4 .$$
The minima stay at $\pm v$, and $-v$ is a minimum for $\epsilon<16/3$. Three numbers $(\lambda,v,\epsilon)$ then **derive** everything the
thin-wall analysis took as input (`field_theory.derived_parameters`):

* the false-vacuum energy $\rho_f$;
* the tension $\sigma=\frac{2\sqrt2}3\sqrt\lambda v^3$;
* the width $w=\sqrt2/\sqrt\lambda v$;
* the scalar masses;
* the flat CdL radius $R_0=3\sigma/\rho_f$;
* the light-wall parameter $\kappa/\chi$;
* the CdL exponent, both thin-wall and from the full gravitational $O(4)$ bounce (`thickwall.cdl_bounce`).

**The chain closes end to end:** couplings → $(\sigma,\rho)$ → $E=(R_{\rm CdL}/R_M)^2$ → $B_{\rm lab}=B_{\rm CdL}\Phi(E)$ (Theorem S3)
(`laboratory_exponent_from_theory`).

| $(\lambda, v/m_P, \epsilon)$ | $E_{\rm vac}$ | $w/R_0$ | $\kappa/\chi$ | $B_{\rm CdL}$ (bounce / thin wall) | $B_{\rm lab}$ for $M=m_P$ |
|---|---|---|---|---|---|
| (1, 0.02, 0.3) | $1.3\times10^{17}$ GeV | 0.037 | 0.30 | 208 425 / 249 543 | $2.1\times10^5$ |
| (4, 0.01, 0.5) | $1.0\times10^{17}$ GeV | 0.062 | 0.12 | 13 073 / 13 475 (for $M=10m_P$: $B_{\rm lab}=2.6\times10^4$) | — |

For these couplings the one-loop correction relative to the tree potential ($\sim\lambda/4\pi^2$) is 0.025 and 0.10. The Planck-suppressed
operators ($\sim8\pi v^2$) contribute 0.010 and 0.0025. The potential is therefore stable under both.

### R5.2 Proposition F1: the frontier needs strong coupling — [Proposition, checked against the bounce to 3–17%]
In the thin-wall regime,
$$B_{\rm CdL}=\frac{27\pi^2\sigma^4}{2\rho^3}\simeq\frac{6745}{\epsilon^3\lambda},\qquad \lambda=\frac{64}{81}\,\frac{(4\rho/\epsilon)^3}{\sigma^4}\ \text{for given }(\sigma,\rho).$$

* The semiclassical frontier of Corollary S7 ($B_{\rm CdL}=10$ at $2.7\times10^{17}$ GeV) needs $\lambda\approx2.6\times10^4$ at $\epsilon=0.3$, or
  $5.6\times10^3$ at $\epsilon=0.5$. There the one-loop correction is 650 and 140 times the tree potential.
  **No perturbative quartic theory reaches the frontier in the thin-wall regime.**
* With $\lambda\le4\pi$ and $\epsilon\le0.5$: $B_{\rm CdL}\ge4.3\times10^3$ and $B_{\rm lab}\ge1095$, so **$P\le10^{-476}$**.
* With $\lambda\le1$ and $\epsilon=0.3$: $B_{\rm lab}\ge8350$, so $P\le10^{-3627}$ (figure, left panel).
* Small $B_{\rm CdL}$ with a perturbative coupling exists only for thick walls ($\epsilon\to16/3$). There the laboratory analysis (S3, S7) has not
  been established; the thick-wall FGG problem (P-T2) is open.

**Consequence for Problem B.** Requiring the material properties to come from a weakly coupled field theory makes the known-physics route far
*less* accessible than the Round 2 frontier suggested. The $10^{-23}$ figure assumed a wall that no perturbative theory produces. Stabilised seeds
(Round 6, G5) inherit this: with $B\ge1095$ even a $10^{41}$ Hz attempt rate gives lifetimes far beyond the age of the universe.

### R5.3 Consistency with observed physics
* **Theory F1.** A hidden-sector scalar at $10^{16}$–$10^{17}$ GeV with Planck-suppressed couplings to the Standard Model is unconstrained by
  collider and fifth-force data. It cannot be excited in a laboratory: preparing its false vacuum needs energy densities of
  $(10^{17}\,{\rm GeV})^4\approx2\times10^{105}$ J/m³.
* **Theory F2 (Route G, the NEC sector).** Galilean genesis,
  $\mathcal L=f^2e^{2\pi}(\partial\pi)^2+(f^3/\Lambda^3)(\partial\pi)^2\Box\pi+(f^3/2\Lambda^3)(\partial\pi)^4$, with $H_0^2\sim\Lambda^3/f$ and
  $h\sim f^2/9M_P^2H_0^2$ (`genesis_scales`; O(1) factors not fixed).
  * For $f=10^{10}$ GeV and $\Lambda=10^5$ GeV, $H_0/\Lambda=0.003$ (EFT-consistent) and $h\approx10^{-54}$ m², so any laboratory seed size is in the
    flat-space genesis regime.
  * The constraints on such a sector are theoretical (positivity, GSL, superluminality; Round 6 C5) and, if it couples to matter, fifth-force
    tests with Vainshtein screening.
  * GW170817 constrains late-time Horndeski operators that change the speed of gravitational waves, not a high-scale genesis sector.

### R5.4 Regime of trust
* Semiclassical gravity needs curvature well below Planckian: $R\gtrsim100\,\ell_P$ and $M\gtrsim m_P$ for the black hole.
* The thin-wall description needs $w/R_0\lesssim0.1$ and perturbative couplings ($\lambda\lesssim4\pi$). Proposition F1 shows these are jointly
  incompatible with $B_{\rm CdL}\lesssim10^3$ at $E_{\rm vac}\sim10^{17}$ GeV.
* Route G's EFT needs $H\ll\Lambda$ throughout the seed's evolution.

**Required result of Round 5 (met for F1, conditional for F2).** Theory F1 predicts the bubble's properties without input by hand, and with them a
laboratory exponent. Theory F2 specifies Route G's NEC-violating phase. Its seed boundary and hand-over remain open (P-R6b).

---

## Round 7 — GENESIS OBSERVABLES

### R7.1 What a detector would record: an explicit quantum field
Take a massless scalar in the s-wave, 1+1 approximation (no angular barrier, no greybody factors). The **laboratory state** is
time-symmetric data at $R_1$ with the field in the conformal vacuum on that slice; no horizon lies on the slice, so the state is regular.

* **Collapse branch.** Every ray reaching $\mathcal I^+$ after the turning-point signal $u_0$ left the wall on the collapse leg, and traces back to the
  initial slice inside the bubble. The map is $u\mapsto U$, with $p'(u)=(dU/d\tau)/(du/d\tau)$ exactly at the exit. The Davies–Fulling flux is
  $\langle T_{uu}\rangle=(1/24\pi)[\tfrac12y'^2-y'']$ with $y=\ln p'$ (`radiation.flux_lab`).
* **Check.** A separate ray-tracing through the bubble for the eternal history (`radiation.flux`) reproduces the Hawking flux $\kappa^2/48\pi$
  ($\kappa=1/4GM$) to 0.02%.

**Result [Proposition-N].** The same three features appear for three different bubbles:

| model | flux at $u_0$ / Hawking | 90% switch-on, in units of $4GM$ | energy relative to immediate thermal emission |
|---|---|---|---|
| $\gamma=0.6$, $E=2.4$ | −0.45 | 3.38 | −1.94 $F_H/\kappa$ |
| $\gamma=0.3$, $E=5$ | −0.13 (min −0.35) | 3.69 | −2.18 $F_H/\kappa$ |
| light wall $\kappa=0.3\chi$ | −0.44 | 3.50 | −2.02 $F_H/\kappa$ |

In words: a negative-energy precursor, a switch-on of Hawking emission over $\approx3.5\times4GM$, and a deficit of about two thermal units.

* **Daughter branch.** The field's state in the newly exterior region after the transition is not fixed by the semiclassical theory. Under the stated
  assumption D1 (the state regular on the Euclidean section), the outgoing flux is thermal from $u_0$ on. **The branches then differ by the
  collapse transient**, about $2F_H/\kappa\approx0.1$ Hawking quanta per event.
* **Late times.** The flux is identical on both branches: Birkhoff fixes the exterior, and the Hawking flux is universal for any state that is
  regular at the horizon.

### R7.2 Correlations, backreaction, measurement
* **Correlations.** In the D1 state the partners of the early Hawking quanta lie across the bifurcation sphere, on the daughter's side. On the
  collapse branch they lie inside the black hole. Either way the exterior reduced state is thermal at late times. Where the purification
  lives is not measurable from outside (no-signalling).
* **Backreaction.** Whatever the transient emits is taken from the hole's mass, which is $M_{\rm ADM}$ on both branches.
* **Measurement.** A detector at infinity sees the light-curve truncation of Round 3 if the wall itself radiates, plus the transient above. Both are
  event-level, history-dependent features, and both are statistical. At $M\sim m_P$ the transient is about 0.1 quanta, so telling the branches
  apart needs $\gtrsim10^2$ events per branch, on top of the $\ge e^{B}$ trials needed to obtain one event.

### R7.3 What is fundamentally inaccessible
* **The daughter's inflation.** Every daughter lies beyond the bridge (G1, G2), and Birkhoff fixes the exterior. No exterior observable can establish
  that the region behind the horizon inflates.
* **Which branch occurred.** At late times this is invisible too. Only the event-level history (light-curve truncation, transient) carries it, and
  only statistically.
* **Visible inflation.** It would be observable, but it means the parent is being consumed (Proposition G3). It is excluded as a goal.

**Required result of Round 7.**
* *Met:* an explicit-field prediction for the collapse branch, the size of the branch difference under a stated state assumption, and the limits on
  what can be inferred.
* *Open:* the daughter branch's state from a full treatment of the transition (P-C4), and 3+1 greybody factors.

---

## Rounds 8–12 — theory-only requirements

### Round 8: experimental tests of the underlying physics
Tests of *ingredients*, none of which creates a daughter:

* **Tunnelling mathematics in controlled systems.** False-vacuum decay by bubble nucleation has been observed in a ferromagnetic superfluid
  (Zenesini et al., Nature Physics 20, 558, 2024). There, nucleation matched instanton theory, and the decay was found to be thermally activated. Analogue systems test the
  bounce formalism of Rounds 1–3, not gravitational creation.
* **Route G's gates.** Each can be tested without making a universe:
  * the existence of a sector with stable NEC violation;
  * forward-scattering positivity in that sector;
  * superluminal propagation around backgrounds;
  * classical violation of the black-hole area theorem.
  Of the five Round 6 gates, (G-a) to (G-d) are these tests; (G-e), preparability, is a calculation (P-R6b).
* **Numerical cross-checks.** An independent Einstein–matter code should reproduce Prop. IV-1, Theorem G1 and Proposition G7.

### Round 9: GENESIS-1 preliminary requirements (derived, not designed)
| item | known-physics route (FGG tunnelling, Theory F1) | Route G (conditional on the gates) |
|---|---|---|
| physical mechanism | quantum tunnelling of a prepared false-vacuum bubble | classical evolution of a self-amplifying NEC-violating seed |
| initial state | false vacuum at $\ge10^{17}$ GeV in a region of $\sim10^{2}$–$10^{3}\,\ell_P$ | a seed of the NEC phase, any size, with a throat-opening boundary layer (G7) |
| energy budget | $\ge m_Pc^2$ concentrated at Planck-adjacent density | the seed's ADM energy, minimised; floor $m_Pc^2\approx2$ GJ |
| probability per attempt | $\le10^{-476}$ (perturbative), $\le10^{-23}$ (strong coupling, outside trust) | ≈ 1 inside the attractor basin |
| control | none (one attempt per bubble) | deterministic formation time $\Delta t(L,h,p_0)$ |
| detection | light-curve truncation; transient ≈ 0.1 quanta per event | the same exterior observables |
| safety | black-hole flash $\ge m_Pc^2$ per event | flash $\ge m_Pc^2$ per event; a seed with the wrong edge consumes the parent (G7), so this is a design-critical theoretical requirement |
| status | not realisable | not specifiable until gates G-a to G-e are resolved |

### Rounds 10–12: subsystems, integration, demonstration
* Nothing in Rounds 10–12 can be specified before Round 9 has a mechanism with resolved gates.
* For the known-physics route the requirements are Planck-adjacent, and in a perturbative theory exponentially suppressed beyond any number of trials.
  There is no subsystem to prototype.
* For Route G, the first subsystem would be a test of the gates (Round 8), not a creation apparatus.
* Replication (12C) would by Round 7 have to rest on event statistics. Even then it could not certify the daughter's inflation.
* **Safety review (Round 11) would have to establish, before anything else, that no visible-bubble branch (G3, G7) is reachable.**
