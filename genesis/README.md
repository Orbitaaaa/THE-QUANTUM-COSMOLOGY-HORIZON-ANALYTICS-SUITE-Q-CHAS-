# Project Baby Universe — Genesis Branch

A research programme on the quantum creation of new spacetime regions. It asks whether an independently
expanding daughter universe can originate from a physically admissible parent spacetime through a consistent
gravitational process.

This directory is self-contained and independent of the Q-CHAS dashboard in the rest of the repository.

## Where to start

* **[DELIVERABLE_1_FOUNDATIONS.md](DELIVERABLE_1_FOUNDATIONS.md)**: the Stage One report (all 14 required items).
* **[RESULTS.md](RESULTS.md)**: verified results, labelled theorem / numerical proposition / reproduction.
* **[PAPER/genesis_foundations.tex](PAPER/genesis_foundations.tex)**: manuscript of the original results.
* **[PHASE_II_VERIFICATION.md](PHASE_II_VERIFICATION.md)**: independent verification (40-digit re-implementation, positivity
  lemma, NG4 critical-limit theorem, SdS Smarr identity, the Euclidean winding question).
* **[PHASE_V_QUANTUM_SYSTEM.md](PHASE_V_QUANTUM_SYSTEM.md)**: the reduced quantum system: two inequivalent channel ends (the daughter
  needs its own clock), sectors, information, topology (preserved; a throat is created), semiclassical validity.
* **[PHASE_VI_OBSERVABILITY.md](PHASE_VI_OBSERVABILITY.md)**: individually inaccessible; laboratory "no free lunch"; characteristic
  daughter-hosting black-hole mass $M_*\approx5.6\,{\rm g}\,(10^{16}{\rm GeV}/E_{\rm vac})^2$ mapped onto PBH windows.
* **[PHASE_IV_FIELD_DYNAMICS.md](PHASE_IV_FIELD_DYNAMICS.md)**: Einstein–scalar evolution of laboratory bubbles (collapse to black holes,
  no anti-trapping: Prop. IV-1) and of both inflating daughter branches, with thin-to-thick and stability studies.
* **[PHASE_III_QUANTUM_TRANSITION.md](PHASE_III_QUANTUM_TRANSITION.md)**: the quantum transition. Why the FGG slicing construction
  fails (III-5); why every formulation hides the FGG ambiguity in one datum (III-2, III-4, III-6); the defensible construction and its
  assumptions; the surviving competitor.

## Stage One in five lines

1. **Classical picture reproduced and sharpened.** Every inflating thin-wall branch lies beyond the Einstein–Rosen
   bridge, is invisible to the parent and has a singular past. Two independent constructions agree on all 15 branches.
2. **Farhi–Guth obstruction stated precisely**, with a map of which assumptions matter and which escape routes exist.
3. **Theorem (NG3):** the Farhi–Guth–Guven transition from laboratory data has **no manifold Euclidean
   instanton**; the "pseudomanifold" is forced.
4. **Theorems (T1–T3):** the canonical exponent is fixed uniquely, has a closed form at $M=0$ (instanton
   action + a Morse sum of horizon areas), and obeys a **first law and a Smarr identity** that express it
   covariantly without the pseudomanifold.
5. **Physics:** laboratory creation is not exponentially easier than the decay of the same false vacuum; for
   GUT-scale parameters and gram-scale bubbles $B\sim10^6$ or more.

## Documents

| file | content |
|---|---|
| [RESEARCH_QUESTIONS.md](RESEARCH_QUESTIONS.md) | problem tree with status; the selected Stage Two problem P\* |
| [LITERATURE.md](LITERATURE.md) | primary sources, verification status, novelty comparison |
| [ASSUMPTIONS.md](ASSUMPTIONS.md) | assumptions register |
| [CLASSICAL_MODEL.md](CLASSICAL_MODEL.md) | action, boundary terms, field equations, admissible data, energy conditions, thick-wall model |
| [JUNCTION_CONDITIONS.md](JUNCTION_CONDITIONS.md) | complete derivation of the Israel conditions, signs, potential |
| [GLOBAL_GEOMETRY.md](GLOBAL_GEOMETRY.md) | branch classification, chart-free itinerary, Penrose diagrams, daughter-universe criteria |
| [NO_GO_THEOREMS.md](NO_GO_THEOREMS.md) | Farhi–Guth, bag-of-gold theorem, exterior embedding obstruction (proof) |
| [TUNNELLING.md](TUNNELLING.md) | FGG reconstruction, branch selection, closed form, first law, Smarr identity, scales |
| [QUANTIZATION.md](QUANTIZATION.md) | reduced phase space, sectors, imaginary area shifts, quantization routes |
| [NEW_APPROACHES.md](NEW_APPROACHES.md) | ranked new approaches |
| [FAILED_APPROACHES.md](FAILED_APPROACHES.md) | dead ends, with reasons |
| [NUMERICAL_METHODS.md](NUMERICAL_METHODS.md) | algorithms and accuracy control |
| [VERIFICATION.md](VERIFICATION.md) | verification dossier and proof gaps |
| [OPEN_PROBLEMS.md](OPEN_PROBLEMS.md) | open problems by phase |

## Phase III in five lines

1. The two FGG prescriptions are two conventions for the boost angle of a flipped wall normal: complex (Sorkin, $i\pi$) or real.
2. The true Hamiltonian is two-sheeted. The tunnelling trajectory touches its branch locus exactly where the Euclidean wall
   folds in Killing time.
3. The FGG slicing (principal) prescription solves no analytic quantum constraint, and $B<0$ in 33/36 models.
4. **Double-flip transitions:** every consistent quantization gives the same exponent $B_{\rm an}$, so the ambiguity is eliminated there.
5. **Single-flip transitions:** two consistent answers remain, $B_{\rm an}$ if the wall's exterior orientation is superselected and
   $B_{\rm sq}=\frac2G\int R\arccos|C|\,dR<B_{\rm an}$ if not. The universal Stokes multiplier is $\sqrt2$.

## Phases IV–VI in brief

* **IV (fields).** The full Einstein–scalar evolution reproduces the thin-wall picture. Laboratory bubbles collapse to black holes and never
  develop an anti-trapped region (Prop. IV-1, valid for any NEC matter). Both inflating daughter branches exist as field-theory solutions,
  but only from data that already lie beyond a bridge, in a white hole, or beyond a cosmological horizon.
* **V (quantum system).** The reduced quantum theory is a scattering problem. The parent channel ends at $t=\infty$; the daughter channel at
  finite $t$, where it needs its own clock. The daughter is invisible to every exterior observable. The transition preserves spatial topology
  and creates a throat.
* **VI (observability).** Individually undetectable. Laboratory creation: best conceivable case $P\le10^{-46}$, at the edge of validity near the
  Planck scale. The only unsuppressed route is cosmological, with a characteristic black-hole mass $M_*\approx5.6\,{\rm g}\,(10^{16}\,{\rm GeV}/E_{\rm vac})^2$.

## Code, data, figures

```
CODE/babyuniverse/   model, symbolic, branches, geometry, plotting, causal, tunnelling, thickwall, quantum, contours,
                     doublenull, stability, quantum_system, observability
CODE/scripts/        run_all.py reproduces every file in DATA/ and FIGURES/
CODE/tests/          pytest suite (symbolic identities, junction, atlas, geometry, tunnelling, thick walls, quantum)
DATA/                JSON results, derivation log
FIGURES/             Penrose diagrams, phase diagram, exponents, Euclidean geometry, thick walls
```

Reproduce (Python ≥ 3.10):

```bash
pip install -r genesis/CODE/requirements.txt
cd genesis/CODE/scripts && python run_all.py      # regenerates DATA/ and FIGURES/ (Stage One + Phases II-VI)
cd .. && python -m pytest tests -q                # 106 tests
```

Units: $G=c=\hbar=1$ unless stated; signature $(-,+,+,+)$.

## Honesty statement

No universe has been created, and nothing here suggests a practical way to create one. The results are
mathematical statements about a specified gravitational model (thin or scalar-field walls, spherical symmetry,
semiclassical quantum theory), with their assumptions listed in ASSUMPTIONS.md and their gaps in
VERIFICATION.md.
