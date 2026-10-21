# Problem tree

Root question: **can an independently expanding region of spacetime arise from a physically admissible
parent through a consistent gravitational process, and can that process be formulated quantum mechanically?**

Status: ✅ answered within the stated model · 🟡 partially answered · 🔴 open · ⛔ answered negatively (obstruction).

```
Q0  Daughter universe from a parent?
├── Q1  Classical GR, asymptotically flat parent
│   ├── Q1.1 Which thin-wall branches inflate and how do they sit globally?        ✅ GLOBAL_GEOMETRY.md
│   ├── Q1.2 Can an inflating branch start from regular laboratory data?            ⛔ Farhi–Guth, NG2
│   ├── Q1.3 Which assumptions does the obstruction need?                            ✅ NO_GO_THEOREMS.md Sec. 1, 7
│   ├── Q1.4 Does a thick (field-theoretic) wall change it?                          ⛔ classically (NEC holds); 🟡 dynamics (Phase IV)
│   └── Q1.5 Does a closed (de Sitter) parent evade it?                              ✅ yes, globally regular solution exists; 🔴 creation from regular data
├── Q2  Semiclassical transition (FGG)
│   ├── Q2.1 Initial/final states, barrier                                           ✅ TUNNELLING.md Sec. 1
│   ├── Q2.2 Is the canonical exponent unambiguous?                                  ✅ T0 (analytic branch selected)
│   ├── Q2.3 Is there a manifold Euclidean instanton?                                ⛔ NG3 (regime A, theorem); NG4 (regime B, numerical)
│   ├── Q2.4 Can the exponent be defined without the pseudomanifold?                 ✅ T3 Smarr identity, integral representation
│   ├── Q2.5 How large is the exponent; how does it compare to vacuum decay?         ✅ T1, T5, T6
│   ├── Q2.6 Do Lorentzian / complex saddles reproduce B? Do they remove the ambiguity or hide it?  🔴 Phase III
│   └── Q2.7 Prefactor, probability interpretation                                   🔴 P-T3
├── Q3  Quantum theory of the transition
│   ├── Q3.1 Reduced phase space, Hamiltonian, sectors                               ✅ QUANTIZATION.md Sec. 2
│   ├── Q3.2 Operator formulation (imaginary area shifts)                            🟡 formulated; eikonal verified
│   ├── Q3.3 Exact solution at M = 0 vs WKB                                          🔴 P-Q2
│   ├── Q3.4 Dirac vs reduced quantization                                           🔴 P-Q1
│   ├── Q3.5 Distinct physical sectors? Superpositions of global geometries?         🟡 sectors = sheets of one constraint
│   └── Q3.6 Unitarity and information across the transition                         🔴 Phase V
├── Q4  Topology and causal structure
│   ├── Q4.1 Does spatial topology change classically?                               ✅ no (R^3 throughout); separation is causal + singular
│   └── Q4.2 Is a non-manifold (degenerate / non-Hausdorff) description needed?      🟡 Euclidean side: yes (NG3)
└── Q5  Observability
    ├── Q5.1 What would an outside observer see?                                     ✅ a black hole of mass M (D3)
    └── Q5.2 Any signature distinguishing models?                                    🔴 Phase VI
```

## The fundamental question selected for Stage Two (task M6)

**Problem P\*.** Construct the transition amplitude $\mathcal A(R_1\to R_2;M)$ of the reduced quantum theory
(QUANTIZATION.md Route A) as an exact solution of the area-shift constraint, at $M=0$ first and then at small
$M$. Determine whether a Lorentzian/complex-saddle evaluation of the gravitational path integral reproduces
(i) the exponent $B$ and (ii) the decomposition of T3. Show that its contour is forced: the real barrier
contour, homotopically trivial with respect to the horizon branch points.

*Assumptions:* A1–A6, A10–A14. *Success criteria:* (a) an exact $M=0$ amplitude whose asymptotics equal
$e^{-B/2\hbar}$ with $B$ from T1; (b) a saddle-point evaluation at $M>0$ whose real part equals T3; (c) a proof that
alternative contours are excluded by a stated physical condition, or a precise statement that they are not.
*Failure modes worth recording:* the Lorentzian contour is not unique (Stokes phenomena), which would mean the
FGG ambiguity is hidden rather than resolved.
