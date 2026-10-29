# Round 6p — independent mathematical audit and a publication-ready paper

**Task.** Audit the mathematics of Rounds 6l–6o independently, and write the results up as a standalone paper.

**Answer: yes.**
* **Audit: 24 of 24 independent checks pass.** Each re-derives a key result by a different route: symbolic algebra, a direct Riemann-tensor
  computation, or an independent numerical method.
* **Paper:** `PAPER/daughter_universes.tex` / `.pdf` is written, with figures, a results table, limitations, the audit appendix and verified
  references.

Code `audit6p.py`; script `scripts/round6p.py` (writes `AUDIT_ROUND6P.md` and `DATA/round6p_audit.json`); tests `tests/test_round6p.py` (3).

## 1. What was audited, and how

| Claim | Independent route | Result |
|---|---|---|
| Hayward μ_ext and T(r₊) (D = 4, 5, 6, 10) | sympy: extremum of μ(r) = r^{D−1}/(r² − α); f'(r₊)/4π with μ(r₊) substituted | exact match |
| Closed-dust inflection a_infl^{D−1} = 2αC/(D−3) | sympy: the derivative's numerator vanishes there, and the second derivative is negative | exact |
| Vacuum threshold Λ_c(C) and its duality with extremality | sympy solve of the low-curvature recollapse condition | exact |
| Hydrodynamic acceleration equation (Round 6n) | sympy derivation from the Misner–Sharp relation, first law and tr equation; GR limit −M/R² − 4πpR | exact identity |
| Kretschmann formula of the static metric (D = 4, 5) | Christoffel symbols and Riemann tensor computed directly for general f(r) | exact |
| Kretschmann formula of closed FRW (D = 4, 5) | direct Riemann computation for general a(t) | exact |
| Dimopoulos–Landsberg radius | derived from the Tangherlini mass with G_D = 1/M_D^{n+2} | relative difference 2×10⁻¹⁶ |
| Head-on trapped fraction, D = 4 | Penrose/Eardley–Giddings value 1/√2 | exact |
| N_acc (Round 6l) | time integration of the interior (DOP853) instead of root-finding and quadrature | 1.9076 vs 1.9073 |
| Curvature bound (Round 6m) | closed-form Kretschmann instead of finite differences; sup over μ = 10⁻³–10⁸ | agreement 3.5×10⁻⁴; sup K l⁴ = 188.0 |
| Cross-section (Round 6l) | Monte Carlo over (x₁, x₂) with the raw NNPDF grid instead of the luminosity table | 4.488(20)×10⁻⁴ vs 4.463×10⁻⁴ pb |
| Event-generator normalisation (Round 6o) | against an independent σ_pp integration | 13.726 vs 13.726 pb |
| Vaidya evaporation (Round 6o) | RK45 instead of LSODA | late mass excess identical to 4 digits |

**Two findings from the audit.**
* The Round 6m finite-difference Kretschmann scalar is noise-limited where K is 10⁶ times below its maximum (the nearly flat exterior).
  This does not affect the bound, which sits at the core.
* The 6l thin-null-shell argument had already been withdrawn in Round 6n (F63). The audit therefore covers the corrected chain.

## 2. The paper

`PAPER/daughter_universes.pdf` ("Daughter spacetimes from regular black holes in higher-dimensional gravity: global structure, collapse
beyond idealisation, collider production and tests") covers:
* model U;
* the global daughter spacetime;
* interiors;
* hydrodynamics beyond idealised collapse;
* collider formation, evaporation and the daughter;
* a results table and limitations;
* the audit appendix.

All references were checked against arXiv and journal records during Rounds 6l–6n. Rounds 7–9 add the predictions, the confrontation with
data and the experimental design.

## 3. Status

| Requirement | Status |
|---|---|
| Every key formula re-derived independently | yes (24/24) |
| Key numbers reproduced by independent methods | yes |
| Standalone paper with figures, limitations and references | yes |
