# Round 8 — model U against current and future data

**Task.** Test the predictions against existing data and project them onto future machines.

**Answer: yes.**
* **Current data.** Model U survives every dataset we can confront it with in its predicted window: LHC monojet and black-hole searches,
  and neutron-star survival.
* **Future machines.** It is decisively testable at FCC-hh. It is not testable at the HL-LHC.
* **A correction.** One was found and made: the M_D convention (F64).

Code `data8.py`; script `scripts/round8.py`; tests `tests/test_round8.py` (2); data `DATA/round8.json`; figure `FIGURES/round8.png`.

## 1. The convention correction [F64]

* **Our convention.** Our horizon radius is the Dimopoulos–Landsberg form, equivalent to G_D = 1/M_D^{n+2} (audited symbolically in 6p).
* **The experiments' convention.** They quote M_D in the PDG (Giudice–Rattazzi–Wells) convention, G_D = (2π)^n/(8πM_D^{n+2}).
* **The conversion.** M_D(DL) = M_D(PDG)·[8π/(2π)^n]^{1/(n+2)} = 0.89, 0.63, 0.50, 0.43, 0.38 × M_D(PDG) for n = 2–6.
* **What was wrong.** Rounds 6k–6l compared the two directly ("M_D ≳ 5–11 TeV allowed").
* **What changes.** Converted, the CMS monojet bounds (arXiv:1712.02345: M_D(PDG) > 9.9, 7.5, 6.3, 5.7, 5.3 TeV for n = 2–6) become
  **M_D(DL) > 8.8, 4.7, 3.2, 2.4, 2.0 TeV**. The physical statements in terms of M_ext (in TeV) are unaffected. The allowed region in
  M_D(DL) is wider than stated before, and every benchmark used (M_D(DL) = 6–12 TeV) is allowed.

## 2. Current data

| Test | Model-U expectation | Data | Verdict |
|---|---|---|---|
| LHC monojet (ATLAS, 139 fb⁻¹: σ×A×ε < 0.3 fb in the tightest inclusive region) | near-threshold holes appear as an initial-state jet + invisible remnant | null | excludes M_ext ≲ 8.0–9.0 TeV with full trapping, or ≲ 4.7–5.5 TeV with the conservative trapped fraction (n = 4–7) |
| LHC multi-object black-hole searches (high S_T, ≥ 3 objects) | only 1–6 % of model-U events pass (GR holes: 75–97 %) | null | these limits **do not apply** to model U; quoting them against model U would be wrong |
| Monojet ADD bounds (P1) | KK-graviton emission is the same in every interior class | M_D(PDG) > 5.3–9.9 TeV | M_D(DL) > 2.0–8.8 TeV; benchmarks allowed |
| Neutron-star / white-dwarf survival (cosmic-ray remnants) | **[Proposition-N 8-1]** a captured remnant accretes ~10¹³–10¹⁴ GeV/s geometrically and re-radiates it, sitting (2–40)×10⁻⁹ above extremality at T ≈ 0.16 GeV | neutron stars exist | **no growth**: consistent and safe (cf. Giddings & Mangano 2008) |

The model-U window, after correction, is M_ext ≈ 5–9 TeV (lower edge, depending on trapping) to ≈ 48 TeV (upper edge at 100 TeV). The
surviving region in (M_D(DL), c_l) for n = 6 is shown in panel (f) of the figure.

## 3. Future data

| Machine | Reach in M_ext (≥ 10 holes, conservative trapping) | Covers the window? |
|---|---|---|
| LHC Run 3 (13.6 TeV, 0.3 ab⁻¹) | 5.7 TeV | no |
| HL-LHC (14 TeV, 3 ab⁻¹) | 6.3 TeV | marginal: only the conservative lower edge |
| HE-LHC (27 TeV, 15 ab⁻¹) | 12.7 TeV | the lower part |
| FCC-hh (85 TeV, 20 ab⁻¹) | 40.6 TeV | most of it |
| 100 TeV (30 ab⁻¹) | 48.2 TeV | essentially all |

The topology is background-poor: an invisible remnant recoiling against a multi-TeV jet, with MET/S_T ≈ 1, and MET far above the
Z(→νν)+jet tail. So a handful of events suffices to establish it (Round 7). Ultra-high-energy cosmic rays and neutrinos produce model-U
holes too rarely to test the window.

## 4. Status

| Requirement | Status |
|---|---|
| Consistent with every current dataset in its window | yes (monojet, black-hole searches, neutron stars) |
| Convention error found and corrected | yes (F64) |
| Testable by planned machines | yes: FCC-hh 85 / 100 TeV cover the window; the HL-LHC does not |
