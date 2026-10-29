# Round 7 — experimentally distinguishable predictions

**Task.** Derive predictions that separate model U from the alternatives: low-scale gravity without black holes, GR black holes, and other
interior classes.

**Answer: yes.**
* **Model U vs a GR hole.** At the same production threshold and cross-section, model U predicts a different event topology. A detector
  could tell them apart with **about 3 events** (5σ, expected log-likelihood).
* **Between the model-U towers.** The two towers differ only in detail: separating them takes ~660 events.
* **The predictions are quantitative and located.** They give the onset energy, the temperature–mass relation, the missing-energy
  topology and the remnant velocity.

Code `predictions7.py`; script `scripts/round7.py`; tests `tests/test_round7.py` (2); data `DATA/round7.json`; figure `FIGURES/round7.png`.

## 1. Method

**Event pipeline.** We take the Round 6o generator and add:
* the hole's rapidity, from the parton momentum fractions (NNPDF2.3 grid);
* hard initial-state radiation: P(p_T > p) = 1 − exp(−A ln²(√ŝ/p)), A ≈ 0.043 (leading-log estimate), with the hole recoiling against it;
* a generic detector: |η| < 4, p_T > 50 GeV; jets, leptons and photons visible; W, Z, H, t visible at 80 %; neutrinos and a neutral
  remnant invisible.

**Discrimination.** We use the binned (multiplicity, MET/S_T) distributions. The expected log-likelihood ratio per event is the
Kullback–Leibler divergence, and n_5σ = 12.5/KL.

**Backgrounds.** These are order-of-magnitude only.
* *QCD 2 → 2:* σ̂ = Kπα_s²/ŝ with K = 10. This gives 3×10⁻⁵ pb above 8 TeV at the LHC, i.e. a handful of events in 138 fb⁻¹, as the
  measured dijet spectra show.
* *Z(→νν)+jet:* taken as 2×10⁻³ of the QCD tail.

## 2. Predictions (100 TeV; n = 6; M_D = 8 TeV; c_l = 0.6, so M_ext = 11.5 TeV; σ ≈ 13.7 pb)

| Observable | Model U (Hayward) | Model U (Dymnikova, c_l = 0.7) | GR hole, same threshold |
|---|---|---|---|
| events with no visible object | 20 % | 22 % | 0 % |
| events that are a single initial-state jet + invisible remnant | 56 % | 63 % | 0 % |
| mean number of visible objects | 1.0 | 0.9 | 2.8 |
| median S_T | 0.63 TeV | 0.45 TeV | 12.0 TeV |
| mean MET / S_T | 0.94 | 0.97 | 0.18 |
| events with MET > 1 TeV | 42 % | 37 % | 58 % |
| remnant mass and velocity | M_ext = 11.5 TeV, ⟨β⟩ = 0.39 | 11.4 TeV, 0.38 | none |
| Hawking temperature vs mass | rises from 0 at M_ext (≈ 1–2 TeV) | similar | ≈ 5 TeV, falling as M^{−1/7} |

The distinguishing predictions:
1. **Onset.** The black-hole rate switches on at M_ext = κ_n c_l^{n+1} M_D, not at a few M_D. Once monojet (KK-graviton) data fix M_D, the
   onset measures the interior class.
2. **Invisible topology.** Near threshold most events look like a monojet with MET ≈ p_T(jet), the remnant recoiling invisibly. This is a
   single-jet + MET signature at multi-TeV MET, unlike the high-multiplicity, balanced events of GR holes.
3. **Temperature–mass relation.** T(M) rises from zero at M_ext (model U) instead of falling (GR). It is measurable in the minority of
   events with several Hawking quanta, from the energy spectra against the reconstructed visible mass plus MET.
4. **Remnants.** They are heavy (≈ M_ext) and slow (β ≈ 0.4).
   * Electrically they are expected to be neutral: Schwinger pair production at r₊ ~ 10⁻²⁰ m discharges any charge, since the field is
     ~10³⁰ V/m ≫ m_e²/e.
   * If colour forces a hadronic dressing, part of them could be charged, as R-hadrons are. With a 30 % charged scenario, 15 % of
     events contain a slow heavy charged track within |η| < 2.5. This fraction is uncertain, and its absence is not decisive.
5. **What is shared with GR.** Species democracy among the emitted quanta (quarks ≈ 50 %, gluons ≈ 16 %, leptons ≈ 8 %, W/Z ≈ 10 %), and
   √s scaling of the rate. Both test that the events are black holes; neither tests the interior.
6. **Low-scale gravity itself (P1).** Monojet / monophoton rates from KK-graviton emission and the high-mass dilepton and diphoton tails.
   These are the same in all interior classes.

**Separation.** Model U against GR at equal threshold: about 3 events (both towers; 4 events at the higher threshold c_l = 0.7). Hayward
against Dymnikova: about 660 events.

## 3. Status

| Requirement | Status |
|---|---|
| Quantitative predictions that differ from GR holes | yes (onset, invisible topology, T(M), remnants) |
| Separation power | 5σ with ~3 events (model U vs GR); ~660 events between towers |
| Backgrounds estimated | order of magnitude (QCD and Z→νν tails); a full detector study is Round 9 / P-R7-a |
