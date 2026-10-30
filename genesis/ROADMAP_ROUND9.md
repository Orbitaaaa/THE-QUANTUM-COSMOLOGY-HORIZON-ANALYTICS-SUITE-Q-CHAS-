# Round 9 — the experimental design

**Task.** Finalise an experiment that would find model-U black holes if they exist in the allowed window, establish their class, and say
what each outcome means for the interior and the daughter.

**Answer: yes.** We specify a design with these parts:
* **Machine and channel.** FCC-hh at 85 TeV (100 TeV option), 20–30 ab⁻¹. The discovery channel is a monojet plus an invisible remnant.
* **Reach.** For the n = 6, M_D(DL) = 8 TeV benchmark it reaches 5σ discovery up to **M_ext ≈ 33 TeV (85 TeV) and ≈ 40 TeV (100 TeV)**.
  The class is established with ~10 selected events up to ≈ 36 / 44 TeV.
* **Charged remnants.** If 30 % of remnants are charged, their mass can be measured directly up to ≈ 39 / 46 TeV.
* **Interpretation.** A six-branch decision tree maps the outcomes onto the interior class.

The backgrounds are order-of-magnitude estimates and the detector is generic. A full simulation is listed as P-R9-a.

Code `design9.py`; script `scripts/round9.py`; tests `tests/test_round9.py` (2); data `DATA/round9.json`; figure `FIGURES/round9.png`.

## 1. Method

**Signal.**
* *Generator.* The Round 7 lab-frame generator: NNPDF2.3 rapidity, leading-log initial-state radiation, generic detector with |η| < 4.
* *Grid.* We scan c_l = 0.58–0.745, which gives M_ext = 9–52 TeV (n = 6, M_D(DL) = 8 TeV, Hayward). We also scan the Dymnikova tower over
  the same M_ext range.
* *Sample size.* 2000–2500 events per point.

**Selection: "monojet + invisible remnant".** We require:
* leading object p_T > p_cut;
* MET > p_cut;
* MET/S_T > 0.8.

At each point, p_cut ∈ {0.5, 1, 2, 3, 5, 8, 12, 16} TeV is chosen to minimise the discovery luminosity. The optimum rises from 3 TeV at
M_ext ≈ 10 TeV to 16 TeV at M_ext ≳ 34 TeV. Signal efficiencies are 1–17 %.

**Background.**
* *Estimate.* Z(→νν)+jet with MET > p_cut, multiplied by 3 to allow for W(→ℓν)+jet with a lost lepton. This is order-of-magnitude only
  (Round 7).
* *Size.* 160 fb at p_cut = 3 TeV and 0.01 fb at 16 TeV (100 TeV).

**Statistics.**
* *Discovery.* The Asimov significance Z = √(2((s+b)ln(1+s/b) − s)) (Cowan, Cranmer, Gross, Vitells 2011). Discovery means Z = 5 and s ≥ 5.
* *Class.* 10 selected events. Round 7 found that 5σ separation from GR needs ~3 events; the margin also measures the onset.
* *Remnant mass.* 10 heavy stable charged tracks with |η| < 2.5, M = p/(βγ).
* *Reach.* Interpolated log-linearly in luminosity between grid points. Because the grid is coarse, quoted reaches carry about ±1–2 TeV.

## 2. Results (n = 6, M_D(DL) = 8 TeV)

| Machine (luminosity) | 5σ discovery reach in M_ext | 10 selected events (class + onset) | 10 charged remnants (if 30 % charged) |
|---|---|---|---|
| HE-LHC 27 TeV (15 ab⁻¹) | 9.7 TeV | 11.5 TeV | 12.1 TeV |
| FCC-hh 85 TeV (20 ab⁻¹), Hayward | **33.5 TeV** | 35.8 TeV | 39.1 TeV |
| FCC-hh 85 TeV, Dymnikova | 33.5 TeV | 36.4 TeV | 38.4 TeV |
| 100 TeV (30 ab⁻¹), Hayward | **40.0 TeV** | 44.2 TeV | 46.3 TeV |
| 100 TeV, Dymnikova | 40.1 TeV | 44.4 TeV | 45.7 TeV |

**Comparison with Round 8.** The raw reach (≥ 10 holes, no selection) was 40.6 TeV at 85 TeV and 48.2 TeV at 100 TeV. The selection and
the background therefore cost 5–8 TeV of reach. The window allowed by current data starts at M_ext ≈ 5–9 TeV (Round 8). So FCC-hh at
85 TeV covers it to ≈ 33 TeV, and 100 TeV to ≈ 40 TeV.

**Signal over background.** After selection S/B > 1 up to M_ext ≈ 31 TeV (85 TeV) and ≈ 36 TeV (100 TeV). Below about 25 TeV,
S/B ≈ 3–15.

**Early running.** The luminosity needed is very small for most of the window:

| M_ext | Discovery at 85 TeV | Discovery at 100 TeV |
|---|---|---|
| 11 TeV | 7 pb⁻¹ | 3 pb⁻¹ |
| 21 TeV | 2 fb⁻¹ | 0.3 fb⁻¹ |
| 29 TeV | 0.3 ab⁻¹ | 0.014 ab⁻¹ |

Most of the window would be decided in the first year of FCC-hh operation.

**Tower independence.** Hayward and Dymnikova towers give the same reach to within 1 TeV at equal M_ext. The design depends on the onset
mass, not on the tower.

## 3. Detector and trigger requirements

| Requirement | Specification | Reason |
|---|---|---|
| Trigger | unprescaled MET > 1 TeV; single jet p_T > 1 TeV; slow-track (late-hit) trigger | the signal is a multi-TeV monojet with MET ≈ p_T(jet) |
| Hermeticity | calorimetry to \|η\| ≈ 6 (FCC-hh reference design) | MET resolution against forward mismeasurement |
| Jet energy scale / MET resolution | 1–2 % at 1–10 TeV / ~5 % at multi-TeV MET | the onset edge in the recoil spectrum |
| Timing | ~30–50 ps layers | a remnant with β = 0.39 arrives **52 ns** late at 10 m |
| Ionisation | dE/dx readout | ~7× MIP at β ≈ 0.4 |
| Mass | M = p/(βγ) from tracking + time of flight | measures M_ext directly, the zero-temperature endpoint |

## 4. Decision tree

1. **Monojet / monophoton excess consistent with KK gravitons.** Low-scale gravity (P1); n and M_D(DL) measured.
2. **No KK excess at FCC-hh.** P1 excluded up to the reach. Then there are no collider black holes of any interior class.
3. **Onset of monojet + invisible-remnant events (MET/S_T → 1) at M_th ≫ M_D.** The regular-core (model-U) class.
   * M_th = M_ext = κ_n c_l^{n+1} M_D measures c_l.
   * By the closure lemma, a remnant means the trapped region never closes. In this class the daughter persists (finite, Round 6m).
4. **Onset of high-multiplicity balanced events at a few M_D, with complete evaporation.** A GR-like exterior. The interior could be
   singular, an FMM-type daughter (unobservable), or a closing trapped region (no daughter); experiment cannot decide between these.
5. **Heavy stable charged tracks with M = p/(βγ) equal to the onset mass.** Direct confirmation that remnants carry M_ext.
6. **KK excess but no black-hole onset up to the reach.** Either model U with c_l above the window (M_ext beyond reach), or no horizon
   formation at these energies.

**What the experiment cannot do.** The daughter itself lies outside J⁻(I⁺) (Round 6m), so no outcome observes it directly. What the
experiment tests is the class of the exterior and the closure of the trapped region. In model U these fix whether the daughter exists.

## 5. Limitations (stated, not hidden)

* **Backgrounds.** The backgrounds are an order-of-magnitude estimate. A factor 10 change in the background shifts the 5σ reach by about
  2–3 TeV at the top of the window, where s ~ b. It does not affect the lower part, where s ≫ b.
* **Detector simulation.** None was done. In particular the FCC-hh pile-up (~10³ interactions per crossing at peak luminosity) is not
  simulated.
* **Charged remnants.** The charged-remnant fraction (30 %) is a scenario, not a prediction. Schwinger discharge favours neutral remnants
  (Round 7).
* **Generator.** As in Round 6o, rotation and greybody factors are omitted.

## 6. Status

| Requirement | Status |
|---|---|
| Machine, luminosity, channel and selection specified | yes (FCC-hh 85 / 100 TeV; monojet + invisible remnant) |
| Discovery, class and remnant-mass reach quantified | yes (33 / 40 TeV; 36 / 44 TeV; 39 / 46 TeV) |
| Detector and trigger requirements | yes (MET and jet triggers, \|η\| ≈ 6, 30–50 ps timing, dE/dx) |
| Interpretation of every outcome | yes (decision tree, six branches) |
| Full detector simulation with pile-up and NLO backgrounds | open (P-R9-a) |
