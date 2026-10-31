# G11-23 — Experimental test protocol

## T0. Incoming inspection

* Visual inspection under a microscope; record the sensor IDs.
* Probe-station IV up to the producer's limit, in 5 V steps with a current compliance of 1 µA. Record the breakdown voltage V_bd.
* CV: gain-layer depletion knee and full depletion.
* **Pass:** V_bd above the planned operating bias with ≥ 20 V margin (AT-1).

## T1. Board-level tests (no sensor)

* Amplifier supply currents; gain and bandwidth with a network analyser, or with pulse injection through the calibration capacitor.
* Baseline noise RMS on the digitiser, with 10⁴ random triggers. **Pass:** ≤ 3 mV (AT-2).
* Channel-to-channel synchronisation with a split pulser. **Pass:** jitter ≤ 5 ps (AT-6).

## T2. Charge collection (Sr-90, single plane)

* Bias scan in 10 V steps up to V_bd − 20 V, 2000 triggers per step.
* Fit a Landau ⊗ Gaussian to the charge integral to get the MPV. Gain = MPV / MPV(PIN reference, or 0.5 fC for 50 µm unity gain).
* **Pass:** gain ≥ 15 at the chosen operating point (AT-3).

## T3. Timing (three-sensor method)

**Data:** ≥ 2000 triple coincidences per configuration.

**Reconstruction:**
1. CFD fraction scan 20–80 % in 10 % steps; choose the minimum of the pairwise widths.
2. Amplitude-walk correction (fit in 1/A).
3. Pairwise widths σ_AB, σ_BC, σ_AC from a Gaussian fit to the central 90 %.
4. Solve σ_A² = (σ_AB² + σ_AC² − σ_BC²)/2, and cyclically for B and C (`timing11.analyze_run`).

**Systematics:**
* Vary the CFD fraction ±10 %, the fit range (80–95 % central), and the amplitude cut (±30 %).
* Repeat with the plane order permuted.

**Pass:** every plane ≤ 45 ps, with total uncertainty ≤ 3 ps (AT-4).

## T4. Stability

* 8 h continuous run. Log temperature, humidity and bias current.
* Resolution in 1 h slices. **Pass:** Δσ ≤ 3 ps and offset drift ≤ 5 ps (AT-5).

## T5. Beam test (V3)

* **Survey and set-up:** survey the baseline (≤ 1 mm); measure cable delays with a pulser, and with a station-swap run if possible.
* **Momentum scan:** 1, 2 and 3 GeV/c, each with ≥ 10⁴ triggers.
* **Analysis:** TOF distribution per momentum; fit the π, K and p peaks; compare the separations with `timing11.beam_tof`.
* **Pass:** p/π at 2 GeV/c separated by ≥ 5σ, and ΔTOF consistent with the prediction within 3σ_stat (AT-7).

## Uncertainty analysis (each reported resolution)

**Statistical:** σ/√(2N) per pairwise width, propagated through the three-sensor solution.

**Systematic budget:**

| Source | Target |
|---|---|
| CFD and walk method | ≤ 3 ps |
| Temperature (from logged correlation) | ≤ 2 ps |
| Offsets and cables | ≤ 2 ps |
| Baseline survey (V3) | ≤ 3 ps at β ≈ 1 |

These are combined in quadrature, with each correlation stated.

## Test report template

1. Configuration: sensors, bias, temperature, boards, digitiser settings, firmware/software versions (git hash).
2. Raw-data location and checksums.
3. Results tables: noise, gain, σ per plane with statistical and systematic uncertainties, stability.
4. Pass/fail against AT-1 … AT-8.
5. Deviations and open issues.
6. Sign-off: operator, reviewer, date.
