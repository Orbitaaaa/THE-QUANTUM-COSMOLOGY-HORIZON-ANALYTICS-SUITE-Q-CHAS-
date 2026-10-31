# G11-17 — GENESIS Timing Demonstrator, Version 1: detector requirements

## Purpose

Show that we can measure the arrival time of single charged particles with silicon LGAD sensors, at the precision that a search for
slow heavy particles needs. The second step is to measure their velocity by time of flight over a 1 m baseline.

The physics motivation is a heavy remnant arriving late. In Round 10, model-U remnants have median β = 0.13 and arrive +32 ns late at
1.5 m. But the instrument is generic: its results stand regardless of which gravitational theory is correct.

## Measurement principle

| Version | Setup | What it measures |
|---|---|---|
| V1 / V2 (bench) | Three LGAD sensors stacked ~12 mm apart, crossed by electrons from a Sr-90 β source | Pairwise time differences; the three-sensor method gives each plane's resolution with no external reference |
| V3 (beam test) | Two LGAD stations 1.000 m apart in a mixed hadron beam, with scintillator trigger and an optional MCP-PMT reference | Time of flight and β for pions, kaons and protons at 0.5–3 GeV/c |

Why the bench and beam versions differ: cosmic rays cannot be used at a 1 m baseline with 1.3 mm pads. The coincidence rate is
~2×10⁻⁵ per day (`timing11.cosmic_coincidence_rate`).

## Performance specification

Model values come from the calibrated LGAD simulation (`timing11`). They are estimates, and the acceptance tests decide.

| Quantity | Requirement | Goal | Basis |
|---|---|---|---|
| Single-plane time resolution (MIP, room temperature) | ≤ 45 ps | ≤ 35 ps | Model: 30 ps at gain 20. Published 1.3 mm HGTD pads: ~40 ps at gain 20, ~27 ps at gain 50 (arXiv:1804.00622) |
| Two-plane TOF resolution | ≤ 64 ps | ≤ 50 ps | √2 × single plane |
| β resolution at β = 0.3, L = 1 m | ≤ 0.2 % | ≤ 0.15 % | σ_β = β² c σ_TOF / L |
| p/π separation at 2 GeV/c, L = 1 m | ≥ 5σ | ≥ 6.8σ | ΔTOF = 341 ps |
| LGAD gain at operating bias | ≥ 15 | 20–30 | Model: resolution degrades below gain ~10 |
| Front-end noise (RMS, at digitiser input) | ≤ 3 mV | ≤ 1.5 mV | Model: 34 ps at 3 mV, 42 ps at 5 mV |
| Analogue bandwidth | 1–1.5 GHz | 1.5 GHz | Model: 2.5 GHz admits more noise; 0.5 GHz slows the edge |
| Sampling rate | ≥ 5 GS/s | 5–10 GS/s | Model: 2 GS/s gives 44 ps |
| Amplitude-walk correction | CFD, offline, 30–60 % fraction scan | — | Standard |
| Clock and channel synchronisation | ≤ 5 ps residual | ≤ 2 ps | Common digitiser clock; calibrated cable delays |
| Event rate, bench | ≥ 5 triggers/s | 15/s | Sr-90 (37 kBq-class) at 10 mm: ~15/s per 1.3 mm pad |
| Statistics per measurement point | ≥ 2000 selected events | — | Width uncertainty σ/√(2N) ≈ 0.8 ps |

## Operating conditions

| Item | Specification |
|---|---|
| Temperature | +20 °C nominal. Optional cooled operation to −20 °C needs dry air or N₂ flush and dew-point monitoring. |
| Bias | Up to ~300 V DC, current-limited (µA range). The actual operating point comes from IV/CV and gain scans of the delivered sensors. |
| Light and electromagnetic environment | Light-tight, grounded enclosure; filtered HV; 50 Ω signal path. |
| Radiation | Sealed Sr-90 source (bench). Beam test: the facility's rules apply. |

## Sensor geometry

* 50 µm active thickness, single 1.3 × 1.3 mm² pads (the ATLAS HGTD / CMS ETL sensor class).
* Bench: 3 sensors on separate carrier boards at 12 mm pitch.
* Beam test: 1–2 sensors per station, plus a scintillator trigger.

## Uncertainty targets

| Source | Target |
|---|---|
| Statistical width uncertainty per point | ≤ 1 ps |
| CFD-fraction and walk systematics | ≤ 3 ps |
| Temperature drift over a run | ≤ 2 ps (logged) |
| Cable and offset calibration | ≤ 2 ps |
| Baseline survey (V3) | ≤ 1 mm, i.e. ≤ 3 ps for β ~ 1 |

## Acceptance tests (pass/fail)

* **AT-1 IV/CV.** Breakdown above the operating bias with ≥ 20 V margin. Gain-layer depletion voltage measured.
* **AT-2 Noise.** RMS noise ≤ 3 mV at the operating point, all channels.
* **AT-3 Charge and gain.** Sr-90 MPV charge measured, giving gain ≥ 15.
* **AT-4 Timing.** Three-sensor method: every plane ≤ 45 ps at ≥ 2000 events, with total uncertainty ≤ 3 ps.
* **AT-5 Stability.** Over 8 h the resolution changes by ≤ 3 ps and the offsets drift by ≤ 5 ps.
* **AT-6 Synchronisation.** Pulser split into all channels gives a channel-to-channel jitter ≤ 5 ps.
* **AT-7 (V3) Time of flight.** Measured ΔTOF(p − π) at 2 GeV/c agrees with the prediction within 3σ_stat, and the separation is ≥ 5σ.
* **AT-8 Safety.** HV interlock and source-handling checks are signed off (G11-24) before first powering.
