# G11-19 — Front-end electronics

## Architecture (follows published LGAD single-channel readout)

```
LGAD pad --(wire bond)--> [stage 1: inverting transimpedance amplifier, single SiGe/RF transistor, ~2 GHz]
        --> [stage 2: broadband RF gain block, e.g. Mini-Circuits GALI-52+ class] --> SMA, 50 ohm --> digitiser (DRS4 / scope)
HV --> [series R + RC filter] --> sensor backside        (guard ring to ground)
```

**Published precedent.** The UCSC single-channel board is a single-transistor common-emitter inverting transimpedance amplifier. One
description gives ~29 dB at 1.9 GHz; another gives a gain of ~10 with ~2 GHz bandwidth. Supplemented by a GALI-52+ stage it gives
≈ 5 kΩ total transimpedance, with good board-to-board uniformity (arXiv:2104.08369).

Our simulation uses R_T = 5 kΩ and a second-order response at 1–1.5 GHz.

## Noise and jitter budget (model `timing11`, 50 µm, gain 20)

| Contribution | Value | Note |
|---|---|---|
| Signal amplitude (MIP, mean) | ~76 mV | 5 kΩ, gain 20 |
| Rise time (10–90 %) | ~0.5 ns | Hole drift time across 50 µm plus bandwidth |
| Jitter, analytic σ_n/(dV/dt) | ~12 ps | At 1.5 mV noise |
| Landau (ionisation-shape) term | ~29 ps | Calibrated to published 50 µm performance |
| Total, model | ~30 ps | Measured HGTD pads: ~40 ps at gain 20, so the model is ~25 % optimistic |

From the design scan:
* **Noise.** 3 mV costs about +4 ps and 5 mV about +12 ps. Hence the requirement noise ≤ 3 mV.
* **Gain.** Below ~10 the jitter term dominates. Hence gain ≥ 15.
* **Bandwidth.** The optimum is about 1–1.5 GHz. At 2.5 GHz more noise is admitted (+4 ps).
* **Sampling.** 5 GS/s is sufficient; 2 GS/s costs about +14 ps.

## Design options

| Option | Description | Pros | Cons / status |
|---|---|---|---|
| A | Use an existing UCSC-type board design | Proven performance | UNRESOLVED: availability of boards or design files must be requested from the authors or groups |
| B | Own two-stage board, following the published topology | Full control; documentable | Needs RF layout expertise and iteration; the transistor part number must be chosen and verified in simulation |
| C | Commercial broadband amplifiers only (two cascaded gain blocks) | Fast to assemble | Usually worse noise for small charges; needs testing |

**Recommendation.** Request option A from the community first. Develop option B in parallel as the documented in-house design.

## Schematic-level design notes (option B)

* **Input.** Short bond wires to a 50 Ω microstrip. A shunt bias resistor to ground at the transistor input sets the DC point.
* **Stage 1.** Common-emitter with a feedback resistor setting the transimpedance (hundreds of Ω). Emitter degeneration for
  stability; RF chokes on the supply.
* **Stage 2.** A packaged broadband gain block. It must be unconditionally stable, have a noise figure of a few dB, and ~20 dB gain
  at 1 GHz.
* **HV.** Series resistor of order 1 MΩ, then an RC low-pass (≥ 500 V rated capacitor) at the sensor backside.
* **Layout.** Continuous ground plane, a 50 Ω controlled-impedance trace to the SMA, and component placement minimising input
  capacitance. A separate HV island with ≥ 2 mm clearance at 300 V (to be confirmed against the board-house rules).
* **Verification before fabrication.** Simulate gain, bandwidth and stability in SPICE with vendor models. The simulation tool and
  models are to be chosen; they were not available in this environment.

## Open items

* Final transistor and gain-block part numbers are UNRESOLVED. GALI-52+ is a published example of the second stage.
* SPICE simulation with vendor models.
* PCB stack-up and board-house rules.
