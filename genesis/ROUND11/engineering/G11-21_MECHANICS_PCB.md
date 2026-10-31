# G11-21 — Mechanical and PCB engineering

Dimensioned layouts are in `FIGURES/round11_timing.png`, panels (d) and (e). This is a concept specification; production drawings
need CAD and a mechanical engineer's review.

## Bench (V1 / V2)

**Enclosure**
* Light-tight aluminium box, inner size ~160 × 90 × 120 mm (indicative), grounded to the digitiser chassis.
* Bulkhead SMA and SHV feedthroughs; a lid microswitch wired to the HV interlock (G11-24).

**Stack**
* Three carrier PCBs (one per sensor) on a vertical rail with 12 mm pitch.
* Alignment to ≤ 0.2 mm laterally, with dowel pins. The 1.3 mm pads must overlap for the β electrons.

**Source holder**
* Sealed Sr-90 source in a holder with a collimator: an aluminium or PMMA plate with a ~1 mm hole, ~5–10 mm above plane A.
* A lockable cover; storage under the radiation-safety rules.

**Carrier PCB**
* Sensor glued with conductive epoxy to the HV pad (backside). Pad bonded to the amplifier input; guard ring bonded to ground.
* Amplifier on the same board (single-channel concept), or a short coax to the amplifier board.

**Thermal**
* Room-temperature operation by default.
* Optional cooled operation: a Peltier stage under the stack, a copper spreader, dry air or N₂ flush, and dew-point and temperature
  sensors (Pt100 or NTC) read out with the DAQ.

## Beam test (V3)

* **Rail.** Two stations on a rigid rail (aluminium profile). Baseline 1.000 m between sensor planes, surveyed to ≤ 1 mm.
* **Each station.** An enclosure like the bench box; a sensor carrier; upstream and downstream plastic-scintillator trigger paddles
  read by PMTs or SiPMs; an optional MCP-PMT reference at station A.
* **Cabling.** Matched-length low-loss coax to one digitiser, with cable delays measured.
* **Mounting.** Adapter to the beam-line table, which depends on the facility.

## PCB engineering notes

* Board stack and impedance control: a 4-layer board with continuous ground and a 50 Ω microstrip to the SMA, following board-house
  rules.
* HV clearance to be set by IPC-2221-type guidelines for ~300 V, and verified by the board house.
* Test points for bias, supply current and a calibration-pulse injection (series capacitor).
* Labelled sensor IDs, and the orientation marked for bonding.

## Shielding and grounding

* A single-point ground at the digitiser.
* Enclosure walls continuous; copper tape on seams if needed.
* HV return path separate from the signal ground, joining at the HV-filter ground.
