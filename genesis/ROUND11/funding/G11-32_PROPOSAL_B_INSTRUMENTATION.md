# G11-32 / Proposal B — Experimental instrumentation: GENESIS Timing Demonstrator V1

**Goal.** Build and characterise a two-plane LGAD timing system:
* ≤ 45 ps per plane, goal 35 ps;
* a bench Sr-90 three-sensor measurement;
* then a 1 m time-of-flight measurement in a mixed hadron beam (p/π separation ≥ 5σ at 2 GeV/c).

The specification is complete in `ROUND11/engineering/` (G11-17 … G11-24).

## Stages and decision points

| Stage | Content | Exit criterion | Variant (G11-22) |
|---|---|---|---|
| S1 | Borrowed scope, 3 sensors, community amplifiers, institutional Sr-90 | AT-1 … AT-4 pass (σ ≤ 45 ps per plane) | V1 minimal |
| S2 | DRS4 readout, multichannel HV, in-house front-end, monitoring | AT-5 stability, AT-6 synchronisation | V2 instrumented |
| S3 | Second station, triggers, beam test at 1 m | AT-7 TOF separation | V3 beam test |

## Why it matters

Heavy, slow, possibly charged particles are predicted by many models, including the model-U remnants (median β = 0.11). They are
identified by time of flight, which needs 30–50 ps timing. The demonstrator trains a team and validates an analysis chain
(`timing11.py`, tested on simulation) on real hardware. It is useful whatever happens to Proposal A.

## Budget

Structure only. Every price is UNRESOLVED until quotes are obtained (see `G11-33_BUDGET.md`). The largest cost items (sensors,
amplifiers, beam time) are expected to come in kind through collaboration. That must be confirmed, not assumed.

## Team and facilities (to be confirmed)

A host lab with a probe station, wire bonder and radiation-source licence. A partner timing-detector group. Beam-time allocation by
proposal.

## Safety

HV interlocks, sealed-source licensing under the host's radiation-safety officer (STUK supervision in Finland), cooling and
condensation interlocks, design review before first power (G11-24).
