# G11-18 — LGAD sensor selection

## What we need

* 50 µm low-gain avalanche diodes with gain 15–30 at room temperature.
* Single pads of ~1.3 × 1.3 mm², wire-bondable to a single-channel amplifier board (G11-19).
* Bias in the few-hundred-volt range.

## Candidate sources (verified to have produced LGADs for timing R&D)

| Producer | Evidence | Typical access route | Status for us |
|---|---|---|---|
| Hamamatsu Photonics (HPK) | HGTD/ETL prototype characterisation; HPK 3.1 LGADs studied with a proton beam, β source and probe station (arXiv:2104.08369) | Production runs for collaborations; small quantities via collaboration partners | UNRESOLVED: whether single-pad samples can be bought by a non-collaboration lab, and at what price/lead time; needs a direct enquiry |
| FBK (Fondazione Bruno Kessler) | LGAD development since RD50; CMS ETL comparisons | Research agreements and shared runs | UNRESOLVED (enquiry) |
| IMB-CNM (Barcelona) | LGAD concept and HGTD layout originated in a CNM RD50 common project (arXiv:1804.00622) | Research collaboration | UNRESOLVED (enquiry) |
| BNL | LGAD development (RD50) | Collaboration | UNRESOLVED |
| IHEP-IME, USTC-IME | Listed as HGTD sensor vendors in a 2025 HGTD talk | Collaboration | UNRESOLVED |

**Practical recommendation.** The realistic route to the first sensors is collaboration: borrow or receive characterised spare pads
from a group in the RD50 / DRD3 community. HIP or a university detector group would be the natural first contact (G11-35). Buying
directly is unlikely to be quick for a small order. I could not verify any catalogue price.

## Compatibility checks before acceptance

* Pad size and capacitance match the amplifier input. A 1.3 mm pad at 50 µm has a few pF; the noise and rise-time model in
  `timing11` assumes this class.
* Bond pads and guard ring are compatible with the carrier PCB. Guard ring grounded; backside contact through conductive epoxy to the
  HV plane.
* Breakdown and operating bias are known from the producer's data, or measured (AT-1).
* Bias supply: ≥ 300 V, current-limited to ≤ 10 µA per channel during scans. A per-channel series resistor of order 1 MΩ, with
  filter capacitors rated ≥ 500 V, is set during board design.

## Unirradiated vs irradiated

The demonstrator uses unirradiated sensors, because the remnant search is a timing question, not a radiation-hardness question.
Irradiation studies are explicitly out of scope for V1–V3.
