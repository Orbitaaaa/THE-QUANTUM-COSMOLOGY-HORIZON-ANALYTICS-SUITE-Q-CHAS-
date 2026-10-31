# G11-33 — Budget (structure; no invented prices)

**Rule.** No price, salary or overhead rate is stated unless verified. Every amount is **UNRESOLVED** until a quotation, a host-
institution rate sheet or a funder rule fixes it. Nothing is marked purchased.

| Line | Proposal A | Proposal B — V1 | V2 | V3 | Source of the amount |
|---|---|---|---|---|---|
| Personnel | 1 postdoc × 36 months + PI time | technician/student support (months to be set) | as V1, plus electronics engineer time | as V2, plus beam-test team | Host institution salary scale: UNRESOLVED |
| Hardware | none | BOM V1 (`../engineering/G11-22_BOM.csv`) | + BOM V2 | + BOM V3 | Vendor quotes: UNRESOLVED |
| In-kind (expected, unconfirmed) | — | sensors, amplifier boards, oscilloscope, Sr-90 source | DAQ loans | trigger scintillators, MCP-PMT, beam time | Partner letters: UNRESOLVED |
| Computing | CSC allocation (application) | workstation (existing) | — | — | Allocation decision |
| Travel | conferences, collaboration visits | — | — | beam-test campaign | Funder per-diem rules |
| Dissemination | open-access fees | — | — | — | Journal fees: UNRESOLVED |
| Contingency | — | 20 % hardware; 30 % in-house PCB | same | same | Planning convention, applied after quotes |
| Overheads | per funder rules (e.g. RCF full-cost model) | same | same | same | Funder rules: confirm in the call text |

## How the amounts get filled

1. Ask the partner lab which items they can lend. This removes lines 1, 2, 4, 7 and 19 of the BOM.
2. Request quotes for the remaining items, with a written quote reference per line.
3. Enter each amount with its quote reference in `G11-22_BOM.csv` (columns `unit_price`, `currency`, `price_source`).
4. Apply contingency and overheads from the funder's rules.
