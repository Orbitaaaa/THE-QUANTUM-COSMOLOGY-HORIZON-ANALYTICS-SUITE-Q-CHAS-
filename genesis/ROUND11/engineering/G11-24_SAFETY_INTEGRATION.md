# G11-24 — Laboratory safety and integration

These requirements must be confirmed with the host institution's safety officers. They are an engineering starting point, not
approval.

## High voltage (sensor bias up to ~300 V DC)

* Use HV supplies with a hardware current limit (≤ 10–100 µA trip) and remote interlock input.
* A lid microswitch on the enclosure, wired to the HV interlock: HV inhibited when the lid is open.
* SHV connectors only on the HV lines, with no exposed conductors. HV cables labelled. A bleeder or discharge path on the filter
  capacitors.
* Operating rule: change sensors or boards only with HV off and discharged. Two-person rule for first power-up.
* PCB clearances sized for the bias (G11-21) and checked at design review.

## Radioactive source (Sr-90 / Y-90 β emitter, sealed)

* A sealed Sr-90 source of tens of kBq is above typical exemption levels. Possession and use need authorisation under national
  rules; in Finland, radiation use is supervised by STUK. The applicable licence or notification, and the responsible person, must
  be confirmed with the host institution's radiation safety officer before acquiring or using a source. Using an existing
  institutional source under its licence is the preferred route.
* Handling:
  * tongs and shielding (PMMA for β, plus a lead layer for bremsstrahlung if required);
  * a locked storage container;
  * a source logbook and periodic leak tests per licence;
  * dosimetry as required.

## Electrical and electronic

* LV supplies current-limited; fused distribution; ESD protection (wrist straps, grounded mats) for sensor handling.
* Grounding plan reviewed. No floating chassis.

## Thermal and cooling (optional)

* Condensation risk below the dew point: dry air or N₂ flush, a dew-point sensor, and an interlock that switches off the cooling
  and bias if humidity rises.
* N₂ in a closed room: assess the oxygen-depletion risk for the flow rate used.

## Beam-test integration (V3)

* Follow the facility's access, radiation and electrical safety training and rules.
* Equipment declarations: HV, sources (none brought to beam), electronics. Cable routing and fire load per facility rules.

## Design-review documentation

**Checklist for the review before first power:**
1. Schematics and layouts, with HV clearances marked.
2. Interlock diagram and test record.
3. Source licence reference and handling procedure.
4. Risk assessment (HV, radiation, cryogenic or N₂, ESD).
5. Operating procedures and emergency contacts.
6. Acceptance-test plan (G11-23).
7. Sign-off by the responsible scientist and the safety officer.
