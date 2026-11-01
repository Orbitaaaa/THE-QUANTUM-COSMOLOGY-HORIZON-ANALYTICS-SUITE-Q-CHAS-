# G11-34 — Funding strategy

Checked by web search on **31 October 2026**.
* **Verified:** the date appears in an official notice or in several consistent listings.
* **Unverified:** sources conflict, or the date is an inference from past rounds.

Every date must be confirmed against the funder's call text before applying.

## Eligibility constraint, which decides everything else

Almost every instrument below requires the applicant to be employed by, or hosted at, a research organisation. The large ones (RCF,
ERC) also require a PhD-holding principal investigator. Unless the project lead already holds such a position, the realistic route is:

1. A host group, such as a timing-detector group at a Finnish university or HIP, whose staff scientist is PI or co-PI.
2. Proposal B (instrumentation) through that group.
3. Proposal A (fundamental physics) as a theory collaboration with a university group.

## Instruments

| Instrument | Fit | Deadline / window | Status | Source |
|---|---|---|---|---|
| Research Council of Finland, Academy Project (autumn 2026 call) | A (theory) and B (instrumentation); Finnish host and PhD PI required | Opens 14 Oct, closes **11 Nov 2026** (marked provisional, "TBC", by RCF); scheme reform announced for 2026 | Verified (provisional) | [RCF advance information](https://www.aka.fi/en/about-the-rcf/whats-new/press-releases/2026/advance-information-on-rcf-call-schedules-for-autumn-2026/) |
| ERC Consolidator Grant 2027 | A, if a PI meets the career-stage window | Opened 24 Sep 2026; closes **12 Jan 2027**, 17:00 Brussels; one proposal per applicant under WP2027; eligibility window reported differently by sources | Deadline verified (several sources); eligibility unverified | [Göttingen internal deadlines](https://uni-goettingen.de/de/document/download/7a617ff94e4354c1b0838fc8d923f5d2.pdf/ERC_CoG_Internal_Deadlines_Call_2027.pdf), [euagenda](https://euagenda.eu/events/2027/01/12/erc-consolidator-grant-erc-2027-cog) |
| ERC Starting Grant 2027 | A | Closed (14 Oct 2026) | Verified closed | (Round 10 check) |
| EIC Pathfinder Open 2027 | B only if framed as breakthrough timing technology. Weak fit; needs a consortium | 2026 round closed 12 May 2026; 2027 date not yet published (WP2027 expected November 2026) | Unverified | [Accelopment](https://accelopment.com/blog/funding/eic-pathfinder-open-2026-limited-funding-high-competition-and-next-steps/) |
| Finnish Cultural Foundation (SKR) | Personal grants, A or B | Single annual round, **15–30 Oct 2026**: closed yesterday. Next round expected October 2027 | Verified | [SKR notice](https://www.sttinfo.fi/tiedote/71844197/kulttuurirahasto-siirtyy-yhteen-hakuaikaan?lang=fi) |
| Magnus Ehrnrooth Foundation | Physics grants (individuals and groups) | Conflicting: 1–31 January (listing checked Oct 2026) vs 15 Oct–30 Nov (2022, 2024 rounds) | **Unverified**: check the foundation portal now | [European Funding Guide](https://www.european-funding-guide.eu/fi/node/74456), [research.fi](https://research.fi/en/results/funding-call/13370) |
| Väisälä Fund (Finnish Academy of Science and Letters) | Physics (doctoral, postdoc pool, project grant for early tenure-track) | 2026 round 14 Aug–14 Sep 2026 (closed); 2027 expected Aug–Sep (inference) | 2026 verified; 2027 unverified | [acadsci.fi project grant](https://acadsci.fi/?p=12649), [HIP news](https://www.hip.fi/?p=99411) |
| Business Finland, Research to Business | Commercialisation by research organisations. Poor fit unless the timing electronics has a market case | Continuous (per earlier check) | Unverified for 2027 | — |
| AIDAinnova successor / DRD3 common funds | Sensor access and test-beam support for B | Unknown | **Unverified** | — |
| Facility beam time (CERN PS/SPS, DESY) | B, campaign V3 | Facility proposal deadlines | Unverified; check each facility | — |

## Recommended sequence

* **Now (Nov 2026):** if a host and PI exist, RCF Academy Project, B first (concrete hardware, 11 Nov 2026 provisional). Otherwise
  start partner discussions (G11-35).
* **Before 12 Jan 2027:** ERC CoG only if an eligible PI wants to lead A. Not recommended for a first proposal from this project.
* **Dec 2026 – Jan 2027:** Magnus Ehrnrooth, if the January window is confirmed.
* **Aug–Oct 2027:** Väisälä and SKR rounds; resubmission to RCF.
* **Parallel:** ask the partner group about sensor loans (RD50/DRD3) and beam-time proposals. These cost no money and remove the
  largest unresolved BOM lines.

## Two separate proposals

* **Proposal A (fundamental physics):** `PROPOSAL_A_PHYSICS.md`.
* **Proposal B (experimental instrumentation):** `G11-32_PROPOSAL_B_INSTRUMENTATION.md`.

They are written to stand alone. B does not depend on A being right: the timing demonstrator is useful for any slow-particle search.

## Update, 1 November 2026 (Round 13)

**Research Council of Finland Academy Project:** re-checked. Opens 14 Oct and closes **11 Nov 2026**, still marked provisional (TBC)
by RCF; aggregators give 23:59 Finnish time. Source:
[RCF advance information](https://www.aka.fi/en/about-the-rcf/whats-new/press-releases/2026/advance-information-on-rcf-call-schedules-for-autumn-2026/).
This is realistic only if a host organisation and an eligible PI are already in place. Otherwise target the next round (autumn 2027)
and use the time for partner agreements.

**Actions for the next 10 days**, if a host exists:
1. Confirm the PI and host (eligibility is limited to Finnish research organisations).
2. Use Proposal B (instrumentation) as the core text, with Proposal A as a theory work package.
3. Attach the reproducibility package (Paper C evidence: `ROUND13/`, `DATA/round13.json`) as preliminary results.

**Detector track:** unchanged, with all prices UNRESOLVED. Next step: request sensor loans and quotes through the partner group.
