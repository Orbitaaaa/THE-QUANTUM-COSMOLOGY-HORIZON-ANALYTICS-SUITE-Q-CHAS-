# Round 13 — robust finite-daughter formation

Main report: [`../ROADMAP_ROUND13.md`](../ROADMAP_ROUND13.md).

| File | Content |
|---|---|
| [`CAUSAL_DOMAIN_PROOF.md`](CAUSAL_DOMAIN_PROOF.md) | Definition 13-0 (daughter region), Theorem 13-1 with proof, Corollary 13-2 (certified cycles), full-cycle window |
| [`COLLIDER_FORMATION.md`](COLLIDER_FORMATION.md) | Collider-hole kinematics, Lemma 13-3 (unbound content), Proposition 13-4, anisotropy estimate |
| [`../PAPER/paperC_daughter_cycle.pdf`](../PAPER/paperC_daughter_cycle.pdf) | Paper C v2 |
| [`../DATA/round13.json`](../DATA/round13.json), [`../FIGURES/round13.png`](../FIGURES/round13.png) | Numerical evidence |

**Reproduce:**

```
cd CODE
python scripts/round13.py          # parallel survey on 4 processes; writes DATA/round13.json and FIGURES/round13.png
python scripts/round13.py --figure-only
python -m pytest tests/test_round13.py -q
```

**Established vs new.** Established inputs: Oppenheimer–Snyder collapse, regular-core towers (Hayward, Dymnikova), the
2D-quotient causal structure, and the Round 6n Lagrangian code. New (to our knowledge): Theorem 13-1 and the certified stage events
(claim C4), the D = 10 fluid survey (C5), and the scope result for collider holes (C6). See
[`../ROUND11/publications/ORIGINALITY_AUDIT.md`](../ROUND11/publications/ORIGINALITY_AUDIT.md).

An external review brief is kept as a doc by the project lead; it is private until shared.
