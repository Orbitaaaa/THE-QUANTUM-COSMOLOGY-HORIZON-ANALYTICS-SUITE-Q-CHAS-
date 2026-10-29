# Round 6o — collider formation, interior evolution and evaporation in one description

**Task.** Describe a collider event in model U as one consistent spacetime, from the parton collision to the daughter's fate, and produce
its observables.

**Answer: yes.** For every event above threshold, one model-U spacetime contains the whole history:
1. **Formation.** A trapped surface forms within ~10⁻²⁸ s.
2. **Interior.** The partons pass through the regular centre about 10⁻²⁸ s later.
3. **Evaporation.** The hole evaporates over 10⁻²⁶–10⁻²⁵ s, approaching a zero-temperature remnant. The trapped region never closes,
   so the daughter is never released back into the parent. The bounced partons reach the boundary of the parent's chart, beyond which
   the daughter region lies, at finite affine parameter (≈ 2 l).
4. **Observables.** An event generator built on the same physics turns this into rates and distributions. The striking result: model-U
   holes made near threshold emit almost nothing (no Hawking quantum in 75–85 % of events), and the remnant carries away 90–97 % of the
   energy.

Code `pipeline6o.py`; script `scripts/round6o.py`; tests `tests/test_round6o.py` (2); data `DATA/round6o.json`; figure `FIGURES/round6o.png`.

## 1. The per-event spacetime

1. **Formation.**
   * √ŝ is drawn from the NNPDF2.3 LO luminosity times F_n πR_S².
   * The impact parameter is drawn with b² uniform below √F_n r_h.
   * The trapped mass is M = y_n √ŝ (head-on area bound, conservative).
   * A horizon forms iff M > M_ext.
2. **Interior.**
   * The partons form an ingoing null shell along v = 0.
   * The centre is regular (f(0) = 1), so the shell passes through r = 0 and continues outward inside the core.
3. **Evaporation [Proposition-N 6o-1].**
   * We use model-U Vaidya: f(v, r) = 1 − r²ψ(μ(v)/r^{D−1}), with dμ/dv = −ε_n r₊²T⁴ (brane Stefan–Boltzmann emission, g* = 106.75;
     ε ≈ 362 for n = 6, c_l = 0.6).
   * We ran 24 cases: n = 4, 6; c_l = 0.6, 0.75; μ₀/μ_ext = 1.5, 3, 10; Hayward and Dymnikova towers.
   * In all of them:
     * **μ − μ_ext ∝ v^{−p}, p = 0.98–1.03**, because T ∝ (μ − μ_ext)^{1/2} near extremality;
     * **the horizon gap r₊ − r₋ ∝ v^{−q}, q = 0.49–0.52**;
     * the horizons never meet, so the trapped region never closes and the closure lemma of Round 6l keeps the daughter;
     * the bounced partons stay inside the core, below r₋(v), and approach the degenerate horizon;
     * their affine parameter converges to **2.00–2.12 l**. Since ∫κ₋ dv diverges, the boundary of the parent's chart is reached in finite
       affine time, and the spacetime continues into the daughter region.
   * A model with complete evaporation would close the trapped region and return the matter to the parent. Model U does not.

**Timeline of one event** (n = 6, M_D = 8 TeV, c_l = 0.6, M = 25 TeV):

| Stage | Value |
|---|---|
| horizon forms (r₊/c) | 7×10⁻²⁹ s |
| partons bounce in the core | 1.1×10⁻²⁸ s |
| evaporation to 1 % above the remnant | 6×10⁻²⁶ s |
| quanta emitted | 3.6 |
| remnant | 11.5 TeV |
| daughter persists | yes |

The other benchmarks (n = 4, 7; Dymnikova tower) give the same ordering.

## 2. The event generator and what it shows

At 100 TeV, n = 6, M_D = 8 TeV, with model-U threshold M_ext = 11.5 TeV (σ ≈ 13.7 pb):

| | model U (Hayward) | model U (Dymnikova) | GR hole, same threshold |
|---|---|---|---|
| mean Hawking quanta | 0.29 | 0.17 | 2.05 |
| events with no quantum at all | 76 % | 85 % | 0 % |
| missing-energy fraction (remnant + ν) | 0.95 | 0.97 | 0.001 |
| peak temperature at 20 TeV | ~1 TeV | ~1.2 TeV | ~5 TeV |

* **Production piles up at threshold.** The cross-section falls steeply above M_ext.
* **Near threshold a model-U hole is nearly extremal.** Its temperature is far below the energy above extremality, so it emits almost
  nothing: the event is a heavy, neutral, stable remnant plus the untrapped collision energy.
* **Consequence for experiments (taken up in Round 7).** Model-U production is largely invisible near threshold. Its handles are:
  * the minority of events with a few hard quanta and large missing energy;
  * recoil against initial-state radiation (monojet with very large missing energy);
  * any charged remnants.

**Limitations, stated.**
* Holes are produced with large spin: the Myers–Perry parameter a* ≈ 3. No rotating model-U solutions are known, so the
  non-rotating treatment and the omission of a spin-down phase are approximations (P-R6n-b).
* Greybody factors are not included.
* The trapped mass is taken b-independent (conservative head-on value).

## 3. Status

| Requirement | Status |
|---|---|
| Formation, interior and evaporation in one spacetime | yes (model-U Vaidya with a regular centre) |
| Daughter survives evaporation | yes: the trapped region never closes, and the daughter's boundary is reached at finite affine parameter |
| Per-event observables from the same physics | yes (event generator) |
| Rotation and greybody factors | approximations (stated) |
