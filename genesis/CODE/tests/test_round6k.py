"""Round 6k: the two doors (Y1 black-hole interiors that become daughters; Y2 black holes at collider energies), the conditional machine,
and natural bridges within reach."""
import warnings

import numpy as np

from babyuniverse import doors6k as D

warnings.simplefilter("ignore")


def test_ecsk_bounce_scales():
    b = D.ecsk_bounce()
    assert 3e58 < b["rho_b_kg_m3"] < 5e58 and b["rho_over_planck"] < 1e-37          # far below Planck density: semiclassical
    m = D.ecsk_min_hole(b["rho_b_kg_m3"])
    assert 3e10 < m["M_min_kg"] < 6e10
    assert 1e17 < D.ecsk_min_hole(1e45)["M_min_kg"] < 5e17


def test_ecsk_interior_bounces_and_inflates_with_vacuum():
    assert abs(D.ecsk_inflation_criterion() - 1.5 * np.sqrt(3)) < 1e-12
    r0 = D.ecsk_interior(1.0, 1e-6)
    assert "cycles" in r0["fate"] and abs(r0["a_min"] / r0["a_min_analytic"] - 1) < 1e-3 and abs(r0["a_max_after"] - 1) < 1e-3
    assert "inflates" in D.ecsk_interior(1.0, 1e-6, 2.0)["fate"]
    assert "cycles" in D.ecsk_interior(1.0, 1e-6, 3.0)["fate"]


def test_limiting_curvature_junction_in_D_dimensions():
    for n in (0, 4, 6):
        j = D.fmm_D(n, 3.0)
        assert j["valid"] and j["r0_over_l"] > 1
        assert abs(j["r0"] - (3.0 ** (n + 1)) ** (1 / (n + 3))) < 1e-12
        assert not D.fmm_D(n, 0.9)["valid"]


def test_conditional_machine_thresholds():
    for n in (4, 5, 6, 7):
        m = D.add_daughter_machine(n, 8.0)
        assert 3.5 < m["M_Y_over_MD"] < 5.0 and m["lifetime_over_interior"] > 5
    assert D.add_daughter_machine(6, 8.0, safety=3.0)["M_Y_over_MD"] > 5e3
    assert D.collider_reach(6, 8.0)["reachable"] and not D.collider_reach(6, 20.0)["reachable"]


def test_regular_core_remnant_signature():
    r = D.regular_core_temperature(6, 1.0, 120)
    T = [x["T"] for x in r["rows"]]
    assert T[0] < 0.05 * max(T) and max(T) > T[-1]                     # T rises from zero at the remnant, peaks, then falls like GR
    assert abs(r["T_large_mu"] / r["gr_T_at_large_mu"] - 1) < 0.1


def test_pbh_partner_and_travel():
    p = D.pbh_partner(5.0)
    assert 80 < p["partner_vac_GeV"] < 90 and 3.0 < p["ew_partner_M_earth"] < 3.6 and 7 < p["travel_years_at_0p001c"] < 9
    assert len(D.y1_theories()) == 5
