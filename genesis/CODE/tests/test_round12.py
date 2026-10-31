"""Round 12: the Cauchy-development daughter (timing law, recollapse and full cycle inside D), the core lemma, and loopholes 1-2."""
import json
import os
import warnings

import numpy as np

from babyuniverse import daughter12 as Dg

warnings.simplefilter("ignore")
DATA = os.path.join(os.path.dirname(__file__), "..", "..", "DATA")


def test_ode_and_quadrature_agree():
    for args in [(10, 1.05, 1.4, "dymnikova"), (10, 1.2, 1.0, "hayward"), (5, 20.0, 0.5, "hayward")]:
        a, b = Dg.cdd(*args), Dg.cdd_fast(*args)
        assert abs(a["eta_cone_minus_eta_b"] - b["cone_budget"]) < 1e-4
        assert abs(a["eta_max_minus_eta_b"] - b["half_period"]) < 1e-4


def test_de_sitter_timing_law():
    r = Dg.cdd(5, 1e4, 0.5)
    assert abs(r["eta_cone_minus_eta_b"] - np.pi / 2) < 1e-3          # the cone is the core's de Sitter horizon
    assert r["untrapped_initially"] and r["a_b"] >= 1.0 - 1e-9
    assert 1.6 < r["a_cone_over_C_root"] < 1.85 and r["efolds_accelerating"] > 3


def test_recollapse_and_full_cycle_inside_D():
    b = Dg.cdd(10, 1.2, 1.0, "hayward")
    assert b["untrapped_initially"] and b["recollapse_begins_inside_D"] and not b["full_cycle_inside_D"]
    c = Dg.cdd(10, 1.05, 1.4, "dymnikova")
    assert c["untrapped_initially"] and c["full_cycle_inside_D"]
    assert c["eta_b2_minus_eta_b"] < c["eta_cone_minus_eta_b"]
    a = Dg.cdd(5, 20.0, 0.5, "hayward")
    assert not a["recollapse_begins_inside_D"] and a["a_cone_over_a_b"] > 5


def test_core_lemma():
    for D in (5, 10):
        for k in ("hayward", "dymnikova"):
            assert Dg.core_monotone(D, 1.05, k)["f_decreasing"]


def test_maps_and_loopholes_from_data():
    d = json.load(open(os.path.join(DATA, "round12.json")))
    assert d["maps"]["D10_dymnikova"]["frac_full_inside"] > 0.01
    assert d["maps"]["D10_hayward"]["frac_full_inside"] == 0.0
    assert d["maps"]["D10_hayward"]["frac_max_inside"] > 0.5 > d["maps"]["D5_hayward"]["frac_max_inside"]
    L1 = d["L1_compact_trapped_region"]
    assert all(x["no_event_horizon"] for x in L1)
    assert L1[2]["blueshift_exponent"] > L1[0]["blueshift_exponent"]       # RSET growth rises with the trapped lifetime
    assert 0.4 < d["L2_attractor"]["kappa_minus_power_q"] < 0.6
    assert d["L2_finite_time_extremality"][0]["blueshift_bounded"] and not d["L2_finite_time_extremality"][2]["blueshift_bounded"]
