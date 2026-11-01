"""Round 13: causal-domain certification, fluid cone tracer validated against the exact dust solution, bound vs unbound content."""
import json
import os
import warnings

import numpy as np

from babyuniverse import robust13 as R
from babyuniverse import daughter12 as Dg

warnings.simplefilter("ignore")
DATA = os.path.join(os.path.dirname(__file__), "..", "..", "DATA")


def test_domain_condition():
    assert R.in_domain(0.5, 0.2, 0.6, 0.3) and not R.in_domain(0.8, 0.2, 0.6, 0.3)


def test_full_cycle_certified():
    c = R.certify_case(10, 1.05, 1.4, "dymnikova")
    ev = c["events"]
    assert ev["second bounce"]["certified_at_centre"] and ev["second bounce"]["centre_margin"] > 1e3 * c["numerical_error"]
    assert ev["maximum"]["chi_in_D"] > 0.8 and 0 < ev["maximum"]["volume_fraction_in_D"] < 0.1
    h = R.certify_case(10, 1.2, 1.0, "hayward")
    assert h["events"]["maximum"]["certified_at_centre"] and not h["events"]["second bounce"]["certified_at_centre"]


def test_fluid_cone_tracer_matches_exact_dust():
    st = R.make_d10_star(1.05, 1.4, "uniform", 0.0, "dymnikova", N=60)
    r = R.run_cycle(st, track_x=(0.025, 0.1), t_max=6)
    ev = Dg.evolve_ball(Dg.cdd(10, 1.05, 1.4, "dymnikova")["C"], 10, "dymnikova")
    tau_cone_exact = ev["tau_b2"] + (Dg.cdd_fast(10, 1.05, 1.4, "dymnikova")["cone_budget"] - (ev["eta_b2"] - ev["eta_b"])) * Dg.turning_points(
        Dg.cdd(10, 1.05, 1.4, "dymnikova")["C"], 10, "dymnikova")[0]
    assert abs(r["t_cone_centre"] / tau_cone_exact - 1) < 0.01
    centre = [e["event"] for e in r["elements"]["0.025"]["events"] if e["in_D"]]
    assert centre[:3] == ["bounce", "maximum", "second bounce"]
    assert "second bounce" not in [e["event"] for e in r["elements"]["0.100"]["events"] if e["in_D"]]


def test_unbound_content_does_not_bounce():
    assert R.unbound_lemma(1.01, 10, 1.05, "dymnikova")["no_turning_point"]
    o = R.open_frw_no_bounce(50.0, 10, "hayward")
    assert o["reaches_zero"] and o["psi_end_alpha"] <= 1.0 + 1e-12 and o["rho_end"] > 1e30


def test_survey_and_collider_from_data():
    path = os.path.join(DATA, "round13.json")
    if not os.path.exists(path):
        return
    d = json.load(open(path))
    s = d["survey"]
    uni = [x for x in s if x["profile"] == "uniform" and x.get("status") == "ok"]
    assert all(x["verdict"] in ("full cycle", "maximum + recollapse") for x in uni)   # homogeneous core: cycle survives pressure, modes
    dym = [x for x in uni if x["kind"] == "dymnikova" and x["mu_over_ext"] == 1.05]
    assert sum(x["verdict"] == "full cycle" for x in dym) >= 0.7 * len(dym)
    peaked = [x for x in s if x["profile"] in ("smooth", "gaussian") and x.get("status") == "ok"]
    assert sum(x["verdict"] == "expands, ends" for x in peaked) > 0.8 * len(peaked)
    assert all(x.get("black_hole") for x in s)
    cc = d["collider_content"]
    assert cc["mu_over_ext_percentiles"][2] < 1.1 and cc["frac_below_1p1"] > 0.6
