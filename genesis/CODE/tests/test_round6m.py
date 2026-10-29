"""Round 6m: the global daughter spacetime of model U."""
import warnings

import numpy as np

from babyuniverse import global6m as G
from babyuniverse import unified6l as U

warnings.simplefilter("ignore")


def test_curvature_is_bounded_independently_of_the_mass():
    c = G.kretschmann_schwarzschild_check()
    assert abs(c["K"] / c["exact"] - 1) < 1e-6
    s = G.curvature_bound_scan(10, 1.0, "hayward", np.geomspace(1e-3, 1e8, 12))
    k = [r["K_max_alpha2"] for r in s["rows"]]
    assert max(k) < 250 and max(k) / min(k) < 1.05
    for C, D in ((1e4, 10), (1e8, 4)):
        fr = G.frw_curvature(C, D)
        assert abs(fr["K_max_alpha2"] / fr["de_sitter_value"] - 1) < 1e-3 and fr["psi_max_alpha"] <= 1.0


def test_ball_reaches_the_new_universe_in_finite_proper_time():
    b = G.ball_history(10, 1.0, 10.0, 0.6)
    assert b["horizon"] and b["R_bounce"] < b["r_minus"] < b["r_plus"] < b["R_max"]
    assert all(np.isfinite(b[k]) and b[k] > 0 for k in ("tau_region_I", "tau_region_II", "tau_core_to_bounce"))
    assert abs(b["tau_region_I"] + b["tau_region_II"] + b["tau_core_to_bounce"] - b["tau_half_cycle"]) < 1e-9
    assert not G.ball_history(10, 1.0, 0.5)["horizon"]
    p = G.penrose_worldline(10, 1.0, 10.0, 0.6, n=400)
    assert p["region_order"] == ["I", "II", "III", "IV", "I'"]


def test_horizons_regular_and_geodesics_complete():
    k = G.kruskal_regularity(10, 1.0, 50.0)
    for name in ("outer", "inner"):
        assert k[name]["finite_nonzero"] and abs(k[name]["limit_ratio"] - 1) < 1e-3
    g = G.geodesic_completeness(10, 1.0, 10.0)
    assert abs(g["f0"] - 1) < 1e-9 and all(x["finite"] for x in g["geodesics"])
    assert any(x["passes_centre"] for x in g["geodesics"]) and any(x["turning_point"] > 0 for x in g["geodesics"])


def test_finite_daughter():
    d = G.daughter_properties(10, 1.0, 1e4, np.pi / 3)
    assert d["forms_horizon"] and d["volume_growth"] > 1e5 and abs(d["a_bounce"] - 1) < 1e-3
    st = G.daughter_properties(4, 1.0, 2e38, np.pi / 3)
    assert st["volume_growth"] > 1e110
    small = G.daughter_properties(10, 1.0, 30.0, 0.3)
    assert small["mu"] < U.extremal_mu(10) and not small["forms_horizon"]
