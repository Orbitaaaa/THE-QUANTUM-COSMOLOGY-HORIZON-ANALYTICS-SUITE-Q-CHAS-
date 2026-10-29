"""Round 6o: collider formation, interior and evaporation in one description; event generator."""
import os
import warnings

import numpy as np

from babyuniverse import pipeline6o as P
from babyuniverse import collider6l as C

warnings.simplefilter("ignore")
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "DATA")


def test_evaporation_never_closes_the_trapped_region():
    v = P.vaidya_history(6, 0.6, 3.0, "hayward", n_pts=500)
    assert v["horizons_never_meet"] and v["shell_below_inner_horizon"]
    assert abs(v["mu_minus_ext_power"] - 1) < 0.05 and abs(v["separation_power"] - 0.5) < 0.05
    assert v["affine_increment_last_decade"] < 1e-3 and 1.0 < v["affine_final"] < 5.0


def test_timeline_and_generator():
    t = P.event_timeline(6, 8.0, 0.6, 25.0, "hayward")
    assert t["t_formation_s"] < t["t_bounce_s"] * 10 and t["t_evaporation_s"] > 10 * t["t_bounce_s"] and t["daughter_persists"]
    L = C.load_lumi(os.path.join(DATA, "round6l_lumi.json"))
    gU = P.Generator(L, 100.0, 6, 8.0, 0.6, "hayward", seed=2).run(300)
    gG = P.Generator(L, 100.0, 6, 8.0, 0.6, "gr", seed=2, M_th_TeV=11.5).run(300)
    assert abs(gU["M_ext_TeV"] - 11.5) < 0.3 and gU["frac_horizon"] == 1.0
    assert np.mean(gU["f_miss"]) > 0.8 and np.mean(gG["f_miss"]) < 0.05
    assert np.mean(gU["N"]) < np.mean(gG["N"])
    assert all(r >= 0.999 * gU["M_ext_TeV"] for r in gU["remnant_TeV"])
