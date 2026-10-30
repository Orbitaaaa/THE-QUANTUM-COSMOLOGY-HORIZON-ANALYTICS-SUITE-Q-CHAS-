"""Round 9: the experimental design."""
import os
import warnings

import numpy as np

from babyuniverse import design9 as X
from babyuniverse import collider6l as C

warnings.simplefilter("ignore")
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "DATA")


def test_statistics_detector_and_decision_tree():
    assert abs(X.z_asimov(100.0, 1e6) - 100.0 / 1e3) < 1e-3          # s << b: Z -> s / sqrt(b)
    assert X.z_asimov(10.0, 1.0) > X.z_asimov(10.0, 10.0) > 0
    d = X.detector_requirements(beta=0.39, path_m=10.0)
    assert abs(d["tof_delay_ns_at_10m"] - 52.2) < 0.5
    assert len(X.decision_tree()) == 6
    rows = [dict(sqrt_s=1.0, reachable=True, L_available_ab=10.0, M_ext_TeV=m, L_discovery_ab=L)
            for m, L in ((10.0, 1.0), (20.0, 100.0))]
    assert abs(X.reach(rows, 1.0) - 15.0) < 1e-9                     # log-linear midpoint


def test_design_point_at_100_TeV():
    L = C.load_lumi(os.path.join(DATA, "round6l_lumi.json"))
    d = X.design_point(L, 100.0, 6, 8.0, 0.66, "hayward", n_events=300, seed=3)
    assert d["reachable"] and 18.0 < d["M_ext_TeV"] < 24.0
    assert np.isfinite(d["L_discovery_ab"]) and d["L_discovery_ab"] < 1.0 < 30.0
    assert d["sigma_sel_fb"] > d["bkg_fb"]                              # signal-dominated after selection
    closed = X.design_point(L, 27.0, 6, 8.0, 0.70, "hayward", n_events=50)
    assert not closed["reachable"]                                      # beyond the kinematic limit
