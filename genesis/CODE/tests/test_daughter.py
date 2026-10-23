"""Roadmap Stage 2: orientation data, Casoratian charges, the regularized daughter edge, the clock fold."""
import warnings

import numpy as np

from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse import daughter as D

warnings.simplefilter("ignore")

MODELS = [ThinWallModel(1.0, 1.0, 0.1545), model_from_dimensionless(0.6, 2.4), model_from_dimensionless(-0.4, 4.2)]


def test_orientation_is_boost_invariant_configuration_data():
    r = D.orientation_is_configuration_data(-1.3, 0.5, np.linspace(-5, 5, 41))
    assert r["invariant"] and r["momentum_sign_changes"]
    assert np.allclose(r["interval"], r["f"])


def test_casoratian_charges():
    for n in (-1, 0, 1, 2):
        c = D.charge_of_wave(0.7, n, 0.05)
        assert c["scaled_spread"] < 1e-12 and c["periodicity_defect"] < 1e-12
        assert abs(c["scaled"][0] - c["expected"]) < 1e-10


def test_daughter_edge_is_removable():
    for m in MODELS[:2]:
        rep = D.edge_report(m)
        Pt = [r["P_reg"] for r in rep]
        assert abs(rep[2]["P_reg"] - rep[3]["P_reg"]) < 1e-4            # continuous across R = 1/chi
        assert abs(rep[2]["P"]) > abs(rep[0]["P"]) + 4 * m.chi**-1 * 0.5 * np.log(10) - 0.5   # bare P ~ (R/2G) ln f_in
        assert max(abs(r["residual"]) for r in rep) < 1e-10
        assert max(r["H_minus_M"] for r in rep) < 1e-10
        assert all(r["dt_dtau"] < 0 and np.isfinite(r["dt_dtau"]) for r in rep)


def test_clock_folds_back_on_every_daughter():
    for m in MODELS:
        c = D.clock_orientation(m)
        assert c["parent_future_directed"] and c["daughter_past_directed"]
        assert D.flux_balance(m, 5.0)["klein"]


def test_global_time_increases_along_parent_and_daughter():
    from babyuniverse.model import model_from_dimensionless
    for m in (model_from_dimensionless(0.6, 2.4), model_from_dimensionless(0.3, 5.0), ThinWallModel(1.0, 1.0, 0.1545)):
        tp = m.turning_points()
        R1, R2, rh = tp[0], tp[-1], m.horizons_out()[0]
        for expanding in (True, False):
            p = D.global_time_orientation(m, 0.5 * (R1 + rh), expanding)
            assert p["global_future_directed"] and not p["killing_clock_reversed"]
        d = D.global_time_orientation(m, 0.5 * (R2 + 1 / m.chi), True)
        assert d["global_future_directed"] and d["killing_clock_reversed"]
        assert d["out"]["Rdot_below_A"] and d["in"]["Rdot_below_A"]
