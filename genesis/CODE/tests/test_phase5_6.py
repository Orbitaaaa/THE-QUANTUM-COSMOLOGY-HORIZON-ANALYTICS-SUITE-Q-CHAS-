"""Phases V and VI: channel structure, critical-mass limit, observability helpers."""
import warnings

import numpy as np

from babyuniverse.model import ThinWallModel
from babyuniverse import quantum_system as Q, observability as O, causal

warnings.simplefilter("ignore")


def test_parent_edge_is_asymptotic_daughter_edge_is_finite():
    m = ThinWallModel(1.0, 1.0, 0.1545)
    tp = Q.asymptotic_time(m, "parent")
    td = Q.asymptotic_time(m, "daughter")
    slope = np.mean(tp["increments"][1:]) / np.log(10)
    assert abs(slope / (2 * m.G * m.M) - 1) < 1e-3          # t ~ -2GM ln(R - 2GM)
    assert abs(td["t"][-1] - td["t"][-2]) < 1e-4 * td["t"][-1]
    P = Q.momentum_near_edge(m, "daughter")
    assert P[-1]["P"] < P[0]["P"] - 4                        # logarithmic divergence


def test_topology_preserved_throat_created():
    m = ThinWallModel(1.0, 1.0, 0.1545)
    R1, R2 = m.turning_points()[0], m.turning_points()[-1]
    sp, sd = causal.time_symmetric_slice(m, R1), causal.time_symmetric_slice(m, R2)
    assert not sp["closed"] and not sd["closed"] and sp["extremal"] == []
    assert any(e.where == "Einstein-Rosen bridge" and e.kind == "min" for e in sd["extremal"])


def test_critical_mass_limit():
    for k in (0.03, 0.01):
        m = ThinWallModel(1.0, k, 1e-9)
        assert abs(m.M_critical * 2 * np.sqrt(1 + k**2) - 1) < 2e-4
    ms = O.masses(1e16, 1e16)
    assert abs(ms["M_cr_g"] / 5.59 - 1) < 0.01


def test_laboratory_bound_reference():
    assert abs(O.B_zero(1e16, 1e16) / 133.24 - 1) < 1e-3
    assert abs(O.B_zero_flat(1e16, 1e16) / 133.24 - 1) < 1e-3


def test_information_and_evaporation():
    info = Q.information_budget(1000.0, 1e-3)
    assert info["log10_p_daughter"] < -434 and info["S_BH_nats"] > 1e10
    ev = Q.evaporation_competition(1000.0, 1e-3)
    assert ev["semiclassical_black_hole"] and ev["max_attempts"] > 1e10


def test_theta_stable_matches_angles_and_survives_thin_barriers():
    from babyuniverse import tunnelling as T
    from babyuniverse.model import model_from_dimensionless
    m = model_from_dimensionless(-0.5, 4.8)
    R1, R2 = m.turning_points()[0], m.turning_points()[-1]
    R = np.linspace(R1, R2, 50)[1:-1]
    assert np.max(np.abs(T.theta_stable(m, R) - T.euclidean_momentum(m, R) * m.G / R)) < 1e-12
    p = O.planck_params(1e14, 0.7e14)                 # theta ~ 1e-8: angle differences lose all digits
    mt = ThinWallModel(1.0, p["k"], 10 * p["chi"])
    from babyuniverse import contours as C
    a, b = T.B_canonical(mt), C.B_sector_blind(mt)
    assert abs(a - b) < 1e-10 * a                      # double flip: B_sq = B_an


def test_best_laboratory_case_small_grid():
    res = O.best_laboratory_case([3e17, 9e17], [0.8, 0.9, 1.0], M_planck=10.0)
    b = res["best"]
    assert b is not None and b["metastable"] and b["wkb_ok"]
    assert 50 < b["B_best"] < 1000
