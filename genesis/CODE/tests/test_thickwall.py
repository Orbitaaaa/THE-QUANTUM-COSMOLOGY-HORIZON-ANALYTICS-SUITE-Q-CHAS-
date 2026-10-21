import numpy as np

from babyuniverse.model import ThinWallModel
from babyuniverse.thickwall import QuarticPotential, thin_wall_comparison, cdl_bounce
from babyuniverse.tunnelling import instanton_M0

RHO_F = 3 / (8 * np.pi)


def test_time_symmetric_data_approach_thin_wall():
    for kap, eps_out in ((0.5, 1), (1.5, -1)):
        pot = QuarticPotential.from_width_tension(0.0025, kap / (4 * np.pi), RHO_F)
        c = thin_wall_comparison(pot, 0.7)
        assert c["eps_out"] == eps_out
        assert abs(c["jump"] / c["jump_thin"] - 1) < 5e-3
        assert abs(c["M_far"] / c["M_thin"] - 1) < 2e-3
        if eps_out < 0:   # bag of gold: the Einstein-Rosen bridge appears at r = 2GM
            assert abs(c["r_min_outside"] / c["bridge_thin"] - 1) < 1e-3


def test_cdl_bounce_thin_wall_limit():
    kap = 0.5
    pot = QuarticPotential.from_width_tension(0.01, kap / (4 * np.pi), RHO_F)
    r = cdl_bounce(pot)
    d = instanton_M0(ThinWallModel(1.0, kap, 0.0))
    assert abs(r["B"] / d["B_CdL_down"] - 1) < 1e-2
    assert abs(r["a_wall"] / d["R0"] - 1) < 1e-2
    assert r["constraint_max"] < 1e-8
