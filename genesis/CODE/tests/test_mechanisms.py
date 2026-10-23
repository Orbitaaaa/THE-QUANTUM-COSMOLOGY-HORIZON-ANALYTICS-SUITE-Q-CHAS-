"""Round 6: the mechanism atlas (budget theorem numbers, thin walls of any matter, charged exteriors, seeds, genesis)."""
import warnings

import numpy as np

from babyuniverse import mechanisms as Mc, tunnelling as T
from babyuniverse.model import ThinWallModel

warnings.simplefilter("ignore")


def test_budget_numbers():
    b = Mc.genesis_budget(0.0)
    assert abs(b["null_force_N"] - 4.815e42) / 4.815e42 < 1e-3
    assert abs(b["fraction_of_max_force"] - 1 / (2 * np.pi)) < 1e-12
    assert abs(Mc.genesis_budget(0.75)["dimensionless"] - 0.25) < 1e-12
    s1, s2 = Mc.stress_scale(1.0), Mc.stress_scale(1e-2)
    assert abs(s2["stress_Pa"] / s1["stress_Pa"] - 1e4) < 1e-6 * 1e4
    assert 20 < s1["scale_GeV"] < 25                                   # ~22 GeV for a 1 m region
    assert Mc.casimir_shortfall(1e-8)["shortfall_orders"] > 50


def test_thin_walls_cannot_cross_the_bridge_classically():
    # any wall matter: at a flip of the exterior orientation V = f_out, so a crossing needs a horizon interior
    for kap in (lambda R: 0.5 * (1 - (0.15 / R) ** 2), lambda R: 0.3 + 0.2 * R, lambda R: 1.0):
        sh = Mc.GeneralShell(lambda R: 1 - R**2, lambda R: 1 - 2 * 0.05 / R, kap)
        for row in sh.bridge_crossing_report(0.16, 3.0):
            assert abs(row["V"] - row["f_out"]) < 1e-9
            assert row["classical_crossing_possible"] == row["inside_horizon"]
    two = Mc.GeneralShell(lambda R: 1 - R**2, lambda R: 1 - 0.1 / R, lambda R: 0.5 * (1 - (0.15 / R) ** 2))
    assert not two.nec_on_wall(0.3)                                    # the negative-dust component violates the NEC


def test_visible_bubbles_need_negative_tension_and_destabilise_minkowski():
    assert Mc.visible_bubble_orientation(1.0, 0.3) == -1 and Mc.visible_bubble_orientation(1.0, -0.3) == 1
    for k in (0.3, 1.0, 3.0):
        up = Mc.upward_nucleation_exponent(1.0, -k)
        assert abs(up["B"] - T.B_canonical(ThinWallModel(1.0, k, 0.0))) < 1e-8 * max(1.0, up["B"])
    assert 550 < Mc.vacuum_survival_exponent() < 570


def test_charged_exterior_route():
    B0 = T.B_canonical(ThinWallModel(1.0, 1.0, 0.1))
    for qf in (0.3, 0.9, 0.99):
        assert abs(Mc.ChargedBubble(1.0, 1.0, 0.1, qf * 0.1).B_lab() / B0 - 1) < 0.03
    r = Mc.ChargedBubble(1.0, 1.0, 0.5, 0.25).classical_route()
    assert r["barrier_free"] and r["flip_inside_horizons"] and r["R_bounce"] < r["r_minus"]
    assert not Mc.ChargedBubble(1.0, 1.0, 0.3, 0.03).classical_route()["barrier_free"]
    mcr = Mc.critical_mass(1.0, 1.0, 0.01, hi=1.5)
    assert abs(mcr - ThinWallModel(1.0, 1.0, 0.1).M_critical) < 2e-3   # continuous with the uncharged threshold
    assert Mc.critical_mass(1.0, 0.3, 0.9, hi=1.5) > Mc.critical_mass(1.0, 0.3, 0.5, hi=1.5)
    t = Mc.classical_threshold_physical(1e16)
    assert 1e-3 < t["M_cr_kg"] < 1e-1


def test_seeds_and_genesis_scaling():
    w = {r["R_lP"]: r for r in Mc.frontier_seed_window()}
    assert w[100]["lifetime_s"] < 1e-17 and 1e-3 < w[300]["lifetime_s"] < 1 and w[1000]["lifetime_s"] > 1e30
    g = Mc.genesis_timeline(1e-2, 1e-8, 1e11)
    assert g["budget_ratio"] > 1 and g["validity_H_t"] < 0.1           # G1 met automatically, flat-space genesis still valid
    assert abs(g["budget_ratio"] - 6 * (1e-4 / 1e-8) ** (1 / 3)) < 1e-6 * g["budget_ratio"]


def test_sharp_edged_genesis_seed_gives_a_visible_bubble():
    for h, rc, M in ((1e-4, 1.0, 1e-3), (1e-2, 0.3, 1e-2), (1e-6, 1.0, 0.1)):
        g = Mc.genesis_interface(h, rc, M)
        assert not g["branches_meet"] and g["antitrapped_reached"]
        assert g["tension_negative_before_antitrapping"] and g["final_kappaR"] < 0
        assert abs(g["HR_when_tension_turns_negative"] / g["HR_predicted"] - 1) < 0.05
