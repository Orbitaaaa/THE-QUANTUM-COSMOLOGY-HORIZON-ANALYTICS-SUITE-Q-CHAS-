import warnings

import numpy as np
import pytest

from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse.branches import classify, bgg_regime
from babyuniverse.geometry import branch_trajectories, regions_along

warnings.simplefilter("ignore")


def comp(lst):
    out = [lst[0]]
    for x in lst[1:]:
        if x != out[-1]:
            out.append(x)
    return out


EXPECTED = {"A": {"bounded": ["IV", "I", "II"], "bounce": ["III"]},
            "B": {"bounded": ["IV", "III", "II"], "bounce": ["III"]},
            "C": {"monotonic-expanding": ["IV", "III"], "monotonic-collapsing": ["III", "II"]}}


@pytest.mark.parametrize("gam", [0.6, -0.5])
@pytest.mark.parametrize("reg", ["A", "B", "C"])
def test_bgg_itineraries_and_chart_agreement(gam, reg):
    gmin, g1 = ThinWallModel.g_min(gam), 2 - 2 * gam
    E = {"A": 1.5 * g1, "B": 0.5 * (gmin + g1), "C": 0.7 * gmin}[reg]
    m = model_from_dimensionless(gam, E)
    assert bgg_regime(m).startswith(reg)
    for b in classify(m):
        assert comp(b.names_out()) == EXPECTED[reg][b.kind]
        tr = branch_trajectories(m, b)
        assert regions_along(tr["out"]) == comp(b.names_out())
        assert regions_along(tr["in"]) == comp(b.names_in())
        if b.flags["future_inflating"]:
            assert b.flags["final_eps_out"] == -1
            assert b.flags["exterior_contains_past_singularity"]
            assert b.flags["daughter_causally_disconnected"]


def test_m0_branch():
    m = ThinWallModel(1.0, 0.5, 0.0, 0.4)
    (b,) = classify(m)
    assert b.kind == "bounce-M0"
    assert comp(b.names_in()) == ["P", "S", "F"]
