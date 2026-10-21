import warnings

import numpy as np
import pytest

from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse import tunnelling as T

warnings.simplefilter("ignore")


@pytest.mark.parametrize("chi,kap,eta", [(1, 0.5, 0), (1, 1.5, 0), (1, 0.5, 0.3), (1, 1.5, 0.5),
                                         (0.3, 1.0, 0.9), (2.0, 3.0, 0.4)])
def test_M0_closed_form_and_morse(chi, kap, eta):
    m = ThinWallModel(chi, kap, 0.0, eta)
    d = T.instanton_M0(m)
    Bn = T.B_canonical(m)
    assert abs(Bn - d["B_closed"]) < 1e-11 * max(1, Bn)
    assert abs(Bn - (d["S_E"] + d["morse"])) < 1e-10 * max(1, Bn)
    assert Bn > 0
    if d["q"] < 0:                       # false vacuum is the large cap: CdL decay exponent
        assert abs(Bn - d["B_CdL_down"]) < 1e-10 * max(1, Bn)


def test_small_tension_precision():
    m = ThinWallModel(1.0, 3.56e-3, 0.0)
    d = T.instanton_M0(m)
    assert abs(d["B_CdL_down"] / T.B_canonical(m) - 1) < 1e-6


@pytest.mark.parametrize("gam,Ef", [(-0.5, 1.6), (0.6, 1.5), (0.3, 3.0), (-0.5, 20.0), ("B", -0.5)])
def test_smarr_identity(gam, Ef):
    if gam == "B":
        g = Ef
        m = model_from_dimensionless(g, 0.5 * (2 - 2 * g + ThinWallModel.g_min(g)))
    else:
        m = model_from_dimensionless(gam, Ef * (2 - 2 * gam))
    s = T.smarr_check(m)
    assert abs(s["residual"]) < 1e-7 * max(1.0, abs(s["B"]))


def test_first_law():
    f = T.first_law_check(model_from_dimensionless(0.6, 1.5 * 0.8))
    assert abs(f["dBt_dM"] / f["minus_dt_W"] - 1) < 1e-6
    assert abs(f["dBt_dLambda"] / f["Vol_over_8piG"] - 1) < 1e-6
    assert abs(f["dBt_dsigma"] / f["A_W"] - 1) < 1e-6


def test_principal_branch_is_unphysical():
    m = model_from_dimensionless(-0.5, 1.6 * 3.0)
    assert T.B_canonical(m, "principal") < 0 < T.B_canonical(m, "analytic")


def test_continuity_at_M0():
    m = model_from_dimensionless(-0.5, 1.6 * 3.0)
    B0 = T.B_canonical(ThinWallModel(m.chi, m.kappa, 0.0))
    Bs = T.B_canonical(ThinWallModel(m.chi, m.kappa, 1e-6 * m.M_critical))
    assert abs(Bs - B0) < 1e-3 * B0


def test_embedding_obstruction_regime_A():
    for gam in (-0.7, 0.0, 0.6):
        m = model_from_dimensionless(gam, 1.5 * (2 - 2 * gam))
        r = T.embedding_report(m, ngrid=81)
        assert r["exterior_flip"]
        assert not r["exterior_embedded"]
        assert r["exterior_multiplicity_range"][1] >= 2


def test_interior_flip_can_be_embedded():
    """Counter-example to a naive interior analogue of the obstruction."""
    m = model_from_dimensionless(0.6, 1.5 * 0.8)
    r = T.embedding_report(m, ngrid=81)
    assert r["interior_flip"] and r["interior_embedded"]
