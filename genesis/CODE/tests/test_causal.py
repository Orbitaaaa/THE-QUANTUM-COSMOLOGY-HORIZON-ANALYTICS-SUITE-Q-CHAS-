import numpy as np

from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse.causal import farhi_guth_audit, nec_wall, time_symmetric_slice, misner_sharp_out


def test_nec_holds_for_positive_tension():
    assert nec_wall(ThinWallModel(1.0, 0.7, 0.1)) >= -1e-15


def test_every_inflating_branch_is_past_singular():
    for gam in (-0.8, 0.0, 0.8):
        gmin, g1 = ThinWallModel.g_min(gam), 2 - 2 * gam
        for E in (3 * g1, 0.5 * (gmin + g1), 0.6 * gmin):
            for row in farhi_guth_audit(model_from_dimensionless(gam, E)):
                if row["future_inflating"]:
                    assert row["interior_reaches_antitrapped"]
                    assert row["exterior_contains_past_singularity"]


def test_fgg_endpoint_slices():
    m = model_from_dimensionless(-0.5, 4.5)
    R1, R2 = m.turning_points()
    s1, s2 = time_symmetric_slice(m, R1), time_symmetric_slice(m, R2)
    assert s1["regular_laboratory_data"] and not s1["extremal"]
    kinds = {(e.kind, e.where) for e in s2["extremal"]}
    assert ("min", "Einstein-Rosen bridge") in kinds and ("max", "wall") in kinds


def test_misner_sharp():
    m = ThinWallModel(1.0, 0.5, 0.2, 0.1)
    r = np.linspace(0.5, 3, 5)
    assert np.allclose(misner_sharp_out(m, r), m.M + m.eta**2 * r**3 / (2 * m.G))
