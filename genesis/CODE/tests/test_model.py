import numpy as np
import pytest

from babyuniverse.model import ThinWallModel, model_from_dimensionless

rng = np.random.default_rng(1)
MODELS = [ThinWallModel(1.0, 0.5, 0.1), ThinWallModel(0.7, 1.3, 0.3, 0.2), ThinWallModel(2.0, 0.3, 0.02, 0.5),
          model_from_dimensionless(0.6, 1.2), model_from_dimensionless(-0.5, 4.5)]


@pytest.mark.parametrize("m", MODELS)
def test_potential_forms_agree(m):
    R = np.geomspace(1e-2, 5.0, 400)
    assert np.allclose(m.V(R), m.V_squared_form(R), rtol=1e-10, atol=1e-10)
    assert np.allclose(m.V(R), m.f_in(R) - m.A_in(R) ** 2, rtol=1e-10, atol=1e-10)
    h = 1e-6
    assert np.allclose((m.V(R + h) - m.V(R - h)) / (2 * h), m.dV(R), rtol=1e-5, atol=1e-6)


@pytest.mark.parametrize("m", MODELS)
def test_unsquared_junction_holds_on_shell(m):
    """Every solution of Rdot^2 + V = 0 satisfies the *unsquared* Israel
    equation with eps = sign(A): squaring introduces no spurious branch."""
    R = np.geomspace(1e-2, 5.0, 400)
    V = m.V(R)
    ok = V <= 0
    Rd = np.sqrt(-V[ok])
    res = m.junction_residual(R[ok], Rd, np.sign(m.A_in(R[ok])), np.sign(m.A_out(R[ok])))
    assert np.max(np.abs(res)) < 1e-9
    # the opposite sign assignment fails wherever A != 0
    bad = m.junction_residual(R[ok], Rd, -np.sign(m.A_in(R[ok])), np.sign(m.A_out(R[ok])))
    assert np.min(np.abs(bad)) > 1e-6 or np.any(np.abs(m.A_in(R[ok])) < 1e-6)


@pytest.mark.parametrize("m", MODELS)
def test_mass_function_inverts_junction(m):
    R = np.geomspace(1e-2, 5.0, 200)
    V = m.V(R)
    ok = V <= 0
    Mrec = m.mass_function(R[ok], np.sqrt(-V[ok]), np.sign(m.A_in(R[ok])))
    assert np.allclose(Mrec, m.M, rtol=1e-9, atol=1e-12)


def test_turning_points_and_critical_mass():
    for gam in (-0.8, -0.2, 0.3, 0.9):
        m = model_from_dimensionless(gam, 3.0 * ThinWallModel.g_min(gam))
        for R in m.turning_points():
            assert abs(m.V(R)) < 1e-10
        mc = m.with_mass(m.M_critical)
        zs = ThinWallModel.z_star(gam) * mc.R_unit
        assert abs(mc.V(zs)) < 1e-9 and abs(mc.dV(zs)) < 1e-8


def test_dimensionless_reduction():
    for m in MODELS:
        z = np.geomspace(0.1, 5, 50)
        R = z * m.R_unit
        assert np.allclose(m.V(R), 1 - ThinWallModel.g(z, m.gamma) / m.E, rtol=1e-10, atol=1e-10)
        assert abs(m.gamma) <= 1


def test_sign_flip_radii():
    m = model_from_dimensionless(0.6, 1.2)
    assert abs(m.A_out(m.flip_radius_out())) < 1e-12
    assert abs(m.A_in(m.flip_radius_in())) < 1e-12
    # the eps_out flip lies inside the forbidden region (regime A) or inside the horizon
    R1, R2 = m.turning_points()
    assert R1 < m.flip_radius_out() < R2
