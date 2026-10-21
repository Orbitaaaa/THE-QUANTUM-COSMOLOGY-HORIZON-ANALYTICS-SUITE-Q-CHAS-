import warnings

import numpy as np
import pytest

from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse import contours as C, tunnelling as T

warnings.simplefilter("ignore")


def test_real_contour_reproduces_B():
    m = model_from_dimensionless(-0.5, 4.8)
    I, _ = C.continue_P(m, C.real_barrier_path(m))
    assert abs(2 * abs(np.imag(I)) - T.B_canonical(m)) < 1e-5
    assert abs(np.real(I)) < 1e-10


@pytest.mark.parametrize("turns", [1, -1, 2])
def test_horizon_monodromy(turns):
    m = model_from_dimensionless(0.6, 1.2)
    for rh in (2 * m.G * m.M, 1 / m.chi):
        r = C.monodromy_shift(m, rh, turns)
        assert abs(abs(r["delta_B"] / r["prediction_per_turn"]) - abs(turns)) < 1e-3
        assert abs(r["delta_Re"]) < 1e-12


def test_cosmological_horizon_monodromy():
    m = ThinWallModel(1.0, 1.0, 0.1, 0.4)
    rc = m.horizons_out()[1]
    r = C.monodromy_shift(m, rc, 1)
    assert abs(abs(r["delta_B"] / r["prediction_per_turn"]) - 1) < 1e-3


def test_sheet_theorem():
    """analytic branch <-> sign-free constraint; principal <-> sign-dressed constraint"""
    for gam in (-0.5, 0.6):
        m = model_from_dimensionless(gam, 1.5 * (2 - 2 * gam))
        s = C.principal_branch_constraint_violation(m)
        assert s["analytic_max_violation"] < 1e-12
        assert s["principal_dressed_violation"] < 1e-12
        assert s["principal_max_violation"] > 0.5
        assert s["analytic_dressed_violation"] > 0.5


@pytest.mark.parametrize("gam,E", [(-0.5, 4.8), (0.6, 1.2)])
def test_branch_locus_at_flip(gam, E):
    """Theorem III-4: the trajectory touches the branch locus of the true
    Hamiltonian at the exterior flip radius (double zero), both sheets give M there."""
    m = model_from_dimensionless(gam, E)
    b = C.branch_locus_check(m)
    assert abs(b["relative_radicand"]) < 1e-12
    assert b["radicand_sign_changes_on_barrier"] == 0
    assert b["H_minus_minus_M"] < 1e-7 and b["H_plus_minus_M"] < 1e-7


def test_radicand_identity_on_trajectory():
    """kappa^2 R^2 - f_in sin^2 theta = kappa^2 R^2 A_out^2 / f_out along the barrier."""
    m = model_from_dimensionless(-0.5, 4.8)
    R1, R2 = m.turning_points()[0], m.turning_points()[-1]
    R = np.linspace(R1, R2, 50)[1:-1]
    th = T.euclidean_momentum(m, R) * m.G / R
    lhs = (m.kappa * R) ** 2 - m.f_in(R) * np.sin(th) ** 2
    rhs = (m.kappa * R) ** 2 * m.A_out(R) ** 2 / m.f_out(R)
    assert np.max(np.abs(lhs - rhs)) < 1e-12


def test_sector_sheets_depend_on_convention():
    """single-flip daughter: H_- in the real-angle convention, H_+ in the complex-angle one"""
    m = model_from_dimensionless(-0.5, 4.8)
    sh = C.sector_hamiltonians(m)
    assert sh["R1"]["artanh"][-1] < 1e-12 and sh["R1"]["sorkin"][-1] < 1e-12
    assert sh["R2"]["eps"] == (1, -1)
    assert sh["R2"]["artanh"][-1] < 1e-12 and sh["R2"]["artanh"][1] > 1e-3
    assert sh["R2"]["sorkin"][1] < 1e-12 and sh["R2"]["sorkin"][-1] > 1e-3


@pytest.mark.parametrize("gam,E", [(-0.5, 4.8), (0.6, 1.2)])
def test_zero_net_period(gam, E):
    m = model_from_dimensionless(gam, E)
    z = C.zero_net_shift(m, m.horizons_out()[0], 1 / m.chi)
    assert z["endpoint_dP"] < 1e-10
    assert abs(abs(z["delta_B"]) / z["period"] - 1) < 1e-5


def test_sector_switch_exponent():
    for gam in (-0.5, 0.0, 0.2):
        m = model_from_dimensionless(gam, 1.5 * (2 - 2 * gam))
        assert C.flip_class(m) == "single"
        Bq, Ban = C.B_sector_switch(m), T.B_canonical(m)
        assert 0 < Bq < Ban
