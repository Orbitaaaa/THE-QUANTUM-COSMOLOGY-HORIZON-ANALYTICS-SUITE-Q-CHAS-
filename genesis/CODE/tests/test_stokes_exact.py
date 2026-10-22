"""Roadmap P-Q2: exact Stokes multiplier of the squared constraint (M = 0, kappa > chi)."""
import numpy as np

from babyuniverse.stokes_exact import SquaredConstraintODE


def test_transport_exponent_is_minus_half_for_every_ordering():
    for lam in (0.0, 0.5, 1.0):
        ode = SquaredConstraintODE(1.0, 1.5, lam)
        assert abs(ode.alpha_from_transport(1) + 0.5) < 1e-5
        assert abs(ode.alpha_from_transport(2) + 0.5) < 1e-5


def test_multiplier_converges_to_sqrt2():
    ode = SquaredConstraintODE(1.0, 1.5, 0.5)
    e1 = abs(ode.multiplier(0.005) / np.sqrt(2) - 1)
    e2 = abs(ode.multiplier(0.0025) / np.sqrt(2) - 1)
    assert e2 < 2.5e-3 and e2 < 0.6 * e1


def test_window_independence():
    ode = SquaredConstraintODE(1.0, 1.5, 1.0)
    s = [ode.multiplier(0.0025, hw) for hw in (0.7, 1.0, 1.3)]
    assert max(s) - min(s) < 1e-4


def test_controls_nu_zero_and_minus_one():
    rq = np.sqrt(1.5**4 - 1)
    ode0 = SquaredConstraintODE(1.0, 1.5, 0.5, beta=rq)          # nu = 0: no switching
    assert abs(ode0.multiplier(0.005)) < 5e-3
    ode1 = SquaredConstraintODE(1.0, 1.5, 0.5, beta=-rq)         # nu = -1: s = sqrt(2 pi)(2 r_q/(a^2 h))^(1/2)
    assert abs(ode1.multiplier(0.005) / ode1.predicted(0.005) - 1) < 0.01
