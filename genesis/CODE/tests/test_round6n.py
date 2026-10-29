"""Round 6n: model U beyond idealised collapse (generalised Misner-Sharp hydrodynamics)."""
import warnings

import numpy as np

from babyuniverse import hydro6n as H

warnings.simplefilter("ignore")


def test_hydro_reproduces_exact_bounce():
    for D in (5, 10):
        r = H.os_check(D, 60, mu_over_ext=20.0)
        assert r["rel_err_R"] < 1e-4 and r["rel_err_t"] < 2e-3


def test_bounce_is_universal_for_any_fluid():
    for D in (5, 10):
        for w in (-0.3, 0.0, 1 / 3, 1.0):
            r = H.frw_bounce_any_fluid(w, 1e4, D)
            assert r["bounces"] and abs(r["a_b"] - 1) < 1e-3 and r["a_max"] > r["a_b"]


def test_inhomogeneous_dust_bounce_law():
    r = H.ltb_bounce_radii(5, "smooth", N=40)
    assert r["bounce_time_monotonic"] and r["spread_of_bounce_times"] > 1.5
    assert all(abs(q["R_b"] / q["R_b_law"] - 1) < 0.03 for q in r["rows"] if q["R_b"] < 0.9)


def test_pressured_collapse_regular_and_energy_conserving():
    r = H.run_to_throat(5, 40, 20.0, "smooth", 0.1, max_steps=3000)
    assert r["first_trapped"] is not None and r["first_bounce"] is not None and r["first_bounce"] > r["first_trapped"]
    assert r["psi_max_alpha"] < 1.0 and r["energy_error"] < 1e-12
    th = H.pressure_threshold(5, 40, 20.0, "smooth", (0.1, 8.0))
    assert th[0]["horizon"] and not th[1]["horizon"]
