"""Problem C: exterior event signature; problem B: assistance frequency."""
import warnings

import numpy as np

from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse import signatures as Sg, bounds as Bd

warnings.simplefilter("ignore")


def test_collapse_tail_and_daughter_cutoff():
    for m in (model_from_dimensionless(0.6, 2.4), model_from_dimensionless(0.3, 5.0)):
        r = Sg.late_photon_fraction(m)
        assert r["N_late_collapse"] > 0 and r["N_late_daughter"] == 0
        assert abs(r["tail_efold"] / r["four_GM"] - 1) < 0.03          # rate ~ exp(-u/4GM)
        R1 = m.turning_points()[0]
        assert abs(r["last_redshift_factor"] - np.sqrt(m.f_out(R1))) < 1e-12


def test_light_curve_is_continuous_at_turning_point():
    m = model_from_dimensionless(0.6, 2.4)
    lc = Sg.light_curve(m, n=120)
    assert abs(lc["u_expanding"][-1] - lc["u_collapse"][0]) < 1e-9
    assert abs(lc["rate_expanding"][-1] - lc["rate_collapse"][0]) < 1e-9
    assert np.all(np.diff(lc["u_collapse"]) > 0) and np.all(np.diff(lc["u_expanding"]) > 0)


def test_assistance_needs_vacuum_scale_frequencies():
    o = Bd.frontier_optimum(10.0, 1.0, 100.0)
    a = Sg.assistance_frequency(o["E_vac_GeV"], o["E_wall_GeV"], 1.0)
    assert 0.1 < a["omega_c_over_Evac"] < 1.0
    assert abs(a["delta_tau_E_lP"] / a["thin_estimate_lP"] - 1) < 0.1
