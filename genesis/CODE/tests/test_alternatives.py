"""Roadmap Stage 4: alternative mechanisms and their costs."""
import warnings

import numpy as np

from babyuniverse.model import ThinWallModel
from babyuniverse import alternatives as Al, tunnelling as T

warnings.simplefilter("ignore")


def test_seeded_model_reduces_to_baseline_and_contains_equator():
    M = 0.3 * ThinWallModel(1.0, 1.0, 1e-9).M_critical
    B0 = T.B_canonical(ThinWallModel(1.0, 1.0, M))
    assert abs(Al.SeededModel(1.0, 1.0, M, 0.0).B() / B0 - 1) < 1e-8
    for f in (0.3, 0.7):
        s = Al.SeededModel(1.0, 1.0, M, f * M)
        R1 = s.turning_points()[0]
        assert s.A_out(R1) > 0                                   # still laboratory-type at R1
        assert 0.7 < s.B() / B0 < 1.0                            # modest reduction only
        assert s.initial_data_contains_equator()


def test_qei_species_bound_scales_with_de_sitter_entropy():
    assert abs(Al.qei_max_antitrapped_radius(1.0) - 0.3455) < 1e-3
    for E in (1e16, 1e17, 1e18):
        r = Al.species_required(E)
        assert abs(r["N_over_SdS"] - 8 / 3) < 1e-9
        assert r["species_cutoff_over_H"] < 1


def test_ds_parent_and_fmm_numbers():
    d = Al.ds_parent_upward_exponent()
    assert 1e122 < d["S_dS_parent"] < 1e123
    f = Al.fmm_junction(1.0)
    assert f["inside_horizon"] and abs(f["M_min"] - 0.9306) < 1e-3
    assert len(Al.mechanism_ledger()) >= 9
