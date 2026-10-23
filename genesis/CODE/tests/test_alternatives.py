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


def test_sds_nonflip_window_needs_comparable_parent_vacuum():
    assert Al.sds_nonflip_masses(0.3, 1e-3) == []
    assert Al.sds_nonflip_masses(0.3, 0.5) == []
    assert len(Al.sds_nonflip_masses(0.3, 0.95)) > 0
    t = Al.sds_nonflip_threshold(0.3, iters=14)
    assert 0.7 < t < 0.9


def test_sds_atlas_worldlines_are_timelike_and_continuous():
    from babyuniverse import sds_atlas as S, branches as Br
    for m, expect in ((ThinWallModel(1.0, 0.3, 0.05, 0.5), ["CP3", "III", "III", "CF3"]),
                      (ThinWallModel(1.0, 0.3, 0.1, 0.95), ["CP", "I", "I", "CF"])):
        for b in Br.classify(m):
            w = S.wall_worldline(m, b, n=150)
            c = S.causal_check(w)
            assert c["timelike_future"] and c["max_join_gap"] < 1e-4
            if b.kind == "bounce":
                assert [p["block"] for p in w] == expect
            else:
                assert [p["block"] for p in w] == ["WH", "I", "I", "BH"]


def test_gvz_scaling():
    assert abs(Al.gvz_fraction(1e-16, 1e20) / Al.gvz_fraction(1e-16, 4e20) - 2) < 1e-12
    r = Al.gvz_lambda_for(1.0, 1e20)
    assert 30 < r["B_inf"] < 45


def test_sds_euclidean_action_is_period_independent():
    from babyuniverse.sds_atlas import SdSExterior, euclidean_section
    for M, eta in ((0.05, 0.5), (0.1, 0.95), (0.01, 1.0)):
        s = SdSExterior(M, eta)
        assert s.kb > s.kc                                               # no period is smooth at both horizons
        ref = -(4 * np.pi * s.r_b**2 + 4 * np.pi * s.r_c**2) / 4
        for beta in (2 * np.pi / s.kb, 2 * np.pi / s.kc, 1.0, 0.0):
            r = euclidean_section(s, beta)
            assert abs(r["I_total"] - ref) < 1e-9 * abs(ref)
            assert r["ricci_max_dev"] < 1e-9 * 12 * eta**2 * 1e3
        assert abs(r["smarr_lhs"] - r["smarr_rhs"]) < 1e-10 * r["smarr_rhs"]
    tiny = euclidean_section(SdSExterior(1e-8, 0.3), 1.0)
    assert abs(tiny["I_total"] + np.pi / 0.3**2) < 1e-5                 # de Sitter limit -3 pi / (G Lambda)
