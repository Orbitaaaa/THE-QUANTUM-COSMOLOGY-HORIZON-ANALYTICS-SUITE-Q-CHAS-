"""Round 6l: model U (D-dimensional Einstein gravity + quasi-topological tower), its interiors, collider production with real parton
distributions, evaporation, multiple daughters from one hole, and experimental discriminants."""
import os
import warnings

import numpy as np

from babyuniverse import unified6l as U
from babyuniverse import collider6l as C

warnings.simplefilter("ignore")
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "DATA")


def test_regular_holes_extremality_and_temperature():
    assert abs(U.extremal_mu(4) - 1.5 * np.sqrt(3)) < 1e-12                     # 4D Hayward: r_s,ext = 3 sqrt3 / 2 l
    for D in (5, 8, 10):
        me = U.extremal_mu(D)
        assert U.static_horizons(0.999 * me, D) is None and U.static_horizons(1.001 * me, D) is not None
        md = U.extremal_mu(D, 1.0, "dymnikova")
        assert U.static_horizons(0.999 * md, D, 1.0, "dymnikova") is None and U.static_horizons(1.001 * md, D, 1.0, "dymnikova") is not None
        t = U.temperature(5 * me, D)
        assert abs(t["T"] / U.hayward_T_analytic(t["r_plus"], D) - 1) < 1e-6
    assert U.temperature(1.0001 * U.extremal_mu(10), 10)["T"] < 0.02 * U.temperature(10 * U.extremal_mu(10), 10)["T"]


def test_closed_dust_bounces_accelerates_and_recollapses():
    for kind in ("hayward", "dymnikova"):
        r = U.closed_dust(1e6, 10, 1.0, kind)
        assert abs(r["a_b"] - 1.0) < 1e-3 and r["a_max"] is not None and 0 < r["N_acc"] < 3
    r = U.closed_dust(1e8, 10)
    assert abs(r["a_infl"] / r["a_infl_exact"] - 1) < 2e-3 and abs(r["N_acc"] - r["N_acc_asymptotic"]) < 0.01
    st = U.closed_dust(2e38, 4)
    assert 29 < st["N_acc"] < 30.5 and abs(st["a_max"] / 2e38 - 1) < 1e-6           # ~30 e-folds inside a stellar hole, then recollapse


def test_vacuum_threshold_matches_einstein_static_limit():
    v = U.vacuum_threshold(1e4, 4)
    assert abs(v["Lam_c"] / (4 / (27 * 1e8)) - 1) < 1e-4
    v10 = U.vacuum_threshold(1e6, 10)
    assert 1.0 < v10["ratio"] < 1.05
    assert "for ever" in U.closed_dust(1e4, 4, 1.0, "hayward", 1.1 * v["Lam_c"])["fate"]


def test_vacuum_collapse_core_is_static_and_lemmas():
    v = U.vacuum_collapse_interior(50 * U.extremal_mu(10), 10)
    assert v["horizon"] and not v["inflating"] and v["core_radius_over_l"] < 1.01
    assert not U.vacuum_collapse_interior(0.8 * U.extremal_mu(10), 10)["horizon"]
    assert U.bounce_untrapped_lemma(500)["frame_identity_max_rel_error"] < 1e-8
    e = U.ecsk_os_mass_deficit(1.0, 1e-6, np.pi / 4)
    assert e["deficit_fraction"] > 0.98 and abs(e["compactness_at_bounce"] - np.sin(np.pi / 4) ** 2) < 1e-6


def test_second_infall_joins_the_existing_daughter_and_saturates_the_bound():
    for dl in (2.0, 50.0):
        r = U.two_epoch_infall(10, 1.0, 20.0, 0.7, dl)
        assert r["meets_inside_core"]
    late = U.two_epoch_infall(10, 1.0, 20.0, 0.7, 50.0)["meeting"]
    assert late["log10_rminus_minus_r"] < -15                                      # meets it at the inner horizon, exponentially close
    c = U.dtr_critical_delay(10, 1.0, 20.0, 1e-3, 0.7)
    assert c["found"] and 5 < c["delay_c"] < 40
    f = U.fmm_patch_separation(1e6)
    assert f["induced_metric_match"] < 1e-9 and abs(f["delta_v_disjoint"] / f["approx_2l_l_over_rs_13"] - 1) < 0.01


def test_parton_luminosity_and_cross_sections():
    assert abs(C.trapped_fraction_headon(0) - 1 / np.sqrt(2)) < 1e-12
    assert 0.55 < C.trapped_fraction_headon(7) < C.trapped_fraction_headon(4) < 0.65
    import json
    tab = json.load(open(os.path.join(DATA, "round6l_lumi.json")))
    for sr in tab["sum_rules"]:
        assert abs(sr["momentum"] - 1) < 0.01 and abs(sr["u_valence"] - 2) < 0.02 and abs(sr["d_valence"] - 1) < 0.02
    L = C.load_lumi(os.path.join(DATA, "round6l_lumi.json"))
    s1 = C.sigma_pp(L, 100.0, 6, 8.0, 38.0, "conservative")
    s2 = C.sigma_pp(L, 100.0, 6, 8.0, 45.0, "conservative")
    assert 2e-4 < s1 < 1e-3 and s2 < s1
    assert C.sigma_pp(L, 100.0, 6, 8.0, 38.0, "optimistic") > 100 * s1
    assert C.sigma_pp(L, 14.0, 6, 8.0, 38.0, "optimistic") == 0.0
    cl = C.c_l_max_reachable(L, 6, 8.0, 100.0, 30.0)
    assert 0.65 < cl < 0.8
    lo = C.c_l_min_lhc(L, 6, 8.0)
    assert 0.5 < lo < cl and 8.0 < C.thresholds(6, 8.0, lo)["model_U_TeV"] < 11.0      # window: M_ext ~ 9-50 TeV


def test_thresholds_evaporation_and_formation():
    t = U.model_u_thresholds(6, 1.0)
    assert abs(t["M_ext_over_MD"] - 51.24) < 0.05 and abs(t["M_fmm_over_MD"] - 4.725) < 0.01
    assert abs(U.model_u_thresholds(6, 0.5)["M_ext_over_MD"] / t["M_ext_over_MD"] - 0.5**7) < 1e-9
    eU = C.evaporation(6, 8.0, 0.7, 15.0, "hayward")
    eG = C.evaporation(6, 8.0, 0.7, 15.0, "gr")
    assert eU["remnant_energy_fraction"] > 0.2 and eU["quanta"] > eG["quanta"] and eU["T_peak_GeV"] < eG["T_peak_GeV"]
    fe = C.formation_vs_evaporation(6, 8.0, 0.7, 2.0, "hayward")
    assert fe["ratio"] > 100
    tc = C.temperature_curves(6, "hayward", 40)
    assert tc["rows"][0]["ratio"] < 0.05 and 2.5 < tc["r_plus_over_l_where_within_10pct"] < 5
    assert len(C.discriminants()) == 7
