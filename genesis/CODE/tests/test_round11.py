"""Round 11: singularity-strength integrals, inner-extremal tower, greybody solver, exterior positivity, collider budget tools, remnant
energy loss and the timing-demonstrator analysis chain."""
import json
import os
import warnings

import numpy as np

from babyuniverse import collider11 as C11
from babyuniverse import greybody11 as G
from babyuniverse import greybody as G0
from babyuniverse import stable11 as ST
from babyuniverse import timing11 as T

warnings.simplefilter("ignore")
DATA = os.path.join(os.path.dirname(__file__), "..", "..", "DATA")


def test_inner_extremal_tower_construction():
    c = ST.ie_checks()
    assert abs(c["mu"] - 16.0) < 1e-12
    r = sorted(c["roots_x"])
    assert max(abs(x - 1.5) for x in r[:3]) < 1e-4 and abs(r[3] - 9.0) < 1e-9      # triple inner root, simple outer root
    assert abs(c["df_dx_at_inner"]) < 1e-12 and c["psi_monotone"] and c["Q_positive"] and c["psi_max"] < 1.0


def test_greybody_solver_limits():
    m = G.Metric(4, mu=2.0, kind="gr")                                  # Schwarzschild, M = 1
    assert abs(m.b_c - np.sqrt(27)) < 1e-6
    for w in (0.02, 0.05):                                              # independent Round 6c code
        assert abs(G.transmission(m, 0, 0, w) / G0.transmission(0, w, 1.0, 0) - 1) < 1e-3
    assert abs(G.transmission(m, 1, 1, 0.3) / G0.transmission(1, 0.3, 1.0, 1) - 1) < 1e-4
    mu = G.Metric(10, mu_over_ext=3.0)                                   # model U: sigma_abs -> horizon area at low frequency
    w = 0.01 / mu.r_plus
    assert abs(G.transmission(mu, 0, 0, w) / (4 * mu.r_plus**2 * w**2) - 1) < 0.01


def test_greybody_tables_physics():
    d = json.load(open(os.path.join(DATA, "round11_greybody.json")))
    h = d["hayward"]
    assert h["1.05"]["spins"]["1"]["number_ratio"] < 0.1                # gauge-boson emission suppressed near extremality
    assert h["1.05"]["spins"]["0"]["number_ratio"] > 1.2
    assert h["10"]["spins"]["1"]["number_ratio"] > 0.9
    for row in (h["3"], d["gr"]["1"]):
        for s in ("0", "0.5", "1"):
            assert 0.8 < row["spins"][s]["number_ratio"] < 1.6


def test_exterior_positivity():
    rows = G.exterior_positivity(D=10, qs=(1.05, 3.0), n=600)
    assert all(v == 0 for r in rows for k, v in r.items() if k.startswith("n_negative"))


def test_physics_results_converged():
    d = json.load(open(os.path.join(DATA, "round11_physics.json")))
    a, b = (d["strength_D5_pulse"][k] for k in ("2000x8000", "4000x16000"))
    assert abs(a["I_T"][-1] / b["I_T"][-1] - 1) < 0.01 and abs(a["r"][-1] - b["r"][-1]) < 1e-4
    assert a["lam"][-1] < 2.0                                            # finite affine distance
    c = [x for x in d["classification"] if x["D"] == 5]
    assert all(x["verdict"] == "singular growth" for x in c)
    ie = d["ie"]["4000x16000"]
    assert max(ie["mu"]) < 16.0 * 1.02 and np.nanmax(np.array(ie["K"], float)) < 120.0


def test_collider_budget_tools():
    s = C11.combine({"a": 0.3, "b": -0.4, "c": 0.4, "d": -0.3})
    assert abs(s["up"] - 0.5) < 1e-12 and abs(s["down"] - 0.5) < 1e-12
    w = C11.ew_weight(np.array([50.0, 1000.0, 2000.0]), 0.0059)
    assert w[0] == 1.0 and abs(w[1] - 0.85) < 0.01 and abs(w[2] - 0.756) < 0.01
    d = json.load(open(os.path.join(DATA, "round11_collider.json")))
    b = d["budget"]
    assert abs(b["nominal"]["discovery_TeV"] - 33.41) < 0.05
    assert "emission: geometric x^2" in b["combined"]["excluded_from_total"]
    assert 0.5 < b["combined"]["total_discovery"]["down"] < 2.0


def test_bethe_and_range():
    dE = C11.bethe_dEdx(np.array([3.0]))                                # near minimum ionisation in iron, ~1.45 MeV cm^2/g (no density effect)
    assert 1.3 < dE[0] < 1.7
    R1 = C11.range_g_cm2(1e3, np.array([0.3]))[0]; R2 = C11.range_g_cm2(2e3, np.array([0.3]))[0]
    assert abs(R2 / R1 - 2) < 1e-6                                       # range scales with mass at fixed velocity


def test_timing_analysis_chain():
    assert np.allclose(T.three_plane(30.0, 40.0, 50.0)["recovered"], (30.0, 40.0, 50.0))
    ev = T.simulate_bench_run(n_events=1200, seed=3)
    r = T.analyze_run(ev)
    assert r["n_selected"] > 600
    assert all(20.0 < p < 45.0 for p in r["plane_ps"])                   # calibrated model: ~30 ps per plane
    tb = T.tof_beta_resolution(35.0, 1.0, betas=(0.3,))[0]
    assert abs(tb["sigma_beta"] - 0.3**2 * 0.299792458 * np.sqrt(2) * 0.035) < 1e-12
