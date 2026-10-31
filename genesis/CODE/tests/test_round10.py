"""Round 10: double-null code, open-problem calculations, statistics and the event pipeline."""
import os
import warnings

import numpy as np

from babyuniverse import doublenull10 as N
from babyuniverse import openprob10 as O
from babyuniverse import stats10 as S
from babyuniverse import pipeline10 as P10
from babyuniverse import lhcdata10 as LD

warnings.simplefilter("ignore")


def test_double_null_static_and_gr_limit():
    # (r, sigma) formulation: the late-v error near the Cauchy horizon converges (it motivates the hybrid formulation)
    e1 = N.static_check(5, 10.0, "hayward", 300, 600, 12.0)["max_rel_mass_error"]
    e2 = N.static_check(5, 10.0, "hayward", 600, 1200, 12.0)["max_rel_mass_error"]
    assert e2 < 0.5 * e1
    res = N.run_hyb(5, 10.0, "hayward", nu=600, nv=1600, v_max=16.0, u_eh=3.5745, n_cluster=200, n_out_u=60, n_out_v=40)
    i = np.argmin(abs(res["u"] - 6.0))
    assert abs(res["r"][i, -1] - res["r_minus"]) < 1e-3          # the static Cauchy horizon is held
    assert abs(res["mu"][i, -1] - 10.0) < 1e-9                     # no flux, mass exact in the hybrid formulation


def test_inner_horizon_curvature_grows_after_a_pulse():
    res = N.run_hyb(5, 10.0, "hayward", nu=1000, nv=4000, v_max=20.0, amp=0.05, u_eh=3.4568, n_cluster=400, n_out_u=120, n_out_v=80)
    i = np.argmin(abs(res["u"] - 6.0))
    K = res["K"][i]; v = res["v"]
    assert np.nanmax(K[v > 14]) > 1e6 * np.nanmax(K[(v > 5) & (v < 8)])
    assert np.nanmax(res["psi"][i]) < 1.0                            # the tower's bound on psi holds


def test_open_problem_calculations():
    assert abs(float(O.rotational_fraction(1.0, 6)) - (1 - 2 ** (-1 / 8))) < 1e-12
    d = O.remnant_discharge(6, 8.0, 0.66)
    assert d["colour_before_hadronisation"] and d["tau_electric_s"] < 1e-20
    assert abs(O.causality_threshold((0.5,))["rows"][0]["M_HS_bound_over_MD"] - 2.0) < 1e-12
    ctrl = O.scalar_from_bounce(10, None, 0.5, 0.0, 1.0, 1.0, t_max=300.0)
    assert ctrl["fate"] == "recollapses"
    hyst = O.scalar_from_bounce(10, None, 0.05, 1e-3, 1.0, 1.0, t_max=3000.0)
    assert hyst["fate"].startswith("lasting") and max(hyst["cycle_maxima"]) > 1.5 * hyst["cycle_maxima"][0]


def test_statistics():
    s = np.array([5.0, 10.0, 3.0]); b = np.array([20.0, 10.0, 2.0])
    m = S.Model((s, s, s), {"z": (b, b, b)}, sig_unc=0.0, lumi_unc=0.0, bkg_unc={"z": 1e-6})
    assert abs(m.Z_asimov() - np.sqrt(np.sum(2 * ((s + b) * np.log(1 + s / b) - s)))) < 1e-6
    r = S.llr_events_for_5sigma([0.5, 0.3, 0.2], [0.2, 0.3, 0.5])
    assert 35 < r["events_5sigma"] < 60


def test_pipeline_pieces():
    det = P10.Detector("fcchh", seed=1)
    assert abs(det.jet_sigma(np.array([1000.0]), np.array([0.0]))[0] - np.hypot(30.0, 0.5 * np.sqrt(1000.0))) < 1e-9
    ev = np.zeros(3, dtype=[(f, "f8") for f in P10.FIELDS])
    ev["met"] = [150, 900, 900]; ev["lead_pt"] = [400, 1000, 1000]; ev["n_jets"] = [1, 1, 1]; ev["dphi_min"] = [2, 2, 0.1]
    assert list(P10.atlas_monojet(ev)) == [False, True, False]
    a = LD.em_arrays()
    assert len(a["obs"]) == 13 and a["obs"][-1] == 207
    out, sig = P10.background("znunu", 13.0, 400, 0, 60, P10.Detector("atlas", seed=2), seed=3)
    assert 0.3 < sig < 1.5 and len(out) > 50
