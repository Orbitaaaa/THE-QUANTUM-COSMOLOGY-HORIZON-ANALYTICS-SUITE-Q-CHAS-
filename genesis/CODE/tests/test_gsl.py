"""Round 6b: genesis and the GSL (dichotomy verdicts, QNEC cost, 2D flips, exact Vaidya flip, growth thresholds, Cauchy route)."""
import warnings

import numpy as np
import sympy as sp

from babyuniverse import gsl as S

warnings.simplefilter("ignore")


def test_product_solutions_and_thresholds_symbolic():
    r, M, Q, L, a = sp.symbols("r M Q Lambda a", positive=True)
    f = 1 - 2 * M / r + Q**2 / r**2 - L * r**2 / 3
    Msol = sp.solve(sp.Eq(f.subs(r, a), 0), M)[0]
    # double root <=> Lambda a^4 - a^2 + Q^2 = 0, and f''(a) = 2 (Q^2/a^4 - Lambda) on that locus
    assert sp.simplify(sp.diff(f, r).subs(r, a).subs(M, Msol) * a**3 + (L * a**4 - a**2 + Q**2)) == 0
    f2 = sp.diff(f, r, 2).subs(r, a).subs(M, Msol)
    assert sp.simplify((f2 - 2 * (Q**2 / a**4 - L)).subs(L, (1 - Q**2 / a**2) / a**2)) == 0
    # numerics: branches, the flat point at Lambda Q^2 = 1/4, no product beyond it
    p = S.product_solutions(1.0, 0.2)
    assert [b["kind"] for b in p["branches"]] == ["AdS2", "dS2"]
    for b in p["branches"]:
        assert abs(1 / b["a"] ** 2 - (0.2 + 1 / b["a"] ** 4)) < 1e-12
    pc = S.product_solutions(1.0, 0.25)
    assert all(b["kind"] == "flat" and abs(b["a"] ** 2 - 2) < 1e-12 for b in pc["branches"])
    assert not S.product_solutions(1.0, 0.2500001)["exists"]
    # growth thresholds: extremal bud 3/16 exactly; uncharged Nariai C_N/2 = 1/(3 sqrt Lambda)
    assert abs(S.bud_threshold() - 3 / 16) < 1e-12
    assert abs(S.switch_on_threshold("conserved") - 0.221212) < 1e-6
    for x in (0.05, 0.15, 0.24):                                       # static throat in a fixed vacuum: C = g(a) < C_N
        a = np.sqrt((1 - np.sqrt(1 - 4 * x)) / (2 * x))
        assert S._g(a, 1.0, x) < S.nariai_parameter(1.0, x, a)
    assert abs(S.nariai_parameter(0.0, 1.0, 1e-6) / 2 - 1 / 3) < 1e-12
    g = S.ks_growth(2.0, 1.0, 0.2, 1.0, 1e3)
    assert g["grows"] and g["accelerating_late"]
    assert not S.ks_growth(2.0, 1.0, 0.18, 1.0, 1e3)["grows"]


def test_qnec_cost_inequality():
    lam = np.linspace(0, 10, 4001)
    Sout = -0.3 * np.tanh(lam - 5) + 0.05 * np.sin(lam)
    eq = S.qnec_raychaudhuri(lam, Sout, 0.0, r0=5.0, rp0=-0.1)
    assert abs(eq["delta_rp"] - eq["bound"]) < 1e-8                  # QNEC saturated: equality
    sl = S.qnec_raychaudhuri(lam, Sout, lambda l: 1e-4 * np.exp(-(l - 3) ** 2), r0=5.0, rp0=-0.1)
    assert sl["delta_rp"] < sl["bound"]                                # slack only adds focusing
    # the cost of a flip in flat space at r = 1 m in Planck units is S_BH * 2/r
    rP = 1 / 1.616255e-35
    assert abs(S.qnec_cost(rP, 2 / rP) / (2 * np.pi * rP) - 1) < 1e-12
    assert abs(S.static_frame_theta(1.0, 1.0) - np.sqrt(2)) < 1e-15


def test_throat_flip_three_ways_and_duration_bound():
    Q, T = 1e3, 1e-10
    th = S.NearExtremalThroat(Q, T)
    assert abs(th.M - (Q + 2 * np.pi**2 * Q**3 * T**2)) / (2 * np.pi**2 * Q**3 * T**2) < 1e-3
    r0 = th.r_at_depth(Q)
    f0 = th.f(r0)
    assert abs(np.sqrt(f0) / th.sqrt_f_approx(Q) - 1) < 1e-4
    for c in (5.0, 20.0, 40.0, 100.0):
        pred = th.flip_condition(Q, c)["flips"]
        cft = S.cft_flip(c, r0, -np.sqrt(f0 / 2), np.sqrt(2) / Q, w=0.002 * Q, lam_max=10 * Q)
        vd = S.vaidya_flip(th.rp, th.rm, c, np.sqrt(f0) / Q, 0.002 * Q / np.sqrt(f0), [r0])
        assert cft["flips"] == pred == (len(vd["antitrapped"]) > 0)
        assert cft["max_dTheta"] < 1e-15                              # quantum focusing along the ray
        assert cft["duration"] < cft["duration_bound"]
    assert S.flip_temperature(40.0, Q) > T > S.flip_temperature(20.0, Q)
    # a bump pulse leaves the ray more focused than it started (quantum interest)
    b = S.cft_flip(100.0, r0, -np.sqrt(f0 / 2), np.sqrt(2) / Q, w=0.05 * Q, lam_max=20 * Q, kind="bump")
    assert b["flips"] and b["theta"][-1] < 2 * (-np.sqrt(f0 / 2)) / r0


def test_dichotomy_verdicts():
    assert S.wall_dichotomy(True, True, "de_sitter", 1.0, 0.0)["allowed"] is False
    assert S.wall_dichotomy(True, True, "flat", 0.0, 0.0)["allowed"] is True
    assert S.wall_dichotomy(True, True, "flat", 0.0, 0.0, neck_persists=False)["allowed"] is False
    assert S.wall_dichotomy(True, False, "de_sitter", 1.0, 0.0)["premise"] == "global hyperbolicity"
    assert S.wall_dichotomy(False, True, "crunch")["allowed"] is True
    assert S.wall_dichotomy(True, True, "de_sitter", 0.5, 1.0)["allowed"] is True
    assert S.wall_dichotomy(True, True, "de_sitter", 2.0, 1.0)["allowed"] is False


def test_cauchy_route_window_and_corona():
    w = S.cauchy_route_window(1.0, 0.3, 0.99, M_lo=0.6, M_hi=1.5, n_scan=31)
    assert w["exists"] and w["M_a"] < w["M_extremal_limit"] < w["M_b"]
    lo, hi = S.flip_between_horizons_bounds(1.0, 0.3, 0.99)
    assert lo <= w["M_a"] + 1e-6 and abs(w["M_b"] - hi) < 1e-5           # upper edge: the outer flip reaches r_-
    assert 0.85 < w["width"] / w["horizon_gap_centre"] <= 1.0           # width <= r_+ - r_- = 2 kappa_- r_-^2
    h = S.electroweak_corona_bound("higgs")
    assert abs(h["K"] - np.pi * S.ALPHA_EM / (2 * 0.129)) < 1e-12 and h["x_max"] < 0.25
    assert not S.electroweak_corona_bound("W")["reaches_decompactification"]
    hd = S.hiding_depth(1 / 32)
    assert abs(hd["depth_approx"] - (1 + np.log(16) / 2)) < 1e-12 and abs(hd["depth_exact"] - 1.9784) < 1e-3
    small = S.hiding_depth(1e-4)
    assert abs(small["depth_exact"] - small["depth_approx"]) < 0.03
    fl = S.antitrapped_floor(1.0)
    assert abs(fl["misner_sharp_kg"] - 6.733e26) / 6.733e26 < 1e-3 and 25 < fl["rho_quarter_GeV"] < 32
