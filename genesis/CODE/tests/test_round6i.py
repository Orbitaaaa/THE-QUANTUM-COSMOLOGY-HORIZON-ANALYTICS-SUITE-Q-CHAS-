"""Round 6i: seven routes to a daughter without a manufactured false-vacuum ball (curvature activation, monopole cores, breathing
bubbles, gauged/vacuum-dependent ghost condensation, resonant tunnelling, classical transitions, four-form membranes), and the bridge
invariant that every classical route meets."""
import warnings
from dataclasses import replace

import numpy as np

from babyuniverse import gauss_bonnet as GB, stable_shells as SS, routes6i as R6
from babyuniverse.model import ThinWallModel
from babyuniverse import tunnelling as T

warnings.simplefilter("ignore")


def test_gauss_bonnet_stress_identities():
    assert GB.verify_divergence_identity()            # div Gamma = -(GB/2) grad F, off shell
    assert GB.null_component_closed_form()            # -Gamma_uu = 4 (2m/r^3) nabla_u nabla_u F + ...
    assert GB.ks_equations()["total_derivative_ok"]   # sqrt(-g) GB = d/dt[8 a'(1 + b'^2)] for Kantowski-Sachs
    assert GB.nariai_gb_invariant() == 8.0 and GB.sds_gb_invariant(0.1, 1.0, 1.0) > 0


def test_weighted_gb_flux_is_pro_collapse_when_F_grows():
    u = np.linspace(0, 40, 40001)
    w = GB.weighted_gb_deficit(0.01 * np.tanh(u / 3) ** 2, u, 3.2)
    assert w["direct"] < 0 and abs(w["direct"] - w["by_parts"]) < 1e-2 * abs(w["direct"])
    w2 = GB.weighted_gb_deficit(-0.01 * np.tanh(u / 3) ** 2, u, 3.2)    # F decreasing: pro-daughter sign
    assert w2["direct"] > 0


def test_scalarised_nariai_throat_never_makes_a_daughter():
    out = [GB.ks_run(GB.coupling_dy, GB.potential_flat(), l, 1e-3, 0.0)["outcome"] for l in (0.05, 0.5, 2.0)]
    assert out == ["black hole", "black hole", "kinetic degeneracy"]
    rev = [GB.ks_run(GB.coupling_reverse_saturating(ps), GB.potential_hilltop(4 * l * 2.0, 0.9), l, 0.0, 1e-3)["outcome"]
           for l, ps in ((2.0, 0.2), (5.0, 0.2), (10.0, 0.4))]
    assert "daughter" not in rev
    assert GB.ks_run(GB.coupling_reverse, GB.potential_hilltop(22.0, 0.9), 0.0, 0.0, 1e-3)["outcome"] == "black hole"


def test_laboratory_curvature_cannot_activate_the_coupling():
    rows = {r["environment"]: r for r in GB.activation_report()}
    assert rows["earth surface"]["log10_shortfall"] > 18 and rows["dense lab mass"]["log10_shortfall"] > 18
    assert rows["neutron star surface"]["log10_shortfall"] < 2
    assert GB.hysteresis_bound()["l_f_over_l_curv_min"] >= 1.0


def test_shell_model_reproduces_thin_wall_code():
    m = ThinWallModel(chi=1.0, kappa=0.5, M=0.05)
    s = SS.Shell(chi=1.0, M=0.05, k_w=0.5)
    R = np.linspace(0.05, 3, 9)
    assert np.allclose(m.V(R), s.V(R))
    tp = m.turning_points()
    assert abs(T.B_canonical(m) - s.B(tp[0], tp[-1])) < 1e-9


def test_breathing_bubbles_need_negative_dust_and_cross_zero_energy():
    sc = SS.breathing_random_scan(800)
    assert sc["wells"] > 0 and sc["wells_crossing_zero_energy"] == sc["wells"]
    assert set(map(tuple, sc["sign_pairs"])) <= {(1, -1), (-1, 1)}
    hid = SS.classify_breathing(SS.Shell(chi=1.0, M=0.05, k_w=1.0, m_d=-0.1))
    eat = SS.classify_breathing(SS.Shell(chi=1.0, M=0.05, k_w=-1.0, m_d=0.1))
    assert hid["outcome"] == "hidden daughter" and eat["outcome"] == "parent-eating bubble"
    assert 0.4 < hid["B"] < 1.0 and hid["crosses_zero_energy"]


def test_monopole_core_thresholds_and_channels():
    assert SS.monopole_static(0.1) is not None and SS.monopole_static(0.3) is None       # eta_c ~ 0.25 m_P (lam = 1, e = 0.5)
    t = SS.monopole_tunnelling(0.01)
    Q = np.sqrt(4 * np.pi) / 0.5
    assert t["M"] < 0.1 * Q                       # own mass: no bridge, the beta_out < 0 branch is disconnected from us
    assert 1.3 < t["B"] / SS.monopole_cdl(1.0) < 1.8
    sh = replace(SS.monopole_shell(0.01), M=1.05 * Q)
    regs = sh.regions(0.01, 3e5)
    assert "".join("A" if r["allowed"] else "F" for r in regs) == "AFA" and sh.B(*regs[-2]["R"]) > 2000
    tn = SS.monopole_tunnelling(0.01, m_d=-0.2 * 4 * np.pi * 0.01 / 0.5)
    assert 0.75 * t["B"] < tn["B"] < t["B"]       # NEC admixture lowers B only modestly


def test_resonance_cannot_remove_the_last_barrier():
    e = SS.resonance_exponent(30.0, 10.0)
    assert e["off"] == 40.0 and e["on"] == 20.0
    assert SS.resonance_exponent(5.0, 10.0)["on"] == 10.0
    found = SS.find_two_well(n=6000)
    assert found and all(SS.resonance_scan(x["shell"])["on"] >= SS.resonance_scan(x["shell"])["B2"] - 1e-12 for x in found[:3])


def test_gauged_and_vacuum_dependent_condensate():
    assert all(R6.gauged_window(c)["empty"] for c in (3.0, 30.0))
    assert R6.coulomb_ratio(0.1, 1.0, 100.0) > 1.0
    v = R6.vacuum_dependent()
    assert 1e-18 < v["r_N_m"] < 3e-18 and 5e8 < v["M_N_kg"] < 2e9
    assert R6.strong_coupling_scale_overall_factor((100 / 1e13) ** 4, 1e13) < 1e-60


def test_collision_pockets_cannot_reach_the_hidden_branch():
    rows = R6.pocket_scan()
    assert rows and all(r["beta_out_R1"] > 0 and r["f_out_R1"] > 0 for r in rows) and all(r["beta_out_R2"] < 0 for r in rows)
    assert R6.pocket_completion(0.2, 0.05, 1.0, 0.5, 0.1)["branch"] == "bound"
    lem = R6.side_change_lemma(n=600)
    assert lem["sign_changes"] > 0 and lem["explained"] == lem["sign_changes"]


def test_membrane_and_bridge_invariant():
    mb = SS.membrane_branches(n=30)
    assert mb["rows"] and all(r["hidden"] for r in mb["rows"]) and mb["M_cr_upper"] is not None
    assert 0.9 < mb["rows"][0]["B"] / SS.cdl_flat(1e-3, 0.05 * 0.075) < 1.2       # no free lunch (R10, R37)
    b = R6.bridge_invariant()
    assert 2.2e65 < b["E_sqrt_rho_SI"] < 2.4e65 and b["log10_gap"] > 39
