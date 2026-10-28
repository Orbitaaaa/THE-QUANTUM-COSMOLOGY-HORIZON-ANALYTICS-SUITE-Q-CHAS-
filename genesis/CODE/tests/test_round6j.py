"""Round 6j: what the bridge invariant bounds (geometry, stored energy, supplied energy), and fifteen loopholes in four groups."""
import warnings

import numpy as np

from babyuniverse import accessibility6j as A, loopholes6j as L, nariai as NA

warnings.simplefilter("ignore")


def test_compactness_identity_and_misner_sharp():
    for E in (1e20, 1e40):
        c = A.compactness_identity(E, 1.0)
        assert abs(c["lhs"] / c["rhs"] - 1) < 1e-12
    c1 = A.compactness_identity(A.C_SI**4 / (2 * A.G_SI), 1.0)
    assert abs(c1["k"] - 1) < 1e-12 and 2.2e65 < c1["lhs"] < 2.4e65
    ms = A.misner_sharp_identity()
    assert ms["d_u"] and ms["d_v"] and ms["orthonormal"]


def test_which_density_peak_never_below_mean():
    rows = A.which_density()
    assert all(r["peak_over_mean"] >= 1 - 1e-3 and r["positive_energy_over_R_half"] >= 1 - 1e-6 for r in rows)
    neg = next(r for r in rows if r["profile"].startswith("negative core"))
    assert neg["positive_energy_over_R_half"] > 1.01             # negative energy raises the positive energy needed


def test_three_bounds_precursors_and_qei():
    b = A.three_bounds(1e-3, E_natural_J=6.0e40)
    assert b["supplied_min_J"] < 1e-2 * b["stored_min_J"]          # supplied energy is not bounded by the invariant
    p = A.precursor_window()
    assert 1e11 < p["M_min_survivor_kg"] < 3e11 and p["delta"] < 1e-8 and p["nariai_decision_time_s"] < 1e-20
    assert A.precursor_window(M_kg=A.M_SUN)["delta"] < 1e-27
    assert all(A.qei_negative_energy(R)["log10_shortfall"] > 6 for R in (1e-30, 1e-18, 1e-3))


def test_penrose_form():
    M = 1e3
    A_horizon = 16 * np.pi * (A.G_SI * M / A.C_SI**2) ** 2
    assert A.penrose_form(A_horizon, M)["satisfies"] and not A.penrose_form(4 * A_horizon, M)["satisfies"]


def test_strong_gravity_core_costs_more_than_it_buys():
    assert abs(A.strong_gravity_threshold() - 0.315) < 0.002
    d = A.dirichlet_check()
    assert abs(d["harmonic"] / d["formula"] - 1) < 1e-3 and d["linear_ramp"] > d["formula"]
    for F in (0.1, 1e-26, 1e-53):
        s = A.strong_gravity_core(F)
        assert s["horizon_class"] and s["core_compactness"] > 1 / s["invariant_gain"] ** 0 and s["core_compactness"] >= 0.75 * np.log(F) ** 2 - 1e-9


def test_pocket_expansion_is_negligible():
    assert A.pocket_expansion(500, 1e-18)["expansion_fraction"] < 1e-15
    assert not A.pocket_expansion(1.5e16, 1e-31)["anti_trapped"]


def test_partner_vacuum_corrected():
    assert abs(NA.partner_vacuum(NA.M_SUN_KG) - 0.329) < 0.002     # rho_f = 1/(72 pi M^2): 9 M^2 Lambda = 1
    assert abs(L.qcd_bag_partner()["M_sun"] - 5.16) < 0.05


def test_accretion_and_bag_star():
    a = L.accretion_requirement(10.0)
    assert a["log10_eddington_shortfall"] > 18 and a["merger_over_required"] < 0.01
    bs = L.bag_star(n=60)
    assert 0.5 < bs["compactness"] < 0.6 and bs["R_over_l_f"] < 0.35


def test_oshita_yokoyama_two_step_exponent():
    r = L.oy_channel(0.15, 1.0, 0.5)
    assert r["exists"] and 0.005 < r["Sigma"] < 0.03
    assert L.oy_expected_daughters(4e-4, 0.05, np.sqrt(5.2e10)) < 1e-100   # GUT vacuum: no daughters from evaporating holes


def test_wineglass_throat_area_law():
    for Lam in (1.0, 0.01):
        f = L.wineglass_family(Lam, np.geomspace(1e-4, 0.01, 6))
        assert abs(f["c_fit"] / (3 * np.pi**3) - 1) < 2e-3
        assert abs(f["S0"] + 12 * np.pi**2 / Lam) < 1e-9
    assert L.wineglass_action(1e-6, 0.1)["S"] - L.wineglass_action(0.0, 0.1)["S"] < 1e-3


def test_fmm_junction_and_extra_dimensions():
    j = L.fmm_junction(10.0)
    assert j["valid"] and abs(j["r0"] - (20.0) ** (1 / 3)) < 1e-12 and abs(j["kretschmann_S"] - 12.0) < 1e-9
    assert not L.fmm_junction(0.4)["valid"]
    b = L.add_black_holes(6, 8.0, 30.0)
    assert b["allowed_by_lhc"] and b["astro_bounds_weak"] and 1 < b["sigma_parton_pb"] < 1e3
    assert abs(L.add_rs_over_fundamental(5.0, 1.0, 6) - 1.0) < 0.05     # collider holes are at the fundamental scale
    assert L.add_semiclassical_mass(6, 3.0) > 5e3 and L.add_semiclassical_mass(4, 3.0) > 5e2
