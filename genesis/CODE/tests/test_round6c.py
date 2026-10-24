"""Round 6c: P-S3 in closed form, the quantum throat (P-R6g, P-R6h, P-A1), and the robustness of the Nariai-throat daughter."""
import warnings

import numpy as np
import sympy as sp

from babyuniverse import crosscheck as X, quantum_throat as QT, throat_evolution as TE, nariai as NA, greybody as GB

warnings.simplefilter("ignore")


def test_s3_coefficient_closed_form():
    assert X.s3_small_E_limit() == sp.Rational(-13, 18)
    # Round 3 multiprecision values at kappa/chi = 0.01 carry O(1e-4) corrections; the closed form is their kappa -> 0 limit
    for E, c3 in ((0.3, -0.6696605441209302), (1.0, -0.5611657273395354), (3.0, -0.35100061535178284), (30.0, -0.02783760208544913)):
        assert abs(X.s3_coefficient(E) - c3) < 1e-4
    assert abs(X.s3_coefficient(1e-4) + 13 / 18) < 1e-4
    assert abs(X.s3_coefficient(1e6) * 1e9 + 16 / np.pi) < 2e-3            # tail -(16/pi) E^(-3/2)


def test_landau_levels_and_strong_qnec():
    sp_ = QT.lll_spectrum(7, a=2.0)
    assert sp_["zero_modes"] == 7 and abs(sp_["gap"] - np.sqrt(8) / 2) < 1e-12
    assert sp_["levels"][0] == (np.sqrt(8) / 2, 9)
    assert QT.lll_central_charge(7, 3) == 21
    assert QT.strong_qnec_curved()["saturated"]


def test_schwarzian_moments():
    T, Q = 0.02, 3.0
    a, b = QT.schwarzian_fluctuation(T, Q), QT.schwarzian_numeric(T, Q)
    assert abs(a["dE_over_E"] - a["dE_over_E_formula"]) < 1e-12
    assert abs(b["E_mean"] / a["E_mean"] - 1) < 1e-3 and abs(b["var"] / (4 * np.pi**2 * Q**3 * T**3 + 1.5 * T**2) - 1) < 1e-2


def test_weighted_anec_bound_is_sharp():
    c, tau = 100.0, 35.0
    bound = QT.weighted_anec_bound(c, tau)
    assert abs(QT.weighted_anec_of_pulse(c, 1 / tau, tau) / bound - 1) < 5e-3           # saturated by A = 1/tau
    for A in (0.5 / tau, 2 / tau, 4 / tau):
        assert QT.weighted_anec_of_pulse(c, A, tau) > bound
    t = QT.quantum_tipping_c(17.32)
    assert 0.85 < t["c_over_rN2"] < 0.95 and t["species_length_over_rN"] > 0.9


def test_nariai_daughter_survives_one_loop_and_needs_nec_violation():
    assert TE.nariai_false_vacuum(Nflux=0.04, c=20.0)["grew"]                          # Bunch-Davies back-reaction, k/r_N^2 = 1.8e-3
    assert not TE.nariai_false_vacuum(Nflux=-0.2)["grew"]                              # positive energy never tips (Farhi-Guth)
    tu, _ = TE.conformal_pulse(100.0, 0.05, 2.0, 1.0, 60.0, 20.0)
    o = TE.nariai_false_vacuum(c=100.0, tu_func=tu)
    assert o["net_energy"] > 0 and o["neg_energy"] > 0.06 and not o["grew"]            # ANEC-respecting CFT pulse does not tip


def test_tipping_law_and_invariant_cost():
    # the weighted deficit separates late pulses and negative/positive pairs exactly as the simulations do
    tau = NA.TAU_OVER_RN * 17.3205
    W = lambda f: NA.weighted_deficit(f, tau)
    g = lambda A, up: (lambda u: -A * np.exp(-((u - up) / 1.0) ** 2))
    th = NA.WTH_OVER_RN * 17.3205
    assert W(g(0.02, 2)) > th and W(g(0.0175, 2)) < th and W(g(0.02, 30)) < th and W(g(0.04, 30)) > th
    pair = lambda D: (lambda u: g(0.04, 2)(u) - g(0.04, 2 + D)(u))
    assert W(pair(20)) < th < W(pair(40))
    d = TE.nariai_false_vacuum(Nflux=0.02)
    assert d["grew"] and d["true_vacuum_untrapped_large_r"] == 0
    assert abs(d["dm_misner_sharp"] / (np.sqrt(d["f0"]) * 0.02 * np.sqrt(np.pi)) - 1) < 5e-3     # Gaussian tail cut at u = 0
    assert abs(NA.invariant_cost(0.032, 1e-8) / d["M_N"] - 4.5e-7) < 0.2e-7
    assert 260 < NA.flat_space_flip(17.3205) / 0.032 < 280


def test_preparation_and_partner_vacuum():
    assert TE.nariai_false_vacuum(Nflux=0.04, u_inner=150.0)["grew"]                 # hole formed inside the false vacuum: grows
    assert not TE.nariai_false_vacuum(Nflux=0.16, u_inner=40.0)["grew"]              # thin layer around a true-vacuum hole: does not
    assert abs(NA.partner_vacuum(NA.M_SUN_KG) - 0.433) < 0.002
    assert abs(NA.partner_mass(1e16) - 3.74e-3) < 0.02e-3
    gc = NA.ghost_condensate_window(17.3205)
    assert abs(gc["need_over_rho_f"] - 4.6e-3) < 0.1e-3 and gc["ceiling_over_need_lH"] > 100


def test_greybody_limits():
    assert abs(GB.transmission(0, 0.005) / (16 * 0.005**2) - 1) < 0.05                 # sigma_abs -> horizon area
    assert GB.transmission(0, 0.6) > 0.999 and GB.transmission(2, 0.1) < 1e-4
    g = [GB.transmission(l, 0.3) for l in (0, 1, 2)]
    assert g[0] > g[1] > g[2] and abs(g[1] - 0.5226) < 5e-3
