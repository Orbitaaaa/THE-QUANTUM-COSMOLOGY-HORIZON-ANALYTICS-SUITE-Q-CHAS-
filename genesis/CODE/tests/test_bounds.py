"""Roadmap Part A: analytic suppression bounds (Lemma S1, Theorems S2-S5, Corollary S6)."""
import warnings

import numpy as np

from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse import bounds as Bd, tunnelling as T, contours as Cc

warnings.simplefilter("ignore")


def test_sine_and_elementary_bounds_hold():
    for gam in (-0.9, -0.3, 0.2, 0.7, 0.95):
        for fac in (1.01, 1.5, 10.0):
            m = model_from_dimensionless(gam, (2 - 2 * gam) * fac)
            Ban, Bsq = T.B_canonical(m), Cc.B_sector_blind(m)
            Bm, Be = Bd.sine_bound(m), Bd.elementary_bound(m)
            assert Ban >= Bsq * (1 - 1e-10) and Bsq >= Bm * (1 - 1e-10) and Bm >= Be * (1 - 1e-10)


def test_elementary_bound_is_sharp_for_thin_barriers():
    m = ThinWallModel(1.0, 0.003, 0.1 * ThinWallModel(1.0, 0.003, 1e-9).M_critical)
    r = Bd.elementary_bound(m) / Bd.thin_barrier_asymptote(m)
    assert 0.97 < r <= 1.0


def test_universal_function_Phi():
    assert abs(Bd.Phi(1e4) - Bd.Phi_large_E(1e4)) < 1e-8
    assert abs(Bd.Phi(1e-3) - Bd.Phi_small_E(1e-3)) < 1e-2
    vals = [Bd.Phi(E) for E in np.geomspace(1e-3, 1e4, 15)]
    assert all(v > 1 for v in vals) and all(a > b for a, b in zip(vals, vals[1:]))


def test_thin_barrier_law_and_its_error_scaling():
    errs = []
    for eps in (0.03, 0.01, 0.003):
        m = ThinWallModel(1.0, eps, 0.5 * ThinWallModel(1.0, eps, 1e-9).M_critical)
        errs.append(abs(T.B_canonical(m) / Bd.thin_barrier_B(m) - 1))
    assert errs[-1] < 1e-5
    assert errs[0] / errs[1] > 5 and errs[1] / errs[2] > 5          # O(eps^2)
    m = ThinWallModel(1.0, 0.003, 0.9 * ThinWallModel(1.0, 0.003, 1e-9).M_critical)
    assert abs(T.B_canonical(m) / Bd.thin_barrier_asymptote(m) - 1) < 1e-4


def test_flip_class_theorem():
    for gam in (-0.4, 0.1, 0.3, 0.6, 0.9):
        for fac in (1.003, 1.05, 1.2, 1.5, 2.0, 5.0, 50.0):
            m = model_from_dimensionless(gam, (2 - 2 * gam) * fac)
            if abs(m.M / max(Bd.M_sd(m), 1e-300) - 1) < 1e-6:
                continue
            assert Bd.flip_class_predicted(m) == Cc.flip_class(m)
    m = model_from_dimensionless(0.4, Bd.E_sd(0.4))
    assert abs(m.M / Bd.M_sd(m) - 1) < 1e-9
    assert abs(T.barrier(m)[1] * m.chi - 1) < 1e-6                  # boundary: R2 = 1/chi


def test_euclidean_time_negative_and_omega_criterion():
    for gam in (-0.9, 0.0, 0.5, 0.95):
        for fac in (1.02, 1.3, 2.0, 20.0, 2000.0):
            E = (2 - 2 * gam) * fac
            assert Bd.delta_tE_dimless(gam, E) < 0
            if Bd.omega_condition(gam, E):
                assert Bd.pointwise_D_min(gam, E) >= 0
        assert Bd.omega_condition(gam, (2 - 2 * gam) * 1.62)


def test_double_flip_exponent_exceeds_cdl():
    for gam in (0.25, 0.5, 0.8):
        m0 = model_from_dimensionless(gam, 10.0)
        B0 = T.B_canonical(ThinWallModel(m0.chi, m0.kappa, 0.0))
        for fac in (1.7, 3.0, 30.0):
            E = (2 - 2 * gam) * fac
            if E <= Bd.E_sd(gam):
                continue
            assert T.B_canonical(model_from_dimensionless(gam, E)) >= B0


def test_laboratory_law_reproduces_phase_vi_optimum():
    lab = Bd.laboratory_law(9e17, 0.8 * 9e17, 10.0)
    assert abs(lab["B_law"] - lab["B_law_identity"]) < 1e-9 * lab["B_law"]
    assert 100 < lab["B_law"] < 112                                  # exact thin-wall value 107 (Phase VI)
    gut = Bd.laboratory_law(1e16, 1e16, 0.5 * 5.59 / 2.176e-5)
    assert abs(gut["B_law"] / 4.4e6 - 1) < 0.02                      # TUNNELLING T6 table
    assert abs(Bd.minimum_exponent(1.0, 1.0, 1.0) - 1.3165) < 1e-3


def test_phi_lower_bound_and_frontier():
    for s in (1e-2, 0.3, 1.0, 3.0, 30.0, 3e3):
        assert Bd.h_even_defect(s) >= -1e-12
    for E in (1e-3, 0.1, 1.0, 5.0):
        assert Bd.Phi(E) >= 8 / (3 * E) + 2 / 3 - 1e-9
    assert abs(Bd.FRONTIER_CONST - 4 / 3 * np.sqrt(np.pi / 2)) < 1e-12
    lab = Bd.laboratory_law(9e17, 0.8 * 9e17, 10.0)
    assert lab["B_law"] >= Bd.frontier_bound(lab["B_flat"], 10.0, lab["R_M_lP"]) * (1 - 1e-9)
    opt = Bd.frontier_optimum(10.0, 1.0, 100.0)
    assert 52 < opt["B_min"] < 54 and opt["kappa_over_chi"] < 0.1


def test_phi_bound_certificate():
    c = Bd.phi_bound_certificate()
    assert c["identity"] and c["cofactor_positive_on_P_ge_1"] and c["c1_roots_ge_1"] == [] and c["f_reduced_matches"]


def test_frontier_against_exact_exponent():
    from babyuniverse import observability as O
    o = Bd.frontier_optimum(10.0, 1.0, 100.0)
    lab = O.laboratory_case(o["E_vac_GeV"], o["E_wall_GeV"], 1.0)
    assert lab["B_an"] >= o["B_min"]
    assert abs(lab["B_an"] - o["B_min"] - 2 / 3 * lab["B_CdL"]) < 0.05 * lab["B_an"]
    assert abs(lab["B_an"] - lab["B_sq"]) < 1e-9 * lab["B_an"]


def test_strip_boundary_limit_and_certificates():
    c = Bd.boundary_certificate()
    assert c["factorization"] and c["h_definition"] and c["barrier_split"]
    for g in (-0.95, -0.5, 0.0, 0.5, 0.9, 0.99):
        lim = Bd.boundary_limit(g)
        assert lim < 0
        J = Bd.delta_tE_dimless(g, (2 - 2 * g) * (1 + 1e-4))
        assert abs(J - lim) < 2e-3 * max(1.0, abs(lim)) + 1e-4
    assert abs(Bd.boundary_limit(0.999)) < 1e-3                     # -> 0- as gamma -> 1
    for g in (-0.9, 0.0, 0.9):
        cert = Bd.strip_lp_certificate(g, (2 - 2 * g) * 1.05)
        assert cert["margin"] > 0 and cert["verified_margin"] > 0
