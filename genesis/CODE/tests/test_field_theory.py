"""Round 5: material properties from a field theory (tilted quartic scalar) and the Route G scales."""
import warnings

import numpy as np

from babyuniverse import field_theory as F

warnings.simplefilter("ignore")


def test_derived_parameters_and_bounce():
    for lam, v, eps in ((1.0, 0.02, 0.3), (4.0, 0.01, 0.5)):
        d = F.derived_parameters(lam, v, eps)
        assert abs(d["B_flat_thin"] / (6745 / (eps**3 * lam)) - 1) < 0.01         # Proposition F1 closed form
        assert d["width_over_R0"] < 0.1 and d["one_loop_relative"] < 0.11
        Bg = F.gravitational_bounce(lam, v, eps)
        assert 0.8 < Bg / d["B_flat_thin"] < 1.0                               # gravity lowers the dS -> flat exponent slightly


def test_frontier_needs_strong_coupling():
    for eps in (0.3, 0.5):
        c = F.couplings_for(2.7e17, 2.17e17, eps)
        assert c["lam"] > 5e3 and not c["perturbative"]
        assert abs(c["lam"] * eps**3 / 690 - 1) < 0.05
    pf = F.perturbative_frontier(4 * np.pi, 0.5)
    assert 1000 < pf["B_lab_min"] < 1200 and pf["log10_P_max"] < -450


def test_end_to_end_chain_and_genesis_scales():
    r = F.laboratory_exponent_from_theory(1.0, 0.02, 0.3)
    assert r["E"] > 10 and abs(r["Phi"] - 1) < 0.05 and r["B_lab"] > 1e5
    g = F.genesis_scales(1e10, 1e5)
    assert g["H0_over_Lambda"] < 0.01 and g["h_m2"] < 1e-40
