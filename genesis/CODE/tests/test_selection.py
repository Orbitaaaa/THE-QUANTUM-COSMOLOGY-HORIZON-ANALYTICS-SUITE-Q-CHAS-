"""Problem A: charge selection rule and its exact single-flip toy."""
import numpy as np

from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse import selection as S


def test_toy_solutions_and_laurent_structure():
    r = S.toy_report(us=(-3, -1, 1, 3))
    for name, d in r.items():
        assert d["residual"] < 1e-30
    rows = {name: d["rows"] for name, d in r.items()}
    # J, Y: q^0 charge in the lab, q^{+1} in the daughter (they grow like e^{+pi u/h})
    for name in ("J", "Y"):
        lab, dau = rows[name][0], rows[name][-1]
        assert abs(abs(lab["K"]) - abs(rows[name][1]["K"])) < 1e-12 * abs(lab["K"])
        assert abs(dau["K_over_q"] - rows[name][2]["K_over_q"]) < 1e-9 * dau["K_over_q"]
    # the daughter-decaying solution carries q^{-1} in the daughter and is non-geometric (q^{-2}) in the lab
    dec = rows["decaying"]
    assert abs(dec[2]["K_times_q"] / dec[3]["K_times_q"] - 1) < 1e-9
    assert dec[0]["K_times_q"] > 1e60 and dec[0]["abs_psi"] > 1e60
    # H1: bounded in the lab, decaying in the daughter, Casoratian identically zero
    h1 = rows["H1"]
    assert h1[0]["abs_psi"] < 1 and h1[-1]["abs_psi"] < 1e-30
    assert max(abs(x["K"]) for x in h1) < 1e-40


def test_toy_amplitude_follows_analytic_angle():
    c = S.toy_amplitude_check()
    assert abs(c["measured_exponent"] - c["analytic_exponent"]) < 0.02 * c["analytic_exponent"]
    assert abs(c["measured_exponent"] - c["sector_blind_exponent"]) > 0.1 * c["analytic_exponent"]


def test_charge_signs_of_physical_daughters():
    s = S.charge_signs(ThinWallModel(1.0, 1.0, 0.1545))
    assert s["n_daughter"] == 1 and s["P_lab"] > 0 and s["P_daughter"] < 0
    for m in (model_from_dimensionless(0.6, 2.4), model_from_dimensionless(0.3, 5.0)):
        s = S.charge_signs(m)
        assert s["n_daughter"] == 0 and s["P_lab"] > 0 and s["P_daughter"] > 0 and s["dRdt_daughter"] < 0
