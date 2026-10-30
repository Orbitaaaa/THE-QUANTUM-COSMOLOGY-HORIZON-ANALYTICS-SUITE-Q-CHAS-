"""Round 8: model U against current and future data."""
import os
import warnings

import numpy as np

from babyuniverse import data8 as E
from babyuniverse import collider6l as C

warnings.simplefilter("ignore")
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "DATA")


def test_convention_and_monojet_translation():
    c = {r["n"]: r["MD_DL_over_MD_PDG"] for r in E.convention()}
    assert abs(c[6] - 0.3770) < 1e-3 and abs(c[2] - 0.8932) < 1e-3
    b = {r["n"]: r["MD_DL_TeV"] for r in E.monojet_bounds()}
    assert 1.9 < b[6] < 2.1 and 3.1 < b[4] < 3.3


def test_lhc_edge_safety_and_future_reach():
    L = C.load_lumi(os.path.join(DATA, "round6l_lumi.json"))
    opt = E.lhc_lower_edge(L, 6, 8.0, "hayward", 0.5, "optimistic")["M_ext_min_TeV"]
    con = E.lhc_lower_edge(L, 6, 8.0, "hayward", 0.5, "conservative")["M_ext_min_TeV"]
    assert 4.0 < con < opt < 10.0
    s = E.remnant_safety(6, 8.0, 0.6, "hayward")
    assert s["equilibrium_excess_fraction"] < 1e-6 and s["equilibrium_T_GeV"] < 1.0
    fr = {(r["machine"], r["N_required"]): r["M_ext_reach_TeV"] for r in E.future_reach(L)}
    assert fr[("HL-LHC 14", 10.0)] < 8.0 < 35.0 < fr[("FCC-hh 85", 10.0)] < fr[("SPPC/FCC 100", 10.0)]
