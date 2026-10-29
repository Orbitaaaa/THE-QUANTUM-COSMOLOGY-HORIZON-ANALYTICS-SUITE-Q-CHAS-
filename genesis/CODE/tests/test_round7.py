"""Round 7: distinguishable predictions (lab-frame events)."""
import os
import warnings

import numpy as np

from babyuniverse import predictions7 as Q
from babyuniverse import collider6l as C
from babyuniverse import pipeline6o as P

warnings.simplefilter("ignore")
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "DATA")


def test_model_U_and_GR_topologies_differ():
    L = C.load_lumi(os.path.join(DATA, "round6l_lumi.json"))
    rU = Q.LabGenerator(P.Generator(L, 100.0, 6, 8.0, 0.6, "hayward", seed=3), seed=3).run(400)
    rG = Q.LabGenerator(P.Generator(L, 100.0, 6, 8.0, 0.6, "gr", seed=3, M_th_TeV=11.5), seed=3).run(400)
    sU, sG = Q.summary(rU), Q.summary(rG)
    assert sU["mean_N_obj"] < sG["mean_N_obj"] and sU["median_ST_TeV"] < 0.3 * sG["median_ST_TeV"]
    assert sU["mean_met_frac"] > 0.7 and sG["mean_met_frac"] < 0.4
    assert 0.1 < sU["mean_beta_remnant"] < 0.8
    assert Q.kl_separation(rU, rG)["events_5sigma"] < 30


def test_background_normalisation_order_of_magnitude():
    L = C.load_lumi(os.path.join(DATA, "round6l_lumi.json"))
    b = Q.qcd_tail_pb(L, 13.0, 8.0)
    assert 1e-6 < b < 1e-3                     # a few events in 138 fb^-1 near m ~ 8 TeV at the LHC, as observed in dijet spectra
    assert Q.qcd_tail_pb(L, 100.0, 20.0) < Q.qcd_tail_pb(L, 100.0, 5.0)
