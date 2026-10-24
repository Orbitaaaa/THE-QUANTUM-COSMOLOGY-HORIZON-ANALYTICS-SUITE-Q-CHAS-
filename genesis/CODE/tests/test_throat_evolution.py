"""Round 6c: semiclassical Einstein-Maxwell-scalar-CFT evolution of the throat, and the obstruction to growth (P-R6f)."""
import warnings

import numpy as np
from scipy.optimize import brentq

from babyuniverse import throat_evolution as TE, obstructions as O

warnings.simplefilter("ignore")
Q = 10.0


def test_reissner_nordstrom_conservation_and_convergence():
    M, r0 = 10.1, 13.0
    f0 = 1 - 2 * M / r0 + Q**2 / r0**2
    dev, R = {}, {}
    for h in (0.2, 0.1):
        a, b, _ = TE.characteristic_data(Q, 0.0, r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, h, int(20 / h), int(20 / h))
        res = TE.evolve(a, b, Q, 0.0, TE.Potential(), h)
        dev[h] = np.nanmax(np.abs(TE.misner_sharp(res) - M))
        R[h] = res["r"]
    assert dev[0.1] < 1e-5 and 3.0 < dev[0.2] / dev[0.1] < 5.0          # second order
    assert np.nanmax(np.abs(R[0.2] - R[0.1][::2, ::2])) < 1e-4


def test_backreacted_qnec_flip_matches_prediction():
    M = 10.0001
    rp, rm = M + np.sqrt(M * M - Q * Q), M - np.sqrt(M * M - Q * Q)
    f0 = 1e-4
    r0 = brentq(lambda r: (r - rp) * (r - rm) / r**2 - f0, rp * (1 + 1e-12), 2 * Q)
    for c, flips in ((50.0, True), (30.0, False)):                    # kick/need = 1.31 and 0.79
        Y, up, wp = 0.1, 5.0, 0.5
        yp = lambda u: 0.5 * Y * (1 + np.tanh((u - up) / wp))
        ypp = lambda u: 0.5 * Y / wp / np.cosh((u - up) / wp) ** 2
        tu = lambda u: c / (24 * np.pi) * (0.5 * yp(u) ** 2 - ypp(u))
        a, b, _ = TE.characteristic_data(Q, c, r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, 0.1, 400, 400, t_u=tu)
        res = TE.evolve(a, b, Q, c, TE.Potential(), 0.1)
        cl = TE.classify(res)
        assert (cl["n_anti"] > 0) == flips
        assert TE.constraint_residual(res, t_u=tu)["Evv_rel"] < 1e-4
        assert 9.0 < cl["r_min"] and cl["r_max"] < 11.0                   # the throat persists


def _fv(r0, ru0, rv0, u2, v2, Lu, Lv, h=0.05):
    x, vphi, eps = 0.5, 0.005, 2.0
    rho_f = x / (8 * np.pi * Q**2)
    lam = 4 * rho_f / (eps * vphi**4)
    w = np.sqrt(2) / (np.sqrt(lam) * vphi)
    phi_v = lambda vv: vphi * np.tanh((np.asarray(vv, float) - v2) / w)
    phi_u = lambda uu: vphi * np.tanh((np.asarray(uu, float) - u2) / w)
    a, b, _ = TE.characteristic_data(Q, 0.0, r0, ru0, rv0, h, int(Lv / h), int(Lu / h), phi_v=phi_v, phi_u=phi_u)
    return TE.evolve(a, b, Q, 0.0, TE.Potential(lam, vphi, rho_f), h, r_inf=1e4)


def test_normality_bound_and_superhubble_growth():
    fa = 0.1
    r0 = Q / (1 - np.sqrt(fa))
    C = r0 * (1 + Q**2 / r0**2 - 0.5 / Q**2 * r0**2 / 3 - fa)
    rc = brentq(lambda r: O.rnds_f(r, C, Q, 0.5 / Q**2), r0, 10 * r0)
    v_c = (rc - r0) / (np.sqrt(fa) / 2)
    ok = TE.classify(_fv(r0, -np.sqrt(fa) / 2, np.sqrt(fa) / 2, 30.0, v_c - 0.6, 6, v_c + 3))["n_anti"]
    bad = TE.classify(_fv(r0, -np.sqrt(fa) / 2, np.sqrt(fa) / 2, 30.0, v_c + 1.0, 6, v_c + 4))["n_anti"]
    assert ok == 0 and bad > 0                                            # Theorem 6c-1, numerically
    big = _fv(15.0, 0.05 / 1.2, 0.3, 50.0, 50.0, 90, 90)
    small = _fv(15.0, 0.05 / 1.2, 0.3, 10.0, 10.0, 50, 50)
    assert np.any(big["end"] == 2) and not np.any(small["end"] == 2)      # super-Hubble inflates to I+


def test_obstruction_numbers():
    for x, top in ((0.3, 0.80), (0.5, 0.55), (1.0, 0.43)):
        assert O.growth_ratio_max(Q, x, 1.01 * Q)["ratio_max"] < top < 1
    s = O.segment_Lmax(Q, 0.5, 1.01 * Q, 1e-4)
    assert abs(s["L_max"] / s["approx"] - 1) < 0.01
    ci = O.corona_identity(3.0)
    assert abs(ci["rho_vac"] / ci["rho_B"] - 1) < 1e-14
    assert O.best_eta_window(Q, 1e-4)["inside_eft_at_eta1"] and not O.best_eta_window(Q, 0.3)["inside_eft_at_eta1"]
    assert O.hiding_growth(0.26)["beyond_mouth"]
    assert abs(O.lifting_focusing(0.5, 2 * np.sqrt(3 * Q**2 / 0.5), Q) - np.sqrt(3 * 0.5)) < 1e-12


def test_nariai_throat_positive_result():
    grow = TE.nariai_false_vacuum(delta=1e-8, Nflux=0.04, h=0.1, Lu=120, Lv=120)
    assert grow["grew"] and 40 < grow["u_growth"] < 60
    assert grow["true_vacuum_untrapped_large_r"] == 0 and grow["true_vacuum_trapped"] > 0     # hidden behind a trapped region
    assert 0.011 < grow["delta_m"] / grow["M_N"] < 0.013                                        # local energy ~1.2% M_N (invariant: sqrt(f0) times that)
    weak = TE.nariai_false_vacuum(delta=1e-8, Nflux=0.01, h=0.1, Lu=120, Lv=120)
    assert not weak["grew"]
    qnec = TE.nariai_false_vacuum(delta=1e-8, c=24 * np.pi * 0.06 / 0.1, Y=0.1, wp=2.0, h=0.1, Lu=120, Lv=120)
    assert not qnec["grew"]                                                                     # QNEC pulse: no growth
