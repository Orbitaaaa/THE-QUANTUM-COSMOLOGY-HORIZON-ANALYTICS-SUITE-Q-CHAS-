"""Round 6d: independent solver, bridge-free data, the complete route, ghost-condensate stability, the solution's generalised entropy, and the
strong-QNEC bound."""
import warnings

import numpy as np

from babyuniverse import (throat_evolution as TE, throat_mol as MOL, ghost_condensate as GC, gsl_solution as GS, preparation as PR,
                          quantum_throat as QT, collapse as CO)

warnings.simplefilter("ignore")
LIGHT = dict(delta=1e-4, eps=4.0, vphi=0.0015, h=0.04)
gauss = lambda A, up=2.0: (lambda u: -A * np.exp(-((np.asarray(u, float) - up) / 1.0) ** 2))


def test_independent_solver_reproduces_the_daughter():
    Q, M, r0 = 10.0, 10.1, 13.0
    f0 = 1 - 2 * M / r0 + Q**2 / r0**2
    err = []
    for h in (0.2, 0.1):
        d = MOL.characteristic_data(r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, h, int(20 / h), int(20 / h), lambda x: 0.0, lambda x: 0.0, n_sub=4)
        err.append(np.nanmax(np.abs(TE.misner_sharp(MOL.evolve(d, h, Q=Q)) - M)))
    assert err[1] < 2e-6 and 3.0 < err[0] / err[1] < 5.0                       # second order, independent scheme
    a = MOL.nariai(flux_u=gauss(0.04))
    assert a["grew"] and a["true_vacuum_untrapped_large_r"] == 0
    assert not MOL.nariai(flux_u=gauss(0.0175))["grew"]


def test_bridge_free_configuration_grows_hidden_and_respects_exterior_gsl():
    o0 = TE.nariai_false_vacuum(**LIGHT)
    assert o0["rv_data"][-1] > 0 and not o0["grew"]                              # exterior normal on the data: no Einstein-Rosen bridge
    o = TE.nariai_false_vacuum(Nflux=0.08, **LIGHT)
    assert o["grew"]
    eh = GS.parent_event_horizon(o)
    assert eh["u_EH"] is not None and eh["deficit"] == 0.0                     # the parent's horizon area never decreases
    dh = GS.daughter_horizons(o, n_points=6)
    assert 0 < dh["max_deficit"] < 0.01 * dh["S_N"]                            # deficit only on daughter-observer horizons, ~ sqrt(delta) S_N


def test_complete_route_from_data_without_a_black_hole():
    E = 7.258540773391722
    meta = PR.build(E, v_shell=4.5, w_shell=0.25, Lv=1500.0, Lu=2.0)[4]
    assert meta["normal_everywhere"] and meta["wall_found"] and 0.005 < meta["delta"] < 0.012
    kw = dict(u_p=1.0, w_p=0.1, v_shell=4.5, w_shell=0.25, Lv=600.0, Lu=30.0)
    res0, cl0, g0, _ = PR.run(E, N_amp=0.0, **kw)
    assert not g0["grew"] and np.sum((cl0["cls"] == -1) & (res0["phi"] < 0)) > 1000     # the shell makes a black hole inside the false vacuum
    res, cl, g, _ = PR.run(E, N_amp=10.0, **kw)
    assert g["grew"] and g["r_max_antifalse"] > 1000
    R = res["r"]
    js = np.where(np.isfinite(R[50]))[0]
    assert R[50, js.max()] < 1.0                                                # daughter-era rows end at the parent's singularity: hidden


def test_ghost_condensate_window_and_modes():
    w = GC.window(1e6)
    assert w["exists"] and 60 < w["width"] < 120
    M = (5 * w["M3rN_min"] / 1e6) ** (1 / 3)
    e0 = -GC.eps0_for_tipping(M, 1e6)
    assert abs(GC.weighted_deficit(e0, M, 1e6) / (GC.WTH_OVER_RN * 1e6) - 1) < 1e-9
    g = GC.mode_growth(M, e0, 1e6, t_end_over_rN=9.6, ks=np.geomspace(1e-4, 3e-4, 6), ls=np.array([0, 1, 5, 20, 80]))
    assert g["max_amplification"] < 4.0                                        # ~1 e-fold over three response times
    e_bad = 100 * e0                                                            # a hundred times deeper NEC phase: unstable
    gb = GC.mode_growth(M, e_bad, 1e6, t_end_over_rN=9.6, ks=np.geomspace(5e-5, 5e-4, 6), ls=np.array([0]))
    assert gb["max_amplification"] > 100


def test_strong_qnec_bound_and_quantum_threshold():
    b = QT.strong_qnec_weighted_bound()
    assert b["matches"]
    t = QT.THRESHOLD_VS_X
    cap = lambda x: x / (48 * np.pi * 3.2 * 17.3205) * 300
    gaps = [t[x] - cap(x) for x in (1.0, 2.0, 3.33, 5.0)]
    assert all(g > 0 for g in gaps) and all(np.diff(gaps) > 0)                 # capacity never catches the threshold


def test_regular_centre_initial_data_are_ordinary():
    pot, info = CO.potential()
    vphi, w = 0.0015, info["w"]
    A, rs, Rw = 0.04317836040165275, 16.643786422454006, 17.85264231414957
    r, s, p, c = CO.cone_data_alpha1(lambda x: vphi * np.tanh((x - Rw) / w), lambda x: A * np.exp(-((x - rs) / 0.08) ** 2), 0.02, 120.0)
    assert np.all(np.diff(r) > 0)                                               # r increases outward along the cone: no bridge
