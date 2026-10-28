"""Round 6f: thermal formation in one GR evolution (it fails), observers that end at singularities, the condensate beyond its IR fluid
with non-spherical modes and the full stability range, the solver cross-check and the continuous route, and problems 4 and 5."""
import warnings

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from babyuniverse import (formation as F, condensate as C, condensate_eft as E, throat_evolution as TE, gsl_solution as GS,
                          accessibility as AC)

warnings.simplefilter("ignore")
RHO_F = 3 / 900 / (8 * np.pi)
V = 0.0015
POT = TE.Potential(4 * RHO_F / (4 * V**4), V, RHO_F)


def _dU(p):
    return POT.lam * p * (p * p - V * V) - RHO_F * (0.75 - 0.75 * (p / V) ** 2) / V


def test_formation_solver_matches_the_homogeneous_ode():
    pl = F.Plasma(RHO_F, V, POT.lam)
    T0 = 1.1 * pl.spinodal() * pl.Tc
    E0 = pl.eos(T0, V)[0] + POT.U(V)
    Tof = lambda p, pd: brentq(lambda T: pl.eos(T, p)[0] - (E0 - 0.5 * pd * pd - POT.U(p)), 1e-6, 10)
    sol = solve_ivp(lambda t, y: [y[1], -pl.Gamma * y[1] - _dU(y[0]) - pl.eos(Tof(y[0], y[1]), y[0])[2]], (0, 1.0), [V, 0.0],
                    max_step=0.0005, rtol=1e-11, atol=1e-16, dense_output=True)
    ph = lambda t: sol.sol(np.asarray(t, float))[0]
    Tt = np.vectorize(lambda t: Tof(*sol.sol(t)))
    err = []
    for h in (0.02, 0.01):
        N = int(1.0 / h)
        A, B, _ = F.corner_data(60.0, -0.5, 0.5, h, N, N, pl, T_v=lambda x: Tt(np.asarray(x) / 2), T_u=lambda x: Tt(np.asarray(x) / 2),
                                phi_v=lambda x: ph(np.asarray(x) / 2), phi_u=lambda x: ph(np.asarray(x) / 2), G=1e-10)
        res = F.evolve_corner(A, B, h, POT, pl, G=1e-10, clip=0)
        k = np.arange(N + 1)
        err.append(np.nanmax(np.abs(res["phi"][k, k] - ph(k * h))) / V)
    assert err[1] < 2e-3 and err[0] / err[1] > 3.0                                     # second order: field, plasma, friction together


def test_regular_centre_keeps_the_vacuum():
    pl = F.Plasma(RHO_F, V, POT.lam)
    h = 0.02
    Z0 = F.centre_cone_data(h, 6.0, pl, phi_r=lambda r: V + 0 * r, T_r=lambda r: pl.T_floor + 0 * r, G=1e-10)
    res = F.evolve_centre(Z0, h, POT, pl, G=1e-10, clip=0)
    n = int(np.sum(np.isfinite(res["r"]).any(axis=1)))
    assert n > 250 and max(abs(res["s"][i, i]) for i in range(1, n)) < 1e-8          # the axis clock fix (Round 6f, F44)


def test_heating_needs_more_than_the_false_vacuum_energy():
    pl = F.Plasma(RHO_F, V, POT.lam, f_H=0.999, x_c=30.0, Gamma_over_m=0.3)
    assert not F.restores_homogeneous(1.5 * RHO_F, pl, POT)[0]                       # stall equilibrium: restoration needs E >~ 2 rho_f
    assert F.restores_homogeneous(2.5 * RHO_F, pl, POT)[0]


def test_observers_ending_at_a_singularity_lose_the_whole_area():
    o = C.route_condensate(lam_f=0.0)                                                  # plain black hole: no daughter at all
    R = o["res"]["r"]
    d = GS.ending_observer_deficits(o, np.isfinite(R) & (R < 0.1 * o["r_N"]), n_points=6)
    assert d["max_deficit"] > 0.5 * np.pi * o["r_N"] ** 2                              # Theorem 6f-2: no GSL for incomplete observers


def test_eft_nonlinear_equals_linear_for_small_seeds():
    mu, al, e0 = 120.0, 1.0, -0.08
    qb = np.sqrt(2 * abs(e0) * mu**2 / al)
    L = 30 * 2 * np.pi / qb
    nl = E.evolve_1p1(mu, al, e0, L=L, N=256, t_end=6.0, seed_amp=1e-7, k_seed=2 * qb)
    li = E.linear_rms_1p1(mu, al, e0, L=L, N=256, t_end=6.0, seed_amp=1e-7, k_seed=2 * qb)
    assert not nl["broke"] and abs(nl["amplification"] / li["amplification"] - 1) < 0.01


def test_growth_exponent_and_universal_window():
    for (mu, e0, sg) in ((30.0, -0.2, 0.0), (30.0, -0.05, 0.05)):
        g = E.growth_exponent(mu, 1.0, e0, sigma=sg, t_end=9.6 if sg == 0 else 20.0)
        qb = np.sqrt(2 * max(abs(e0), sg) * mu**2)
        a = E.linear_rms_1p1(mu, 1.0, e0, sigma=sg, L=30 * 2 * np.pi / qb, N=256, t_end=9.6 if sg == 0 else 20.0)["amplification"]
        assert 0.4 < np.exp(g["G"]) / a < 2.5                                         # exponential growth only, tracks the full ODE
    x = np.linspace(0, 4.5, 60)
    nu, rr = -1.8e-3 * np.exp(-x), np.ones_like(x)
    g6 = [E.growth_exponent_history(m, 1e6, 1.0, x, nu, rr)["G"] for m in (3e-4, 3e-3, 3e-2)]
    g20 = [E.growth_exponent_history(m, 1e20, 1.0, x, nu, rr)["G"] for m in (3e-4, 3e-3, 3e-2)]
    assert np.allclose(g6, g20, rtol=1e-2)                                             # the window depends on M^3 r_N only (to the k-grid sampling)
    assert g6[0] > g6[1] < g6[2]                                                       # gradient below, Jeans above


def test_nonspherical_modes_in_the_regular_region():
    o = C.nariai_condensate(-1e-5, frame="corner", h=0.08)
    rN = o["r_N"]
    sc = E.lmode_amplification(o["res"], ells=(0, 2, 10), pot=o["pot"], r_cut=0.5 * rN)
    rw = E.lmode_amplification(o["res"], ells=(2, 5, 20), kappa=1.0, pot=o["pot"], r_cut=0.5 * rN, power=3)
    assert max(sc.values()) < 1.5 and max(rw.values()) < 1.5                         # no growth outside the crushing region


def test_two_solvers_agree_on_the_route():
    o = C.route_condensate(lam_f=1e-6, Lv=300.0)
    ra = o["res"]
    rN, r0, w_shell = o["r_N"], 15.0, 0.25
    rv0 = np.sqrt(1 - 3 / 900 * r0**2 / 3) / 2
    w = np.sqrt(2) / (np.sqrt(POT.lam) * V)
    shell = lambda v: 7.258540773391722 / (np.sqrt(np.pi) * w_shell) * np.exp(-((np.asarray(v, float) - 4.5) / w_shell) ** 2)
    pl = F.Plasma(RHO_F, V, POT.lam, g_star=1e-12)
    A, B, _ = F.corner_data(r0, -rv0, rv0, 0.04, 7500, 750, pl, phi_v=lambda v: V * np.tanh((np.asarray(v, float) - o["v_wall"]) / w),
                            phi_u=lambda u: -V + 0 * np.asarray(u, float), extra_Nv=shell)
    rb = F.evolve_corner(A, B, 0.04, POT, pl, lam_f=1e-6, clip=0, r_inf=1e4)
    both = np.isfinite(ra["r"]) & np.isfinite(rb["r"]) & (ra["r"] > 0.5 * rN) & (ra["r"] < 3 * rN)
    assert o["grew"] and TE.inflated(rb, TE.classify(rb), rN)["grew"]
    assert np.median(np.abs(rb["r"][both] / ra["r"][both] - 1)) < 1e-4


def test_exit_must_come_after_the_throat_forms():
    late = C.route_condensate(lam_f=1e-6, eps=4.0, clip=1, dphi0=1.5e-9, r_inf=1e40)   # spinodal at T = 44, before the throat forms
    assert not late["grew"]
    ok = C.route_condensate(lam_f=1e-6, eps=2.0, E=7.233019253534085, clip=1, dphi0=1.5e-9, r_inf=1e40)   # spinodal at T = 221
    assert ok["normal_on_data"] and ok["grew"] and C.hiding_report(ok)["hidden"]


def test_accessibility_invariant_and_gaps():
    inv = AC.requirement_invariant()
    for rN in (1e5, 1e30, 1e60):
        s = AC.physical_scales(rN)
        assert abs(s["E_N_J"] * np.sqrt(s["rho_f_J_m3"]) / inv - 1) < 1e-10        # Theorem 6f-7: independent of the vacuum scale
    assert AC.gap_report()["best_technology_gap_log10"] > 38
    assert AC.r_N_min_from_condensate_bound(100.0)["M_N_Msun"] > 1e8                   # an observationally allowed condensate: >= 1e8 M_sun
    assert AC.condensate_needed_for_mass(1e6) > 1e12                                    # a laboratory-mass throat needs M >~ 1e13 GeV
