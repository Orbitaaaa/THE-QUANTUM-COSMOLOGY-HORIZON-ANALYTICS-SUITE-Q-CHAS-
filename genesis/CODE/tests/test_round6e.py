"""Round 6e: the ghost condensate as a dynamical source, its self-generated NEC phase, the generalised entropy of eternal and transient
daughters, and the creation of the false-vacuum ball by supercooling."""
import warnings

import numpy as np

from babyuniverse import condensate as C, throat_evolution as TE, ghost_condensate as GC, gsl_solution as GS, fireball as FB

warnings.simplefilter("ignore")


def test_condensate_solver_reproduces_and_converges():
    M, r0 = 1.0, 6.0
    f0 = 1 - 2 * M / r0
    z = lambda x: 0 * np.asarray(x, float)
    A, B, _ = C.condensate_data(0.0, 0.0, r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, 0.1, 60, 60, z, z, frame="corner")
    a, b, _ = TE.characteristic_data(0.0, 0.0, r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, 0.1, 60, 60)
    assert np.nanmax(np.abs(C.evolve(A, B, 0.0, 0.0, TE.Potential(), 0.1)["r"] - TE.evolve(a, b, 0.0, 0.0, TE.Potential(), 0.1)["r"])) == 0.0
    lump = lambda x: 2e-3 * np.exp(-((np.asarray(x, float) - 2.0) / 0.8) ** 2)
    err = []
    for h in (0.1, 0.05):
        N = int(6 / h)
        A, B, _ = C.condensate_data(0.0, 0.0, r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, h, N, N, lump, lump, frame="corner")
        err.append(C.constraint_residual(C.evolve(A, B, 0.0, 0.0, TE.Potential(), h))["Evv"])
    assert 3.0 < err[0] / err[1] < 5.0                                         # second order with a back-reacting condensate


def test_clock_and_tilt_source():
    pot = TE.Potential(0.0, 1.0, 0.0)
    z = lambda x: 0 * np.asarray(x, float)
    m1 = lambda x: -1 + 0 * np.asarray(x, float)
    A, B, _ = C.condensate_data(0.0, 0.0, 5.0, -0.5, 0.5, 0.1, 40, 40, z, z, phi_v=m1, phi_u=m1, pot=pot, frame="corner")
    res = C.evolve(A, B, 0.0, 0.0, pot, 0.1, lam_f=1e-6)
    n, T = C.stress(res)["n"], res["T"]
    assert abs(n[-1, -1] / (-1e-6 * T[-1, -1]) - 1) < 1e-3                      # div(n u) = -lam_f in the false vacuum
    res0 = C.evolve(A, B, 0.0, 0.0, pot, 0.1)
    assert np.nanmax(np.abs(res0["D"])) == 0.0                                  # Theorem 6e-1: no charge in, no charge out


def test_condensate_tips_the_throat_by_itself():
    o = C.nariai_condensate(-1e-5, frame="corner", h=0.08)
    assert o["grew"] and o["true_vacuum_untrapped_large_r"] == 0               # its own NEC phase, no prescribed flux; hidden
    assert not C.nariai_condensate(-4e-6, frame="corner", h=0.08)["grew"]
    seeded = C.nariai_condensate(-1e-5, frame="corner", h=0.08, dphi0=1.5e-9)   # 1e-6 v off the false vacuum: the implicit potential keeps
    assert seeded["grew"] and seeded["n_big"] == o["n_big"]                    # the stiff scalar stable where m h e^s / 2 > 1 (F42)
    dh, eh = GS.daughter_horizons(o, n_points=6), GS.parent_event_horizon(o)
    assert eh["deficit"] == 0.0 and dh["max_deficit"] > 0                      # eternal daughter: local GSL deficit only on daughter horizons
    assert C.hiding_report(o)["hidden"]                                          # no daughter row reaches the outer edge


def test_self_generated_nec_phase():
    o = C.nariai_condensate(0.0, frame="corner", h=0.08, lam_f=4e-7)
    assert o["grew"] and o["true_vacuum_untrapped_large_r"] == 0               # zero charge on the data: the tilt makes it
    assert not C.nariai_condensate(0.0, frame="corner", h=0.08)["grew"]
    good = C.nariai_condensate(0.0, frame="corner", h=0.08, lam_f=4e-7, Lu=40.0, Lv=40.0)
    bad = C.nariai_condensate(0.0, frame="corner", h=0.08, lam_f=4e-7, lam_s=0.0, Lu=40.0, Lv=40.0)
    g = C.constraint_residual(good["res"], r_lo=0.5 * good["r_N"], r_hi=3 * good["r_N"], s_max=5.0)
    b = C.constraint_residual(bad["res"], r_lo=0.5 * bad["r_N"], r_hi=3 * bad["r_N"], s_max=5.0)
    assert g["Evv_rel"] < 1e-5 and b["Evv_rel"] > 0.1                          # only the conserving combination satisfies the constraints


def test_stability_with_a_measured_history():
    rN = 1e6
    x = np.linspace(0, 3.6, 50)
    nu = -1.8e-3 * np.ones_like(x)
    M = (2e-3 / rN) ** (1 / 3)
    eps_t = lambda t: np.interp(t / rN, x, nu, right=nu[-1] * np.exp(-(t / rN - 3.6))) / rN**2 / (2 * M**4)
    g = GC.mode_growth(M, 0.0, rN, t_end_over_rN=9.6, eps_t=eps_t, ks=np.geomspace(1e-2 / rN, 3e-4, 12), ls=np.array([0, 2, 10, 40]))
    assert g["max_amplification"] < 4.0


def test_transient_exit_of_the_false_vacuum():
    vphi, eps = 0.0015, 4.0
    rho_f = 3 / 900 / (8 * np.pi)
    pot = TE.Potential(4 * rho_f / (eps * vphi**4), vphi, rho_f)
    z = lambda x: 0 * np.asarray(x, float)
    fv = lambda x: -vphi + 1.5e-9 + 0 * np.asarray(x, float)
    A, B, _ = C.condensate_data(0.0, 0.0, 20.0, -0.5, 0.5, 0.02, 600, 600, z, z, phi_v=fv, phi_u=fv, pot=pot, frame="corner")
    lam = 1e-5                                                                  # T_exit = rho_f (16/(3 eps) - 1) / lam = 4.4
    r1 = C.evolve(A, B, 0.0, 0.0, pot, 0.02, lam_f=lam, clip=1, phi_max=10 * vphi)
    assert np.nanmax(r1["phi"]) > 0.5 * vphi                                   # the tilt destabilises the false vacuum: graceful exit
    r0 = C.evolve(A, B, 0.0, 0.0, pot, 0.02, clip=1, phi_max=10 * vphi)
    assert np.nanmax(np.abs(r0["phi"] + vphi)) < 1e-6                           # without it the false vacuum stays


def test_supercooling_creates_the_ball():
    fb = FB.Fireball(eps=4.0, alpha_pt=3.0)
    s = FB.static_hot_ball(fb, T_over_Tc=1.4)
    assert s["normal"] and s["x_end"] < -0.999 and s["kinetic_end_over_rho_f"] < 1e-6 and s["adiabaticity"] < 0.01
    weak = FB.static_hot_ball(FB.Fireball(eps=4.0, alpha_pt=1.5), T_over_Tc=1.4)
    assert not weak["normal"]                                                   # weak transitions: the hot ball would be trapped
    fb1 = FB.Fireball(eps=4.0, alpha_pt=1.0)
    rho_i = fb1.thermo(fb1.v, 0.0, 1.8 * fb1.Tc)[0]
    ad = 0.95 * 0.8 / 0.6
    sol = fb1.run(1.8 * fb1.Tc, np.sqrt((1 + ad**2) / (8 * np.pi / 3 * rho_i)), t_max=400.0)
    sm = fb1.summary(sol, sin_chi_b=0.6)
    assert sm["turnaround"] and sm["normal_throughout"] and sm["x_end"] < -0.99 and sm["rad_over_rho_f"] > 1.0   # confined: Einstein-static bound
