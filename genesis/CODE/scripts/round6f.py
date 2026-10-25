"""Round 6f: (1) the thermal formation of the false-vacuum ball in one continuous GR evolution; (2) the generalised entropy of observers
whose histories end at singularities; (3) the condensate beyond its IR fluid, non-spherical perturbations and the full stability range;
(4) the cross-check of the solvers and the continuous route from the ball to the hidden daughter; problems 4 and 5 (does nature support
the route; is it technologically accessible).  Writes DATA/round6f.json and FIGURES/round6f.png; `--figure-only` redraws the figure."""
import copy
import gc
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from _common import FIG, DATA, dump
from babyuniverse import (condensate as C, throat_evolution as TE, formation as F, condensate_eft as E, gsl_solution as GS,
                          accessibility as AC)

warnings.simplefilter("ignore")
RHO_F = 3 / 900 / (8 * np.pi)
VPHI = 0.0015
RN = np.sqrt(300.0)


def pot_for(eps):
    return TE.Potential(4 * RHO_F / (eps * VPHI**4), VPHI, RHO_F)


# ------------------------------------------------------------------ item 1
def threshold_job(a):
    eps, Gm, fH, xc = a
    pot = pot_for(eps)
    pl = F.Plasma(RHO_F, VPHI, pot.lam, f_H=fH, x_c=xc, Gamma_over_m=Gm)
    return dict(eps=eps, Gamma_over_m=Gm, f_H=fH, x_c=xc, E_min_over_rho_f=F.restoration_threshold(pl, pot))


def centre_job(a):
    """A hot ball at fac x E_min (largest cone-normal radius), evolved from a regular centre with gravity."""
    eps, Gm, fH, xc, Emin, fac = a
    pot = pot_for(eps)
    pl = F.Plasma(RHO_F, VPHI, pot.lam, f_H=fH, x_c=xc, Gamma_over_m=Gm)
    Ef = fac * Emin
    Th = brentq(lambda T: pl.eos(T, VPHI)[0] - Ef * RHO_F, pl.Tc * 0.3, pl.Tc * 6)
    R = F.largest_normal_ball(pl, Th, np.arange(15.0, 22.0, 0.25))
    out = dict(eps=eps, Gamma_over_m=Gm, E_over_rho_f=Ef, R_normal_max=R, R_target=0.595 * 30)
    if R is None:
        return out
    h = 0.04
    ball = lambda r: pl.T_floor + (Th - pl.T_floor) * 0.5 * (1 - np.tanh((np.asarray(r, float) - R) / 0.5))
    Z0 = F.centre_cone_data(h, 140.0, pl, phi_r=lambda r: VPHI + 0 * np.asarray(r, float), T_r=ball)
    res = F.evolve_centre(Z0, h, pot, pl, n_rays=int(70 / h))
    out["history"] = F.centre_history(res, (5, 10, 15, 20, 30, 40, 50, 60, 70), VPHI)
    out["restoration_time_homogeneous"] = F.restores_homogeneous(Ef * RHO_F, pl, pot)[1]
    return out


def verification_item1():
    """PDE vs the homogeneous ODE (G -> 0), the vacuum centre stage, the centre-to-corner hand-off, constraint convergence."""
    v = {}
    pot = pot_for(4.0)
    U, dU = pot.U, (lambda p: pot.lam * p * (p * p - VPHI**2) - RHO_F * (0.75 - 0.75 * (p / VPHI) ** 2) / VPHI)
    pl = F.Plasma(RHO_F, VPHI, pot.lam)
    T0 = 1.1 * pl.spinodal() * pl.Tc
    ode = {}
    for h in (0.02, 0.01, 0.005):
        E0 = pl.eos(T0, VPHI)[0] + U(VPHI)
        Tof = lambda p, pd: brentq(lambda T: pl.eos(T, p)[0] - (E0 - 0.5 * pd * pd - U(p)), 1e-6, 10)
        sol = solve_ivp(lambda t, y: [y[1], -pl.Gamma * y[1] - dU(y[0]) - pl.eos(Tof(y[0], y[1]), y[0])[2]], (0, 2.0), [VPHI, 0.0],
                        max_step=0.0005, rtol=1e-11, atol=1e-16, dense_output=True)
        ph = lambda t: sol.sol(np.asarray(t, float))[0]
        Tt = np.vectorize(lambda t: Tof(*sol.sol(t)))
        N = int(2.0 / h)
        A, B, _ = F.corner_data(60.0, -0.5, 0.5, h, N, N, pl, T_v=lambda x: Tt(np.asarray(x) / 2), T_u=lambda x: Tt(np.asarray(x) / 2),
                                phi_v=lambda x: ph(np.asarray(x) / 2), phi_u=lambda x: ph(np.asarray(x) / 2), G=1e-10)
        res = F.evolve_corner(A, B, h, pot, pl, G=1e-10, clip=0)
        k = np.arange(N + 1)
        ode[f"h={h}"] = dict(phi_err_over_v=float(np.nanmax(np.abs(res["phi"][k, k] - ph(k * h))) / VPHI),
                             T_err_rel=float(np.nanmax(np.abs(np.exp(res["lnTf"][k, k]) / Tt(k * h) - 1))))
    v["pde_vs_homogeneous_ode"] = ode
    h = 0.02
    Z0 = F.centre_cone_data(h, 8.0, pl, phi_r=lambda r: VPHI + 0 * r, T_r=lambda r: pl.T_floor + 0 * r, G=1e-10)
    res = F.evolve_centre(Z0, h, pot, pl, G=1e-10, clip=0)
    n = int(np.sum(np.isfinite(res["r"]).any(axis=1)))
    v["vacuum_centre"] = dict(rows=n, max_abs_s_axis=float(max(abs(res["s"][i, i]) for i in range(1, min(n, 390)))))
    plc = F.Plasma(RHO_F, VPHI, pot.lam, f_H=0.98, x_c=8.0)
    Th = 1.02 * plc.spinodal() * plc.Tc
    ball = lambda r: plc.T_floor + (Th - plc.T_floor) * 0.5 * (1 - np.tanh((np.asarray(r, float) - 19.0) / 1.0))
    ov = []
    for h in (0.04, 0.02):
        Z0 = F.centre_cone_data(h, 80.0, plc, phi_r=lambda r: VPHI + 0 * np.asarray(r, float), T_r=ball)
        resA = F.evolve_centre(Z0, h, pot, plc, n_rays=int(30 / h), lam_f=1e-6)
        i1, j1 = int(2.0 / h), int(32.0 / h)
        A, B = F.handoff(resA, i1, j1, Lu=26.0, Lv=40.0)
        resB = F.evolve_corner(A, B, h, pot, plc, lam_f=1e-6)
        ov.append(dict(h=h, **F.overlap_check(resA, resB, i1, j1)))
        del resA, resB
        gc.collect()
    v["handoff_overlap"] = ov
    return v


# ------------------------------------------------------------------ item 2
def gsl_job(tag):
    if tag == "plain_black_hole":
        o = C.route_condensate(lam_f=0.0)
    elif tag == "transient_daughter":
        o = C.nariai_condensate(0.0, frame="corner", h=0.08, lam_f=2e-7, Lu=520.0, Lv=520.0, clip=1, dphi0=1.5e-9)
    else:
        o = C.nariai_condensate(0.0, frame="corner", h=0.08)
    R = o["res"]["r"]
    d = GS.ending_observer_deficits(o, np.isfinite(R) & (R < 0.1 * o["r_N"]), n_points=10)
    return dict(case=tag, points=d.get("points"), max_deficit=d.get("max_deficit"), S_N=float(np.pi * o["r_N"] ** 2))


# ------------------------------------------------------------------ item 3
def history(tag):
    n0, lf = (-6e-6, 0.0) if tag == "initial_charge" else (0.0, 4e-7)
    o = C.nariai_condensate(n0, frame="corner", h=0.04, lam_f=lf)
    res, rN = o["res"], o["r_N"]
    n = C.stress(res)["n"]
    k = np.arange(min(n.shape))
    T, nd, rd = res["T"][k, k], n[k, k], res["r"][k, k]
    ok = np.isfinite(T) & np.isfinite(nd)
    return tag, T[ok] / rN, nd[ok] * rN**2, rd[ok] / rN


GRID = np.geomspace(1e-5, 1.0, 61)


def window_job(a):
    tag, x, nu, rr, rN, al = a
    G = [E.growth_exponent_history(m, rN, al, x, nu, rr, decay_after=(tag == "initial_charge"))["G"] for m in GRID]
    ok = [float(m) for m, g in zip(GRID, G) if g < 1]
    return dict(history=tag, r_N=rN, alpha=al, lo=min(ok) if ok else None, hi=max(ok) if ok else None, G=[float(g) for g in G])


def lmode_job(tag):
    if tag == "band_self_tipping":
        o = C.nariai_condensate(-1e-5, frame="corner", h=0.08)
    elif tag == "band_control_no_tip":
        o = C.nariai_condensate(-4e-6, frame="corner", h=0.08)
    elif tag == "route_self_made_push":
        o = C.route_condensate(lam_f=1e-6)
    elif tag == "route_control_no_push":
        o = C.route_condensate(lam_f=0.0)
    else:
        o = C.nariai_condensate(0.0, frame="corner", h=0.08, lam_f=2e-7, Lu=400.0, Lv=400.0, clip=1, dphi0=1.5e-9)
    rN = o["r_N"]
    return dict(case=tag, grew=bool(o["grew"]),
                scalar_all=E.lmode_amplification(o["res"], ells=(0, 1, 2, 3, 5, 10, 20), pot=o["pot"]),
                scalar_regular=E.lmode_amplification(o["res"], ells=(0, 1, 2, 3, 5, 10, 20), pot=o["pot"], r_cut=0.5 * rN),
                rw_weyl_regular=E.lmode_amplification(o["res"], ells=(2, 3, 5, 10, 20), kappa=1.0, pot=o["pot"], r_cut=0.5 * rN, power=3))


def eft_job(a):
    M3rN, mode, seed_amp = a
    rN, al, K = 1e6, 1.0, 7.5e-4
    mu = (M3rN * rN**2) ** (1 / 3)
    eps0 = -K * rN**2 / mu**4
    J = 8 * np.pi * al * (mu / rN) ** 2
    sig = abs(eps0) if mode == "self" else 0.0
    qb = np.sqrt(2 * (abs(eps0) + J) * mu**2 / al)
    L = 30 * 2 * np.pi / qb
    nl = E.evolve_1p1(mu, al, eps0, J=J, sigma=sig, L=L, N=512, t_end=9.6, seed_amp=seed_amp, k_seed=2 * qb)
    li = E.linear_rms_1p1(mu, al, eps0, J=J, sigma=sig, L=L, N=512, t_end=9.6, seed_amp=seed_amp, k_seed=2 * qb)
    return dict(M3rN=M3rN, mode=mode, seed=seed_amp, mu=float(mu), eps0=float(eps0), dX_over_eps=float((seed_amp * qb) ** 2 / abs(eps0)),
                nonlinear=nl["amplification"], linear=li["amplification"], caustic=bool(nl["broke"]))


# ------------------------------------------------------------------ item 4
def xcheck_job(lam_f):
    """The route evolved by the formation solver (negligible plasma) and by condensate.evolve, same data."""
    o = C.route_condensate(lam_f=lam_f)
    ra = o["res"]
    # rebuild the same data for the formation solver
    E_, l_H, r0, eps, w_shell, v_shell = 7.258540773391722, 30.0, 15.0, 4.0, 0.25, 4.5
    Lam = 3 / l_H**2
    lam = 4 * RHO_F / (eps * VPHI**4)
    w = np.sqrt(2) / (np.sqrt(lam) * VPHI)
    rv0 = np.sqrt(1 - Lam * r0**2 / 3) / 2
    shell = lambda v: E_ / (np.sqrt(np.pi) * w_shell) * np.exp(-((np.asarray(v, float) - v_shell) / w_shell) ** 2)
    h, Nv, Nu = 0.04, int(600 / 0.04), int(30 / 0.04)
    vw = o["v_wall"]
    pl = F.Plasma(RHO_F, VPHI, lam, g_star=1e-12)
    A2, B2, _ = F.corner_data(r0, -rv0, rv0, h, Nv, Nu, pl, phi_v=lambda v: VPHI * np.tanh((np.asarray(v, float) - vw) / w),
                              phi_u=lambda u: -VPHI + 0 * np.asarray(u, float), extra_Nv=shell)
    rb = F.evolve_corner(A2, B2, h, o["pot"], pl, lam_f=lam_f, clip=0, r_inf=1e4)
    gb = TE.inflated(rb, TE.classify(rb), RN)
    both = np.isfinite(ra["r"]) & np.isfinite(rb["r"]) & (ra["r"] > 0.5 * RN) & (ra["r"] < 3 * RN)
    return dict(lam_f=lam_f, grew=(bool(o["grew"]), bool(gb["grew"])), n_big=(int(o["n_big"]), int(gb["n_big"])),
                r_rel_median=float(np.median(np.abs(rb["r"][both] / ra["r"][both] - 1))),
                r_rel_p95=float(np.percentile(np.abs(rb["r"][both] / ra["r"][both] - 1), 95)),
                phi_p95_over_v=float(np.percentile(np.abs(rb["phi"][both] - ra["phi"][both]), 95) / VPHI),
                D_p95_rel=float(np.percentile(np.abs(rb["D"][both] - ra["D"][both]), 95) / max(np.nanmax(np.abs(ra["D"][both])), 1e-300)))


E_SHELL = {4.0: 7.258540773391722, 2.0: 7.233019253534085}        # shell energy that puts the data ray on the marginal surface


def continuous_job(a):
    """One evolution from the route's ball (zero charge) through the self-generated NEC phase, tipping and the hidden daughter."""
    eps, lam_f, dphi0 = a
    o = C.route_condensate(lam_f=lam_f, eps=eps, E=E_SHELL[eps], clip=1, dphi0=dphi0, Lv=600.0, r_inf=1e40)
    res = o["res"]
    R, P, T, S = res["r"], res["phi"], res["T"], res["s"]
    ok = np.isfinite(R)
    dau = ok & (R > 3 * RN) & (np.abs(P / VPHI + 1) < 0.5)
    T_exit = RHO_F * (16 / (3 * eps) - 1) / lam_f if lam_f > 0 else np.inf
    # the stiffness of the tachyonic false vacuum after the spinodal: h^2 e^{2s} |U''| / 16 at the last daughter cells
    lam = 4 * RHO_F / (eps * VPHI**4)
    late = dau & (T > T_exit)
    m2 = 1.5 * lam_f * np.maximum(T - T_exit, 0) / VPHI**2 if lam_f > 0 else 0 * T
    stiff = (res["h"] ** 2 * np.exp(2 * S) * m2 / 16)[late] if late.any() else np.array([0.0])
    hid = C.hiding_report(o)
    return dict(eps=eps, E_shell=E_SHELL[eps], lam_f=lam_f, dphi0=dphi0, grew=bool(o["grew"]), n_big=int(o["n_big"]), T_exit_spinodal=float(T_exit),
                T_daughter_p99=float(np.percentile(T[dau], 99)) if dau.any() else None,
                daughter_cells=int(dau.sum()),
                daughter_cells_after_spinodal=int(late.sum()), stiffness_median_after_spinodal=float(np.median(stiff)),
                normal_on_data=bool(o["normal_on_data"]), hidden=bool(hid["hidden"]))


# ------------------------------------------------------------------ run
if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round6f.json"))
else:
    out = {}
    i1 = {}
    jobs = [(eps, 0.3, 0.999, 30.0) for eps in (1.0, 2.0, 3.0, 4.0, 5.0)] + [(4.0, Gm, 0.999, 30.0) for Gm in (0.03, 1.0, 3.0)] + \
           [(5.2, 0.03, 0.999, 30.0), (4.0, 0.3, 0.95, 5.0)]
    with ProcessPoolExecutor(4) as ex:
        i1["restoration_thresholds"] = list(ex.map(threshold_job, jobs))
        th = {(r["eps"], r["Gamma_over_m"], r["f_H"], r["x_c"]): r["E_min_over_rho_f"] for r in i1["restoration_thresholds"]}
        cj = [(4.0, 0.3, 0.999, 30.0, th[(4.0, 0.3, 0.999, 30.0)], f) for f in (1.03, 1.15)] + \
             [(5.2, 0.03, 0.999, 30.0, th[(5.2, 0.03, 0.999, 30.0)], f) for f in (1.03, 1.1)]
        i1["centre_runs"] = list(ex.map(centre_job, cj))
    i1["route_target"] = dict(r_wall_over_lH=0.595, ADM_mass=8.93, M_N=float(RN / 3), note="Round 6e, R117: the ball the route starts from")
    i1["verification"] = verification_item1()
    out["item1"] = i1
    gc.collect()

    with ProcessPoolExecutor(3) as ex:
        out["item2"] = dict(ending_observers=list(ex.map(gsl_job, ["plain_black_hole", "near_nariai_collapse", "transient_daughter"])))
    gc.collect()

    i3 = {}
    with ProcessPoolExecutor(2) as ex:
        H = {t: (x, nu, rr) for t, x, nu, rr in ex.map(history, ["initial_charge", "self_generated"])}
    i3["histories"] = {t: dict(x_max=float(x.max()), int_abs_nu=float(np.trapezoid(np.abs(nu), x))) for t, (x, nu, rr) in H.items()}
    with ProcessPoolExecutor(4) as ex:
        i3["windows"] = list(ex.map(window_job, [(t, *H[t], rN, al) for t in H for rN in (1e6, 1e20, 1e40) for al in (0.1, 1.0, 10.0)]))
        i3["eft_nonlinear"] = list(ex.map(eft_job, [(m, mode, s) for m in (3.8e-4, 3.2e-2) for mode in ("dilute", "self")
                                                     for s in (4e-8, 1e-5, 1e-4, 1e-3)]))
    i3["grid"] = GRID.tolist()
    with ProcessPoolExecutor(3) as ex:
        i3["lmodes"] = list(ex.map(lmode_job, ["band_control_no_tip", "band_self_tipping", "route_control_no_push", "route_self_made_push",
                                               "transient_daughter"]))
    out["item3"] = i3
    gc.collect()

    i4 = dict(solver_crosscheck=xcheck_job(1e-6))
    gc.collect()
    with ProcessPoolExecutor(2) as ex:
        i4["continuous"] = list(ex.map(continuous_job, [(4.0, 1e-6, 1.5e-9), (2.0, 1e-6, 1.5e-9), (4.0, 1e-6, 0.0), (2.0, 0.0, 1.5e-9)]))
    out["item4"] = i4
    gc.collect()

    w1 = [w for w in i3["windows"] if w["alpha"] == 1.0 and w["r_N"] == 1e6]
    win = {w["history"]: (w["lo"], w["hi"]) for w in w1}
    rep = {}
    for tag, wv in win.items():
        g = AC.gap_report(window=wv)
        rep[tag] = dict(window=wv, r_N_min=g["r_N_min"], condensate_for_lab=g["condensate_for_lab"], lightest_throat=g["lightest_throat"],
                        condensate_at_lightest=AC.condensate_mass_range(g["lightest_throat"]["r_N"], wv))
    g = AC.gap_report(window=win["initial_charge"])
    out["problems45"] = dict(invariant=g["invariant"], capabilities=g["capabilities"], best_technology_gap_log10=g["best_technology_gap_log10"],
                             by_history=rep, jeans_bound_code_conventions_GeV=AC.jeans_bound_GeV(),
                             literature_bounds_GeV=dict(nonlinear=100.0, linear=1e-2))
    dump("round6f.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(1, 4, figsize=(19, 4.6))
i1, i3, i4, p45 = out["item1"], out["item3"], out["item4"], out["problems45"]
# (a) the centre-stage histories: false-vacuum radius against time, trapped points marked
for k, run in enumerate(i1["centre_runs"]):
    if "history" not in run:
        continue
    hs = [x for x in run["history"] if x["computed"]]
    col = [C1, C2, C3, C4][k % 4]
    ax[0].plot([x["t"] for x in hs], [x["r_false_vacuum"] for x in hs], "o-", color=col, ms=5,
               label=f"eps={run['eps']:g}, E={run['E_over_rho_f']:.2f} rho_f, R={run['R_normal_max']:g}")
    tr = [x for x in hs if x["trapped"]]
    ax[0].plot([x["t"] for x in tr], [x["r_false_vacuum"] for x in tr], "x", color=INK, ms=7, mew=1.5)
ax[0].axhline(i1["route_target"]["r_wall_over_lH"] * 30, color=INK, lw=0.8, ls="--")
ax[0].text(1, i1["route_target"]["r_wall_over_lH"] * 30 + 0.4, "ball the route needs (r = 17.9)", fontsize=7, color=INK)
ax[0].set_xlabel("t = (u + v) / 2"); ax[0].set_ylabel("radius still in the false vacuum")
ax[0].set_title("(a) heating: the false vacuum shrinks or is trapped (x = trapped)", fontsize=9)
ax[0].legend(frameon=False, fontsize=6.5, loc="upper right")
# (b) stability: growth exponent against M^3 r_N, alpha = 1, both histories
g = np.array(i3["grid"])
for w, col, lab in ((next(x for x in i3["windows"] if x["history"] == "initial_charge" and x["alpha"] == 1.0 and x["r_N"] == 1e6), C1, "initial charge"),
                    (next(x for x in i3["windows"] if x["history"] == "self_generated" and x["alpha"] == 1.0 and x["r_N"] == 1e6), C2, "self-generated")):
    ax[1].plot(g, np.maximum(w["G"], 1e-3), color=col, label=f"{lab}: window {w['lo']:.1e} to {w['hi']:.1e}")
ax[1].axhline(1, color=INK, lw=0.8, ls="--")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_ylim(1e-3, 1e3)
ax[1].set_xlabel("M^3 r_N (Planck units; identical for r_N = 1e6, 1e20, 1e40)"); ax[1].set_ylabel("largest growth exponent (k and l modes)")
ax[1].set_title("(b) full linear stability range (alpha = 1; edges scale as alpha^-1/2)", fontsize=9)
ax[1].legend(frameon=False, fontsize=7)
# (c) non-spherical perturbations in the regular region (r > r_N / 2): tipping backgrounds against their controls
cases = i3["lmodes"]
ells = ["2", "3", "5", "10", "20"]
xs = np.arange(len(ells))
for k, c in enumerate(cases):
    col = [MUTED, C1, "#a9a8a2", C2, C3][k]
    ax[2].plot(xs, [c["rw_weyl_regular"][l] for l in ells], "o-", color=col, ms=5, lw=1.5, label=c["case"].replace("_", " "))
ax[2].set_xticks(xs); ax[2].set_xticklabels(ells)
ax[2].set_xlabel("l (Regge-Wheeler proxy, curvature amplitude)"); ax[2].set_ylabel("largest amplification for r > r_N / 2")
ax[2].set_ylim(0.8, 5)
ax[2].set_title("(c) non-spherical modes: tipping adds no growth over the controls", fontsize=9)
ax[2].legend(frameon=False, fontsize=6.5)
# (d) problems 4 and 5: energy against energy density
caps = p45["capabilities"]
Eg = np.geomspace(1e-6, 1e70, 50)
inv = p45["invariant"]
ax[3].plot(Eg, (inv / Eg) ** 2, color=INK, lw=1.2, ls="--", label="any near-Nariai throat (Theorem 6f-7)")
for c in caps:
    col = C1 if c["kind"] == "technology" else C3
    ax[3].plot(c["E_J"], c["rho_J_m3"], "o" if c["kind"] == "technology" else "s", color=col, ms=7)
    ax[3].annotate(c["name"].split(" (")[0], (c["E_J"], c["rho_J_m3"]), fontsize=6, color=INK, xytext=(4, 3), textcoords="offset points")
bh = p45["by_history"]["initial_charge"]
for key, mk in (("M<=100 GeV", "^"), ("M<=0.01 GeV", "v")):
    s = bh["r_N_min"][key]
    ax[3].plot(s["E_N_J"], s["rho_f_J_m3"], mk, color=C2, ms=8)
    ax[3].annotate(f"smallest throat, condensate {key.replace('M<=', 'M <= ')}", (s["E_N_J"], s["rho_f_J_m3"]), fontsize=6, color=INK,
                   xytext=(-150, -12), textcoords="offset points")
s = bh["lightest_throat"]
ax[3].plot(s["E_N_J"], s["rho_f_J_m3"], "D", color=C4, ms=7)
ax[3].annotate("throat at 1.5e16 GeV (needs M ~ 1e16 GeV)", (s["E_N_J"], s["rho_f_J_m3"]), fontsize=6, color=INK, xytext=(5, -2), textcoords="offset points")
ax[3].set_xscale("log"); ax[3].set_yscale("log")
ax[3].set_xlim(1e-6, 1e70); ax[3].set_ylim(1e-10, 1e110)
ax[3].set_xlabel("energy delivered (J)"); ax[3].set_ylabel("energy density reached (J / m^3)")
ax[3].set_title(f"(d) best technology is 10^{p45['best_technology_gap_log10']:.0f} short; compact-object astrophysics is not", fontsize=9)
ax[3].legend(frameon=False, fontsize=7, loc="lower left")
fig.tight_layout()
fig.savefig(f"{FIG}/round6f.png", dpi=150, facecolor="#fcfcfb")
print("done")
