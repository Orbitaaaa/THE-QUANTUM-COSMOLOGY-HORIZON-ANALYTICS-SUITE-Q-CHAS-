"""Round 6e (problems 1-3): the ghost condensate as a dynamical source (its own NEC phase, located, stable, also self-generated through a
false-vacuum tilt), the generalised entropy of the complete transition (eternal versus transient daughters), and the creation of the
false-vacuum ball from the true vacuum by supercooling."""
import gc
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse import (condensate as C, throat_evolution as TE, ghost_condensate as GC, gsl_solution as GS, preparation as PR,
                          fireball as FB)

warnings.simplefilter("ignore")
RHO_F = 3 / 900 / (8 * np.pi)


def bisect(fun, lo, hi, n):
    """fun(x) True = grows; geometric bisection on |x| in [lo, hi]."""
    for _ in range(n):
        mid = np.sqrt(lo * hi)
        if fun(mid):
            hi = mid
        else:
            lo = mid
    return [float(lo), float(hi)]


def thr_fill(kw):
    kw = dict(kw)
    lo, hi, n = kw.pop("lo"), kw.pop("hi"), kw.pop("n", 6)
    return dict(kw=kw, n_th=bisect(lambda x: C.nariai_condensate(-x, frame="corner", **kw)["grew"], lo, hi, n))


def thr_tilt(h):
    return dict(h=h, lam_th=bisect(lambda x: C.nariai_condensate(0.0, frame="corner", h=h, lam_f=x)["grew"], 1e-8, 4e-7, 7))


def thr_route(_):
    return bisect(lambda x: C.route_condensate(lam_f=x)["grew"], 3e-7, 1e-6, 4)


def tilt_diag(a):
    tag, lf, ls, h, L = a
    o = C.nariai_condensate(0.0, frame="corner", h=h, lam_f=lf, lam_s=ls, Lu=L, Lv=L)
    cr = C.constraint_residual(o["res"], r_lo=0.5 * o["r_N"], r_hi=3 * o["r_N"], s_max=5.0)
    return dict(tag=tag, lam_f=lf, lam_s=ls, h=h, L=L, grew=bool(o["grew"]), Evv_rel=cr["Evv_rel"], Euu_rel=cr["Euu_rel"])


def history(tag):
    n0, lf = (-6e-6, 0.0) if tag == "initial_charge" else (0.0, 4e-7)
    o = C.nariai_condensate(n0, frame="corner", h=0.04, lam_f=lf)
    res, rN = o["res"], o["r_N"]
    n = C.stress(res)["n"]
    k = np.arange(min(n.shape))
    T, nd, rd = res["T"][k, k], n[k, k], res["r"][k, k]
    ok = np.isfinite(T) & np.isfinite(nd)
    return tag, T[ok] / rN, nd[ok] * rN**2, rd[ok] / rN


def modes(a):
    tag, x, nu, rr, f = a
    rN = 1e6
    xm = x.max()
    if tag == "initial_charge":
        nu_f = lambda X: np.where(X <= xm, np.interp(X, x, nu), nu[-1] * np.exp(-(X - xm)))
        r_f = lambda X: rN * np.where(X <= xm, np.interp(X, x, rr), rr[-1])
    else:
        nu_f = lambda X: np.where(X <= xm, np.interp(X, x, nu), nu[-1])
        r_f = lambda X: rN * np.where(X <= xm, np.interp(X, x, rr), rr[-1] * np.exp((X - xm) / np.sqrt(3)))
    xs = np.linspace(0, 9.6, 2001)
    I = float(np.trapezoid(np.abs(nu_f(xs)), xs))
    M3rN = I / 4 * f
    M = (M3rN / rN) ** (1 / 3)
    eps_t = lambda t: nu_f(t / rN) / rN**2 / (2 * M**4)
    kmax = np.sqrt(np.abs(nu_f(np.linspace(0, 9.6, 200))).max() / rN**2 / (2 * M**4) * 2 * M * M) * 3
    g = GC.mode_growth(M, 0.0, rN, t_end_over_rN=9.6, eps_t=eps_t, r_t=lambda t: float(r_f(t / rN)), ks=np.geomspace(1e-2 / rN, kmax, 30),
                       ls=np.unique(np.concatenate([np.arange(6), np.geomspace(6, max(7, kmax * rN), 15).astype(int)])))
    return dict(tag=tag, I_3tau=I, M3rN=float(M3rN), factor=f, window_lo=I / 4, window_hi=GC.window(rN)["M3rN_max"], **g)


def gsl_case(job):
    tag, val = job
    if tag == "fill":
        o = C.nariai_condensate(val, frame="corner", h=0.04)
    elif tag == "tilt":
        o = C.nariai_condensate(0.0, frame="corner", h=0.04, lam_f=val)
    elif tag == "route_tilt":
        o = C.route_condensate(lam_f=val)
    else:
        res, cl, g, mt = PR.run(7.258540773391722, N_amp=val, u_p=1.0, w_p=0.1, v_shell=4.5, w_shell=0.25, Lv=600.0, Lu=30.0)
        o = dict(res=res, cls=cl["cls"], r_N=mt["r_N"], grew=g["grew"])
    dh = GS.daughter_horizons(o, n_points=12)
    eh = GS.parent_event_horizon(o)
    return dict(case=f"{tag}:{val}", grew=bool(o["grew"]), parent_deficit=eh.get("deficit"), max_daughter_deficit=dh.get("max_deficit"),
                S_N=dh.get("S_N"))


def transient(seed):
    """Proposition-N 6e-7: lam_f = 2e-7, eps = 4 (T_exit = rho_f / (3 lam_f) = 221); seed = 1e-6 v off the false vacuum lets the
    destabilised false vacuum roll; seed = 0 is the control that cannot exit."""
    h, L = 0.08, 520.0
    o = C.nariai_condensate(0.0, frame="corner", h=h, lam_f=2e-7, Lu=L, Lv=L, clip=1, dphi0=seed)
    res, Cc = o["res"], o["cls"]
    R, P, T = res["r"], res["phi"], res["T"]
    dau = np.isfinite(R) & (R > 5 * o["r_N"]) & (P < 0) & (Cc == 1)
    windows = {}
    for lo, hi in ((0, 200), (200, 320), (320, 420), (420, 520)):
        sl = slice(int(lo / h), int(hi / h))
        fin = np.isfinite(R[sl])
        windows[f"{lo}-{hi}"] = dict(daughter_points=int(dau[sl].sum()), r_max=float(np.nanmax(np.where(fin, R[sl], np.nan))),
                                     max_phi_over_v=float(np.nanmax(np.where(fin, P[sl], np.nan)) / 0.0015))
    band = np.zeros_like(P, dtype=bool)
    band[:, : int(25 / h)] = True                                       # inside the old wall position (v2 = 32.9): band and daughter
    first_exit = np.argwhere(band & np.isfinite(P) & (P / 0.0015 > -0.5))
    eh = GS.parent_event_horizon(o)
    return dict(seed=seed, grew=bool(o["grew"]), T_exit_pred=RHO_F / (3 * 2e-7), windows=windows, parent_deficit=eh.get("deficit"),
                hiding=C.hiding_report(o),
                T_at_first_exit=float(T[tuple(first_exit[np.argmin(T[tuple(first_exit.T)])])]) if len(first_exit) else None)


if "--figure-only" in sys.argv:                  # redraw from DATA/round6e.json without recomputing
    out = json.load(open(f"{FIG}/../DATA/round6e.json"))
    p1, p2, p3 = out["problem1"], out["problem2"], out["problem3"]
else:
    out = {}
    # ------------------------------------------------ problem 1: verification of the dynamical condensate
    ver = dict(schwarzschild=[], tilt=[])
    M_, r0_ = 1.0, 6.0
    f0_ = 1 - 2 * M_ / r0_
    lump = lambda x: 2e-3 * np.exp(-((np.asarray(x, float) - 2.0) / 0.8) ** 2)
    for h in (0.1, 0.05, 0.025):
        N = int(6 / h)
        A, B, _ = C.condensate_data(0.0, 0.0, r0_, -np.sqrt(f0_) / 2, np.sqrt(f0_) / 2, h, N, N, lump, lump, frame="corner")
        cr = C.constraint_residual(C.evolve(A, B, 0.0, 0.0, TE.Potential(), h))
        ver["schwarzschild"].append(dict(h=h, Evv=cr["Evv"], Euu=cr["Euu"]))
    z = lambda x: 0 * np.asarray(x, float)
    A, B, _ = C.condensate_data(0.0, 0.0, r0_, -np.sqrt(f0_) / 2, np.sqrt(f0_) / 2, 0.05, 200, 200, z, z, frame="corner")
    a, b, _ = TE.characteristic_data(0.0, 0.0, r0_, -np.sqrt(f0_) / 2, np.sqrt(f0_) / 2, 0.05, 200, 200)
    ver["zero_condensate_reproduces_diamond_solver"] = float(np.nanmax(np.abs(C.evolve(A, B, 0.0, 0.0, TE.Potential(), 0.05)["r"]
                                                                              - TE.evolve(a, b, 0.0, 0.0, TE.Potential(), 0.05)["r"])))
    with ProcessPoolExecutor(4) as ex:
        ver["tilt"] = list(ex.map(tilt_diag, [("consistent", 4e-7, 4e-7, h, 40.0) for h in (0.08, 0.04, 0.02)] +
                                  [("drift_only", 4e-7, 0.0, 0.08, 40.0), ("drift_only_full", 4e-7, 0.0, 0.08, 200.0),
                                   ("source_only", 0.0, 4e-7, 0.08, 200.0)]))
    out["verification"] = ver

    # ------------------------------------------------ problem 1: thresholds, location, self-generation, route
    p1 = {}
    jobs = [dict(h=0.08, lo=5e-6, hi=1e-5), dict(h=0.04, lo=5e-6, hi=1e-5),
            dict(h=0.08, profile="lump", lump_width=4.0, lo=1e-6, hi=1e-3), dict(h=0.08, profile="lump", lump_width=8.0, lo=1e-6, hi=1e-3),
            dict(h=0.08, profile="lump", lump_width=16.0, lo=1e-6, hi=1e-3), dict(h=0.08, v_end=20.0, lo=1e-6, hi=1e-3),
            dict(h=0.08, u_end_frac=0.3, lo=1e-6, hi=1e-3), dict(h=0.08, u_end_frac=0.9, lo=1e-6, hi=1e-3)]
    with ProcessPoolExecutor(4) as ex:
        p1["initial_charge_thresholds"] = list(ex.map(thr_fill, jobs))
        p1["tilt_thresholds"] = list(ex.map(thr_tilt, [0.08, 0.04]))
        p1["route_tilt_threshold"] = ex.submit(thr_route, 0).result()
        hist = list(ex.map(history, ["initial_charge", "self_generated"]))
        p1["stability"] = list(ex.map(modes, [(t, x, nu, rr, f) for (t, x, nu, rr) in hist for f in (0.5, 1.5, 5.0)]))
    p1["route_runs"] = []
    for lf in (0.0, 3e-7, 1e-6, 3e-6):
        o = C.route_condensate(lam_f=lf)
        R, P, Cc = o["res"]["r"], o["res"]["phi"], o["cls"]
        p1["route_runs"].append(dict(lam_f=lf, grew=bool(o["grew"]), n_big=o["n_big"], normal_on_data=o["normal_on_data"], **C.hiding_report(o)))
        del o
        gc.collect()
    for kk in ("h0.02_n4.7e-6", "h0.02_n5.1e-6"):
        nval = -float(kk.split("_n")[1])
        o = C.nariai_condensate(nval, frame="corner", h=0.02)
        p1[kk] = dict(grew=bool(o["grew"]), n_big=o["n_big"])
        del o
        gc.collect()
    p1["CLNS_ceiling_nu"] = 1 / (4 * np.pi)
    out["problem1"] = p1

    # ------------------------------------------------ problem 2: generalised entropy, eternal versus transient
    with ProcessPoolExecutor(3) as ex:
        p2 = dict(eternal=list(ex.map(gsl_case, [("fill", -6e-6), ("tilt", 1.3e-7), ("tilt", 4e-7), ("route_tilt", 1e-6), ("route_tilt", 3e-6),
                                                   ("route_push", 10.0)])))
        p2["transient"] = list(ex.map(transient, [1.5e-9, 0.0]))
    out["problem2"] = p2

    # ------------------------------------------------ problem 3: creating the ball by supercooling
    p3 = dict(static_ball=[], confined=[])
    for al in (2.0, 2.5, 3.0, 4.0):
        fb = FB.Fireball(eps=4.0, alpha_pt=al)
        p3["static_ball"].append(dict(alpha_pt=al, **FB.static_hot_ball(fb, T_over_Tc=1.4)))
    fb = FB.Fireball(eps=4.0, alpha_pt=1.0)
    xs = np.linspace(-0.7, 1.3, 4000) * fb.v
    for rT in np.linspace(1.0, 3.0, 401):
        Fp = fb.cT * (rT * fb.Tc) ** 2 * (xs + fb.v) + fb.dU(xs)
        if not np.any((Fp[:-1] < 0) & (Fp[1:] >= 0)):
            p3["spinodal_T1_over_Tc_eps4"] = float(rT)
            break
    for al, Ti, sb in ((1.0, 1.8, 0.6), (2.0, 1.8, 0.75), (0.5, 1.6, 0.6)):
        fb = FB.Fireball(eps=4.0, alpha_pt=al)
        rho_i = fb.thermo(fb.v, 0.0, Ti * fb.Tc)[0]
        ad = 0.95 * np.sqrt(1 - sb * sb) / sb
        a_i = np.sqrt((1 + ad**2) / (8 * np.pi / 3 * rho_i))
        p3["confined"].append(dict(alpha_pt=al, Ti_over_Tc=Ti, **fb.summary(fb.run(Ti * fb.Tc, a_i, t_max=400.0), sin_chi_b=sb)))
    meta = PR.build(7.258540773391722, v_shell=4.5, w_shell=0.25, Lv=1500.0, Lu=2.0)
    a, b, info, pot, mt = meta
    m_ray = PR._ray_mass(a, b, info, 0.04, 1.0, 0.0, 0.0, pot)[0]
    p3["route_data"] = dict(normal_everywhere=mt["normal_everywhere"], r_wall_over_lH=float(np.interp(mt["v_wall"], np.arange(len(a[0])) * 0.04, a[0]) / 30),
                            ADM_mass=float(np.nanmedian(m_ray[-200:])), M_N=mt["M_N"], delta=mt["delta"])
    out["problem3"] = p3
    dump("round6e.json", out)

# ------------------------------------------------ figure
C1, C2, C3, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(1, 4, figsize=(19, 4.4))
lumps = [r for r in p1["initial_charge_thresholds"] if r["kw"].get("profile") == "lump"]
wd = [r["kw"]["lump_width"] for r in lumps]
nth = [np.sqrt(r["n_th"][0] * r["n_th"][1]) * 300 for r in lumps]
ax[0].plot(wd, nth, "o-", color=C2, ms=6, label="threshold |n| r_N^2 (lump at the corner)")
ax[0].plot(wd, [nth[0] * wd[0] / w for w in wd], "--", color=MUTED, lw=1, label="charge x width = const")
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel("width of the NEC-phase lump (u, v units; r_N = 17.3)")
ax[0].set_title("(a) the condensate tips the throat by itself", fontsize=9)
ax[0].legend(frameon=False, fontsize=8)
o = C.nariai_condensate(0.0, frame="corner", h=0.08, lam_f=4e-7)
Cm, Pm = o["cls"], o["res"]["phi"]
cmap = matplotlib.colors.ListedColormap(["#9cc3ee", "#f4f3ef", "#f2a07e"])
ax[1].imshow(Cm, extent=[0, Cm.shape[1] * 0.08, Cm.shape[0] * 0.08, 0], aspect="auto", cmap=cmap, vmin=-1, vmax=1)
ax[1].contour(np.arange(Cm.shape[1]) * 0.08, np.arange(Cm.shape[0]) * 0.08, Pm, levels=[0], colors=INK, linewidths=1)
ax[1].set_xlabel("v"); ax[1].set_ylabel("u")
ax[1].set_title("(b) self-generated: zero charge on the data, tilt makes it (orange = anti-trapped)", fontsize=8)
names = [e["case"] for e in p2["eternal"]]
ax[2].barh(range(len(names)), [e["max_daughter_deficit"] or 0 for e in p2["eternal"]], color=C2, label="daughter-observer deficit (eternal daughter)")
ax[2].barh(range(len(names)), [e["parent_deficit"] or 0 for e in p2["eternal"]], color=C1, label="parent event horizon")
ax[2].set_yticks(range(len(names))); ax[2].set_yticklabels(names, fontsize=7)
ax[2].set_xlabel("largest decrease of A/4G (nats); S_N = 942")
ax[2].set_title("(c) eternal daughters violate the local GSL (Theorem 6e-6)", fontsize=9)
ax[2].legend(frameon=False, fontsize=7)
fb = FB.Fireball(eps=4.0, alpha_pt=3.0)
from scipy.integrate import solve_ivp
TA = 1.4 * fb.Tc
R = 0.6 * 30
f = lambda t, y: [y[1], -fb.Gamma * y[1] - (fb.cT * (TA * np.exp(-max(t - 6, 0) / R)) ** 2 * (y[0] + fb.v) + fb.dU(y[0]))]
s = solve_ivp(f, (0, 6 + 4 * R), [fb.v, 0.0], method="DOP853", rtol=1e-10, atol=1e-14, max_step=0.01)
ax[3].plot(s.t, s.y[0] / fb.v, color=C3, lw=1.2)
ax[3].axhline(-1, color=INK, lw=0.8); ax[3].axhline(1, color=INK, lw=0.8, ls=":")
ax[3].axvline(6, color=MUTED, lw=0.8, ls="--")
ax[3].set_xscale("symlog", linthresh=1)
ax[3].set_xlim(0, 6 + 4 * R)
ax[3].set_xlabel("time (light-crossing of the ball = 18)"); ax[3].set_ylabel("phi / v")
ax[3].set_title("(d) ball created: heating moves phi to the false vacuum;\nit stays there while the radiation leaves", fontsize=8)
fig.tight_layout()
fig.savefig(f"{FIG}/round6e.png", dpi=150, facecolor="#fcfcfb")
print("done")
