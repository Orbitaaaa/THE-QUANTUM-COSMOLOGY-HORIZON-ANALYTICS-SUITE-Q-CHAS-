"""Round 6d: independent verification, an asymptotically connected (bridge-free) configuration, the complete route from data without a black
hole, ghost-condensate stability in the Nariai geometry, the generalised entropy of the solution, and the quantum-field threshold versus c."""
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse import (throat_evolution as TE, throat_mol as MOL, ghost_condensate as GC, gsl_solution as GS, preparation as PR,
                          quantum_throat as QT)

warnings.simplefilter("ignore")
out = {}
gauss = lambda A, up=2.0: (lambda u: -A * np.exp(-((np.asarray(u, float) - up) / 1.0) ** 2))
LIGHT = dict(delta=1e-4, eps=4.0, vphi=0.0015, h=0.04)

# ------------------------------------------------ 1. independent solver (method of lines) against the diamond scheme
ver = dict(rn=[], nariai=[])
Q, M, r0 = 10.0, 10.1, 13.0
f0 = 1 - 2 * M / r0 + Q**2 / r0**2
for h in (0.2, 0.1, 0.05):
    d = MOL.characteristic_data(r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, h, int(20 / h), int(20 / h), lambda x: 0.0, lambda x: 0.0, n_sub=4)
    ver["rn"].append(dict(h=h, max_mass_error=float(np.nanmax(np.abs(TE.misner_sharp(MOL.evolve(d, h, Q=Q)) - M)))))
for A in (0.0175, 0.02, 0.025, 0.03, 0.04):
    a = MOL.nariai(flux_u=gauss(A))
    b = TE.nariai_false_vacuum(Nflux=A)
    ver["nariai"].append(dict(Nflux=A, W=float(A * np.sqrt(np.pi)), mol_grew=bool(a["grew"]), mol_hidden=a["true_vacuum_untrapped_large_r"] == 0,
                              diamond_grew=bool(b["grew"]), diamond_hidden=b["true_vacuum_untrapped_large_r"] == 0))
out["independent_verification"] = ver

# ------------------------------------------------ 2. bridge-free configuration (light wall: exterior normal on the data)
lw = dict(runs=[], inner=[])
for N in (0.0, 0.04, 0.06, 0.08, 0.16):
    o = TE.nariai_false_vacuum(Nflux=N, **LIGHT)
    rv = o["rv_data"]
    lw["runs"].append(dict(Nflux=N, grew=bool(o["grew"]), n_big=o["n_big"], r_max=o["r_max"], exterior_rv_on_data=float(rv[-1]),
                           parent_EH=GS.parent_event_horizon(o), daughter=GS.daughter_horizons(o, n_points=8).get("max_deficit")))
for ui in (35.0, 60.0, 100.0, 150.0):
    o = TE.nariai_false_vacuum(Nflux=0.08, u_inner=ui, **LIGHT)
    lw["inner"].append(dict(u_inner=ui, grew=bool(o["grew"]), u_hole=o["u_hole"]))
iso = TE.nariai_false_vacuum(tu_func=gauss(0.08), Nv_func=lambda v: -0.08 * np.exp(-((np.asarray(v, float) - 2) / 1.0) ** 2), **LIGHT)
lw["isotropic"] = dict(grew=bool(iso["grew"]), n_big=iso["n_big"], parent_EH=GS.parent_event_horizon(iso),
                       daughter=GS.daughter_horizons(iso, n_points=8).get("max_deficit"))
out["bridge_free"] = lw

# ------------------------------------------------ 3. the complete route from data without a black hole
E = 7.258540773391722          # PR.tune(v_shell=4.5, w_shell=0.25, Lv=1500, Lu=2): largest E with normal, bridge-free data (delta ~ 0.009)
meta = PR.build(E, v_shell=4.5, w_shell=0.25, Lv=1500.0, Lu=2.0)[4]
assert meta["normal_everywhere"] and meta["wall_found"]
route = dict(E=E, meta={k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in meta.items()}, runs=[])
route_maps = {}
for N in (0.0, 3.0, 6.0, 10.0):
    res, cl, g, mt = PR.run(E, N_amp=N, u_p=1.0, w_p=0.1, v_shell=4.5, w_shell=0.25, Lv=600.0, Lu=30.0)
    C, P, R = cl["cls"], res["phi"], res["r"]
    tvn = (P > 0) & (C == 0) & np.isfinite(R)
    rows_tvn = np.where(tvn.any(axis=1))[0]
    ends_singular = []
    for i in (50, 80):
        js = np.where(np.isfinite(R[i]))[0]
        ends_singular.append(float(R[i, js.max()]))
    route["runs"].append(dict(N=N, grew=bool(g["grew"]), n_big=g["n_big"], r_max=g["r_max_antifalse"], data_normal=mt["normal_everywhere"],
                              last_u_with_normal_TV=float(rows_tvn.max() * 0.04) if len(rows_tvn) else None,
                              r_at_end_of_rows_u2_u3p2=ends_singular, parent_EH=GS.parent_event_horizon(dict(res=res, cls=C, r_N=mt["r_N"]))))
    if N in (0.0, 10.0):
        route_maps[N] = (C, P)
out["route"] = route

# ------------------------------------------------ 4. ghost condensate in the Nariai geometry
gc = dict(window={f"rN=1e{e}": GC.window(10.0**e) for e in (4, 6, 8)}, modes=[])
for rN in (1e4, 1e6):
    w = GC.window(rN)
    for f in (1.5, 5.0, 20.0, 60.0):
        M3rN = w["M3rN_min"] * f
        Mgc = (M3rN / rN) ** (1 / 3)
        e0 = -GC.eps0_for_tipping(Mgc, rN)
        gc["modes"].append(dict(r_N=rN, M3rN=M3rN, M=Mgc, eps0=e0, M_rN=Mgc * rN, **GC.mode_growth(Mgc, e0, rN, t_end_over_rN=9.6)))
out["ghost_condensate"] = gc

# ------------------------------------------------ 5. quantum fields: strong-QNEC bound and the threshold versus c
qf = dict(bound=QT.strong_qnec_weighted_bound(1.0, 1.0)["bound"].__str__(), table=[])
for c, lo, hi in ((300.0, 0.04, 0.1), (600.0, 0.04, 0.08), (1000.0, 0.08, 0.16), (1500.0, 0.08, 0.32)):
    for _ in range(5):
        N = 0.5 * (lo + hi)
        if TE.nariai_false_vacuum(c=c, Nflux=N)["grew"]:
            hi = N
        else:
            lo = N
    qf["table"].append(dict(c=c, x=c / 300.0, W_th=(float(lo * np.sqrt(np.pi)), float(hi * np.sqrt(np.pi))),
                            capacity=float(c / (48 * np.pi * 3.2 * 17.3205))))
out["quantum_fields"] = qf
dump("round6d.json", out)

# ------------------------------------------------ figure
C1, C2, C3, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))
Cm, Pm = route_maps[10.0]
cmap = matplotlib.colors.ListedColormap(["#9cc3ee", "#f4f3ef", "#f2a07e"])
h = 0.04
ax[0].imshow(Cm, extent=[0, Cm.shape[1] * h, Cm.shape[0] * h, 0], aspect="auto", cmap=cmap, vmin=-1, vmax=1)
ax[0].contour(np.arange(Cm.shape[1]) * h, np.arange(Cm.shape[0]) * h, Pm, levels=[0], colors=INK, linewidths=1)
ax[0].set_xscale("symlog", linthresh=10)
ax[0].set_xlabel("v")
ax[0].set_ylabel("u")
ax[0].set_title("(a) complete route: no hole in the data -> hole in false vacuum -> push -> hidden daughter", fontsize=8)
ax[0].text(1, 20, "inflating false vacuum\n(daughter)", fontsize=8)
ax[0].text(300, 2.5, "trapped:\nparent black hole", fontsize=8)
xs = [r["x"] for r in qf["table"]]
ax[1].fill_between(xs, [r["W_th"][0] for r in qf["table"]], [r["W_th"][1] for r in qf["table"]], color=C2, alpha=0.5)
ax[1].plot([0] + xs, [0.032] + [np.mean(r["W_th"]) for r in qf["table"]], "o-", color=C2, ms=5, label="threshold W_th (same c fields' back-reaction)")
ax[1].plot([0] + xs, [0] + [r["capacity"] for r in qf["table"]], "s--", color=C1, ms=5, label="most a CFT can supply, c / 48 pi tau")
ax[1].set_xlabel("x = c / r_N^2")
ax[1].set_ylabel("weighted deficit W (local, l_H = 30)")
ax[1].set_title("(b) quantum fields never catch up (strong-QNEC bound)", fontsize=9)
ax[1].legend(frameon=False, fontsize=8)
for rN, col in ((1e4, C1), (1e6, C3)):
    rows = [m for m in gc["modes"] if m["r_N"] == rN]
    ax[2].plot([m["M3rN"] for m in rows], [m["max_amplification"] for m in rows], "o-", color=col, ms=6, mec="#fcfcfb", label=f"r_N = {rN:.0e} l_P")
w6 = gc["window"]["rN=1e6"]
ax[2].axvspan(w6["M3rN_min"], w6["M3rN_max"], color=C3, alpha=0.12)
ax[2].axhline(np.e, color=INK, lw=1)
ax[2].set_xscale("log")
ax[2].set_xlabel("M^3 r_N  (ghost-condensate scale, Planck units)")
ax[2].set_ylabel("largest perturbation growth over 3 tau")
ax[2].set_title("(c) ghost condensate stable while it tips the throat", fontsize=9)
ax[2].legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(f"{FIG}/round6d.png", dpi=150, facecolor="#fcfcfb")
print("route:", route["runs"])
