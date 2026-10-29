"""Round 6m: the global daughter spacetime of model U (maximal extension of Oppenheimer-Snyder collapse, curvature bound, regularity of
every horizon, geodesic completeness, the daughter region and the finite closed daughter).  Writes DATA/round6m.json and
FIGURES/round6m.png; `--figure-only` redraws the figure."""
import json
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import global6m as G
from babyuniverse import unified6l as U

warnings.simplefilter("ignore")

if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round6m.json"))
else:
    out = {}
    out["schwarzschild_check"] = G.kretschmann_schwarzschild_check()
    out["curvature_bound"] = [G.curvature_bound_scan(D, 1.0, k, np.geomspace(1e-3, 1e8, 23)) for D in (4, 5, 10) for k in ("hayward", "dymnikova")]
    out["frw_curvature"] = [G.frw_curvature(c, D) for D in (4, 10) for c in (1e2, 1e4, 1e8)]
    out["ball_history"] = [G.ball_history(D, 1.0, m, E) for D in (5, 10) for m in (2.0, 10.0, 100.0) for E in (0.3, 0.6, 0.9)]
    out["ball_no_horizon"] = G.ball_history(10, 1.0, 0.5)
    out["penrose"] = G.penrose_worldline(10, 1.0, 10.0, 0.6)
    out["completeness"] = [G.geodesic_completeness(D, 1.0, m) for D in (5, 10) for m in (2.0, 50.0)]
    out["kruskal"] = [dict(D=D, mu_over_ext=m, **G.kruskal_regularity(D, 1.0, m)) for D in (5, 10) for m in (2.0, 50.0)]
    out["daughter"] = [G.daughter_properties(D, 1.0, c, chi) for D in (4, 10) for c in (1e3, 1e6, 1e10) for chi in (0.5, np.pi / 3, 2.0)]
    out["daughter_4D_stellar"] = G.daughter_properties(4, 1.0, 2e38, np.pi / 3)
    # history of the ball's areal radius through the regions (for the figure)
    mu = 10 * U.extremal_mu(10)
    f = lambda r: float(U.static_f(r, mu, 10))
    E = 0.6
    rm, rp = U.static_horizons(mu, 10)
    Rmax = __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: f(r) - E**2, rp * (1 + 1e-12), rp * 1e4)
    Rb = __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: f(r) - E**2, 1e-8, rm * (1 - 1e-12))
    th = np.linspace(0, np.pi, 1200)
    R = Rb + (Rmax - Rb) * np.cos(th / 2) ** 2 * np.ones_like(th)
    R = np.concatenate([Rb + (Rmax - Rb) * np.cos(np.linspace(0, np.pi / 2, 600)) ** 2, (Rb + (Rmax - Rb) * np.cos(np.linspace(0, np.pi / 2, 600)) ** 2)[::-1]])
    tau = [0.0]
    for k in range(1, len(R)):
        rmid = 0.5 * (R[k] + R[k - 1])
        tau.append(tau[-1] + abs(R[k] - R[k - 1]) / np.sqrt(max(E**2 - f(rmid), 1e-12)))
    out["radius_history"] = dict(tau=tau, R=R.tolist(), r_plus=rp, r_minus=rm)
    dump("round6m.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig = plt.figure(figsize=(17, 9.6))
gs = fig.add_gridspec(2, 3)
ax0 = fig.add_subplot(gs[:, 0])
ax = [ax0] + [fig.add_subplot(gs[i, j]) for i in range(2) for j in (1, 2)]
# (a) Penrose diagram
p = out["penrose"]
names = {"I": "I  parent\nuniverse", "II": "II  black hole", "III": "III  de Sitter\ncore", "IV": "IV  white hole", "I'": "I'  new universe\n(daughter)"}
for k, c in p["centres"].items():
    c = np.array(c)
    d = np.array([[0, -1], [1, 0], [0, 1], [-1, 0], [0, -1]]) + c
    ax0.plot(d[:, 0], d[:, 1], color=MUTED, lw=0.8)
    ax0.text(c[0] + 0.15, c[1] - 0.1, names[k], fontsize=8, color=INK)
x, t = np.array(p["x"]), np.array(p["t"])
ax0.plot(x, t, ".", ms=1.5, color=C2)
ax0.fill_betweenx(t, np.minimum(x, x) - 0.0, x, where=None, color=C2, alpha=0.0)
ax0.annotate("ball surface: collapse, bounce at f = E^2,\nwhite-hole explosion into the new universe", xy=(x[len(x) // 2], t[len(t) // 2]),
             xytext=(-2.95, 4.6), fontsize=7, arrowprops=dict(arrowstyle="->", color=MUTED))
ax0.text(-2.9, 2.5, "closed FRW interior\n(left of the surface):\nthe finite daughter,\nbounces at a = l", fontsize=7, color=C2)
ax0.set_aspect("equal"); ax0.set_xticks([]); ax0.set_yticks([]); ax0.grid(False)
for sp in ax0.spines.values():
    sp.set_visible(False)
ax0.set_title("(a) global structure of model-U collapse (D = 10, mu = 10 mu_ext):\nevery horizon is regular; the daughter lies outside J^-(I+ of the parent)", fontsize=8)
# (b) curvature bound
for cb, col, ls in zip(out["curvature_bound"], (C1, C1, C3, C3, C2, C2), ("-", "--") * 3):
    rows = cb["rows"]
    ax[1].semilogx([r["mu"] for r in rows], [r["K_max_alpha2"] for r in rows], ls=ls, color=col, label=f"D = {cb['D']}, {cb['kind']}")
ax[1].set_yscale("log")
ax[1].set_xlabel("mass parameter mu (units of l^{D-3})"); ax[1].set_ylabel("max Kretschmann x l^4")
ax[1].set_title("(b) a solution-independent curvature bound: max K l^4 is the same for\nholes of every mass (vacuum exterior + core)", fontsize=8)
ax[1].legend(frameon=False, fontsize=7, ncol=2)
# (c) radius history
rh = out["radius_history"]
ax[2].plot(rh["tau"], rh["R"], color=C2)
ax[2].axhline(rh["r_plus"], color=INK, ls="--", lw=0.9); ax[2].axhline(rh["r_minus"], color=INK, ls=":", lw=0.9)
ax[2].text(0.05, rh["r_plus"] * 1.02, "r_+", fontsize=7); ax[2].text(0.05, rh["r_minus"] * 1.02, "r_-", fontsize=7)
ax[2].set_xlabel("proper time of the ball surface (units of l)"); ax[2].set_ylabel("areal radius R / l")
bh = [b for b in out["ball_history"] if b.get("horizon") and b["D"] == 10 and b["mu_over_ext"] == 10.0 and b["E"] == 0.6][0]
ax[2].set_title(f"(c) the surface goes I -> II -> III -> IV -> I' in finite proper time\n(parent horizon to daughter horizon: {bh['tau_parent_horizon_to_daughter_horizon']:.2f} l)", fontsize=8)
# (d) Kruskal regularity
k = [q for q in out["kruskal"] if q["D"] == 10 and q["mu_over_ext"] == 50.0][0]
for name, col in (("outer", C1), ("inner", C3)):
    rows = k[name]["rows"]
    ax[3].semilogx([r["eps"] for r in rows], [r["F_inside"] for r in rows], "o-", color=col, label=f"{name} horizon, inside")
    ax[3].semilogx([r["eps"] for r in rows], [r["F_outside"] for r in rows], "s--", color=col, label=f"{name} horizon, outside")
ax[3].invert_xaxis()
ax[3].set_xlabel("distance to the horizon |r - r_h| / r_h"); ax[3].set_ylabel("Kruskal metric function f e^{-2 kappa r*} / kappa^2")
ax[3].set_title("(d) both horizons are regular: the Kruskal metric function has the same\nfinite non-zero limit from both sides (analytic extension)", fontsize=8)
ax[3].legend(frameon=False, fontsize=7)
# (e) finite daughter
for D, col in ((4, C1), (10, C2)):
    rows = [d for d in out["daughter"] if abs(d["chi_b"] - np.pi / 3) < 1e-9 and d["a_max"] and round(np.log10(d["C"])) in (3, 6, 10) and (d["C"] in (1e3, 1e6, 1e10))]
    rows = [d for d in rows if abs(np.log(d["volume_growth"]) - (np.log(d["a_max"] / d["a_bounce"]) * (D - 1))) < 1e-6]
    ax[4].loglog([d["C"] for d in rows], [d["volume_growth"] for d in rows], "o-", color=col, label=f"D = {D}")
st = out["daughter_4D_stellar"]
ax[4].set_xlabel("ball mass parameter C (units of l^{D-3})"); ax[4].set_ylabel("daughter volume growth  V_max / V_bounce")
ax[4].set_title(f"(e) the finite daughter: the closed interior re-expands from a = l to a_max\n(4D stellar ball: volume x 10^{np.log10(st['volume_growth']):.0f})", fontsize=8)
ax[4].legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(f"{FIG}/round6m.png", dpi=150, facecolor="#fcfcfb")
print("done")
