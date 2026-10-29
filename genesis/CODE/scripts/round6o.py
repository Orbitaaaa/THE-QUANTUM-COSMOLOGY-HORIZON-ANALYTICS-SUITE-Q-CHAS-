"""Round 6o: collider formation, interior evolution and evaporation in one description (model U): the evaporating Vaidya geometry, the
fate of the bounced partons and of the daughter, the per-event timeline, and an event generator for the observables.  Writes
DATA/round6o.json and FIGURES/round6o.png; `--figure-only` redraws the figure."""
import json
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import pipeline6o as P
from babyuniverse import collider6l as C

warnings.simplefilter("ignore")

if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round6o.json"))
else:
    L = C.load_lumi(f"{DATA}/round6l_lumi.json")
    out = {}
    out["vaidya"] = [P.vaidya_history(n, cl, m, k, n_pts=900) for n in (4, 6) for cl in (0.6, 0.75) for m in (1.5, 3.0, 10.0)
                     for k in ("hayward", "dymnikova")]
    out["timeline"] = [P.event_timeline(n, md, cl, M, k) for (n, md, cl, M, k) in
                       ((6, 8.0, 0.6, 25.0, "hayward"), (6, 8.0, 0.6, 40.0, "hayward"), (4, 6.0, 0.6, 25.0, "hayward"),
                        (6, 8.0, 0.7, 25.0, "dymnikova"), (7, 10.0, 0.6, 40.0, "hayward"))]
    gens = {}
    for name, kind, cl, mth in (("U-hayward", "hayward", 0.6, None), ("U-dymnikova", "dymnikova", 0.7, None), ("GR", "gr", 0.6, 11.5)):
        g = P.Generator(L, 100.0, 6, 8.0, cl, kind, seed=11, M_th_TeV=mth)
        gens[name] = g.run(6000)
    out["generator"] = gens
    dump("round6o.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
ax = ax.ravel()
v = [q for q in out["vaidya"] if q["n"] == 6 and q["c_l"] == 0.6 and q["mu0_over_ext"] == 3.0 and q["kind"] == "hayward"][0]
vv = np.array(v["v"]) / v["t_scale"]
ax[0].semilogx(vv[1:], v["r_plus"][1:], color=C1, label="outer trapping horizon r_+(v)")
ax[0].semilogx(vv[1:], v["r_minus"][1:], color=C3, label="inner trapping horizon r_-(v)")
ax[0].semilogx(vv[1:], v["r_shell"][1:], color=C2, ls="--", label="bounced partons (outgoing from the regular centre)")
ax[0].set_xlabel("advanced time v / evaporation time scale"); ax[0].set_ylabel("areal radius r / l")
ax[0].set_title("(a) one event (n = 6, mu0 = 3 mu_ext): the hole evaporates toward its remnant, the horizons\napproach but never meet, the bounced partons stay inside the core", fontsize=8)
ax[0].legend(frameon=False, fontsize=7)
for q, col in zip([q for q in out["vaidya"] if q["n"] == 6 and q["c_l"] == 0.6 and q["kind"] == "hayward"], (C1, C2, C4)):
    x = np.array(q["v"][1:]) / q["t_scale"]
    ax[1].loglog(x, np.array(q["mu_over_ext"][1:]) - 1, color=col, label=f"mu0 = {q['mu0_over_ext']} mu_ext: mu - mu_ext ~ v^-{q['mu_minus_ext_power']:.2f}")
    sep = np.array(q["r_plus"][1:]) - np.array(q["r_minus"][1:])
    ax[1].loglog(x, sep, color=col, ls=":", lw=1.2)
ax[1].set_xlabel("advanced time v / evaporation time scale"); ax[1].set_ylabel("mu/mu_ext - 1 (solid);  r_+ - r_- (dotted)")
ax[1].set_title("(b) T -> 0 at the remnant: mass excess ~ 1/v, horizon gap ~ v^-1/2 -> the trapped region\nnever closes, so the daughter is never released back into the parent", fontsize=8)
ax[1].legend(frameon=False, fontsize=7)
for q, col in zip([q for q in out["vaidya"] if q["n"] == 6 and q["c_l"] == 0.6 and q["kind"] == "hayward"], (C1, C2, C4)):
    x = np.array(q["v"][1:]) / q["t_scale"]
    ax[2].semilogx(x, q["affine"][1:], color=col, label=f"mu0 = {q['mu0_over_ext']} mu_ext: affine length -> {q['affine_final']:.2f} l")
ax[2].set_xlabel("advanced time v / evaporation time scale"); ax[2].set_ylabel("affine parameter of the bounced partons")
ax[2].set_title("(c) the partons reach the boundary of the parent's chart (v -> infinity) at FINITE affine\nparameter: the spacetime continues into the daughter region", fontsize=8)
ax[2].legend(frameon=False, fontsize=7)
G = out["generator"]
bins = np.linspace(10, 60, 26)
for name, col in (("U-hayward", C1), ("U-dymnikova", C3), ("GR", C2)):
    ax[3].hist(G[name]["M_TeV"], bins=bins, histtype="step", lw=2, color=col, label=f"{name}: sigma = {G[name]['sigma_pb']:.1f} pb")
ax[3].set_yscale("log"); ax[3].set_xlabel("black-hole mass (TeV)"); ax[3].set_ylabel("events (6000 generated)")
ax[3].set_title("(d) 100 TeV, n = 6, M_D = 8 TeV (NNPDF2.3 LO): production piles up at threshold\n(model U: M_ext = 11.5 TeV)", fontsize=8)
ax[3].legend(frameon=False, fontsize=7)
for name, col in (("U-hayward", C1), ("U-dymnikova", C3), ("GR", C2)):
    M = np.array(G[name]["M_TeV"]); N = np.array(G[name]["N"])
    edges = np.array([11.5, 13, 15, 18, 22, 28, 36, 50])
    mids, mean = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (M >= a) & (M < b)
        if sel.sum() > 3:
            mids.append(0.5 * (a + b)); mean.append(N[sel].mean())
    ax[4].plot(mids, mean, "o-", color=col, label=name)
ax[4].set_xlabel("black-hole mass (TeV)"); ax[4].set_ylabel("mean number of Hawking quanta")
ax[4].set_title("(e) multiplicity: model-U holes near threshold emit almost nothing (T -> 0 at M_ext);\nGR holes (T ~ 5 TeV) decay at once", fontsize=8)
ax[4].legend(frameon=False, fontsize=7)
for name, col in (("U-hayward", C1), ("U-dymnikova", C3), ("GR", C2)):
    M = np.array(G[name]["M_TeV"]); f = np.array(G[name]["f_miss"])
    ax[5].scatter(M, f, s=3, alpha=0.4, color=col, label=name)
mm = np.linspace(11.5, 60, 100)
ax[5].plot(mm, 11.5 / mm, color=INK, lw=1, ls="--", label="M_ext / M (remnant)")
tl = out["timeline"][0]
ax[5].set_xlabel("black-hole mass (TeV)"); ax[5].set_ylabel("missing-energy fraction (remnant + neutrinos)")
ax[5].set_title(f"(f) the remnant carries M_ext/M of the hole: formation {tl['t_formation_s']:.0e} s, bounce {tl['t_bounce_s']:.0e} s,\nevaporation {tl['t_evaporation_s']:.0e} s (25 TeV hole), daughter persists: {tl['daughter_persists']}", fontsize=8)
ax[5].legend(frameon=False, fontsize=7, markerscale=4)
fig.tight_layout()
fig.savefig(f"{FIG}/round6o.png", dpi=150, facecolor="#fcfcfb")
print("done")
