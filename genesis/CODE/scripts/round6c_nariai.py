"""Round 6c (positive result): the Nariai throat.  A near-Nariai Schwarzschild-de Sitter hole inside a false-vacuum bubble, pushed by a
small outgoing negative-energy pulse, grows a hidden inflating daughter; QNEC-bounded CFT pulses do not, because of their back-reaction."""
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse import throat_evolution as TE

warnings.simplefilter("ignore")
out = {}
strip = lambda d: {k: v for k, v in d.items() if k not in ("res", "cls")}

# threshold in delta m versus Nariai tuning (classical null fluid, width 1)
thr = []
for delta in (1e-4, 1e-6, 1e-8):
    for N in (0.01, 0.02, 0.03, 0.04):
        thr.append(dict(delta=delta, N=N, **strip(TE.nariai_false_vacuum(delta=delta, Nflux=N))))
out["threshold_scan"] = thr

# convergence and hiding of a growing case
conv = {}
for h in (0.1, 0.05):
    conv[f"h={h}"] = strip(TE.nariai_false_vacuum(delta=1e-8, Nflux=0.04, h=h, Lu=200, Lv=200))
out["growing_case_convergence"] = conv

# robustness: wall position and tuning
rob = []
for wf in (0.9, 0.95, 0.995):
    for delta in (1e-6, 1e-10):
        rob.append(dict(wall_frac=wf, delta=delta, **strip(TE.nariai_false_vacuum(delta=delta, Nflux=0.04, wall_frac=wf))))
out["robustness"] = rob

# QNEC-saturating pulses: back-reaction of the same species pushes the region towards collapse
qn = []
for Y in (0.01, 0.02, 0.06, 0.1):
    for dm in (0.03, 0.06):
        c = 24 * np.pi * dm / Y
        d = TE.nariai_false_vacuum(delta=1e-8, c=c, Y=Y, wp=0.2 / Y)
        qn.append(dict(Y=Y, c=float(c), k_over_r2=float(c / (12 * np.pi * d["r_N"] ** 2)), **strip(d)))
for c in (50.0, 200.0):
    d = TE.nariai_false_vacuum(delta=1e-8, c=c, Y=0.0)
    qn.append(dict(Y=0.0, c=c, k_over_r2=float(c / (12 * np.pi * d["r_N"] ** 2)), **strip(d)))
out["qnec_pulses"] = qn
dump("round6c_nariai.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED,
                     "ytick.color": MUTED})
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))
d = TE.nariai_false_vacuum(delta=1e-8, Nflux=0.04, h=0.08)
R, P, C = d["res"]["r"], d["res"]["phi"], d["cls"]
h = 0.08
ext = [0, R.shape[1] * h, R.shape[0] * h, 0]
cmap = matplotlib.colors.ListedColormap(["#9cc3ee", "#f4f3ef", "#f2a07e"])
ax[0].imshow(C, extent=ext, aspect="auto", cmap=cmap, vmin=-1, vmax=1)
ax[0].contour(np.arange(R.shape[1]) * h, np.arange(R.shape[0]) * h, P, levels=[0], colors=INK, linewidths=1)
ax[0].set_xlim(0, 60)
ax[0].set_ylim(120, 0)
ax[0].text(4, 100, "anti-trapped false vacuum:\ninflates to I+ (r > 10^4)", fontsize=8, color=INK)
ax[0].text(36, 30, "trapped true vacuum:\nthe parent's black hole\n(neck; daughter hidden)", fontsize=8, color=INK)
ax[0].set_xlabel("v")
ax[0].set_ylabel("u")
ax[0].set_title("(a) Nariai throat + 1.2% mass pulse: hidden growing daughter", fontsize=9)
for delta, col in ((1e-4, C1), (1e-6, C2), (1e-8, C3)):
    rows = [t for t in thr if t["delta"] == delta]
    ax[1].plot([t["delta_m"] / t["M_N"] * 100 for t in rows], [t["r_max"] / t["r_N"] if t["r_max"] else 1 for t in rows], "o-",
               color=col, ms=6, mec="#fcfcfb", label=f"delta = {delta:g}")
ax[1].axhline(5, color=INK, lw=1)
ax[1].text(0.32, 5.4, "growth criterion (r > 5 r_N, coherent region)", fontsize=8, color=INK)
ax[1].set_yscale("log")
ax[1].set_xlabel("pulse delta m / M_N  (%)")
ax[1].set_ylabel("largest r / r_N in the anti-trapped false vacuum")
ax[1].set_title("(b) Threshold 0.6-0.9% of the hole mass (flat space needs O(1))", fontsize=9)
ax[1].legend(frameon=False, fontsize=8, loc="upper left")
qq = [q for q in qn if q["Y"] > 0]
ax[2].scatter([q["k_over_r2"] for q in qq], [q["r_max"] / q["r_N"] if q["r_max"] else 1 for q in qq], color=C2, s=40,
              edgecolor="#fcfcfb", label="QNEC-saturating pulses (delta m = 0.5-1% M_N)")
ax[2].axhline(5, color=INK, lw=1)
ax[2].set_xscale("log")
ax[2].set_yscale("log")
ax[2].set_ylim(0.8, 20)
ax[2].set_xlabel("c / (12 pi r_N^2)  (strength of the same species' back-reaction)")
ax[2].set_ylabel("largest r / r_N in the anti-trapped false vacuum")
ax[2].set_title("(c) QNEC pulses fail: their species' back-reaction wins", fontsize=9)
ax[2].legend(frameon=False, fontsize=8, loc="upper left")
fig.tight_layout()
fig.savefig(f"{FIG}/round6c_nariai.png", dpi=150, facecolor="#fcfcfb")
print("wrote FIGURES/round6c_nariai.png")
