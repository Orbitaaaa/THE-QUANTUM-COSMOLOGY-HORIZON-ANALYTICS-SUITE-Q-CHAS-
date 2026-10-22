"""Stage 1 task 10: conformal diagrams of thin-wall worldlines in Schwarzschild-de Sitter exteriors (sds_atlas.py)."""
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse.model import ThinWallModel
from babyuniverse import sds_atlas as S, branches as Br

warnings.simplefilter("ignore")
C1, C2, MUTED, INK = "#2a78d6", "#eb6834", "#8a8984", "#0b0b0b"
CASES = [("usual: daughter beyond the bridge (eps_out(R2) = -1)", ThinWallModel(1.0, 0.3, 0.05, 0.5)),
         ("non-flip window: near-side daughter (eps_out(R2) = +1)", ThinWallModel(1.0, 0.3, 0.1, 0.95))]


def outlines(ax):
    diamonds = [0, -2, 2, -4]
    for c in diamonds:
        ax.plot([c - 1, c, c + 1, c, c - 1], [0, 1, 0, -1, 0], color=MUTED, lw=0.8)
    for xl, xr, T in ((-2, 0, 1), (-2, 0, -1), (0, 2, 1), (0, 2, -1), (-4, -2, 1), (-4, -2, -1)):
        ax.plot([xl, xr], [T, T], color=MUTED, lw=0.8, ls=(0, (2, 2)))
    ax.text(0, 0.05, "I", ha="center", color=MUTED); ax.text(-2, 0.05, "III", ha="center", color=MUTED)
    ax.text(-1, 0.7, "BH", ha="center", color=MUTED); ax.text(-1, -0.75, "WH", ha="center", color=MUTED)
    ax.text(1, 0.7, "cosmo\nfuture", ha="center", color=MUTED, fontsize=7); ax.text(1, -0.8, "cosmo\npast", ha="center", color=MUTED, fontsize=7)
    ax.text(-3, 0.7, "cosmo\nfuture", ha="center", color=MUTED, fontsize=7); ax.text(-3, -0.8, "cosmo\npast", ha="center", color=MUTED, fontsize=7)


fig, axes = plt.subplots(1, 2, figsize=(11, 3.3))
summary = []
for ax, (title, m) in zip(axes, CASES):
    outlines(ax)
    for b in Br.classify(m):
        if b.kind == "static-unstable":
            continue
        w = S.wall_worldline(m, b)
        chk = S.causal_check(w)
        col = C1 if b.kind == "bounded" else C2
        for k, p in enumerate(w):
            ax.plot(p["x"], p["T"], color=col, lw=2.0, label=b.kind if k == 0 else None)
        summary.append(dict(case=title, branch=b.kind, blocks=[p["block"] for p in w], matched=[p.get("matched") for p in w], **chk))
    ax.set_title(title, loc="left", fontsize=9)
    ax.set_aspect("equal"); ax.set_xlim(-4.6, 2.6); ax.set_ylim(-1.35, 1.35)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.legend(fontsize=8, frameon=False, loc="lower left")
fig.suptitle("Thin-wall worldlines in Schwarzschild-de Sitter exteriors (conformal atlas; dashed: r=0 and I$^\\pm$, schematic)", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "sds_atlas.png"), dpi=150)
dump("sds_atlas.json", summary)
for s in summary:
    print(s["case"][:20], s["branch"], s["blocks"], s["matched"], s["timelike_future"], f"{s['max_join_gap']:.1e}")
