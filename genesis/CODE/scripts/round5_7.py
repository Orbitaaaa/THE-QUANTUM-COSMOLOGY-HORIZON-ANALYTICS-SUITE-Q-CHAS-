"""Rounds 5 and 7: material properties from a field theory, and radiation from an explicit quantum field."""
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse import field_theory as F, radiation as Rd
from babyuniverse.model import model_from_dimensionless, ThinWallModel

warnings.simplefilter("ignore")
out = {}
out["frontier_couplings"] = {f"eps={e}": F.couplings_for(2.7e17, 2.17e17, e) for e in (0.2, 0.3, 0.5)}
out["perturbative_frontier"] = {f"lam={l:.3g},eps={e}": F.perturbative_frontier(l, e) for l in (1.0, 4 * np.pi) for e in (0.3, 0.5)}
out["theories"] = [dict(lam=l, v=v, eps=e, **F.derived_parameters(l, v, e), chain=F.laboratory_exponent_from_theory(l, v, e))
                   for l, v, e in ((1.0, 0.02, 0.3), (4.0, 0.01, 0.5))]
out["genesis_scales"] = {f"f={f:g},Lambda={L:g}": F.genesis_scales(f, L) for f, L in ((1e10, 1e5), (1e3, 1e-3))}
models = {"gamma=0.6,E=2.4": model_from_dimensionless(0.6, 2.4), "gamma=0.3,E=5": model_from_dimensionless(0.3, 5.0),
          "light wall kappa=0.3": ThinWallModel(1.0, 0.3, 0.1)}
fl = {k: Rd.flux_lab(m) for k, m in models.items()}
out["radiation"] = {k: {kk: vv for kk, vv in f.items() if not isinstance(vv, np.ndarray)} for k, f in fl.items()}
dump("round5_7.json", out)

C1, C2, C3, C4, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0})
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.1))
lam = np.logspace(-1, 4.6, 200)
for e, col in ((0.3, C1), (0.5, C2)):
    ax[0].loglog(lam, [F.perturbative_frontier(l, e)["B_lab_min"] for l in lam], color=col, label=rf"frontier from the field theory, $\epsilon={e}$")
ax[0].axhline(52.8, color=MUTED, lw=1, ls="--")
ax[0].text(0.12, 60, "frontier assumed in Round 2 ($B_{\\rm CdL}=10$): $B=52.8$", fontsize=7.5, color=MUTED)
ax[0].axvspan(4 * np.pi, 4e4, color="#f1f0ec", zorder=0)
ax[0].text(30, 4e3, "one loop $\\gtrsim$ tree\n(non-perturbative)", fontsize=8, color=MUTED)
ax[0].set_xlabel(r"quartic coupling $\lambda$"); ax[0].set_ylabel(r"minimum laboratory exponent $B_{\rm lab}$")
ax[0].set_title("Round 5: what a real field theory allows", loc="left", fontsize=10)
ax[0].legend(fontsize=7.5, frameon=False, loc="lower left")
for (k, f), col in zip(fl.items(), (C1, C2, C3)):
    ax[1].plot(f["u"] * f["kappa"], f["F"] / f["F_hawking"], color=col, label=f"collapse branch, {k}")
ax[1].plot([-0.5, 0, 0, 7], [np.nan, np.nan, 1, 1], color=C4, lw=3, ls="--", label="daughter branch (state D1, assumed)")
ax[1].axhline(0, color=MUTED, lw=0.8)
ax[1].set_xlim(-0.2, 7); ax[1].set_ylim(-0.6, 1.2)
ax[1].set_xlabel(r"$(u-u_0)/4GM$ after the turning-point signal"); ax[1].set_ylabel(r"$\langle T_{uu}\rangle$ / Hawking flux")
ax[1].set_title("Round 7: radiation from an explicit 1+1 field", loc="left", fontsize=10)
ax[1].legend(fontsize=7.5, frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "round5_7.png"), dpi=150)
print({k: (round(v["switch_on_in_4GM"], 2), round(v["excess_in_thermal_units"], 2)) for k, v in out["radiation"].items()})
