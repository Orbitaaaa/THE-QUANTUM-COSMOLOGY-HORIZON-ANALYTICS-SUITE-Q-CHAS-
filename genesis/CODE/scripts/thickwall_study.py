"""Thick scalar-field walls: convergence of time-symmetric data to the thin-wall
junction (including the bag-of-gold case) and the O(4) CdL bounce versus the
closed-form thin-wall action."""
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse.model import ThinWallModel
from babyuniverse.thickwall import QuarticPotential, thin_wall_comparison, cdl_bounce
from babyuniverse.tunnelling import instanton_M0, B_canonical

warnings.simplefilter("ignore")
G = 1.0
chi = 1.0
rho_f = 3 * chi**2 / (8 * np.pi * G)
conv = {}
fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
for (kap, lw, label, col) in ((0.5, 0.7, "normal (eps_out=+1)", "C0"), (1.5, 0.7, "bag of gold (eps_out=-1)", "C3")):
    sigma = kap / (4 * np.pi * G)
    rows = []
    for w in (0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125):
        pot = QuarticPotential.from_width_tension(w, sigma, rho_f)
        c = thin_wall_comparison(pot, lw, G)
        rows.append({k: v for k, v in c.items() if k != "data"})
        if w == 0.01:
            d = c["data"]
            axs[0].plot(d["l"], d["r"], color=col, lw=1.6, label=label)
            axs[0].axvline(lw, color=col, lw=0.5, ls=":")
            axs[0].axhline(c["bridge_thin"], color=col, lw=0.6, ls="--")
    conv[label] = rows
    ws = np.array([r["w"] for r in rows])
    err = np.array([abs(r["jump"] / r["jump_thin"] - 1) for r in rows])
    errM = np.array([abs(r["M_far"] / r["M_thin"] - 1) for r in rows])
    axs[1].loglog(ws, err, "o-", color=col, label=f"{label}: junction")
    axs[1].loglog(ws, errM, "s--", color=col, label=f"{label}: mass formula")
axs[1].loglog([1e-3, 4e-2], [1e-3 / 0.0025 * 2.5e-3 * 0.4, 4e-2 / 0.0025 * 2.5e-3 * 0.4], "k:", label="O(w)")
axs[0].set_xlabel("proper radial distance l")
axs[0].set_ylabel("area radius r(l)")
axs[0].set_xlim(0, 3.0)
axs[0].legend(fontsize=7)
axs[0].set_title("time-symmetric scalar-field data (w=0.01); dashed: 2GM", fontsize=9)
axs[1].set_xlabel("wall width w")
axs[1].set_ylabel("relative deviation from thin wall")
axs[1].legend(fontsize=6)
axs[1].set_title("thin-wall limit of the Israel condition", fontsize=9)

cdl = {}
for kap, col in ((0.5, "C0"), (1.5, "C3")):
    sigma = kap / (4 * np.pi * G)
    m = ThinWallModel(chi, kap, 0.0, 0.0)
    d = instanton_M0(m)
    rows = []
    for w in (0.04, 0.02, 0.01, 0.005):
        pot = QuarticPotential.from_width_tension(w, sigma, rho_f)
        r = cdl_bounce(pot, G)
        rows.append(dict(w=w, B_numeric=r["B"], a_wall=r["a_wall"], constraint_max=r["constraint_max"],
                         B_thin=d["B_CdL_down"], R0=d["R0"], rel=abs(r["B"] / d["B_CdL_down"] - 1)))
        if w == 0.01:
            ic = r["i_cut"]
            xx = np.concatenate([[0.0, r["x_start"]], r["x"][: ic + 1]])
            ph = np.concatenate([[1.0, r["phi"][0] / pot.v], r["phi"][: ic + 1] / pot.v])
            aa = np.concatenate([[0.0, r["x_start"]], r["a"][: ic + 1]])
            axs[2].plot(xx, ph, color=col, lw=1.4, label=f"kappa={kap}: phi/v")
            axs[2].plot(xx, aa, color=col, lw=1.0, ls="--", label=f"kappa={kap}: a(xi)")
            axs[2].axvline(r["xi_pole"], color=col, lw=0.6, ls=":")
    cdl[f"kappa={kap}"] = dict(rows=rows, B_canonical_M0=B_canonical(m), p=d["p"], q=d["q"])
    print("CdL kappa", kap, [(r["w"], round(r["B_numeric"], 6), round(r["rel"], 5)) for r in rows],
          "thin", d["B_CdL_down"], "B_can(M=0)", B_canonical(m))
axs[2].set_xlabel(r"$\xi$")
axs[2].legend(fontsize=7)
axs[2].set_title("O(4) CdL bounces with gravity (w=0.01); dotted: far pole\n"
                 "(numerics cut at closest approach to the false vacuum, completed analytically)", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "thick_walls.png"), dpi=150)
plt.close(fig)
dump("thickwall_convergence.json", dict(time_symmetric=conv, cdl=cdl))
