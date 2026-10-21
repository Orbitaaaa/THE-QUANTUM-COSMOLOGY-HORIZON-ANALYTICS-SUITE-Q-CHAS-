"""Farhi-Guth audit of every branch in the atlas, and embedding diagrams of
the initial (laboratory) and final (inflating) time-symmetric slices of the
FGG transition."""
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse.causal import farhi_guth_audit, embedding_profile, nec_wall, time_symmetric_slice

warnings.simplefilter("ignore")
rows = []
for gam in np.linspace(-0.9, 0.9, 7):
    gmin, g1 = ThinWallModel.g_min(gam), 2 - 2 * gam
    for E in (3 * g1, 0.5 * (gmin + g1), 0.6 * gmin):
        m = model_from_dimensionless(gam, E)
        for r in farhi_guth_audit(m):
            r.update(gamma=gam, E=E, nec_wall_min=nec_wall(m))
            rows.append(r)
infl = [r for r in rows if r["future_inflating"]]
summary = dict(n_branches=len(rows), n_inflating=len(infl),
               inflating_with_past_singular_exterior=sum(bool(r["exterior_contains_past_singularity"]) for r in infl),
               inflating_reaching_antitrapped_interior=sum(bool(r["interior_reaches_antitrapped"]) for r in infl),
               min_NEC_wall=min(r["nec_wall_min"] for r in rows))
print(summary)
dump("farhi_guth_audit.json", dict(summary=summary, rows=rows))

# embedding diagrams for the FGG pair (regime A)
m = model_from_dimensionless(-0.5, 4.5)
R1, R2 = m.turning_points()
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, R, title in ((axs[0], R1, "initial slice (bounded branch at R1): regular laboratory data"),
                     (axs[1], R2, "final slice (bounce at R2): bag of gold + Einstein-Rosen bridge")):
    pr = embedding_profile(m, R, r_far=3.5 * max(R2, 2 * m.M))
    ri, zi = pr["interior"]
    re, ze = pr["exterior"]
    for sgn in (1, -1):
        ax.plot(sgn * ri, zi, color="#d35400", lw=2)
        ax.plot(sgn * re, ze, color="#2471a3", lw=2)
    ax.axhline(zi[-1], color="#c0392b", lw=0.6, ls=":")
    ts = time_symmetric_slice(m, R)
    txt = "; ".join(f"{e.kind} sphere at {e.where} (A={e.area:.2f})" for e in ts["extremal"]) or "no extremal spheres"
    ax.set_title(title + "\n" + txt, fontsize=8)
    ax.set_xlabel("r (embedding radius)")
    ax.set_ylabel("z (embedding height, increases along the slice)")
    ax.set_aspect("equal", adjustable="datalim")
fig.suptitle("Embedding diagrams of the K=0 slices, orange: false vacuum (dS cap), blue: exterior (Flamm)",
             fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fgg_slices_embedding.png"), dpi=150)
plt.close(fig)
