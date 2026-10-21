"""Gravitational branch atlas: Penrose diagrams for every regime and a
phase diagram of the Schwarzschild-exterior model in the (gamma, E) plane."""
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse.branches import classify, bgg_regime
from babyuniverse.geometry import branch_trajectories, regions_along
from babyuniverse.plotting import branch_figure

warnings.simplefilter("ignore")


def comp(lst):
    out = [lst[0]] if lst else []
    for x in lst[1:]:
        if x != out[-1]:
            out.append(x)
    return out


cases = []
for gam in (0.6, -0.5):
    gmin = ThinWallModel.g_min(gam)
    g1 = 2 - 2 * gam
    for tag, E in (("A", 1.5 * g1), ("B", 0.5 * (gmin + g1)), ("C", 0.7 * gmin)):
        cases.append((f"gamma{gam:+.1f}_{tag}", model_from_dimensionless(gam, E)))
cases.append(("M0_flat_exterior", ThinWallModel(chi=1.0, kappa=0.5, M=0.0, eta=0.0)))
cases.append(("M0_dS_exterior", ThinWallModel(chi=1.0, kappa=0.5, M=0.0, eta=0.4)))
cases.append(("M0_dS_exterior_large_tension", ThinWallModel(chi=1.0, kappa=1.5, M=0.0, eta=0.4)))

records = []
for name, m in cases:
    for b in classify(m):
        rec = dict(case=name, regime=bgg_regime(m) if m.eta == 0 else "dS exterior",
                   params=m.describe(), branch=b.kind,
                   itinerary_out=comp(b.names_out()), itinerary_in=comp(b.names_in()),
                   flags={k: v for k, v in b.flags.items()})
        try:
            tr = branch_trajectories(m, b)
            rec["chart_regions_out"] = regions_along(tr["out"])
            rec["chart_regions_in"] = regions_along(tr["in"])
            rec["chart_agrees"] = (rec["chart_regions_out"] == rec["itinerary_out"]
                                   if tr["out"].chart == "Kruskal" else None)
            rec["chart_agrees_in"] = rec["chart_regions_in"] == rec["itinerary_in"]
            fn = os.path.join(FIG, f"penrose_{name}_{b.kind}.png")
            branch_figure(m, b, tr, fn, f"{name}: {rec['regime']}")
            rec["figure"] = os.path.basename(fn)
        except Exception as exc:  # documented, never silently dropped
            rec["chart_error"] = repr(exc)
        records.append(rec)
        print(name, b.kind, rec["itinerary_out"], rec["itinerary_in"], rec.get("chart_agrees"),
              rec.get("chart_agrees_in"))
dump("branch_atlas_cases.json", records)

# ---------------------------------------------------------------- phase diagram
gams = np.linspace(-0.95, 0.95, 191)
logE = np.linspace(-1.0, 1.5, 251)
code = np.zeros((len(logE), len(gams)))
flip_in_barrier = np.zeros_like(code)
eps_in_R2 = np.zeros_like(code)
for j, gam in enumerate(gams):
    gmin, g1 = ThinWallModel.g_min(gam), 2 - 2 * gam
    for i, le in enumerate(logE):
        E = 10**le * gmin
        if E < gmin:
            code[i, j] = 3
            continue
        code[i, j] = 1 if E > g1 else 2
        m = model_from_dimensionless(gam, E)
        tps = m.turning_points()
        if len(tps) == 2:
            fi = m.flip_radius_in()
            flip_in_barrier[i, j] = 1.0 if (fi and tps[0] < fi < tps[1]) else 0.0
            eps_in_R2[i, j] = np.sign(m.A_in(tps[1]))
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
cm = matplotlib.colors.ListedColormap(["#aed6f1", "#f9e79f", "#f5b7b1"])
ax[0].pcolormesh(gams, logE, code, cmap=cm, shading="auto", vmin=0.5, vmax=3.5)
ax[0].plot(gams, np.log10((2 - 2 * gams) / np.array([ThinWallModel.g_min(g) for g in gams])), "k-", lw=1)
ax[0].axhline(0, color="k", lw=1)
ax[0].set_xlabel(r"$\gamma=(\chi^2-\kappa^2)/(\chi^2+\kappa^2)$")
ax[0].set_ylabel(r"$\log_{10}(E/g_{\min})$   ($E\propto M^{-2/3}$)")
ax[0].text(0.0, 1.0, "A: bounded branch visible in region I\n+ inflating bounce in region III", ha="center", fontsize=8)
ax[0].text(-0.6, 0.03, "B: bounded branch hidden (III)", fontsize=8)
ax[0].text(-0.3, -0.5, "C: no barrier, monotonic\n(white hole -> III)", fontsize=8)
ax[0].set_title("Schwarzschild exterior: regimes", fontsize=10)
masked = np.ma.masked_where(code == 3, flip_in_barrier + 2 * (eps_in_R2 < 0))
cm2 = matplotlib.colors.ListedColormap(["#d5f5e3", "#fdebd0", "#d2b4de", "#a9cce3"])
ax[1].pcolormesh(gams, logE, masked, cmap=cm2, shading="auto", vmin=-0.5, vmax=3.5)
ax[1].axhline(0, color="k", lw=1)
ax[1].set_xlabel(r"$\gamma$")
ax[1].set_title(r"interior orientation at the inflating turning point $R_2$", fontsize=10)
ax[1].text(-0.8, 0.8, r"$\epsilon_{in}(R_2)=+1$", fontsize=9)
ax[1].text(0.35, 0.8, r"$\epsilon_{in}(R_2)=-1$" "\n" r"($\epsilon_{in}$ flips in barrier)", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "phase_diagram_schwarzschild.png"), dpi=150)
plt.close(fig)
dump("phase_diagram_schwarzschild.json", dict(gamma=gams, log10_E_over_gmin=logE, regime_code=code,
                                              flip_in_in_barrier=flip_in_barrier, eps_in_R2=eps_in_R2,
                                              legend={"1": "A", "2": "B", "3": "C"}))
