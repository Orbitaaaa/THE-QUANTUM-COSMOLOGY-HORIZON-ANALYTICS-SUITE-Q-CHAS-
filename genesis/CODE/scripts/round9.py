"""Round 9: the experimental design (machine, luminosity, selections, detector requirements, decision tree).  Writes DATA/round9.json and
FIGURES/round9.png."""
import json
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import design9 as X
from babyuniverse import collider6l as C

warnings.simplefilter("ignore")

if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round9.json"))
    out["reach"] = X.reach_table(out["requirements"])
    out["reach_dymnikova"] = X.reach_table(out["requirements_dymnikova"])
    dump("round9.json", out)
else:
    L = C.load_lumi(f"{DATA}/round6l_lumi.json")
    out = {}
    out["requirements"] = X.luminosity_requirements(L, c_ls=np.linspace(0.58, 0.745, 12), n_events=2500)
    out["requirements_dymnikova"] = X.luminosity_requirements(L, machines=((85.0, 20.0), (100.0, 30.0)), kind="dymnikova",
                                                              c_ls=np.linspace(0.68, 0.87, 9), n_events=2000)
    out["reach"] = X.reach_table(out["requirements"])
    out["reach_dymnikova"] = X.reach_table(out["requirements_dymnikova"])
    out["detector"] = X.detector_requirements()
    out["decision_tree"] = X.decision_tree()
    dump("round9.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
ax = ax.ravel()
req = out["requirements"]
cols = {27.0: C4, 85.0: C1, 100.0: C2}
avail = {27.0: 15.0, 85.0: 20.0, 100.0: 30.0}
for rs, col in cols.items():
    rows = [r for r in req if r["sqrt_s"] == rs and r.get("reachable") and np.isfinite(r.get("L_discovery_ab", np.inf))]
    ax[0].semilogy([r["M_ext_TeV"] for r in rows], [r["L_discovery_ab"] for r in rows], "o-", color=col, label=f"{rs:.0f} TeV: 5 sigma discovery (monojet + remnant)")
    ax[0].axhline(avail[rs], color=col, ls=":", lw=1)
ax[0].set_ylim(1e-5, 1e3); ax[0].set_xlabel("model-U threshold M_ext (TeV)"); ax[0].set_ylabel("integrated luminosity needed (ab^-1)")
ax[0].set_title("(a) discovery luminosity in the monojet + invisible-remnant channel (dotted: planned\nluminosity); n = 6, M_D = 8 TeV, Hayward tower; order-of-magnitude backgrounds", fontsize=8)
ax[0].legend(frameon=False, fontsize=7)
for rs, col in ((85.0, C1), (100.0, C2)):
    rows = [r for r in req if r["sqrt_s"] == rs and r.get("reachable")]
    ax[1].semilogy([r["M_ext_TeV"] for r in rows], [r["L_class_ab"] for r in rows], "s--", color=col, label=f"{rs:.0f} TeV: 10 selected events (class + onset)")
    ax[1].semilogy([r["M_ext_TeV"] for r in rows], [r["L_remnant_mass_ab"] for r in rows], "^:", color=col, label=f"{rs:.0f} TeV: 10 charged remnants (if 30 % charged)")
    ax[1].axhline(avail[rs], color=col, ls=":", lw=1)
ax[1].set_ylim(1e-5, 1e3); ax[1].set_xlabel("model-U threshold M_ext (TeV)"); ax[1].set_ylabel("integrated luminosity needed (ab^-1)")
ax[1].set_title("(b) establishing the class (separation from GR needs ~3 events) and measuring the\nremnant mass directly with slow charged tracks, M = p/(beta gamma)", fontsize=8)
ax[1].legend(frameon=False, fontsize=7)
for rs, col in ((85.0, C1), (100.0, C2)):
    rows = [r for r in req if r["sqrt_s"] == rs and r.get("reachable") and r.get("bkg_fb", 0) > 0 and np.isfinite(r["L_discovery_ab"])]
    ax[2].semilogy([r["M_ext_TeV"] for r in rows], [r["sigma_sel_fb"] / r["bkg_fb"] for r in rows], "o-", color=col, label=f"{rs:.0f} TeV: S/B after selection")
ax[2].axhline(1, color=INK, lw=0.7)
ax[2].set_xlabel("model-U threshold M_ext (TeV)"); ax[2].set_ylabel("signal / background")
ax[2].set_title("(c) signal over background after the optimised selection (lead pT, MET > cut,\nMET/S_T > 0.8): signal-dominated for most of the window", fontsize=8)
ax[2].legend(frameon=False, fontsize=7)
b = np.linspace(0.05, 0.95, 100)
ax[3].semilogy(b, 10 / 0.2998 * (1 / b - 1), color=C4, label="time-of-flight delay at 10 m (ns)")
ax[3].semilogy(b, 1 / b**2, color=C3, label="dE/dx in units of a minimum-ionising particle")
ax[3].axvspan(0.2, 0.6, color=C4, alpha=0.12, label="predicted remnant velocities")
ax[3].set_xlabel("remnant velocity beta"); ax[3].set_ylabel("signature size")
ax[3].set_title("(d) detector requirements for charged remnants: ~50 ns delays and ~7x MIP ionisation\n(trivially separable with ~30-50 ps timing layers)", fontsize=8)
ax[3].legend(frameon=False, fontsize=7)
ax[4].axis("off")
import textwrap
y = 1.0
for d in out["decision_tree"]:
    obs = textwrap.fill(f"{d['step']}. {d['observation']}", 70)
    imp = textwrap.fill("-> " + d["implies"], 78)
    ax[4].text(0.0, y, obs, fontsize=6.8, color=INK, weight="bold", transform=ax[4].transAxes, va="top")
    y -= 0.045 * (obs.count("\n") + 1)
    ax[4].text(0.03, y, imp, fontsize=6.5, color=MUTED, transform=ax[4].transAxes, va="top")
    y -= 0.042 * (imp.count("\n") + 1) + 0.02
ax[4].set_title("(e) decision tree: what each outcome means for the interior class and the daughter", fontsize=8)
ax[5].axis("off")
rc = {r["sqrt_s"]: r for r in out["reach"]}
spec = [("machine", "FCC-hh, 85 TeV (100 TeV option), 20-30 ab^-1"),
        ("first step", "monojet / monophoton: n and M_D (P1)"),
        ("discovery channel", "monojet + invisible remnant, MET/S_T > 0.8"),
        ("5 sigma reach (M_ext)", "%.0f TeV at 85 TeV; %.0f TeV at 100 TeV" % (rc[85.0]["discovery_TeV"], rc[100.0]["discovery_TeV"])),
        ("10 selected events", "%.0f TeV at 85 TeV; %.0f TeV at 100 TeV" % (rc[85.0]["class_TeV"], rc[100.0]["class_TeV"])),
        ("10 charged tracks", "%.0f TeV at 85 TeV; %.0f TeV at 100 TeV" % (rc[85.0]["remnant_mass_TeV"], rc[100.0]["remnant_mass_TeV"])),
        ("class test", "~3 events separate model U from GR at 5 sigma"),
        ("triggers", "MET > 1 TeV; jet pT > 1 TeV; slow-track trigger"),
        ("detector", "hermetic to |eta| ~ 6; 30-50 ps timing; dE/dx"),
        ("earlier option", "HE-LHC 27 TeV: 5 sigma to M_ext ~ %.0f TeV" % rc[27.0]["discovery_TeV"])]
tb = ax[5].table(cellText=[[a, b] for a, b in spec], colLabels=["design element", "specification"], loc="upper center", cellLoc="left",
                 colWidths=[0.3, 0.7])
tb.auto_set_font_size(False); tb.set_fontsize(6.8); tb.scale(1, 1.55)
ax[5].set_title("(f) the experimental design", fontsize=8)
fig.subplots_adjust(left=0.05, right=0.98, top=0.93, bottom=0.06, wspace=0.28, hspace=0.38)
fig.savefig(f"{FIG}/round9.png", dpi=150, facecolor="#fcfcfb")
print("done")
