"""Round 7: experimentally distinguishable predictions of model U (lab-frame events with initial-state radiation and a generic detector;
model U vs a GR hole at the same threshold; backgrounds to order of magnitude).  Writes DATA/round7.json and FIGURES/round7.png."""
import json
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import predictions7 as Q
from babyuniverse import collider6l as C
from babyuniverse import pipeline6o as P

warnings.simplefilter("ignore")

if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round7.json"))
else:
    L = C.load_lumi(f"{DATA}/round6l_lumi.json")
    out = {}
    out["benchmarks"] = Q.predictions_table(L, 4000)
    # temperature-mass relation from events with at least two quanta (BH rest frame, generator level)
    tm = {}
    for name, kind, cl, mth in (("U-hayward", "hayward", 0.6, None), ("GR", "gr", 0.6, 11.5)):
        g = P.Generator(L, 100.0, 6, 8.0, cl, kind, seed=31, M_th_TeV=mth)
        r = g.run(8000)
        M = np.array(r["M_TeV"]); T0 = np.array(r["T0_GeV"]); mE = np.array(r["mean_E_GeV"]); N = np.array(r["N"])
        tm[name] = dict(M=M.tolist(), T0=T0.tolist(), meanE=mE.tolist(), N=N.tolist())
    out["temperature_mass"] = tm
    # onset: cross-section above threshold vs threshold, 85 and 100 TeV, with the model-U threshold for c_l values
    Mth = np.linspace(8, 70, 63)
    out["onset"] = {str(rs): [C.sigma_pp(L, rs, 6, 8.0, m, "conservative") for m in Mth] for rs in (85.0, 100.0)}
    out["onset_M"] = Mth.tolist()
    out["backgrounds"] = dict(ST=[(s, Q.qcd_tail_pb(L, 100.0, s)) for s in (2, 5, 10, 15, 20, 30)],
                              MET=[(m, Q.znunu_met_tail_pb(L, 100.0, m)) for m in (0.5, 1, 2, 3, 5)],
                              lhc_dijet_8TeV_pb=Q.qcd_tail_pb(L, 13.0, 8.0))
    # remnant charge scenario: heavy stable charged tracks if 30 % of remnants end up charged after hadronisation
    g = P.Generator(L, 100.0, 6, 8.0, 0.6, "hayward", seed=41)
    out["hscp"] = Q.LabGenerator(g, f_charged=0.3, seed=41).run(3000)
    dump("round7.json", {k: v for k, v in out.items()})

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
ax = ax.ravel()
B = out["benchmarks"]
cols = {"U-hayward-0.6": C1, "U-dymnikova-0.7": C3, "GR-11.5": C2}
for name, col in cols.items():
    r = B[name]["raw"]
    ax[0].hist(np.array(r["ST"]) / 1e3, bins=np.linspace(0, 40, 41), histtype="step", lw=2, color=col, density=True, label=name)
ax[0].set_yscale("log"); ax[0].set_xlabel("S_T, scalar sum of visible pT (TeV)"); ax[0].set_ylabel("fraction of events per TeV")
ax[0].set_title("(a) visible activity: GR holes put ~M into many objects; model-U holes near\nthreshold are nearly invisible (remnant carries ~95 %)", fontsize=8)
ax[0].legend(frameon=False, fontsize=7)
for name, col in cols.items():
    r = B[name]["raw"]
    f = np.array(r["met_frac"]); st = np.array(r["ST"])
    ax[1].hist(f[st > 0], bins=np.linspace(0, 1, 26), histtype="step", lw=2, color=col, density=True, label=name)
ax[1].set_xlabel("MET / S_T"); ax[1].set_ylabel("fraction of events")
ax[1].set_title("(b) the model-U signature: a hard initial-state jet balanced by an invisible\nremnant (MET/S_T -> 1); GR holes are balanced (MET/S_T ~ 0.15)", fontsize=8)
ax[1].legend(frameon=False, fontsize=7)
for name, col in cols.items():
    r = B[name]["raw"]
    ax[2].hist(np.minimum(r["N_obj"], 10), bins=np.arange(-0.5, 11.5), histtype="step", lw=2, color=col, density=True, label=name)
ax[2].set_xlabel("number of visible objects (pT > 50 GeV, |eta| < 4)"); ax[2].set_ylabel("fraction of events")
sep = out["benchmarks"]["separation"]
ax[2].set_title(f"(c) multiplicity; model U vs GR separable at 5 sigma with ~{sep['U-hayward-0.6 vs GR-11.5']['events_5sigma']:.0f} events,\n"
                f"Hayward vs Dymnikova tower with ~{sep['U-hayward-0.6 vs U-dymnikova-0.7']['events_5sigma']:.0f}", fontsize=8)
ax[2].legend(frameon=False, fontsize=7)
tm = out["temperature_mass"]
for name, col in (("U-hayward", C1), ("GR", C2)):
    M = np.array(tm[name]["M"]); T = np.array(tm[name]["T0"])
    ax[3].scatter(M, T, s=2, alpha=0.3, color=col, label=f"{name}: initial Hawking temperature")
ax[3].set_xlabel("black-hole mass (TeV)"); ax[3].set_ylabel("Hawking temperature (GeV)")
ax[3].set_yscale("log")
ax[3].set_title("(d) the temperature-mass relation: model U is cold near threshold (T -> 0 at\nM_ext) and rises with mass; GR is hot and falls as M^-1/7", fontsize=8)
ax[3].legend(frameon=False, fontsize=7, markerscale=4)
M = np.array(out["onset_M"])
for rs, col in (("85.0", C1), ("100.0", C2)):
    ax[4].semilogy(M, np.maximum(out["onset"][rs], 1e-12), color=col, label=f"{float(rs):.0f} TeV: sigma(M_BH > M_th)")
for s, b in out["backgrounds"]["ST"]:
    ax[4].plot(s, max(b, 1e-12), "s", color=MUTED)
ax[4].plot([], [], "s", color=MUTED, label="QCD S_T tail (order of magnitude)")
ax[4].set_ylim(1e-7, 1e3)
ax[4].set_xlabel("threshold / S_T (TeV)"); ax[4].set_ylabel("cross-section (pb)")
ax[4].set_title("(e) onset: the rate switches on at M_ext = kappa_n c_l^{n+1} M_D, a point measurement\nof the interior class once M_D is known from monojets", fontsize=8)
ax[4].legend(frameon=False, fontsize=7)
h = out["hscp"]
br = np.array(h["beta_remnant"]); hs = np.array(h["hscp"]) > 0
ax[5].hist(br[np.array(h["remnant"]) > 0], bins=np.linspace(0, 1, 26), histtype="step", lw=2, color=C4, label="all remnants")
ax[5].hist(br[hs], bins=np.linspace(0, 1, 26), histtype="stepfilled", alpha=0.4, color=C4, label="charged (30 % scenario), |eta| < 2.5")
ax[5].set_xlabel("remnant velocity beta"); ax[5].set_ylabel("events")
ax[5].set_title("(f) remnants are slow (beta ~ 0.4) and heavy (~M_ext): if some are charged after\nhadronisation, they leave slow heavy tracks (time of flight, dE/dx)", fontsize=8)
ax[5].legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(f"{FIG}/round7.png", dpi=150, facecolor="#fcfcfb")
print("done")
