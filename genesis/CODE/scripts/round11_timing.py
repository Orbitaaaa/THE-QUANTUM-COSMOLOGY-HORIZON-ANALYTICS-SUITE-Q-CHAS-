"""Round 11 engineering (G11-17 .. G11-23): GENESIS Timing Demonstrator - design scan, three-sensor method on a simulated bench run,
time-of-flight and beta resolution, rates, and dimensioned layouts of the bench (V1/V2) and beam-test (V3) setups.
Writes DATA/round11_timing.json and FIGURES/round11_timing.png."""
import json
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

from _common import FIG, DATA, dump
from babyuniverse import timing11 as T

warnings.simplefilter("ignore")

if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round11_timing.json"))
else:
    out = {}
    out["design_scan"] = T.design_scan(400)
    ev = T.simulate_bench_run(n_events=2000)
    out["bench_run"] = T.analyze_run(ev)
    out["bench_truth_ps"] = [T.plane_resolution(n=500, gain=g, noise_mV=1.5, fs_GSps=5.0)["sigma_ps"] for g in (20.0, 15.0, 25.0)]
    out["tof_beta"] = T.tof_beta_resolution(35.0) + T.tof_beta_resolution(45.0)
    out["beam_tof"] = T.beam_tof(sigma_tof_ps=50.0) + T.beam_tof(sigma_tof_ps=64.0)
    out["rates"] = dict(cosmic_per_day_1m=T.cosmic_coincidence_rate(1.3, 1.0), cosmic_per_day_1cm=T.cosmic_coincidence_rate(1.3, 0.01),
                        sr90_per_s_10mm_37kBq=T.sr90_rate(37.0, 10.0))
    out["jitter_check"] = dict(analytic_ps=T.jitter_estimate(1.5, 76.0, 0.5), note="sigma_n / (dV/dt) with 10-90 % rise ~0.5 ns, amplitude 76 mV")
    dump("round11_timing.json", out)

C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 1.8, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
ax = ax.ravel()
ds = out["design_scan"]
labs = [r["case"] for r in ds]; vals = [r["sigma_ps"] for r in ds]
cols = [C4 if "Landau" in l else (C1 if "gain" in l else (C2 if "noise" in l else (C3 if "bandwidth" in l else MUTED))) for l in labs]
ax[0].barh(range(len(labs)), vals, color=cols)
ax[0].axvline(45, color=INK, ls="--", lw=1); ax[0].text(45.5, 0.2, "requirement 45 ps", fontsize=7)
ax[0].set_yticks(range(len(labs))); ax[0].set_yticklabels(labs, fontsize=7); ax[0].set_xlabel("per-plane time resolution (ps)")
ax[0].set_title("(a) G11-19 front-end design scan (calibrated LGAD model, 50 um, 1.3 mm pad):\ngain >= 15, noise <= 3 mV, 1-1.5 GHz, >= 5 GS/s keep sigma ~30-35 ps", fontsize=8)
br = out["bench_run"]
x = np.arange(3)
ax[1].bar(x - 0.2, out["bench_truth_ps"], 0.4, color=MUTED, label="true single-plane resolution (MC)")
ax[1].bar(x + 0.2, br["plane_ps"], 0.4, color=C1, label="recovered by the three-sensor method")
ax[1].set_xticks(x); ax[1].set_xticklabels(["plane A (G=20)", "plane B (G=15)", "plane C (G=25)"])
ax[1].set_ylabel("ps"); ax[1].set_ylim(0, 42); ax[1].legend(frameon=False, fontsize=7, loc="upper center", ncol=2)
ax[1].set_title("(b) G11-20/23 acceptance-test analysis on a simulated Sr-90 stack run (%d events):\npairwise widths %.0f/%.0f/%.0f ps -> planes recovered within ~2 ps (stat. %.1f ps)" % (br["n_selected"], *br["pairwise_ps"], br["stat_unc_ps"]), fontsize=8)
bt = out["beam_tof"][:4]
p = [r["p_GeV"] for r in bt]
ax[2].semilogy(p, [r["sep_p_pi_sigma"] for r in bt], "o-", color=C2, label="p / pi separation (sigma), L = 1 m, 50 ps")
ax[2].semilogy(p, [r["sep_K_pi_sigma"] for r in bt], "s-", color=C1, label="K / pi separation (sigma)")
bt2 = out["beam_tof"][4:]
ax[2].semilogy(p, [r["sep_p_pi_sigma"] for r in bt2], "o:", color=C2, label="p / pi, 64 ps (conservative)")
ax[2].axhline(3, color=INK, lw=0.8)
ax[2].set_xlabel("beam momentum (GeV/c)"); ax[2].set_ylabel("separation (sigma)")
ax[2].set_title("(c) V3 beam test: slow-particle proxies (protons, kaons) separated from pions by time\nof flight over 1 m - the measurement the remnant search needs, at accessible energies", fontsize=8)
ax[2].legend(frameon=False, fontsize=7)
# (d) bench layout (side view)
a3 = ax[3]; a3.set_aspect("equal"); a3.axis("off")
a3.add_patch(Rectangle((0, 0), 160, 90, fill=False, ec=INK, lw=1.5)); a3.text(2, 92, "light-tight, grounded box (~160 x 90 mm inner, schematic)", fontsize=7)
a3.add_patch(Rectangle((70, 70), 20, 8, color=C2)); a3.text(92, 72, "Sr-90 sealed source + collimator", fontsize=7)
for k, y in enumerate((50, 38, 26)):
    a3.add_patch(Rectangle((40, y), 80, 3, color=C3)); a3.add_patch(Rectangle((78.7, y + 3), 2.6, 0.8, color=C4))
    a3.text(122, y, "carrier PCB %s + LGAD (1.3 mm pad)" % "ABC"[k], fontsize=6.5)
a3.add_patch(FancyArrowPatch((80, 70), (80, 22), arrowstyle="-|>", mutation_scale=10, color=C2))
a3.annotate("", xy=(30, 50), xytext=(30, 26), arrowprops=dict(arrowstyle="<->", color=MUTED)); a3.text(5, 37, "12 mm\npitch", fontsize=6.5)
a3.add_patch(Rectangle((40, 6), 80, 8, color="#d9d8d4")); a3.text(42, 8, "SMA to amplifiers -> DRS4 (5 GS/s) or scope; HV filter", fontsize=6.3)
a3.set_xlim(-5, 200); a3.set_ylim(0, 100)
a3.set_title("(d) G11-21 V1/V2 bench layout (side view, dimensions indicative)", fontsize=8)
# (e) beam-test layout (top view)
a4 = ax[4]; a4.axis("off"); a4.set_xlim(-0.1, 1.45); a4.set_ylim(-0.25, 0.3)
a4.plot([-0.05, 1.4], [0, 0], color=MUTED, ls="--"); a4.text(1.3, 0.02, "beam", fontsize=7)
for xpos, lab in ((0.0, "station A: scint. trigger + LGAD + MCP-PMT ref."), (1.0, "station B: LGAD + scint.")):
    a4.add_patch(Rectangle((xpos - 0.02, -0.06), 0.04, 0.12, color=C3)); a4.text(xpos - 0.05, 0.1 if xpos == 0 else -0.16, lab, fontsize=7)
a4.annotate("", xy=(1.0, -0.09), xytext=(0.0, -0.09), arrowprops=dict(arrowstyle="<->", color=INK)); a4.text(0.45, -0.07, "L = 1.000 m (surveyed)", fontsize=7)
a4.text(0.0, 0.22, "rail with surveyed stations; common clock; DRS4/CAEN digitiser; cables matched and calibrated", fontsize=7)
a4.set_title("(e) G11-21 V3 beam-test layout (top view)", fontsize=8)
# (f) beta resolution
tb = out["tof_beta"]
b35 = tb[:5]; b45 = tb[5:]
ax[5].plot([r["beta"] for r in b35], [r["sigma_beta"] for r in b35], "o-", color=C1, label="35 ps per plane")
ax[5].plot([r["beta"] for r in b45], [r["sigma_beta"] for r in b45], "s--", color=C2, label="45 ps per plane")
ax[5].set_xlabel("particle velocity beta"); ax[5].set_ylabel("sigma_beta (L = 1 m)")
ax[5].set_title("(f) beta resolution over a 1 m baseline: sigma_beta ~ beta^2 c sigma_TOF / L;\nbeta ~ 0.3 to ~0.1-0.2 %%; cosmics at 1 m: %.0e/day -> use a beam" % out["rates"]["cosmic_per_day_1m"], fontsize=8)
ax[5].legend(frameon=False, fontsize=7)
fig.subplots_adjust(left=0.1, right=0.98, top=0.93, bottom=0.06, wspace=0.35, hspace=0.4)
fig.savefig(f"{FIG}/round11_timing.png", dpi=150, facecolor="#fcfcfb")
print("done")
