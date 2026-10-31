"""The machine: FCC-hh ring, the detector with the model-U-specific elements, and a model-U event.  Writes FIGURES/machine.png.
Schematic; dimensions follow the FCC-hh reference detector (FCC CDR, Eur. Phys. J. ST 228, 755 (2019)) and the Delphes FCC-hh card
(tracker R = 1.5 m, half-length 5 m, B = 4 T, calorimetry to |eta| = 6); ring after the FCC Feasibility Study (2025, ~91 km)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Wedge, FancyArrowPatch, Arc

from _common import FIG

INK, MUTED, BG = "#1d1d1b", "#6b6a66", "#fcfcfb"
TRK, TIM, ECAL, HCAL, SOL, MUON, FWD = "#9fc5e8", "#8e5bd6", "#7fc8a9", "#f4b183", "#b7b7b7", "#d5a6bd", "#e6c86e"
plt.rcParams.update({"font.size": 9})
fig = plt.figure(figsize=(17, 10.5), facecolor=BG)

# ------------------------------------------------------------------ (a) the ring
ax = fig.add_axes([0.02, 0.50, 0.30, 0.46]); ax.set_aspect("equal"); ax.axis("off")
R_fcc = 90.7 / (2 * np.pi); R_lhc = 26.7 / (2 * np.pi)
t = np.linspace(0, 2 * np.pi, 400)
ax.fill(R_fcc * np.cos(t) * 1.0 + 0, R_fcc * np.sin(t) * 0.92, color="#eef3fb")
ax.plot(R_fcc * np.cos(t), R_fcc * np.sin(t) * 0.92, color="#2a78d6", lw=4, label="FCC-hh tunnel, ~91 km")
ax.plot(R_lhc * np.cos(t) - 7.0, R_lhc * np.sin(t) - 5.5, color=MUTED, lw=2, label="LHC, 26.7 km")
for ang, lab in ((np.pi * 0.5, "IP: model-U experiment"), (np.pi * 1.5, "IP: second experiment")):
    x, y = R_fcc * np.cos(ang), R_fcc * np.sin(ang) * 0.92
    ax.plot(x, y, "o", ms=12, color="#eb6834", mec="w", mew=2)
    ax.text(x, y + (1.6 if y > 0 else -1.9), lab, ha="center", fontsize=8, color=INK)
ax.text(1.5, 2.5, "~14 T Nb3Sn dipoles\nsqrt(s) = 85 TeV (100 TeV option)\n20-30 ab^-1, up to ~1000 collisions / crossing",
        ha="center", fontsize=8.5, color=INK)
ax.text(-7.0, -5.5, "LHC", ha="center", va="center", fontsize=8, color=MUTED)
ax.plot([-R_fcc - 0.5, -R_fcc + 4.5], [-R_fcc * 0.92 - 1.2] * 2, color=INK, lw=2)
ax.text(-R_fcc + 2, -R_fcc * 0.92 - 2.2, "5 km", ha="center", fontsize=8)
ax.legend(loc="lower right", frameon=False, fontsize=8)
ax.set_xlim(-R_fcc - 2, R_fcc + 2); ax.set_ylim(-R_fcc - 3.5, R_fcc + 2.5)
ax.set_title("(a) the collider: a ~91 km ring at CERN (schematic, not to geographic scale)", fontsize=10, loc="left")

# ------------------------------------------------------------------ (b) detector quarter section (r-z)
ax = fig.add_axes([0.34, 0.50, 0.64, 0.46]); ax.set_aspect("equal")
def box(z0, z1, r0, r1, col, lab=None, alpha=1.0):
    for s in (1, -1):
        ax.add_patch(Rectangle((min(s * z0, s * z1), r0), abs(z1 - z0), r1 - r0, color=col, alpha=alpha, lw=0))
    if lab:
        ax.text(0.5 * (z0 + z1), 0.5 * (r0 + r1), lab, ha="center", va="center", fontsize=7.5, color=INK)
box(0, 5.0, 0.03, 1.5, TRK, "tracker  |eta| < 6\nB = 4 T")
box(0, 5.0, 1.5, 1.62, TIM)
box(0, 5.3, 1.7, 2.6, ECAL, "ECAL 10%/sqrtE")
box(0, 5.3, 2.8, 5.0, HCAL, "HCAL 50%/sqrtE + 3%")
box(0, 10.0, 5.0, 5.4, SOL, "solenoid 4 T")
box(0, 10.0, 6.0, 9.5, MUON, "muon system + multi-crossing late-hit trigger (remnants arrive ~0.2 us late)", alpha=0.85)
box(5.4, 7.0, 0.3, 5.0, HCAL, "endcap\ncalo", alpha=0.85)
box(15.5, 17.5, 0.06, 1.6, FWD, "forward calo\n4 < |eta| < 6")
ax.annotate("timing layer, 30-50 ps\n(median remnant: +34 ns at 1.5 m)", xy=(3.5, 1.56), xytext=(6.5, 0.4), fontsize=8, color=TIM,
            arrowprops=dict(arrowstyle="->", color=TIM))
ax.axhline(0, color=INK, lw=1); ax.text(17.8, 0.15, "beam", fontsize=8)
for eta in (1, 2, 4):
    th = 2 * np.arctan(np.exp(-eta))
    ax.plot([0, 18 * np.cos(th)], [0, 18 * np.sin(th)], color=MUTED, lw=0.6, ls=":")
    ax.text(18 * np.cos(th) * (0.97 if eta < 4 else 0.995), min(18 * np.sin(th), 9.6) + 0.15, "eta=%d" % eta, fontsize=7, color=MUTED)
ax.plot(0, 0, "*", ms=12, color="#eb6834")
ax.set_xlim(-18, 18); ax.set_ylim(-0.2, 10.2)
ax.set_xlabel("z (m)"); ax.set_ylabel("r (m)")
ax.set_title("(b) the detector (upper half, r-z; schematic after the FCC-hh reference design): hermetic to |eta| = 6 for MET,\n"
             "a precision timing layer and dE/dx for slow heavy tracks, unprescaled MET > 1 TeV and jet > 1 TeV triggers", fontsize=10, loc="left")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# ------------------------------------------------------------------ (c) transverse event display
ax = fig.add_axes([0.02, 0.03, 0.34, 0.44]); ax.set_aspect("equal"); ax.axis("off")
for r0, r1, col in ((0.03, 1.5, TRK), (1.5, 1.62, TIM), (1.7, 2.6, ECAL), (2.8, 5.0, HCAL), (5.0, 5.4, SOL), (6.0, 9.5, MUON)):
    ax.add_patch(Wedge((0, 0), r1, 0, 360, width=r1 - r0, color=col, alpha=0.55, lw=0))
phi_j = 0.35
# ISR jet: tracks + calorimeter tower
rng = np.random.default_rng(4)
for k in range(14):
    dphi = rng.normal(0, 0.05); c = rng.uniform(-0.25, 0.25)
    r = np.linspace(0, 1.5, 40); ph = phi_j + dphi + c * r / 1.5 * 0.2
    ax.plot(r * np.cos(ph), r * np.sin(ph), color="#2a78d6", lw=0.8)
ax.add_patch(Wedge((0, 0), 5.0, np.degrees(phi_j) - 6, np.degrees(phi_j) + 6, width=3.3, color="#eb6834", alpha=0.85))
ax.text(5.6 * np.cos(phi_j), 5.6 * np.sin(phi_j) + 0.3, "ISR jet\npT = 1.2 TeV", fontsize=8.5, color=INK)
# a soft colour-discharge jet
phi_d = 2.4
ax.add_patch(Wedge((0, 0), 3.6, np.degrees(phi_d) - 4, np.degrees(phi_d) + 4, width=1.9, color="#eb6834", alpha=0.45))
ax.text(4.0 * np.cos(phi_d) - 1.6, 4.0 * np.sin(phi_d), "discharge jet\n~0.5 TeV", fontsize=7.5, color=MUTED)
# MET arrow (invisible remnant)
phi_m = phi_j + np.pi + 0.15
ax.add_patch(FancyArrowPatch((0, 0), (8.6 * np.cos(phi_m), 8.6 * np.sin(phi_m)), arrowstyle="-|>", mutation_scale=22, lw=2.4,
                             color="#8e5bd6", ls="--"))
ax.text(-10.3, -9.9, "MET = 1.3 TeV: the invisible remnant\n(M ~ M_ext, beta ~ 0.02-0.5, median 0.13)", fontsize=8.5, color="#8e5bd6",
        bbox=dict(facecolor=BG, edgecolor="none", alpha=0.85))
# optional charged slow remnant
r = np.linspace(0, 9.2, 200); ph = phi_m + 0.02 * r
ax.plot(r * np.cos(ph), r * np.sin(ph), color=INK, lw=1.2, ls=":")
ax.plot(1.56 * np.cos(phi_m + 0.03), 1.56 * np.sin(phi_m + 0.03), "o", ms=6, color=TIM)
ax.text(-9.6, 8.6, "dotted: the same remnant if charged (alternative scenario):\nlate timing hit (median +34 ns at 1.5 m), dE/dx >~ 10x MIP,\nM = p/(beta gamma); the slowest may stop in the calorimeter",
        fontsize=7.5, color=INK)
ax.set_xlim(-10, 10); ax.set_ylim(-10, 10)
ax.set_title("(c) a model-U event near threshold (transverse view, representative)", fontsize=10, loc="left")

# ------------------------------------------------------------------ (d) what each part measures
ax = fig.add_axes([0.40, 0.03, 0.58, 0.42]); ax.axis("off")
rows = [("calorimetry to |eta| = 6", "MET from the invisible remnant; onset of the black-hole rate at M_ext", "class (regular core) via the onset"),
        ("MET > 1 TeV, jet > 1 TeV triggers", "monojet + invisible remnant topology (MET/S_T, multiplicity)", "separation from GR holes (1-13 selected events)"),
        ("tracker + timing (30-50 ps)", "slow heavy tracks, beta ~ 0.02-0.5; pile-up rejection", "remnant mass M = p/(beta gamma), if charged"),
        ("dE/dx, late-hit muon trigger", "highly ionising, out-of-time particles (> 25 ns)", "confirms a stable remnant"),
        ("low-MET monojet / monophoton", "KK-graviton emission", "fixes n and M_D (low-scale gravity)"),
        ("cannot be measured", "the daughter region itself (outside the parent's causal past)", "inferred only through the class")]
tb = ax.table(cellText=rows, colLabels=["detector element", "what it records", "what it establishes"], loc="upper center", cellLoc="left",
              colWidths=[0.22, 0.44, 0.34])
tb.auto_set_font_size(False); tb.set_fontsize(8.2); tb.scale(1, 2.0)
for (i, j), c in tb.get_celld().items():
    c.set_edgecolor("#d9d8d4")
    if i == 0:
        c.set_facecolor("#efeeea"); c.set_text_props(weight="bold")
    if i == len(rows):
        c.set_facecolor("#fbecec")
ax.set_title("(d) what each part of the machine measures - and the one thing no machine can", fontsize=10, loc="left")
fig.savefig(f"{FIG}/machine.png", dpi=150, facecolor=BG)
print("done")
