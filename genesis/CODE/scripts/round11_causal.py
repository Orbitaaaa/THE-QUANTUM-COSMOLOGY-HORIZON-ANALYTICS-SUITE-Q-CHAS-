"""Round 11 (G11-04): schematic conformal diagrams of the three continuation cases -> FIGURES/round11_causal.png.
(a) unperturbed model U (Round 6m: daughter IV u I' beyond r_-); (b) generic perturbed tower in D = 5: weak null singularity on the Cauchy
horizon (G11-01), continuation C^0 only and not fixed by the field equations; (c) inner-extremal tower (G11-02): rays cross r_- at finite v
with bounded curvature in the evolved range.  Schematic: not computed conformal compactifications."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _common import FIG

INK, C1, C2, C3, MUTED = "#0b0b0b", "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"


def diamond(ax, cx, cy, s, **kw):
    ax.plot([cx - s, cx, cx + s, cx, cx - s], [cy, cy + s, cy, cy - s, cy], **kw)


def base(ax, title):
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.25, 2.6)
    ax.set_title(title, fontsize=9)
    # parent exterior (region I) and black-hole region II
    ax.plot([0, 1, 0], [-1, 0, 1], color=INK, lw=1.2)                       # I+ / I- of the parent
    ax.text(0.62, 0.55, "I+", fontsize=8); ax.text(0.62, -0.62, "I-", fontsize=8)
    ax.text(0.45, -0.05, "I (parent)", fontsize=8)
    ax.plot([0, -1], [1, 0], color=INK, lw=0.8, ls=":")
    ax.fill([0, 1, 0, -1], [-1, 0, 1, 0], color="#f3f2ee", zorder=0)
    ax.plot([-0.5, 0.0], [0.5, 1.0], color=MUTED, lw=0.6)
    ax.text(-0.25, 0.35, "II (trapped)", fontsize=8, ha="center")
    # collapsing matter / pulse
    ax.annotate("", xy=(-0.2, 0.8), xytext=(0.7, -0.1), arrowprops=dict(arrowstyle="->", color=C1, lw=1.0))
    ax.text(0.35, 0.45, "pulse", color=C1, fontsize=7, rotation=-45)


fig, ax = plt.subplots(1, 3, figsize=(15, 5.6))
# (a)
base(ax[0], "(a) unperturbed model U (Round 6m)\ndaughter beyond r_-, outside J-(I+) of the parent")
ax[0].plot([0, -1, 0, 1, 0], [1, 2, 3, 2, 1], color=INK, lw=0.8)
ax[0].plot([0, -0.6], [1, 1.6], color=C3, lw=2); ax[0].plot([0, 0.6], [1, 1.6], color=C3, lw=2)
ax[0].text(-0.75, 1.35, "r = r_-", color=C3, fontsize=8); ax[0].text(0.35, 1.25, "r = r_-", color=C3, fontsize=8)
ax[0].text(0, 1.75, "III/IV regular core\n(bounce)", ha="center", fontsize=8)
ax[0].text(0, 2.3, "I' daughter\n(own null infinity)", ha="center", fontsize=8, color=C3)
# (b)
base(ax[1], "(b) generic tower, D = 5 (G11-01): weak null singularity\nTipler-weak, r finite; extension C^0 only -> not determined")
t = np.linspace(0, 1, 200)
ax[1].plot(-0.6 * t + 0.02 * np.sin(60 * t), 1 + 0.6 * t + 0.02 * np.sin(60 * t), color=C2, lw=1.6)
ax[1].plot(0.6 * t - 0.02 * np.sin(60 * t), 1 + 0.6 * t + 0.02 * np.sin(60 * t), color=C2, lw=1.6)
ax[1].text(-1.05, 1.62, "CH: mass inflation\nK x l^4 > 1e10", color=C2, fontsize=7)
ax[1].plot([-0.6, 0, 0.6], [1.6, 2.2, 1.6], color=MUTED, lw=0.8, ls="--")
ax[1].text(0, 1.75, "?", fontsize=22, ha="center", color=MUTED)
ax[1].text(0, 2.3, "daughter not fixed by\nthe classical equations", ha="center", fontsize=8, color=MUTED)
# (c)
base(ax[2], "(c) inner-extremal tower (G11-02): kappa_- = 0\nrays cross r_- at v ~ 27 l, K <= 106 / l^4 (v <= 36, converged)")
ax[2].plot([0, -1, 0, 1, 0], [1, 2, 3, 2, 1], color=INK, lw=0.8)
ax[2].plot([0, -0.6], [1, 1.6], color=C3, lw=2, ls="-."); ax[2].plot([0, 0.6], [1, 1.6], color=C3, lw=2, ls="-.")
ax[2].annotate("", xy=(-0.05, 1.9), xytext=(0.5, 0.6), arrowprops=dict(arrowstyle="->", color=C1, lw=1.0))
ax[2].text(0, 2.05, "daughter reached by\nthe field equations\n(tuned mass; detuning\nand semiclassical caveats)", ha="center", fontsize=7.5, color=C3)
fig.tight_layout()
fig.savefig(f"{FIG}/round11_causal.png", dpi=140)
print("wrote FIGURES/round11_causal.png")
