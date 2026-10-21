"""Figures: effective potentials and Penrose diagrams built from integrated
trajectories (never from the potential alone).

Kruskal charts are compactified with U~ = arctan U, V~ = arctan V; in these
coordinates the Schwarzschild singularities UV = 1 are the straight lines
Y = +-pi/4 exactly (arctan U + arctan 1/U = pi/2).  de Sitter is drawn in
global conformal coordinates (psi, T) in [0, pi] x (-pi/2, pi/2).
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .model import ThinWallModel  # noqa: E402

WALL = "#c0392b"
HOR = "#7f8c8d"
SING = "#2c3e50"
SHADE_IN = "#f5cba7"
SHADE_OUT = "#aed6f1"


def plot_potential(model: ThinWallModel, ax, branches=None, Rmax=None):
    tps = model.turning_points()
    hz = model.horizons_out()
    Rmax = Rmax or 2.2 * max(tps + hz + ([1 / model.chi] if model.chi else []))
    R = np.linspace(Rmax * 1e-3, Rmax, 3000)
    V = model.V(R)
    ax.plot(R, V, color="k", lw=1.4, label="V(R)")
    ax.axhline(0, color="k", lw=0.6)
    lo = np.nanpercentile(V, 5)
    ax.set_ylim(max(lo, -3.0), max(1.2, np.nanmax(V) * 1.2 if np.nanmax(V) > 0 else 0.5))
    for h in hz:
        ax.axvline(h, color=HOR, ls="--", lw=0.8)
        ax.text(h, ax.get_ylim()[1] * 0.92, " $r_h$", color=HOR, fontsize=8)
    if model.chi:
        ax.axvline(1 / model.chi, color="#b9770e", ls=":", lw=0.9)
        ax.text(1 / model.chi, ax.get_ylim()[1] * 0.80, r" $1/\chi$", color="#b9770e", fontsize=8)
    for Rf, lab in ((model.flip_radius_out(), r"$\epsilon_{out}$ flip"),
                    (model.flip_radius_in(), r"$\epsilon_{in}$ flip")):
        if Rf and Rf < Rmax:
            ax.axvline(Rf, color="#8e44ad", ls="-.", lw=0.8)
            ax.text(Rf, ax.get_ylim()[0] * 0.85, " " + lab, color="#8e44ad", fontsize=7, rotation=90)
    ax.fill_between(R, V, 0, where=V > 0, color="#f9e79f", alpha=0.6, label="forbidden")
    for t in tps:
        ax.plot([t], [0], "o", color=WALL, ms=4)
    ax.set_xlabel("R")
    ax.set_ylabel("V(R)")
    ax.legend(fontsize=7, loc="lower left")


def _wall_x_of_y(X, Y):
    o = np.argsort(Y)
    Ys, Xs = Y[o], X[o]
    keep = np.concatenate([[True], np.diff(Ys) > 1e-12])
    return Xs[keep], Ys[keep]


def plot_kruskal(traj, ax, title="exterior (Kruskal)"):
    s = np.linspace(-np.pi / 2, np.pi / 2, 400)
    # boundaries: horizons U=0, V=0 ; singularities Y = +-pi/4 ; scri
    ax.plot([-np.pi / 4, np.pi / 4], [-np.pi / 4, np.pi / 4], color=HOR, lw=0.8, ls="--")
    ax.plot([-np.pi / 4, np.pi / 4], [np.pi / 4, -np.pi / 4], color=HOR, lw=0.8, ls="--")
    ax.plot([-np.pi / 4, np.pi / 4], [np.pi / 4, np.pi / 4], color=SING, lw=2.2)
    ax.plot([-np.pi / 4, np.pi / 4], [-np.pi / 4, -np.pi / 4], color=SING, lw=2.2)
    for sx in (1, -1):
        ax.plot([sx * np.pi / 4, sx * np.pi / 2], [np.pi / 4, 0], color="k", lw=0.8)
        ax.plot([sx * np.pi / 4, sx * np.pi / 2], [-np.pi / 4, 0], color="k", lw=0.8)
    # shading of the exterior side: X > X_w(Y)
    Xw, Yw = _wall_x_of_y(traj.X, traj.Y)
    gx, gy = np.meshgrid(np.linspace(-np.pi / 2, np.pi / 2, 500), np.linspace(-np.pi / 4, np.pi / 4, 300))
    Ut, Vt = gy - gx, gy + gx
    inside = (np.abs(Ut) < np.pi / 2) & (np.abs(Vt) < np.pi / 2) & (np.abs(gy) < np.pi / 4)
    xw = np.interp(gy, Yw, Xw)
    ext = inside & (gx > xw)
    ax.contourf(gx, gy, ext.astype(float), levels=[0.5, 1.5], colors=[SHADE_OUT], alpha=0.8)
    ax.plot(traj.X, traj.Y, color=WALL, lw=2)
    ax.text(np.pi / 2 * 0.62, 0, "I", fontsize=11, ha="center")
    ax.text(-np.pi / 2 * 0.62, 0, "III", fontsize=11, ha="center")
    ax.text(0, np.pi / 4 * 0.6, "II", fontsize=11, ha="center")
    ax.text(0, -np.pi / 4 * 0.6, "IV", fontsize=11, ha="center")
    ax.set_xlim(-np.pi / 2 - 0.05, np.pi / 2 + 0.05)
    ax.set_ylim(-np.pi / 2 * 0.6, np.pi / 2 * 0.6)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=9)
    del s


def plot_ds(traj, ax, interior=True, title="interior (de Sitter)"):
    ax.plot([0, np.pi, np.pi, 0, 0], [-np.pi / 2, -np.pi / 2, np.pi / 2, np.pi / 2, -np.pi / 2],
            color="k", lw=0.8)
    ax.plot([0, np.pi], [-np.pi / 2, np.pi / 2], color=HOR, lw=0.8, ls="--")
    ax.plot([0, np.pi], [np.pi / 2, -np.pi / 2], color=HOR, lw=0.8, ls="--")
    Xw, Yw = _wall_x_of_y(traj.X, traj.Y)
    gx, gy = np.meshgrid(np.linspace(0, np.pi, 400), np.linspace(-np.pi / 2, np.pi / 2, 400))
    xw = np.interp(gy, Yw, Xw)
    reg = (gx < xw) if interior else (gx > xw)
    # before the wall exists / after it ends at a pole the extension handles it
    ax.contourf(gx, gy, reg.astype(float), levels=[0.5, 1.5],
                colors=[SHADE_IN if interior else SHADE_OUT], alpha=0.8)
    ax.plot(traj.X, traj.Y, color=WALL, lw=2)
    ax.text(0.25, 0, "N", fontsize=11)
    ax.text(np.pi - 0.45, 0, "S", fontsize=11)
    ax.text(np.pi / 2 - 0.08, 0.95, "F", fontsize=11)
    ax.text(np.pi / 2 - 0.08, -1.1, "P", fontsize=11)
    ax.set_xlim(-0.05, np.pi + 0.05)
    ax.set_ylim(-np.pi / 2 - 0.05, np.pi / 2 + 0.05)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=9)


def branch_figure(model, branch, trajs, path, suptitle=""):
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), gridspec_kw=dict(width_ratios=[1.2, 1, 1.25]))
    plot_potential(model, axs[0])
    lo = min(trajs["out"].R.min(), trajs["in"].R.min())
    hi = max(trajs["out"].R.max(), trajs["in"].R.max())
    axs[0].axvspan(lo, min(hi, axs[0].get_xlim()[1]), color=WALL, alpha=0.08)
    axs[0].set_title(f"branch: {branch.kind}", fontsize=9)
    plot_ds(trajs["in"], axs[1], interior=True)
    if trajs["out"].chart == "Kruskal":
        plot_kruskal(trajs["out"], axs[2])
    elif trajs["out"].chart == "dS":
        plot_ds(trajs["out"], axs[2], interior=False, title="exterior (de Sitter)")
    fig.suptitle(suptitle, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
