"""Round 3, Stage 1: the A/B end of the T5 strip, NG4 at both ends, the independent check of S3, Euclidean SdS."""
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse import bounds as Bd, tunnelling as T, crosscheck as X
from babyuniverse.model import ThinWallModel
from babyuniverse.sds_atlas import SdSExterior, euclidean_section

warnings.simplefilter("ignore")
out = {}

# ------------------------------------------------ T5 strip (regime A, E_AB < E < omega threshold)
out["boundary_certificate"] = Bd.boundary_certificate()
out["boundary_limit"] = {f"{g:+.2f}": Bd.boundary_limit(g) for g in (-0.99, -0.95, -0.5, 0.0, 0.5, 0.9, 0.99)}
gams = np.linspace(-0.98, 0.98, 50)
eps = np.concatenate([np.logspace(-7, -1, 31), np.linspace(0.12, 0.62, 11)])
Jmax, arg = -np.inf, None
for g in gams:
    for e in eps:
        J = Bd.delta_tE_dimless(g, (2 - 2 * g) * (1 + e))
        if J > Jmax:
            Jmax, arg = J, (g, e)
out["strip_scan"] = dict(points=len(gams) * len(eps), max_J=Jmax, at_gamma=arg[0], at_eps=arg[1])
out["lp_certificates"] = {f"{g:+.1f}": {f"{f:.2f}": Bd.strip_lp_certificate(g, (2 - 2 * g) * f) for f in (1.01, 1.02, 1.05, 1.3, 1.6)}
                          for g in (-0.9, -0.5, 0.0, 0.5, 0.9)}

# ------------------------------------------------ NG4 (regime B)
viol, minr, pts = 0, np.inf, 0
for g in np.linspace(-0.98, 0.98, 60):
    b, gmin = 2 - 2 * g, ThinWallModel.g_min(g)
    x = (1 - np.cos(np.linspace(0, np.pi, 82)[1:-1])) / 2
    th = np.array([Bd.theta_D_over_pi(g, gmin + (b - gmin) * xi) for xi in x])
    viol += int((np.diff(th) <= 0).sum())
    minr = min(minr, th.min())
    pts += len(x)
out["ng4"] = dict(points=pts, non_increasing_steps=viol, min_theta_over_pi=minr,
                  ends={f"{g:+.2f}": dict(critical=T.theta_D_critical(g), ab=Bd.ng4_ab_limit(g)) for g in (-0.9, -0.5, 0.0, 0.5, 0.9, 0.99)})

# ------------------------------------------------ S3, independent
s3 = []
for e in (0.1, 0.03, 0.01):
    EAB = 4 * e**2 / (1 + e**2)
    for E in sorted({1.5 * EAB, 10 * EAB, 0.03, 0.3, 1.0, 3.0, 30.0, 100.0}):
        if E > EAB:
            s3.append(X.s3_check(e, E, dps=25))
out["s3_crosscheck"] = s3
out["s3_small_E_coefficient"] = X.s3_check(0.001, 1.5 * 4e-6 / (1 + 1e-6), dps=35)["rel_dev_over_eps2"]

# ------------------------------------------------ Euclidean SdS
sds = []
for M, eta in ((0.05, 0.5), (0.1, 0.95), (0.01, 1.0)):
    s = SdSExterior(M, eta)
    for beta in (2 * np.pi / s.kb, 2 * np.pi / s.kc, 1.0):
        sds.append(dict(M=M, eta=eta, **euclidean_section(s, beta)))
out["sds_euclidean"] = sds
dump("round3_stage1.json", out)

# ------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0})
fig, ax = plt.subplots(1, 3, figsize=(13, 4.0))
ee = np.logspace(-7, np.log10(0.6), 60)
for g, col in ((-0.9, C1), (0.0, C2), (0.9, C3)):
    ax[0].semilogx(ee, [Bd.delta_tE_dimless(g, (2 - 2 * g) * (1 + e)) for e in ee], color=col, label=rf"$\gamma={g:+.1f}$")
    ax[0].plot([ee[0]], [Bd.boundary_limit(g)], "o", ms=7, mfc="none", mec=col)
ax[0].axhline(0, color=MUTED, lw=0.8)
ax[0].set_xlabel(r"$E/E_{AB}-1$"); ax[0].set_ylabel(r"$J\propto\Delta t^W_E$")
ax[0].set_title("T5 strip: J < 0, limit $J_0+\\pi$ (circles)", loc="left", fontsize=10)
ax[0].legend(fontsize=8, frameon=False, loc="lower left")

for g, col in ((-0.9, C1), (0.0, C2), (0.9, C3)):
    b, gmin = 2 - 2 * g, ThinWallModel.g_min(g)
    x = (1 - np.cos(np.linspace(0, np.pi, 60)[1:-1])) / 2
    ax[1].plot(x, [Bd.theta_D_over_pi(g, gmin + (b - gmin) * xi) for xi in x], color=col, label=rf"$\gamma={g:+.1f}$")
    ax[1].plot([0], [T.theta_D_critical(g)], "s", ms=7, mfc="none", mec=col)
    ax[1].plot([1], [Bd.ng4_ab_limit(g)], "o", ms=7, mfc="none", mec=col)
ax[1].axhline(1, color=MUTED, lw=0.8)
ax[1].set_xlabel(r"$(E-g_{\min})/(E_{AB}-g_{\min})$ (regime B)"); ax[1].set_ylabel(r"$|\theta_D|/\pi$")
ax[1].set_title("NG4: both ends proved (markers), interior monotone", loc="left", fontsize=10)
ax[1].legend(fontsize=8, frameon=False, loc="center right")

for e, col, mk in ((0.1, C1, "o"), (0.03, C2, "s"), (0.01, C4, "^")):
    rows = [r for r in s3 if r["eps"] == e]
    ax[2].semilogx([r["E"] for r in rows], [r["rel_dev_over_eps2"] for r in rows], mk, ms=7, mfc="none", mec=col,
                   label=rf"$\kappa/\chi={e}$")
EE = np.logspace(1.3, 2.2, 20)
ax[2].semilogx(EE, -5.0 * EE**-1.5, color=MUTED, ls="--", lw=1.2, label=r"$-5E^{-3/2}$ (large $E$)")
ax[2].axhline(-13 / 18, color=C3, lw=1.2, label=r"$-13/18$")
ax[2].set_xlabel(r"$E=(R_{\rm CdL}/R_M)^2$"); ax[2].set_ylabel(r"$(B_{\rm an}/B_0\Phi-1)/(\kappa/\chi)^2$")
ax[2].set_title("S3, independent check: an $O(\\kappa^2/\\chi^2)$ correction", loc="left", fontsize=10)
ax[2].legend(fontsize=8, frameon=False, loc="center right")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "round3_stage1.png"), dpi=150)
print("strip:", out["strip_scan"])
print("ng4:", {k: v for k, v in out["ng4"].items() if k != "ends"})
print("s3 small-E coefficient:", out["s3_small_E_coefficient"])
