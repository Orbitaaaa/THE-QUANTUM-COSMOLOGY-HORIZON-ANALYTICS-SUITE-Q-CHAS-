"""Tunnelling analysis: canonical exponent, branch prescriptions, M -> 0
limit, first law, Smarr identity, Euclidean wall geometry, multiplicity."""
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse import tunnelling as T
from babyuniverse.branches import bgg_regime

warnings.simplefilter("ignore")

# ----------------------------------------------------- B(M) for several gamma
out = {}
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
for gam, col in ((-0.5, "C0"), (0.0, "C1"), (0.6, "C3")):
    base = model_from_dimensionless(gam, 10.0)       # fixes chi, kappa
    chi, kap = base.chi, base.kappa
    Mcr = base.M_critical
    Ms = Mcr * np.geomspace(1e-4, 0.999, 60)
    Ban, Bpr, Bt, regime = [], [], [], []
    for M in Ms:
        m = ThinWallModel(chi, kap, M)
        Ban.append(T.B_canonical(m, "analytic"))
        Bpr.append(T.B_canonical(m, "principal"))
        Bt.append(Ban[-1] - T.endpoint_area_term(m))
        regime.append(bgg_regime(m)[0])
    m0 = ThinWallModel(chi, kap, 0.0)
    d0 = T.instanton_M0(m0)
    out[f"gamma={gam}"] = dict(chi=chi, kappa=kap, M_critical=Mcr, M=Ms, B_analytic=Ban,
                               B_principal=Bpr, B_tilde=Bt, regime=regime,
                               M0=dict(S_E=d0["S_E"], morse=d0["morse"], B_closed=d0["B_closed"],
                                       B_numeric=T.B_canonical(m0)))
    axs[0].plot(Ms / Mcr, Ban, color=col, lw=1.8, label=fr"analytic branch, $\gamma={gam}$")
    axs[0].plot(Ms / Mcr, Bpr, color=col, lw=1.0, ls="--", label=fr"principal branch, $\gamma={gam}$")
    axs[0].plot([0], [d0["B_closed"]], "o", color=col, ms=5)
    regs = np.array(regime)
    iB = np.where(regs == "B")[0]
    if len(iB):
        axs[0].axvspan(Ms[iB[0]] / Mcr, 1.0, color=col, alpha=0.06)
axs[0].axhline(0, color="k", lw=0.6)
axs[0].set_xscale("log")
axs[0].set_xlabel(r"$M/M_{\rm cr}$")
axs[0].set_ylabel(r"$B$  (units $G=\hbar=1$, $\kappa=1$)")
axs[0].set_title("canonical WKB exponent; dots: closed form at M = 0", fontsize=9)
axs[0].legend(fontsize=6.5)

# Smarr residuals over a grid
rows = []
for gam in np.linspace(-0.85, 0.85, 9):
    gmin, g1 = ThinWallModel.g_min(gam), 2 - 2 * gam
    for Efac in (1.02, 1.3, 2.0, 5.0, 20.0):
        E = Efac * g1 if Efac > 1.05 else 0.5 * (gmin + g1)
        m = model_from_dimensionless(gam, E)
        s = T.smarr_check(m)
        rep = T.embedding_report(m)
        rows.append(dict(gamma=gam, E=E, regime=bgg_regime(m)[0], M=m.M, chi=m.chi, kappa=m.kappa,
                         B=s["B"], B_principal=T.B_canonical(m, "principal"), B_tilde=s["B_tilde"],
                         covariant=s["covariant"], residual=s["residual"],
                         endpoint=s["endpoint"], **{k: v for k, v in s["pieces"].items()},
                         eps_out=rep["eps_out"], eps_in=rep["eps_in"],
                         exterior_embedded=rep.get("exterior_embedded"),
                         interior_embedded=rep.get("interior_embedded"),
                         theta_D=rep.get("theta_D"), theta_max=rep.get("theta_max")))
res = np.array([abs(r["residual"]) / max(1.0, abs(r["B"])) for r in rows])
print("Smarr identity: max relative residual over", len(rows), "models =", res.max())
dump("smarr_identity_grid.json", rows)
fl = [T.first_law_check(model_from_dimensionless(g, Ef * (2 - 2 * g))) for g, Ef in
      ((-0.5, 1.6), (0.0, 3.0), (0.6, 1.5), (0.3, 3.0))]
dump("first_law_checks.json", fl)
for f in fl:
    print("first law:", {k: round(float(v), 8) for k, v in f.items()})

# M = 0 closed forms over a grid (incl. de Sitter exteriors)
m0rows = []
for chi in (0.5, 1.0, 2.0):
    for kap in (0.2, 0.5, 1.0, 1.5, 3.0):
        for eta in (0.0, 0.3 * chi, 0.8 * chi, 1.3 * chi):
            m = ThinWallModel(chi, kap, 0.0, eta)
            d = T.instanton_M0(m)
            Bn = T.B_canonical(m)
            m0rows.append(dict(chi=chi, kappa=kap, eta=eta, p=d["p"], q=d["q"], S_E=d["S_E"],
                               morse=d["morse"], B_closed=d["B_closed"], B_numeric=Bn,
                               B_CdL_down=d["B_CdL_down"], B_up=d["B_up"],
                               err_closed=abs(Bn - d["B_closed"]), err_morse=abs(Bn - d["S_E"] - d["morse"])))
print("M=0: max |B_num - B_closed| =", max(r["err_closed"] for r in m0rows),
      " max |B_num - (S_E + Morse)| =", max(r["err_morse"] for r in m0rows))
dump("M0_closed_form_grid.json", m0rows)

# Smarr bar chart for one model
m = model_from_dimensionless(-0.5, 1.6 * 3.0)
s = T.smarr_check(m)
p = s["pieces"]
labels = [r"$-\frac{\Lambda_{in}}{8\pi G}V_{in}$", r"$-\frac{\sigma}{2}A_W$", r"$-\frac{M}{2}\Delta t_W$",
          r"$\sum$", r"$B-\frac{1}{4G}\sum\pm A$"]
vals = [p["S_in"], p["S_wall"], p["S_ext"], s["covariant"], s["B_tilde"]]
axs[1].bar(range(5), vals, color=["#f5cba7", "#c0392b", "#aed6f1", "#7f8c8d", "#2c3e50"])
axs[1].set_xticks(range(5))
axs[1].set_xticklabels(labels, fontsize=9)
axs[1].axhline(0, color="k", lw=0.6)
axs[1].set_title(f"Smarr identity, gamma=-0.5, M={m.M:.3f}: residual {s['residual']:.1e}", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "tunnelling_exponent.png"), dpi=150)
plt.close(fig)
dump("B_of_M.json", out)

# ------------------------------------------------ Euclidean wall geometry figure
fig = plt.figure(figsize=(13, 4.4))
for k, (gam, Ef, title) in enumerate(((-0.5, 1.6, r"$\gamma=-0.5$ (regime A)"),
                                      (0.6, 1.5, r"$\gamma=0.6$ (regime A)"),
                                      (-0.5, 0.5, "regime B (hidden bounded branch)"))):
    if Ef > 1:
        m = model_from_dimensionless(gam, Ef * (2 - 2 * gam))
    else:
        g1, gmin = 2 - 2 * gam, ThinWallModel.g_min(gam)
        m = model_from_dimensionless(gam, 0.5 * (g1 + gmin))
    ew = T.euclidean_wall(m, n=500)
    ax = fig.add_subplot(1, 3, k + 1, projection="polar")
    rh = 2 * m.G * m.M
    rr = np.sqrt(ew.R / rh - 1)
    eo = np.sign(m.A_out(ew.R))
    ax.plot(ew.theta_out[eo > 0], rr[eo > 0], color="#c0392b", lw=2, label=r"$\epsilon_{out}=+1$")
    ax.plot(ew.theta_out[eo < 0], rr[eo < 0], color="#8e44ad", lw=2, label=r"$\epsilon_{out}=-1$")
    ax.plot([0, 0], [rr[0], rr.max() * 1.6], color="k", lw=1)
    thD = ew.theta_out[-1]
    if eo[-1] < 0:
        ax.plot([thD, thD], [0, rr[-1]], color="#2471a3", lw=1)
        ax.plot([thD + np.pi] * 2, [0, rr.max() * 1.6], color="#2471a3", lw=1)
    ax.plot([0], [0], "ko", ms=3)
    ax.set_title(title + f"\n theta_D={thD:.2f}, theta_max={ew.theta_out.max():.2f}", fontsize=8)
    ax.set_rticks([])
    ax.legend(fontsize=6, loc="lower left")
fig.suptitle("Euclidean wall in the exterior (r, t_E) cigar: angle = t_E/4GM, radius ~ sqrt(r/2GM - 1); "
             "black: initial slice, blue: final slice through the horizon", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "euclidean_wall_exterior.png"), dpi=150)
plt.close(fig)

# multiplicity map
m = model_from_dimensionless(-0.5, 1.6 * 3.0)
GX, GY, mult, (xs, ys) = T.winding_multiplicity(m, ngrid=201, rho_max=1.6)
fig, ax = plt.subplots(figsize=(5.2, 5))
im = ax.pcolormesh(GX, GY, mult, cmap="RdBu_r", vmin=-2, vmax=2, shading="auto")
ax.plot(xs, ys, "k-", lw=0.8)
ax.set_aspect("equal")
ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-1.6, 1.6)
plt.colorbar(im, ax=ax, label="covering multiplicity of the exterior region")
ax.set_title("FGG exterior Euclidean region: multiplicity\n(values outside {0,1} = not a manifold region)",
             fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "pseudomanifold_multiplicity.png"), dpi=150)
plt.close(fig)
vals, counts = np.unique(mult, return_counts=True)
dump("pseudomanifold_multiplicity.json", dict(values=vals, counts=counts,
                                              params=m.describe(), report=T.embedding_report(m)))
print("multiplicity values:", dict(zip(vals.tolist(), counts.tolist())))

# ------------------------------------- laboratory regime A: range of B / B(0)
gam_grid = np.linspace(-0.95, 0.95, 39)
labrows = []
for gam in gam_grid:
    g1 = 2 - 2 * gam
    mA = model_from_dimensionless(gam, g1 * (1 + 1e-7))     # A/B boundary: R1 = 2GM
    B0 = T.B_canonical(ThinWallModel(mA.chi, mA.kappa, 0.0))
    Ms = np.geomspace(1e-4, 1.0, 50) * mA.M
    Bs = np.array([T.B_canonical(ThinWallModel(mA.chi, mA.kappa, M)) for M in Ms])
    labrows.append(dict(gamma=gam, M_AB_over_Mcr=mA.M / mA.M_critical, B0=B0, B_at_AB=Bs[-1],
                        min_ratio=float(Bs.min() / B0), max_ratio=float(Bs.max() / B0),
                        argmin_M_over_MAB=float(Ms[np.argmin(Bs)] / mA.M),
                        R1_over_2GM_at_AB=mA.turning_points()[0] / (2 * mA.M)))
dump("laboratory_regime_bounds.json", labrows)
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(gam_grid, [r["min_ratio"] for r in labrows], "C0-", label=r"$\min_{A} B / B(M\to0)$")
ax.plot(gam_grid, [r["max_ratio"] for r in labrows], "C3-", label=r"$\max_{A} B / B(M\to0)$")
ax.axhline(1, color="k", lw=0.6)
ax.set_yscale("log")
ax.set_xlabel(r"$\gamma$")
ax.set_title("FGG exponent over laboratory initial data (regime A)\n"
             r"$B(M\to0)=B_{\rm CdL}$ (decay of the false vacuum) when $\epsilon_{in}(R_2)=-1$", fontsize=9)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "laboratory_regime_bounds.png"), dpi=150)
plt.close(fig)
print("regime A: min over gamma of min_ratio =", min(r["min_ratio"] for r in labrows))
