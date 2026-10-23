"""Round 6: the mechanism atlas -- budget theorem numbers, charged-exterior route, stabilised seeds, genesis scaling."""
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse import mechanisms as Mc, tunnelling as T
from babyuniverse.model import ThinWallModel

warnings.simplefilter("ignore")
out = {}

# ------------------------------------------------ G1: the classical budget
out["budget"] = dict(null_force_N=Mc.nec_force_SI(), **{f"C0={c}": Mc.genesis_budget(c) for c in (0.0, 0.5, 0.9)})
out["stress_scale"] = {f"L={L:g}m": Mc.stress_scale(L) for L in (1.0, 1e-2, 1e-6, 1e-15, 100 * Mc.L_P_SI)}
out["casimir"] = {f"a={a:g}m": Mc.casimir_shortfall(a) for a in (1e-6, 1e-8, 1e-9)}

# ------------------------------------------------ G2/G3: thin walls of any matter
gs = Mc.GeneralShell(lambda R: 1 - R**2, lambda R: 1 - 2 * 0.05 / R, lambda R: 0.5 * (1 - (0.15 / R) ** 2))
out["two_component_wall"] = dict(nec_at_R0p3=gs.nec_on_wall(0.3), flips=gs.bridge_crossing_report(0.15 * 1.0001, 3.0))
out["upward_nucleation"] = [dict(kappa=k, **Mc.upward_nucleation_exponent(1.0, -k), B_positive_tension=T.B_canonical(ThinWallModel(1.0, k, 0.0)))
                            for k in (0.3, 1.0, 3.0)]
out["vacuum_survival_exponent"] = Mc.vacuum_survival_exponent()

# ------------------------------------------------ G4: charged exteriors
rows = []
for k in (0.3, 1.0):
    for M in (0.05, 0.1, 0.2):
        B0 = T.B_canonical(ThinWallModel(1.0, k, M))
        for qf in (0.3, 0.6, 0.9, 0.99):
            rows.append(dict(kappa=k, M=M, Q_over_M=qf, B=Mc.ChargedBubble(1.0, k, M, qf * M).B_lab(), B_Q0=B0))
out["charged_B"] = rows
out["charged_threshold"] = [dict(kappa=k, Q_over_M=q, M_cr=Mc.critical_mass(1.0, k, q, hi=1.5), M_cr_FGG=ThinWallModel(1.0, k, 0.1).M_critical)
                            for k, q in ((0.3, 0.01), (0.3, 0.5), (0.3, 0.9), (1.0, 0.01), (1.0, 0.5), (1.0, 0.7))]
out["charged_route_example"] = Mc.ChargedBubble(1.0, 1.0, 0.5, 0.25).classical_route()
out["threshold_physical"] = {f"{E:.0e}GeV": Mc.classical_threshold_physical(E) for E in (1.22e19, 1e17, 1e16, 1e10, 1e3)}

# ------------------------------------------------ G5, G6
out["frontier_seeds"] = Mc.frontier_seed_window()
out["genesis"] = [dict(L_m=L, h_m2=h, **Mc.genesis_timeline(L, h, 1e11)) for L, h in ((1.0, 1e-4), (1e-2, 1e-8), (1e-6, 1e-16))]
dump("round6_atlas.json", out)

# ------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0})
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.1))

L = np.logspace(np.log10(Mc.L_P_SI), 1, 200)
need = Mc.C_SI**4 / (8 * np.pi * Mc.G_SI * L**2)
casimir = np.pi**2 * Mc.HBAR_SI * Mc.C_SI / (720 * L**4)
ax[0].loglog(L, need, color=C2, label=r"required NEC violation $c^4/8\pi G L^2$ (Theorem G1)")
ax[0].loglog(L, casimir, color=C1, label=r"Casimir $|\rho|=\pi^2\hbar c/720L^4$ (QEI-limited)")
ax[0].axhline(1e11, color=MUTED, lw=1, ls="--")
ax[0].text(3e-30, 3e11, "strongest materials (positive stress)", fontsize=7.5, color=MUTED)
ax[0].axvline(Mc.L_P_SI, color=MUTED, lw=0.8)
ax[0].text(2.2e-35, 1e40, r"$\ell_P$", fontsize=8, color=MUTED)
ax[0].set_xlabel("size of the region $L$ (m)"); ax[0].set_ylabel("stress (Pa)")
ax[0].set_title("Classical creation: the NEC budget", loc="left", fontsize=10)
ax[0].legend(fontsize=7.5, frameon=False, loc="upper right")

E = np.logspace(2, np.log10(1.22e19), 200)
Mcr = np.array([Mc.classical_threshold_physical(e)["M_cr_kg"] for e in E])
ax[1].loglog(E, Mcr, color=C3, label=r"charged route: $M\geq M_{\rm cr}\approx0.16\,m_P^3/\sqrt{\rho}$ (beyond a Cauchy horizon)")
ax[1].loglog([2.7e17], [2.18e-8], "o", ms=9, mfc="none", mec=C2, mew=2)
ax[1].annotate(r"tunnelling frontier, $P\leq10^{-23}$", xy=(2.7e17, 2.18e-8), xytext=(1e6, 1e-12), fontsize=8, color=C2,
               arrowprops=dict(arrowstyle="->", color=C2))
ax[1].axhline(1.0, color=MUTED, lw=1, ls="--")
ax[1].text(2e2, 3.0, "1 kg", fontsize=7.5, color=MUTED)
ax[1].axhline(5.97e24, color=MUTED, lw=1, ls=":")
ax[1].text(2e2, 1.5e25, "Earth mass", fontsize=7.5, color=MUTED)
ax[1].set_xlabel(r"false-vacuum energy scale $\rho^{1/4}$ (GeV)"); ax[1].set_ylabel("required mass of the prepared object (kg)")
ax[1].set_title("Routes without exotic matter", loc="left", fontsize=10)
ax[1].legend(fontsize=7.5, frameon=False, loc="upper right")
ax[1].set_ylim(1e-14, 1e32)

Rs = np.linspace(80, 420, 120)
from babyuniverse.bounds import FRONTIER_CONST
life = [Mc.seed_lifetime(FRONTIER_CONST * np.sqrt(10 * R), R)["lifetime_s"] for R in Rs]
ax[2].semilogy(Rs, life, color=C4, label=r"stabilised seed at the frontier, $\tau=e^{B}/\omega$")
ax[2].axhspan(1e-9, 3e7, color="#f1f0ec", zorder=0)
ax[2].text(90, 1e3, "laboratory window\n(1 ns - 1 yr)", fontsize=8, color=MUTED)
ax[2].text(215, 1e-19, "an unstabilised bubble gets one attempt:\n$P\\leq e^{-53}$ at $R=100\\,\\ell_P$", fontsize=7.5, color=C2)
ax[2].set_xlabel(r"seed radius $R$ ($\ell_P$), $B_{\rm CdL}=10$, $M=m_P$"); ax[2].set_ylabel("lifetime until conversion (s)")
ax[2].set_title("Stabilised seeds: suppression is not the bottleneck", loc="left", fontsize=10)
ax[2].legend(fontsize=7.5, frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "round6_atlas.png"), dpi=150)
print("budget N:", out["budget"]["null_force_N"])
print("thresholds:", {k: round(v["M_cr"], 4) for k, v in zip(range(6), out["charged_threshold"])})
print("seeds:", [(r["R_lP"], round(r["B"], 1), f"{r['lifetime_s']:.2e}") for r in out["frontier_seeds"]])
