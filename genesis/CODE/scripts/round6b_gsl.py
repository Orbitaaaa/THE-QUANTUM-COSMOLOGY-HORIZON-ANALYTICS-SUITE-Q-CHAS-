"""Round 6b: genesis and the generalized second law -- dichotomy, QNEC cost, throat flips, growth, Cauchy-horizon route."""
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse import gsl as S

warnings.simplefilter("ignore")
out = {}

# Theorem GSL-3 verdicts
cases = {"eternal dS daughter, globally hyperbolic": dict(observer_complete=True, final_phase="de_sitter", Lambda_final=1.0),
         "disconnected daughter (neck closes)": dict(observer_complete=True, final_phase="flat", neck_persists=False),
         "transient connected daughter": dict(observer_complete=True, final_phase="flat"),
         "crunching daughter": dict(observer_complete=False, final_phase="crunch"),
         "daughter beyond a Cauchy horizon": dict(observer_complete=True, globally_hyperbolic=False, final_phase="de_sitter",
                                                  Lambda_final=1.0),
         "dS parent, Lambda_d < Lambda_p": dict(observer_complete=True, final_phase="de_sitter", Lambda_final=0.5, Lambda_parent=1.0),
         "dS parent, Lambda_d > Lambda_p": dict(observer_complete=True, final_phase="de_sitter", Lambda_final=2.0, Lambda_parent=1.0)}
out["dichotomy"] = {k: S.wall_dichotomy(**v) for k, v in cases.items()}

# Theorem GSL-4: the cost of a flip for a 1 m laboratory sphere, and its throat discount
rP = 1.0 / 1.616255e-35
out["qnec_cost_1m_flat"] = dict(S_BH=float(np.pi * rP**2), cost_per_planck_length=S.qnec_cost(rP, 2 / rP),
                                cost_over_S_BH_times_2_over_r=1.0)
out["antitrapped_floor"] = {f"R={R:g} m": S.antitrapped_floor(R) for R in (1e-13, 1e-3, 1.0, 1e3)}

# Proposition GSL-6/7: flip window in a near-extremal throat (Planck-unit model and three-way check)
Q, T = 1e3, 1e-10
th = S.NearExtremalThroat(Q, T)
r0 = th.r_at_depth(Q)
f0 = th.f(r0)
checks = []
for c in (5.0, 20.0, 40.0, 100.0):
    cf = S.cft_flip(c, r0, -np.sqrt(f0 / 2), np.sqrt(2) / Q, w=0.002 * Q, lam_max=10 * Q)
    vd = S.vaidya_flip(th.rp, th.rm, c, np.sqrt(f0) / Q, 0.002 * Q / np.sqrt(f0), [r0])
    checks.append(dict(c=c, predicted=th.flip_condition(Q, c)["flips"], cft=cf["flips"], vaidya=bool(len(vd["antitrapped"])),
                       kick_over_need=cf["kick_bound"] / np.sqrt(f0 / 2), dm_over_budget=vd["delta_m_peak"] / vd["budget"][0],
                       duration=cf["duration"], duration_bound=cf["duration_bound"], max_dTheta=cf["max_dTheta"]))
out["throat_model"] = dict(Q=Q, T=T, T_over_Tbrk=T / th.T_brk, sqrt_f_at_rho_Q=float(np.sqrt(f0)), checks=checks)
out["throat_examples_SI"] = {f"Q={Qm:g} m, c={c:g}": S.throat_window_SI(Qm, c) for Qm, c in ((1e-13, 1e21), (1e-13, 1e3), (1e-3, 1e3))}
out["magnetic_charge"] = {f"q={q:g}": S.magnetic_charge_radius(q) for q in (1e20, 1.06e22)}

# Proposition GSL-8: growth thresholds and the electroweak corona
out["growth"] = dict(bud_threshold=S.bud_threshold(), switch_on_conserved=S.switch_on_threshold("conserved"), decompactification=0.25,
                     nariai_M_sqrtLambda=S.nariai_parameter(0.0, 1.0, 1e-6) / 2,
                     bud_density={f"Q={Qm:g} m": S.bud_density(Qm) for Qm in (1e-13, 1e-3, 1.0)})
out["electroweak_corona"] = {k: S.electroweak_corona_bound(k) for k in ("higgs", "W")}
out["hiding_depth"] = {f"dM/Q={d:g}": S.hiding_depth(d) for d in (1 / 32, 1e-2, 1e-4, 1e-8)}

# Proposition GSL-9: Cauchy-horizon route near extremality
qs = (0.9, 0.97, 0.99, 0.997, 0.999)
win = {kap: [dict(q=q, **S.cauchy_route_window(1.0, kap, q, M_lo=0.3, M_hi=2.2, n_scan=77)) for q in qs] for kap in (0.3, 1.0)}
out["cauchy_route"] = {f"kappa/chi={k}": v for k, v in win.items()}
dump("round6b_gsl.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.2))

# (a) product branches
a = ax[0, 0]
x = np.linspace(1e-4, 0.25, 600)
aA = (1 - np.sqrt(1 - 4 * x)) / (2 * x)
aN = (1 + np.sqrt(1 - 4 * x)) / (2 * x)
a.plot(x, aA, color=C1, label=r"AdS$_2\times S^2$ (throat)")
a.plot(x, aN, color=C2, label=r"dS$_2\times S^2$ (charged Nariai)")
for xv, lab in ((3 / 16, "bud grows\n(3/16)"), (0.25, "no static\nthroat (1/4)")):
    a.axvline(xv, color=MUTED, lw=1, ls="--")
    a.text(xv - 0.004, 9.0, lab, ha="right", va="top", color=MUTED, fontsize=8)
for k, mk in (("higgs", "o"), ("W", "s")):
    xm = out["electroweak_corona"][k]["x_max"]
    a.plot([xm], [(1 - np.sqrt(1 - 4 * xm)) / (2 * xm)], mk, color=C3, ms=8, mec="#fcfcfb", mew=2)
a.annotate("electroweak corona, self-consistent max\n(circle: $eB = m_H^2$; square: $eB = m_W^2$ with full\n$\\rho_{EW}$, an overestimate)",
           xy=(out["electroweak_corona"]["higgs"]["x_max"], 1.08), xytext=(0.012, 3.6),
           color=INK, fontsize=8, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
a.set_ylim(0.8, 10)
a.set_xlabel(r"$\Lambda Q^2$")
a.set_ylabel(r"$a^2/Q^2$ (sphere radius$^2$)")
a.set_title("(a) Einstein-Maxwell-$\\Lambda$ products: throat and its decompactification", fontsize=9)
a.legend(frameon=False, fontsize=8, loc="upper left")

# (b) flip window
b = ax[0, 1]
rho = np.linspace(0.05, 3.0, 300)
for c, col in ((1e2, C1), (1e4, C2), (1e6, C3)):
    tf = c / 12 / np.sinh(rho)
    b.semilogy(rho, tf, color=col, label=f"$T_{{flip}}$, c = $10^{{{int(np.log10(c))}}}$")
b.axhline(1.0, color=INK, lw=1.2)
b.text(0.1, 0.55, r"$T = T_{brk}$: Schwarzian regime below (semiclassics fails)", ha="left", fontsize=8, color=INK)
b.text(0.45, 3.0, "window, c = $10^2$", fontsize=8, color=INK)
b.fill_between(rho, 1.0, 1e2 / 12 / np.sinh(rho), where=1e2 / 12 / np.sinh(rho) > 1, color=C1, alpha=0.12, lw=0)
b.set_xlabel(r"proper distance $\rho/Q$ outward from the horizon (AdS$_2$ region $\rho \lesssim$ few $Q$)")
b.set_ylabel(r"$T / T_{brk}$")
b.set_ylim(0.3, 1e6)
b.set_title("(b) QNEC-compatible anti-trapping window: $T_{brk} < T < (c/12)T_{brk}/\\sinh(\\rho/Q)$", fontsize=9)
b.legend(frameon=False, fontsize=8, loc="upper right")

# (c) theta along the ray
cax = ax[1, 0]
for c, col in ((20.0, C1), (40.0, C2), (100.0, C3)):
    cf = S.cft_flip(c, r0, -np.sqrt(f0 / 2), np.sqrt(2) / Q, w=0.002 * Q, lam_max=3 * Q)
    cax.plot(cf["lam"] / Q, cf["theta"] / abs(cf["theta"][0]), color=col, label=f"c = {c:g}")
cax.axhline(0, color=INK, lw=1)
cax.text(2.95, 0.15, "anti-trapped above", ha="right", fontsize=8, color=INK)
cax.set_xlabel(r"affine parameter $\lambda/Q$ along the ingoing ray")
cax.set_ylabel(r"$\theta_-/|\theta_-(0)|$")
cax.set_title(f"(c) Pulse of scale Q at depth Q, T = {T / th.T_brk:.2f} $T_{{brk}}$: flip iff c > 12 $T\\sinh(1)/T_{{brk}}$ = "
              f"{12 * T / th.T_brk * np.sinh(1):.0f}", fontsize=9)
cax.legend(frameon=False, fontsize=8, loc="upper right")

# (d) Cauchy-route window
d = ax[1, 1]
for kap, col in ((0.3, C1), (1.0, C2)):
    rows = win[kap]
    om = np.array([1 - r["q"] for r in rows])
    d.fill_between(om, [r["M_a"] for r in rows], [r["M_b"] for r in rows], color=col, alpha=0.25, lw=0)
    d.plot(om, [r["M_a"] for r in rows], color=col, lw=1.5)
    d.plot(om, [r["M_b"] for r in rows], color=col, lw=1.5, label=f"$\\kappa/\\chi$ = {kap}: barrier-free window")
    d.axhline(1 / np.sqrt(1 + kap**2), color=col, lw=1, ls=":")
d.plot([], [], color=MUTED, lw=1, ls=":", label=r"extremal limit $M = 1/\sqrt{\chi^2+\kappa^2}$")
d.set_xscale("log")
d.set_xlabel(r"$1 - Q/M$")
d.set_ylabel(r"$M\chi$")
d.set_title("(d) Classical route through the Cauchy horizon: window width $\\leq r_+ - r_- = 2\\kappa_- r_-^2$", fontsize=9)
d.legend(frameon=False, fontsize=8, loc="upper left")
fig.tight_layout()
fig.savefig(f"{FIG}/round6b_gsl.png", dpi=150, facecolor="#fcfcfb")
print("wrote FIGURES/round6b_gsl.png")
