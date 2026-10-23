"""Round 3 (problems A, B, C): selection rule and its exact toy, exterior event signature, assistance frequency."""
import os
import warnings

import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse import selection as S, signatures as Sg, bounds as Bd

warnings.simplefilter("ignore")
out = {}

# ------------------------------------------------ A: toy
out["toy"] = S.toy_report(us=(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3))
out["toy_amplitude"] = [dict(h=h, a=a, **S.toy_amplitude_check(h=h, a=a)) for h in (0.1, 0.05) for a in (0.3, 0.5, 0.8)]
out["charge_signs"] = {name: S.charge_signs(m) for name, m in
                       (("single_chi=kappa=1", ThinWallModel(1.0, 1.0, 0.1545)),
                        ("double_gamma=0.6", model_from_dimensionless(0.6, 2.4)),
                        ("double_gamma=0.3", model_from_dimensionless(0.3, 5.0)),
                        ("double_light_wall", ThinWallModel(1.0, 0.05, 0.3 * ThinWallModel(1.0, 0.05, 1e-9).M_critical)))}

# ------------------------------------------------ C: light curves
mC = model_from_dimensionless(0.6, 2.4)
lc = Sg.light_curve(mC)
out["late_photons"] = {name: Sg.late_photon_fraction(m) for name, m in
                       (("gamma=0.6,E=2.4", mC), ("gamma=0.3,E=5", model_from_dimensionless(0.3, 5.0)))}

# ------------------------------------------------ B: assistance
o = Bd.frontier_optimum(10.0, 1.0, 100.0)
out["assistance"] = dict(frontier=Sg.assistance_frequency(o["E_vac_GeV"], o["E_wall_GeV"], 1.0),
                         gut=Sg.assistance_frequency(1e16, 1e16, 10.0))
dump("round3_abc.json", out)

# ------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0})
fig, ax = plt.subplots(1, 3, figsize=(13, 4.0))
lam, h, x, sols = S._toy(1.0, 0.1, 30)
us = np.linspace(-2.5, 2.5, 161)
for name, col, lab in (("H1", C1, r"$H^{(1)}_\nu$: geometric in lab, bounded in daughter (K$\equiv$0)"),
                       ("J", C2, r"$J_\nu$: grows $e^{+\pi u/h}$ in daughter"),
                       ("decaying", C3, r"$e^{-i\pi\nu}J_{-\nu}$: non-geometric in lab")):
    if name == "H1":   # exact reflection identity H1_{-nu} = e^{i pi nu} H1_nu: |psi(u)| = e^{-pi u/h} |psi(-u)| (avoids cancellation)
        ys = [float(abs(sols[name](mp.mpf(u)))) if u <= 0 else float(mp.e ** (-mp.pi * u / h) * abs(sols[name](-mp.mpf(u)))) for u in us]
    else:
        ys = [float(abs(sols[name](mp.mpf(u)))) for u in us]
    ax[0].semilogy(us, ys, color=col, label=lab)
ax[0].axvspan(-1, 1, color="#f1f0ec", zorder=0)
ax[0].text(-2.4, 1e-60, "laboratory\n(C>1)", color=MUTED, fontsize=8); ax[0].text(-0.45, 1e-60, "barrier", color=MUTED, fontsize=8)
ax[0].text(1.3, 1e-60, "daughter\n(C<-1)", color=MUTED, fontsize=8)
ax[0].set_ylim(1e-70, 1e45)
ax[0].set_xlabel(r"$u$ (toy: $C=-u$, $h=0.1$)"); ax[0].set_ylabel(r"$|\psi|$")
ax[0].set_title("A. Exact single-flip toy: the admissible solution", loc="left", fontsize=10)
ax[0].legend(fontsize=7, frameon=False, loc="upper left")

aa = np.linspace(0.05, 0.95, 19)
for hh, mk in ((0.1, "o"), (0.05, "s")):
    meas = [S.toy_amplitude_check(h=hh, a=a)["measured_exponent"] * hh for a in aa]
    ax[1].plot(aa, meas, mk, ms=8, mfc="none", mec=INK, label=f"exact $H^{{(1)}}$, h={hh}")
ax[1].plot(aa, np.pi * aa, color=C1, label=r"analytic angle: $\pi a$ ($B_{\rm an}$)")
ax[1].plot(aa, 2 * (aa * np.arccos(aa) - np.sqrt(1 - aa**2) + 1), color=C2, ls="--", label=r"sector-blind: $\int\arccos|C|$ ($B_{\rm sq}$)")
ax[1].set_xlabel(r"half-width $a$ of the barrier segment"); ax[1].set_ylabel(r"$h\,\ln|\psi(-a)/\psi(a)|$")
ax[1].set_title("A. Decay across the coalescence point", loc="left", fontsize=10)
ax[1].legend(fontsize=7.5, frameon=False)

uu = np.concatenate([lc["u_expanding"], lc["u_collapse"]])
rr = np.concatenate([lc["rate_expanding"], lc["rate_collapse"]])
u0 = lc["u_collapse"][0]
ax[2].semilogy(uu - u0, rr, color=C2, label="collapse branch (exponential tail)")
ax[2].semilogy(lc["u_expanding"] - u0, lc["rate_expanding"], color=C1, lw=3.5, label="daughter branch: emission stops at $u_0$")
ax[2].annotate("no signal after $u_0$\n(last redshift $\\sqrt{f(R_1)}$)", xy=(0.02, lc["rate_expanding"][-1]), xytext=(0.6, 3e-3),
               fontsize=8, color=C1, arrowprops=dict(arrowstyle="->", color=C1))
ax[2].axvline(0, color=MUTED, lw=0.8)
ax[2].set_xlim(-0.6, 3.0); ax[2].set_ylim(1e-8, 2)
ax[2].set_xlabel(r"arrival time $u-u_0$ (units $1/\chi$)"); ax[2].set_ylabel(r"photon rate $d\tau/du$ (per unit emission rate)")
ax[2].set_title("C. Exterior event signature", loc="left", fontsize=10)
ax[2].legend(fontsize=8, frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "round3_abc.png"), dpi=150)
print("toy amplitude:", [(r["h"], r["a"], round(r["measured_exponent"], 6), round(r["analytic_exponent"], 6)) for r in out["toy_amplitude"]])
print("late photons:", out["late_photons"])
print("assistance:", out["assistance"])
