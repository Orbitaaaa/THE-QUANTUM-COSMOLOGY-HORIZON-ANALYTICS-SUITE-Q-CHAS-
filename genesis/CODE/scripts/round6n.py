"""Round 6n: model U beyond idealised collapse: perfect-fluid hydrodynamics (pressure, inhomogeneity, shocks) with the generalised
Misner-Sharp equations, the universality of the bounce, inhomogeneous dust, and the inner-horizon phase.  Writes DATA/round6n.json and
FIGURES/round6n.png; `--figure-only` redraws the figure."""
import json
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import hydro6n as H
from babyuniverse import unified6l as U

warnings.simplefilter("ignore")


def os_track(D, N=80, mu_over_ext=20.0):
    st = H.make_conservative(D, N, mu_over_ext, "uniform", 0.0)
    psi0 = float(H.psi_and_derivative(st.mu[-1] / st.R[-1] ** (D - 1))[0])
    amax = 1 / np.sqrt(psi0)
    sinchi = st.R[-1] / amax
    C = st.mu[-1] / st.R[-1] ** (D - 1) * amax ** (D - 1)
    an = U.closed_dust(C, D)
    T = 2.0 * an["tau_half_cycle"]
    ts, Rs = [0.0], [st.R[-1]]
    while st.t < T:
        st.step()
        ts.append(st.t); Rs.append(st.R[-1])
    # analytic: integrate adot^2 = a^2 psi - 1 from a_max (contracting), bounce, re-expansion
    a = np.concatenate([np.linspace(amax, an["a_b"], 2000), np.linspace(an["a_b"], amax, 2000)[1:]])
    F = np.maximum(U._adot2(a, C, D, 1.0, "hayward", 0.0), 1e-300)
    tau = np.concatenate([[0], np.cumsum(np.abs(np.diff(a)) / np.sqrt(0.5 * (F[1:] + F[:-1])))])
    return dict(D=D, t=ts[::max(1, len(ts) // 800)], R=Rs[::max(1, len(Rs) // 800)], t_an=tau[::5].tolist(), R_an=(a * sinchi)[::5].tolist(),
                check=H.os_check(D, N, mu_over_ext=mu_over_ext))


if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round6n.json"))
else:
    out = {}
    out["os"] = [os_track(5), os_track(10)]
    out["universal_bounce"] = [H.frw_bounce_any_fluid(w, c, D) for D in (5, 10) for w in (-0.3, 0.0, 0.25, 1 / 3, 0.5, 1.0) for c in (1e2, 1e4, 1e8)]
    out["ltb"] = [H.ltb_bounce_radii(5, p) for p in ("smooth", "cored", "gaussian")]
    out["throat_runs"] = [dict(**{k: v for k, v in H.run_to_throat(5, 50, mu, p, e).items() if k != "profile_end"})
                          for (p, e, mu) in (("smooth", 0.1, 20.0), ("smooth", 0.1, 60.0), ("cored", 0.1, 20.0), ("smooth", 0.3, 20.0),
                                             ("gaussian", 0.2, 20.0), ("smooth", 0.0, 20.0))]
    out["throat_profile"] = H.run_to_throat(5, 50, 20.0, "smooth", 0.1)["profile_end"]
    out["trajectories"] = H.trajectories(5, 50, 20.0, "smooth", 0.1)
    out["pressure_threshold"] = H.pressure_threshold()
    dump("round6n.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
ax = ax.ravel()
# (a) OS validation
for o, col in zip(out["os"], (C1, C2)):
    ax[0].plot(o["t_an"], o["R_an"], color=col, lw=1, ls="--")
    ax[0].plot(o["t"], o["R"], color=col, lw=2.5, alpha=0.6, label=f"D = {o['D']}: hydro (solid) vs exact (dashed), bounce error {o['check']['rel_err_R']:.0e}")
ax[0].set_xlabel("proper time of the surface (units of l)"); ax[0].set_ylabel("surface areal radius R / l")
ax[0].set_title("(a) the hydrodynamics code against the exact Oppenheimer-Snyder bounce of model U\n(energy conserved to round-off)", fontsize=8)
ax[0].legend(frameon=False, fontsize=7)
# (b) universality of the bounce
for D, col in ((5, C1), (10, C2)):
    rows = [r for r in out["universal_bounce"] if r["D"] == D and r["C"] == 1e4]
    ax[1].plot([r["w"] for r in rows], [r["a_b"] for r in rows], "o-", color=col, label=f"bounce a_b / l, D = {D}")
    ax[1].plot([r["w"] for r in rows], [r["a_max"] for r in rows], "s--", color=col, alpha=0.6, label=f"turnaround a_max / l, D = {D}")
ax[1].set_yscale("log"); ax[1].axhline(1, color=INK, lw=0.7)
ax[1].set_xlabel("equation of state w = p / rho"); ax[1].set_ylabel("scale factor (units of l)")
ax[1].set_title("(b) every closed universe bounces at a_b >= l, whatever the fluid (psi <= 1/alpha):\nProposition 6n-1", fontsize=8)
ax[1].legend(frameon=False, fontsize=7)
# (c) inhomogeneous dust (LTB)
for r, col in zip(out["ltb"], (C1, C3, C2)):
    rows = r["rows"]
    ax[2].plot([q["R0"] for q in rows], [q["tau_bounce"] for q in rows], color=col, label=f"{r['profile']}: bounce time of each shell")
ax2 = ax[2].twinx()
for r, col in zip(out["ltb"], (C1, C3, C2)):
    rows = r["rows"]
    ax2.plot([q["R0"] for q in rows], [q["R_b"] / max(q["R_b_law"], 1e-9) for q in rows], ls=":", color=col)
ax2.set_ylim(0.95, 1.1); ax2.set_ylabel("R_bounce / sqrt(k alpha)  (dotted)", color=MUTED); ax2.grid(False)
ax[2].set_xlabel("initial shell radius R_0 / l"); ax[2].set_ylabel("proper time to bounce (units of l)")
ax[2].set_title("(c) inhomogeneous dust: every shell bounces at R_b = sqrt(k alpha) (within 2 %),\ninner shells first, so the bounced core meets the envelope", fontsize=8)
ax[2].legend(frameon=False, fontsize=7, loc="upper left")
# (d) spacetime diagram of the pressured run
tr = out["trajectories"]
t = np.array(tr["t"]); R = np.array(tr["R"]); S = np.array(tr["status"])
cmap = {0: C3, 1: C1, 2: C2}
for j in range(R.shape[1]):
    for code in (0, 1, 2):
        m = S[:, j] == code
        ax[3].scatter(t[m], R[m, j], s=1.5, color=cmap[code])
ax[3].set_yscale("log"); ax[3].axhline(1.0, color=INK, lw=0.7, ls=":")
from matplotlib.lines import Line2D
handles = [Line2D([], [], marker="o", ls="", color=c, ms=4, label=l) for c, l in
           ((C3, "untrapped"), (C1, "trapped (black-hole region)"), (C2, "anti-trapped (white-hole region)"))]
ax[3].set_xlabel("proper time of the surface (units of l)"); ax[3].set_ylabel("shell areal radius R / l")
ax[3].set_title("(d) pressured, inhomogeneous collapse (D = 5): horizon, untrapped core bounce, pile-up at the\ninner horizon (R ~ l); squares: last slice, where a shell has entered the white-hole region", fontsize=8)
prof = out["throat_profile"]
ax[3].scatter([t[-1] * 1.004] * len(prof["R"]), prof["R"], s=8, marker="s",
              c=[C3 if g >= 0 else (C1 if u < 0 else C2) for g, u in zip(prof["gradR2"], prof["U"])], zorder=5)
ax[3].legend(handles=handles, frameon=False, fontsize=7, loc="lower left")
# (e) diagnostics
ax[4].semilogy(t, tr["rho_c"], color=C1, label="central density")
ax[4].semilogy(t, np.array(tr["psi_max"]), color=C2, label="max psi x alpha (bounded by 1)")
ax[4].semilogy(t, np.maximum(np.array(tr["energy"]), 1e-17), color=MUTED, ls=":", label="total-energy error")
ax[4].set_xlabel("proper time of the surface (units of l)")
ax[4].set_title("(e) curvature stays below the bound and energy is conserved to round-off\nthrough horizon formation, bounce and the inner-horizon collision", fontsize=8)
ax[4].legend(frameon=False, fontsize=7)
# (f) outcomes table
rows = out["throat_runs"]
labels = [f"{r['profile']}, eps0={r['eps0']}, mu={r['D'] and round(r.get('rho_c_max', 0), 0) and ''}" for r in rows]
mus = [20, 60, 20, 20, 20, 20]
cell = [[r["profile"], f"{r['eps0']}", f"{m}", "yes" if r["first_trapped"] else "no", "yes" if r["first_bounce"] else "no",
         "yes" if r["first_anti_trapped"] else "-", f"{r['psi_max_alpha']:.3f}", "throat" if r["end"] == "throat" else "shock"]
        for r, m in zip(rows, mus)]
thr = out["pressure_threshold"]
ax[5].axis("off")
tb = ax[5].table(cellText=cell, colLabels=["profile", "eps0", "mu/mu_ext", "horizon", "bounce", "white-hole\nentry", "max psi l^2", "slicing\nends at"],
                 loc="upper center", cellLoc="center")
tb.auto_set_font_size(False); tb.set_fontsize(7); tb.scale(1, 1.5)
ax[5].text(0.0, 0.18, "pressure threshold (smooth profile, mu = 20 mu_ext):  " + ", ".join(f"eps0={r['eps0']}: {'BH' if r['horizon'] else 'no BH'}" for r in thr),
           fontsize=7, transform=ax[5].transAxes)
ax[5].text(0.0, 0.08, "Beyond the inner-horizon collision the comoving (Gamma > 0) slicing ends: at a throat (anti-trapped next to trapped,\nGamma -> 0) or at a shock; which one depends on resolution. A characteristic code is needed there (P-R6n-a).",
           fontsize=7, transform=ax[5].transAxes)
ax[5].set_title("(f) outcomes of the pressured runs (D = 5)", fontsize=8)
fig.tight_layout()
fig.savefig(f"{FIG}/round6n.png", dpi=150, facecolor="#fcfcfb")
print("done")
