"""Round 6c, the answers table: the universal tipping law, invariant costs, quantum fields (designed pulses and vacuum polarisation),
the ghost-condensate window, what the initial geometry needs, and the 3+1 greybody factors of the Round 7 flux (P-R7)."""
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse import throat_evolution as TE, nariai as NA, quantum_throat as QT, greybody as GB

warnings.simplefilter("ignore")
out = {}
strip = lambda d: {k: v for k, v in d.items() if k not in ("res", "cls")}
gauss = lambda A, up: (lambda u: -A * np.exp(-((u - up) / 1.0) ** 2))
dipole = lambda D: (lambda u: -0.04 * np.exp(-((u - 2) / 1.0) ** 2) + 0.04 * np.exp(-((u - 2 - D) / 1.0) ** 2))
tophat = lambda a, T: (lambda u: -a * 0.25 * (1 + np.tanh((u - 2) / 1.0)) * (1 - np.tanh((u - 2 - T) / 1.0)))

# ------------------------------------------------ 1. the universal tipping law (delta = 1e-8, l_H = 30)
cases = [("gauss", A, 2, gauss(A, 2)) for A in (0.0025, 0.005, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.04, 0.08)]
cases += [("late", A, up, gauss(A, up)) for up, A in ((10, 0.02), (30, 0.02), (60, 0.02), (30, 0.04), (60, 0.06))]
cases += [("dipole", 0.04, D, dipole(D)) for D in (5, 10, 20, 40, 80)]
cases += [("tophat", a, T, tophat(a, T)) for T, a in ((10, 0.004), (10, 0.006), (40, 0.0015), (40, 0.003), (100, 0.0006), (100, 0.0012))]
rows = []
for kind, a, b, f in cases:
    L = 400.0 if (kind == "late" and b == 60) else 200.0
    d = TE.nariai_false_vacuum(tu_func=f, Lu=L, Lv=L)
    rows.append(dict(kind=kind, amplitude=a, parameter=b, grew=bool(d["grew"]), hidden=d["true_vacuum_untrapped_large_r"] == 0,
                     neg=d["neg_energy"], net=d["net_energy"], dm_misner_sharp=d["dm_misner_sharp"]))
fit = NA.tipping_law_fit([c[3] for c in cases], [r["grew"] for r in rows])
out["tipping_law"] = dict(runs=rows, fit=fit, r_N=17.3205)

# ------------------------------------------------ 2. threshold versus tuning, local and invariant
thr = []
for delta in (1e-10, 1e-8, 1e-6, 1e-4, 1e-3):
    for A in (0.005, 0.01, 0.02, 0.04, 0.08):
        d = TE.nariai_false_vacuum(delta=delta, Nflux=A)
        thr.append(dict(delta=delta, W_local=float(A * np.sqrt(np.pi)), dm_misner_sharp=d["dm_misner_sharp"], M_N=d["M_N"],
                        grew=bool(d["grew"]), hidden=d["true_vacuum_untrapped_large_r"] == 0))
out["threshold_vs_delta"] = thr
out["flat_space_flip_over_W"] = NA.flat_space_flip(17.3205) / (0.032)

# ------------------------------------------------ 3. quantum fields
q = dict(bound=QT.quantum_tipping_c(17.3205), one_loop=[], designed=[], weighted_anec=[])
for c in (30.0, 100.0, 300.0, 1000.0):
    q["one_loop"].append(dict(c=c, **{k: v for k, v in strip(TE.nariai_false_vacuum(c=c)).items() if k in ("grew", "r_max", "n_big")}))
for c, N in ((20.0, 0.04), (60.0, 0.04), (150.0, 0.04), (300.0, 0.04), (300.0, 0.1)):
    d = TE.nariai_false_vacuum(c=c, Nflux=N)
    q["one_loop"].append(dict(c=c, Nflux=N, grew=bool(d["grew"]), hidden=d["true_vacuum_untrapped_large_r"] == 0))
for c, A, w, L, W2 in ((100, 0.05, 1, 60, 20), (100, 0.1, 1, 60, 20), (300, 0.03, 1, 150, 30), (300, 0.015, 1, 150, 30),
                       (500, 0.03, 10, 150, 30)):
    tu, _ = TE.conformal_pulse(c, A, 2.0, w, L, W2)
    d = TE.nariai_false_vacuum(c=c, tu_func=tu)
    q["designed"].append(dict(c=c, A=A, w=w, L=L, W2=W2, grew=bool(d["grew"]), neg=d["neg_energy"], net=d["net_energy"],
                              W=NA.weighted_deficit(tu, NA.TAU_OVER_RN * 17.3205)))
for x in (0.25, 0.5, 1.0, 2.0, 4.0):
    q["weighted_anec"].append(dict(A_tau=x, ratio_to_bound=QT.weighted_anec_of_pulse(100.0, x / 35.0, 35.0) / QT.weighted_anec_bound(100.0, 35.0)))
out["quantum_fields"] = q

# ------------------------------------------------ 4. classical field: the ghost-condensate window
out["ghost_condensate"] = {f"T/tau={t}": NA.ghost_condensate_window(17.3205, t) for t in (0.5, 1.0, 3.0, np.inf)}

# ------------------------------------------------ 5. what the initial geometry needs
prep = []
for ui in (None, 6, 12, 25, 40, 60, 100, 150, 195):
    for A in ((0.04, 0.16) if ui == 60 else (0.04,)):
        d = TE.nariai_false_vacuum(Nflux=A, u_inner=ui)
        prep.append(dict(u_inner=ui, Nflux=A, grew=bool(d["grew"]), hidden=d["true_vacuum_untrapped_large_r"] == 0, u_hole=d["u_hole"]))
for wf in (0.5, 0.9, 0.995):
    d = TE.nariai_false_vacuum(Nflux=0.04, wall_frac=wf)
    prep.append(dict(wall_frac=wf, Nflux=0.04, grew=bool(d["grew"]), hidden=d["true_vacuum_untrapped_large_r"] == 0))
out["preparation"] = prep
out["partner_vacuum"] = dict(M_sun_GeV=NA.partner_vacuum(NA.M_SUN_KG), GUT_mass_kg=NA.partner_mass(1e16),
                             EW_mass_earths=NA.partner_mass(100.0) / 5.972e24, QCD_02_mass_suns=NA.partner_mass(0.2) / NA.M_SUN_KG)

# ------------------------------------------------ 6. P-R7: greybody factors and the 3+1 branch difference
gb = dict(low=[(w, GB.transmission(0, w) / (16 * w * w)) for w in (0.005, 0.01, 0.02)],
          table={l: [(w, GB.transmission(l, w)) for w in (0.1, 0.2, 0.3, 0.4, 0.6)] for l in (0, 1, 2)})
rates = GB.hawking_rates()
gb["rates"] = {k: v for k, v in rates.items() if k != "per_l"}
gb["per_l"] = rates["per_l"]
gb["branch_difference"] = {f"D={D}": GB.branch_difference(D=D, rates=rates) for D in (1.94, 2.02, 2.18)}
out["greybody"] = gb
dump("round6c_answers.json", out)

# ------------------------------------------------ figure
C1, C2, C3, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED,
                     "ytick.color": MUTED})
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))
tau = NA.TAU_OVER_RN * 17.3205
W = np.array([NA.weighted_deficit(c[3], tau) for c in cases])
g = np.array([r["grew"] for r in rows])
kinds = [r["kind"] for r in rows]
for k, col, mk in (("gauss", C1, "o"), ("late", C2, "s"), ("dipole", C3, "^"), ("tophat", INK, "D")):
    sel = np.array([kk == k for kk in kinds])
    ax[0].scatter(W[sel], g[sel] + 0.04 * np.random.default_rng(1).standard_normal(sel.sum()), color=col, marker=mk, s=36,
                  edgecolor="#fcfcfb", label=k)
ax[0].axvline(NA.WTH_OVER_RN * 17.3205, color=INK, lw=1)
ax[0].set_xscale("log")
ax[0].set_xlim(3e-3, 0.3)
ax[0].set_yticks([0, 1], ["no daughter", "daughter grows"])
ax[0].set_xlabel("weighted NEC deficit W = -int N e^{-u/tau} du  (tau = 3.2 r_N)")
ax[0].set_title("(a) one law for 25 runs: grows iff W > 1.85e-3 r_N/G", fontsize=9)
ax[0].legend(frameon=False, fontsize=8, loc="center left")
for delta, col in ((1e-8, C3), (1e-6, C2), (1e-4, C1)):
    rr = [t for t in thr if t["delta"] == delta]
    ax[1].plot([t["dm_misner_sharp"] / t["M_N"] for t in rr], [1 if t["grew"] else 0 for t in rr], "o-", color=col, ms=6,
               mec="#fcfcfb", label=f"delta = {delta:g}")
ax[1].set_xscale("log")
ax[1].set_yticks([0, 1], ["no", "yes"])
ax[1].set_xlabel("invariant (Misner-Sharp) pulse energy / M_N")
ax[1].set_title("(b) invariant cost ~ 1e-2 sqrt(delta) M_N", fontsize=9)
ax[1].legend(frameon=False, fontsize=8, loc="center right")
gbl = gb["table"]
for l, col in ((0, C1), (1, C2), (2, C3)):
    ax[2].plot([w for w, _ in gbl[l]], [G for _, G in gbl[l]], "o-", color=col, ms=5, mec="#fcfcfb", label=f"l = {l}")
ax[2].set_xlabel("M omega")
ax[2].set_ylabel("Gamma_l (massless scalar)")
ax[2].set_title("(c) P-R7: Regge-Wheeler greybody factors", fontsize=9)
ax[2].legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(f"{FIG}/round6c_answers.png", dpi=150, facecolor="#fcfcfb")
print("tipping fit:", fit["tau_range"], "rates:", gb["rates"], "branch:", gb["branch_difference"])
