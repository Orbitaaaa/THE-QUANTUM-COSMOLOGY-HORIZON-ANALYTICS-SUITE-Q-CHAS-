"""Round 6i: the seven candidate routes to a daughter that does not need a manufactured false-vacuum ball (curvature activation,
monopole cores, breathing bubbles, gauged/vacuum-dependent ghost condensation, resonant tunnelling, classical transitions, four-form
membranes).  Writes DATA/round6i.json and FIGURES/round6i.png; `--figure-only` redraws the figure."""
import json
import sys
import warnings
from dataclasses import replace

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import gauss_bonnet as GB, stable_shells as SS, routes6i as R6, accessibility as AC

warnings.simplefilter("ignore")


def monopole_section():
    out = dict(critical=[SS.monopole_critical_eta(l, e, 0.0, np.linspace(0.05, 1.0, 96)) for l, e in ((1.0, 0.5), (0.1, 0.5), (1.0, 1.0))])
    etas = [0.003, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.24]
    out["static"] = [dict(eta=x, **(SS.monopole_static(x) or {})) for x in etas]
    out["own_mass_channel"] = [dict(eta=x, **SS.monopole_tunnelling(x)) for x in etas]
    Mf = 4 * np.pi * 0.01 / 0.5
    out["nec_admixture"] = [dict(m_d_over_M=f, **SS.monopole_tunnelling(0.01, 1.0, 0.5, f * Mf)) for f in (0.0, -0.1, -0.2, -0.24)]
    con = []
    for eta in (0.003, 0.01):
        Q = np.sqrt(4 * np.pi) / 0.5
        for MQ in (1.05, 1.5, 2.0, 4.0, 10.0):
            sh = replace(SS.monopole_shell(eta, 1.0, 0.5), M=MQ * Q)
            regs = sh.regions(0.01, 3e5)
            bar = regs[-2]
            con.append(dict(eta=eta, M_over_Q=MQ, pattern="".join("A" if r["allowed"] else "F" for r in regs), barrier=bar["R"],
                            B=float(sh.B(*bar["R"])), beta_out_inflating=regs[-1]["beta_out"]))
    out["connected_channel"] = con
    st = SS.monopole_static(0.05)
    added = []
    for fac in (1.0, 2.0, 4.0, 6.0, 10.0):
        sh = replace(SS.monopole_shell(0.05), M=st["M"] * fac)
        regs = sh.regions(0.05, 5000)
        added.append(dict(M_over_M0=fac, horizon=bool(sh.M > np.sqrt(sh.Q2)), barrier=[r["R"] for r in regs if not r["allowed"]][-1],
                          pattern="".join("A" if r["allowed"] else "F" for r in regs)))
    out["added_mass"] = added
    return out


if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round6i.json"))
else:
    out = {}
    # ---------------------------------------------------------------- route 1
    eq = GB.ks_equations()
    u = np.linspace(0, 40, 40001)
    out["gb_identities"] = dict(divergence=GB.verify_divergence_identity(), null_closed_form=GB.null_component_closed_form(),
                                ks_total_derivative=eq["total_derivative_ok"],
                                weighted_example=GB.weighted_gb_deficit(0.01 * np.tanh(u / 3) ** 2, u, 3.2))
    out["scalarisation"] = GB.scalarisation_scan()
    out["scalarisation_paths"] = [GB.ks_run(GB.coupling_dy, GB.potential_flat(), l, 1e-3, 0.0, path=True) for l in (0.05, 0.5, 2.0)]
    rev = GB.reverse_scan()
    out["reverse_counts"] = {k: sum(r["outcome"] == k for r in rev) for k in ("daughter", "black hole", "kinetic degeneracy", "undecided")}
    out["reverse_runs"] = rev
    out["reverse_paths"] = [GB.ks_run(GB.coupling_reverse, GB.potential_hilltop(4 * 10 * 2.0, 0.9), 10.0, 0.0, 1e-3, path=True),
                            GB.ks_run(GB.coupling_reverse_saturating(0.2), GB.potential_hilltop(4 * 5 * 2.0, 0.9), 5.0, 0.0, 1e-3, path=True)]
    out["control_no_gb"] = GB.ks_run(GB.coupling_reverse, GB.potential_hilltop(22.0, 0.9), 0.0, 0.0, 1e-3)
    out["activation"] = GB.activation_report()
    out["hysteresis"] = GB.hysteresis_bound()
    # ---------------------------------------------------------------- routes 2, 3, 5, 7
    out["monopole"] = monopole_section()
    out["breathing_scan"] = SS.breathing_random_scan(6000)
    out["breathing_family"] = SS.breathing_barrier_vs_mass(1.0, -0.1, n=25)
    ex = SS.Shell(chi=1.0, M=0.05, k_w=1.0, m_d=-0.1)
    Rg = np.geomspace(0.02, 2.0, 600)
    out["breathing_example"] = dict(R=Rg.tolist(), V=np.clip(ex.V(Rg), -3, 2).tolist(), R0=float(np.sqrt(0.1)),
                                    regions=ex.regions(1e-3, 20.0), classify=SS.classify_breathing(ex))
    tw = SS.find_two_well()
    out["two_well"] = [dict(shell=str(x["shell"]), **SS.resonance_scan(x["shell"])) for x in tw]
    if tw:
        out["resonant_masses_example"] = SS.resonant_masses(tw[0]["shell"], 10.0)
    out["membrane"] = SS.membrane_branches()
    # ---------------------------------------------------------------- routes 4, 6, ledger
    out["gauged_window"] = [R6.gauged_window(c) for c in (3.0, 30.0)]
    out["vacuum_dependent"] = R6.vacuum_dependent()
    out["overall_factor_cutoff_GeV"] = R6.strong_coupling_scale_overall_factor((100 / 1e13) ** 4, 1e13)
    out["pocket_scan"] = R6.pocket_scan()
    out["dtr_example"] = R6.dtr_region_c(10.0, 1.0, 0.5, 0.03)
    out["ledger"] = R6.ledger()
    out["bridge_invariant"] = R6.bridge_invariant()
    out["side_change_lemma"] = R6.side_change_lemma()
    mb0 = out["membrane"]["rows"][0]
    out["no_free_lunch"] = dict(monopole_B_own_mass=out["monopole"]["own_mass_channel"][1]["B"], monopole_cdl=SS.monopole_cdl(1.0),
                                membrane_B_M0=mb0["B"], membrane_cdl=SS.cdl_flat(1e-3, 0.05 * (0.05 + 0.025)))
    dump("round6i.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.4))
ax = ax.ravel()
# (a) the homogeneous Nariai throat with the Gauss-Bonnet scalar
for p, col, lab in zip(out["scalarisation_paths"], (C1, C3, C4), ("lam^2 = 0.05", "lam^2 = 0.5", "lam^2 = 2")):
    ax[0].plot(p["t"], p["b"], color=col, label=f"scalarising, {lab} r_N^2: {p['outcome']}")
for p, col, lab in zip(out["reverse_paths"], (C2, MUTED), ("f = -phi^2/2, lam^2 = 10", "saturating f, lam^2 = 5")):
    ax[0].plot(p["t"], p["b"], color=col, ls="--", label=f"reversed sign, {lab}: {p['outcome']}")
ax[0].axhline(1.0, color=INK, lw=0.6)
ax[0].set_ylim(0, 1.6)
ax[0].set_xlabel("proper time / r_N"); ax[0].set_ylabel("2-sphere radius b / r_N")
c = out["reverse_counts"]
ax[0].set_title(f"(a) Nariai throat + Gauss-Bonnet scalar (nonlinear, homogeneous)\nreversed-sign scan: {c['daughter']} daughters, "
                f"{c['black hole']} black holes, {c['kinetic degeneracy']} degenerate", fontsize=8)
ax[0].legend(frameon=False, fontsize=6.5, loc="lower left")
# (b) breathing bubble potential
be = out["breathing_example"]
ax[1].plot(be["R"], be["V"], color=C1)
ax[1].axhline(0, color=INK, lw=0.6)
ax[1].axvline(be["R0"], color=C2, ls="--", lw=1)
ax[1].text(be["R0"] * 1.05, 1.2, "net surface energy = 0\n(shell moves at c)", color=C2, fontsize=7)
ax[1].set_xscale("log"); ax[1].set_ylim(-3, 1.8)
ax[1].set_xlabel("shell radius / false-vacuum Hubble length"); ax[1].set_ylabel("V(R)   (Rdot^2 = -V)")
bs = out["breathing_scan"]
ax[1].set_title(f"(b) breathing bubble (wall k_w = 1, dust m_d = -0.1, M = 0.05): well, barrier, hidden daughter\n"
                f"random scan: {bs['wells_crossing_zero_energy']} of {bs['wells']} wells cross zero net energy", fontsize=8)
# (c) monopole tunnelling exponents
mo = out["monopole"]
own = [r for r in mo["own_mass_channel"] if r.get("B")]
ax[2].plot([r["eta"] for r in own], [r["B"] for r in own], "o-", color=C3, label="own mass (M << Q): ends disconnected from our universe")
for eta, col in ((0.003, C2), (0.01, C1)):
    rows = [r for r in mo["connected_channel"] if r["eta"] == eta]
    ax[2].plot([eta] * len(rows), [r["B"] for r in rows], "s", color=col, ms=6,
               label=f"connected (M = 1.05-10 Q, a magnetic black hole), eta = {eta}")
for cr in mo["critical"]:
    lo, hi = cr["eta_c_between"]
    ax[2].axvline(hi, color=MUTED, lw=0.8, ls=":")
ax[2].text(0.27, 3e3, "eta_c\n(0.26-0.40)", fontsize=7, color=MUTED)
ax[2].set_xscale("log"); ax[2].set_yscale("log")
ax[2].set_xlabel("symmetry-breaking scale eta / m_P"); ax[2].set_ylabel("tunnelling exponent B")
ax[2].set_title("(c) monopole core (thin shell, Sakai et al. 2006): stable below eta_c;\na daughter attached to us needs M > Q first", fontsize=8)
ax[2].legend(frameon=False, fontsize=6.5, loc="lower left")
# (d) quantum barriers: breathing family and membrane, in false-vacuum Hubble units
bf = out["breathing_family"]
ax[3].plot([r["M"] for r in bf], [r["B"] for r in bf], "o-", color=C1, ms=3, label="breathing bubble (B / (l_f/l_P)^2)")
mb = out["membrane"]["rows"]
ax[3].plot([r["M"] * out["membrane"]["chi"] for r in mb], [r["B"] * out["membrane"]["chi"] ** 2 for r in mb], "-", color=C4,
           label="four-form membrane (Farhi-Guth wall)")
ax[3].set_xscale("log"); ax[3].set_yscale("log")
ax[3].set_xlabel("mass x chi (false-vacuum Hubble units)"); ax[3].set_ylabel("B / (l_f / l_P)^2")
ax[3].set_title("(d) shape factors of the tunnelling exponent; B = shape x (l_f/l_P)^2\n(GUT vacuum: (l_f/l_P)^2 ~ 1e10)", fontsize=8)
ax[3].legend(frameon=False, fontsize=7)
# (e) the gravitational invariant: every classical route and humanity's capabilities
inv = AC.requirement_invariant()
vd = out["vacuum_dependent"]
scen = [("universal condensate floor\n(2.5e9 M_sun)", AC.r_N_min_from_condensate_bound(100.0)),
        ("vacuum-dependent condensate\n(9e8 kg, Prop. 6i-9)", None)]
labels, vals = [], []
for name, sc in scen:
    if sc is None:
        E, rho = vd["M_N_kg"] * AC.C_LIGHT**2, AC.RHO_P / (8 * np.pi * (vd["r_N_m"] / AC.L_P) ** 2)
    else:
        E, rho = sc["E_N_J"], sc["rho_f_J_m3"]
    labels.append(name); vals.append(np.log10(E * np.sqrt(rho)))
for row in AC.accessibility_table():
    labels.append(row["name"]); vals.append(np.log10(row["fom"]))
cols = [C2, C2] + [C1 if r["kind"] == "technology" else C3 for r in AC.accessibility_table()]
ax[4].barh(range(len(labels)), vals, color=cols)
ax[4].axvline(np.log10(inv), color=INK, lw=1, ls="--")
ax[4].text(np.log10(inv) - 1, len(labels) - 0.6, "invariant\n(Theorem 6f-7)", ha="right", fontsize=7)
ax[4].set_yticks(range(len(labels))); ax[4].set_yticklabels(labels, fontsize=6.5)
ax[4].set_xlabel("log10 E sqrt(rho)   (J^{3/2} m^{-3/2})")
ax[4].set_title("(e) a smaller background does not shrink the gravitational invariant:\nevery near-Nariai route needs E sqrt(rho) ~ 1e65", fontsize=8)
# (f) the ledger: state / trigger / completion
L = out["ledger"]
grid = np.array([[1, 0, 0], [1, 0, 1], [0, 0, 1], [1, 2, 2], [1, 0, 1], [1, 0, 0], [1, 0, 1]], float)   # 2 achieved, 1 conditional, 0 blocked
cm = matplotlib.colors.ListedColormap(["#f2a07e", "#f6dc8a", "#9fd8bd"])
ax[5].imshow(grid, cmap=cm, vmin=0, vmax=2, aspect="auto")
ax[5].set_xticks(range(3)); ax[5].set_xticklabels(["create the state", "trigger", "complete the daughter"], fontsize=7)
ax[5].set_yticks(range(len(L))); ax[5].set_yticklabels([r["route"] for r in L], fontsize=7)
ax[5].grid(False)
txt = [["horizon-class\ncurvature", "pro-collapse", "degenerates"], ["natural relic?\nM>Q needs m_P", "e^-B, B>3e3", "after tunnel"],
       ["NEC-violating\ndust", "e^-B", "after tunnel"], ["9e8 kg ball\n(conditional)", "condensate", "simulated"],
       ["two wells exist", "removes inner\nbarrier only", "as tunnelling"], ["pocket forms", "collision", "normal parent:\nno"],
       ["flux pocket", "e^-B", "after tunnel"]]
for i in range(len(L)):
    for j in range(3):
        ax[5].text(j, i, txt[i][j], ha="center", va="center", fontsize=6.3, color=INK)
ax[5].set_title("(f) the seven routes: green done, yellow conditional or natural only, red blocked", fontsize=8)
fig.tight_layout()
fig.savefig(f"{FIG}/round6i.png", dpi=150, facecolor="#fcfcfb")
print("done")
