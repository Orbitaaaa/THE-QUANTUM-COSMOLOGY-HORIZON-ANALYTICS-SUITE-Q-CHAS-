"""Round 6l: unified quantum-gravity formation.  Model U (D = 4 + n Einstein gravity + a quasi-topological tower: P1 and P2 in one action),
its interiors, collider production with real parton distributions, evaporation, multiple daughters from one hole, and what experiments can
distinguish.  Writes DATA/round6l.json (and, if a Pythia 8 PDF grid is available, refreshes DATA/round6l_lumi*.json) and FIGURES/round6l.png;
`--figure-only` redraws the figure."""
import json
import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import unified6l as U
from babyuniverse import collider6l as C

warnings.simplefilter("ignore")
LUMI = f"{DATA}/round6l_lumi.json"
LUMI31 = f"{DATA}/round6l_lumi_nnpdf31.json"


def refresh_lumi():
    base = "/usr/local/share/Pythia8/pdfdata/"
    for f, out in (("NNPDF23_lo_as_0130_qed_0000.dat", LUMI), ("NNPDF31_lo_as_0118_0000.dat", LUMI31)):
        if os.path.exists(base + f):
            p = C.PDF(base + f)
            tab = C.lumi_table(p)
            tab["sum_rules"] = [p.sum_rules(100.0), p.sum_rules(1e4)]
            json.dump(tab, open(out, "w"))


def dust_path(C_, D, n_pts=1500):
    r = U.closed_dust(C_, D)
    a_b, a_max = r["a_b"], r["a_max"]
    a = a_b + (a_max - a_b) * np.sin(np.linspace(1e-6, np.pi / 2, n_pts)) ** 2
    F = np.maximum(U._adot2(a, C_, D, 1.0, "hayward", 0.0), 1e-300)
    tau = np.concatenate([[0], np.cumsum(np.diff(a) / np.sqrt(0.5 * (F[1:] + F[:-1])))])
    return dict(tau=tau.tolist(), a=a.tolist(), N_acc=r["N_acc"], a_infl=r["a_infl"])


if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round6l.json"))
else:
    refresh_lumi()
    L = C.load_lumi(LUMI)
    L31 = C.load_lumi(LUMI31)
    out = {}
    # ---------------- P1: model U
    out["thresholds"] = [U.model_u_thresholds(n, 1.0, k) for n in (1, 2, 3, 4, 5, 6, 7) for k in ("hayward", "dymnikova")]
    out["add_scales"] = [U.add_scales(n, md) for n in (4, 5, 6, 7) for md in (6.0, 8.0, 12.0)]
    out["static_example"] = {str(x): dict(r=np.geomspace(0.05, 6, 300).tolist(),
                                          f=U.static_f(np.geomspace(0.05, 6, 300), x * U.extremal_mu(10), 10).tolist()) for x in (0.5, 1.0, 3.0, 20.0)}
    out["closed_dust"] = [U.closed_dust(c, D) for D in (4, 5, 10) for c in np.geomspace(30, 1e12, 12)]
    out["closed_dust_4D_stellar"] = U.closed_dust(2e38, 4)
    out["closed_dust_4D_smbh"] = U.closed_dust(2e47, 4)
    out["dust_paths"] = {"1e4": dust_path(1e4, 10), "1e10": dust_path(1e10, 10)}
    out["universality"] = [dict(kind=k, **{q: U.closed_dust(1e6, 10, 1.0, k)[q] for q in ("a_b", "a_max", "N_acc")}) for k in ("hayward", "dymnikova")]
    out["vacuum_threshold"] = [U.vacuum_threshold(c, D) for D in (4, 10) for c in (30.0, 1e2, 1e3, 1e4, 1e6)]
    out["vacuum_collapse"] = [U.vacuum_collapse_interior(x * U.extremal_mu(10), 10) for x in (0.8, 1.5, 5.0, 50.0)]
    out["temperature_curves"] = [C.temperature_curves(n, k) for n in (4, 6, 7) for k in ("hayward", "dymnikova")]
    # ---------------- P2: collider
    out["pdf_sum_rules"] = json.load(open(LUMI))["sum_rules"]
    out["trapped_fraction"] = {n: C.trapped_fraction_headon(n) for n in range(8)}
    Mth = np.linspace(5, 95, 91)
    out["sigma_curves"] = {name: dict(sqrt_s=rs, L_ab=lab, M_th=Mth.tolist(),
                                      cons=[C.sigma_pp(L, rs, 6, 8.0, m, "conservative") for m in Mth],
                                      opt=[C.sigma_pp(L, rs, 6, 8.0, m, "optimistic") for m in Mth],
                                      cons_RS=[C.sigma_pp(L, rs, 6, 8.0, m, "conservative", "RS") for m in Mth],
                                      cons_nnpdf31=[C.sigma_pp(L31, rs, 6, 8.0, m, "conservative") for m in Mth])
                           for name, (rs, lab) in C.MACHINES.items()}
    scan = []
    for n in (4, 5, 6, 7):
        for md in (6.0, 8.0, 10.0, 12.0, 15.0):
            for cl in (0.5, 0.6, 0.7, 0.8, 1.0):
                for k in ("hayward", "dymnikova"):
                    th = C.thresholds(n, md, cl, k)
                    row = dict(n=n, M_D_TeV=md, c_l=cl, kind=k, M_ext_TeV=th["model_U_TeV"], M_fmm_TeV=th["fmm_TeV"])
                    for name in ("FCC-hh 85", "SPPC/FCC 100"):
                        rs, lab = C.MACHINES[name]
                        for tr in ("conservative", "optimistic"):
                            sg = C.sigma_pp(L, rs, n, md, th["model_U_TeV"], tr)
                            row[f"{name}|{tr}|sigma_pb"] = sg
                            row[f"{name}|{tr}|events"] = C.events(sg, lab)
                    scan.append(row)
    out["scan"] = scan
    out["c_l_max"] = [dict(n=n, M_D_TeV=md, kind=k, machine=name, trapped=tr,
                           c_l_max=C.c_l_max_reachable(L, n, md, C.MACHINES[name][0], C.MACHINES[name][1], 10.0, k, tr))
                      for n in (4, 5, 6, 7) for md in (5.0, 6.0, 8.0, 10.0, 12.0, 15.0) for k in ("hayward", "dymnikova")
                      for name in ("FCC-hh 85", "SPPC/FCC 100") for tr in ("conservative", "optimistic")]
    out["c_l_min_lhc"] = [dict(n=n, M_D_TeV=md, kind=k, c_l_min=C.c_l_min_lhc(L, n, md, kind=k),
                               M_ext_at_min_TeV=C.thresholds(n, md, C.c_l_min_lhc(L, n, md, kind=k), k)["model_U_TeV"])
                          for n in (4, 5, 6, 7) for md in (5.0, 6.0, 8.0, 10.0, 12.0, 15.0) for k in ("hayward", "dymnikova")]
    out["evaporation"] = [C.evaporation(n, 8.0, 0.7, 15.0, k) for n in (4, 6, 7) for k in ("hayward", "dymnikova", "gr")]
    out["formation_vs_evaporation"] = [C.formation_vs_evaporation(n, md, cl, x, k) for n in (4, 6, 7) for md in (6.0, 10.0)
                                       for cl in (0.6, 0.8) for x in (1.5, 3.0) for k in ("hayward", "dymnikova")]
    # ---------------- P3: daughters
    out["bounce_lemma"] = U.bounce_untrapped_lemma()
    out["ecsk_deficit"] = [U.ecsk_os_mass_deficit(1.0, s, c) for s in (1e-9, 1e-6, 1e-3) for c in (0.3, np.pi / 4, 1.2)]
    out["two_epoch"] = [U.two_epoch_infall(D, 1.0, mue, E, dl) for (D, mue, E) in ((10, 20.0, 0.7), (5, 10.0, 0.5), (10, 3.0, 0.9))
                        for dl in (2.0, 10.0, 50.0, 500.0)]
    out["dtr_critical"] = [U.dtr_critical_delay(D, 1.0, mue, q, E) for (D, mue, E) in ((10, 20.0, 0.7), (5, 10.0, 0.5))
                           for q in (1e-1, 1e-3, 1e-6, 1e-12)]
    mu1 = 20 * U.extremal_mu(10)
    dtr_rows = []
    for dl in np.geomspace(1.0, 60.0, 40):
        m = U.two_epoch_infall(10, 1.0, 20.0, 0.7, dl)["meeting"]
        c = U.dtr_crossing(mu1, 1e-3 * mu1, m["r"], 10, log10_gap=m["log10_rminus_minus_r"] if m["log10_rminus_minus_r"] < -6 else None)
        dtr_rows.append(dict(delay=float(dl), **c))
    out["dtr_scan"] = dtr_rows
    out["fmm_patches"] = [U.fmm_patch_separation(rs) for rs in (3.0, 30.0, 1e3, 1e6, 1e38)]
    # ---------------- P4
    out["discriminants"] = C.discriminants()
    dump("round6l.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
ax = ax.ravel()
# (a) model U statics and T(M)
for key, col in (("0.5", MUTED), ("1.0", C4), ("3.0", C1), ("20.0", C2)):
    s = out["static_example"][key]
    ax[0].plot(s["r"], s["f"], color=col, label=f"mu = {key} mu_ext")
ax[0].axhline(0, color=INK, lw=0.7)
ax[0].plot(np.linspace(0, 1.4, 50), 1 - np.linspace(0, 1.4, 50) ** 2, ls=":", color=INK, lw=1, label="de Sitter core 1 - r^2/l^2")
ax[0].set_ylim(-1.2, 1.05); ax[0].set_xlim(0, 6)
ax[0].set_xlabel("areal radius r / l"); ax[0].set_ylabel("f(r)")
ax[0].set_title("(a) model U (D = 10, Hayward tower): regular holes with a static de Sitter core;\nhorizons only above mu_ext, where T -> 0 (a remnant)", fontsize=8)
ax[0].legend(frameon=False, fontsize=7, loc="lower right")
# (b) closed-dust interiors: bounce, transient acceleration, recollapse
for key, col in (("1e4", C1), ("1e10", C2)):
    p = out["dust_paths"][key]
    ax[1].semilogy(p["tau"], p["a"], color=col, label=f"C = {key}: N_acc = {p['N_acc']:.2f} e-folds")
    ax[1].axhline(p["a_infl"], color=col, ls=":", lw=1)
ax[1].set_xlabel("proper time from the bounce (units of l)"); ax[1].set_ylabel("scale factor a / l")
st = out["closed_dust_4D_stellar"]
ax[1].set_title(f"(b) matter interiors bounce at a = l, accelerate for N = ((D-3)/(D-1)) ln(a_max/l) e-folds,\n"
                f"then decelerate and recollapse (4D stellar hole: N = {st['N_acc']:.1f}); lasting inflation needs vacuum energy", fontsize=8)
ax[1].legend(frameon=False, fontsize=7)
# (c) collider cross-sections
cols = {"HL-LHC 14": MUTED, "HE-LHC 27": C4, "FCC-hh 85": C1, "SPPC/FCC 100": C2}
for name, col in cols.items():
    s = out["sigma_curves"][name]
    m = np.array(s["M_th"]); cons = np.array(s["cons"]); opt = np.array(s["opt"])
    ax[2].fill_between(m, np.maximum(cons, 1e-12), np.maximum(opt, 1e-12), color=col, alpha=0.25)
    ax[2].semilogy(m, np.where(cons > 0, cons, np.nan), color=col, label=f"{name} TeV ({s['L_ab']} ab^-1)")
    ax[2].axhline(10 / (s["L_ab"] * 1e6), color=col, ls=":", lw=0.8)
for cl, ls in ((0.6, "--"), (0.7, "-."), (0.8, ":")):
    th = C.thresholds(6, 8.0, cl)["model_U_TeV"]
    ax[2].axvline(th, color=INK, ls=ls, lw=0.9)
    ax[2].text(th + 0.6, 1e-8, f"M_ext (c_l = {cl})", rotation=90, fontsize=6)
ax[2].set_ylim(1e-9, 1e4)
ax[2].set_xlabel("black-hole mass threshold (TeV)"); ax[2].set_ylabel("sigma(pp -> BH, M > threshold) (pb)")
ax[2].set_title("(c) production with NNPDF2.3 LO (n = 6, M_D = 8 TeV): band = trapped mass\n0.59-1 x sqrt(s_hat); dotted = 10 events; black = model-U daughter thresholds", fontsize=8)
ax[2].legend(frameon=False, fontsize=7, loc="upper right")
# (d) parameter space: largest c_l giving >= 10 daughters at 100 TeV
for n, col in ((4, C1), (5, C3), (6, C2), (7, C4)):
    for k, ls in (("hayward", "-"), ("dymnikova", "--")):
        rows = [r for r in out["c_l_max"] if r["n"] == n and r["kind"] == k and r["machine"] == "SPPC/FCC 100" and r["trapped"] == "conservative"]
        lo = [r for r in out["c_l_min_lhc"] if r["n"] == n and r["kind"] == k]
        ax[3].plot([r["M_D_TeV"] for r in rows], [r["c_l_max"] for r in rows], ls=ls, color=col, marker="o", ms=3, label=f"n = {n}, {k}")
        ax[3].plot([r["M_D_TeV"] for r in lo], [r["c_l_min"] for r in lo], ls=ls, color=col, lw=1, alpha=0.6)
ax[3].axhline(1.0, color=INK, lw=0.7)
ax[3].set_xlabel("fundamental scale M_D (TeV)"); ax[3].set_ylabel("largest limiting length c_l = l M_D")
ax[3].text(5.2, 0.585, "thin lines: LHC lower edge (no black holes above M_ext seen at 13 TeV)", fontsize=6.5)
ax[3].set_title("(d) the window for the unified model: thick = largest l M_D giving >= 10 daughter-bearing\nholes at 100 TeV (30 ab^-1, conservative); thin = smallest allowed by the LHC", fontsize=8)
ax[3].legend(frameon=False, fontsize=6, ncol=2)
# (e) multiple daughters: DTR demand vs the curvature bound
d = out["dtr_scan"]
dl = np.array([r["delay"] for r in d]); req = np.array([r["f_F_required"] for r in d]); fmin = np.array([r["f_F_min_model_U"] for r in d])
ax[4].semilogx(dl, np.sign(req) * np.log10(1 + np.abs(req)), color=C1, label="f in the crossing region demanded by DTR")
ax[4].semilogx(dl, np.sign(fmin) * np.log10(1 + np.abs(fmin)), color=C2, ls="--", label="lowest f model U allows (psi = 1/alpha)")
dc = [r for r in out["dtr_critical"] if r["D"] == 10 and abs(r["mu2_over_mu1"] - 1e-3) < 1e-9][0]["delay_c"]
ax[4].axvline(dc, color=INK, ls=":", lw=1)
ax[4].text(dc * 1.05, -4, f"delay_c = {dc:.1f} l", fontsize=7)
ax[4].set_xlabel("advanced-time delay of the second infall (units of l)"); ax[4].set_ylabel("sign(f) log10(1 + |f|)")
ax[4].set_title("(e) a second infall meets the bounced first collapse inside the same core; after\n~7-30 l the crossing hits the curvature bound (locally de Sitter): no thin-shell solution", fontsize=8)
ax[4].legend(frameon=False, fontsize=7, loc="lower left")
# (f) status
prem = [("P1 gravity strong at M_D (n >= 4)", "allowed; tested by monojet + BH onset at 85-100 TeV", 1),
        ("P2 singularity resolved (model U)", "one consistent pure-gravity action with P1", 2),
        ("Daughter region in model U", "yes: bounce -> white hole -> new universe", 2),
        ("Inflating daughter", "transient only (N ~ ln); lasting needs vacuum energy", 0),
        ("Threshold reachable at 100 TeV", "iff M_ext ~ 9-50 TeV: a window in l M_D", 1),
        ("Forms before evaporating", "interior 1e-28 s; lifetime 1e-27-1e-24 s", 2),
        ("Second daughter from one hole", "no: joins core; late infall saturates the bound", 1),
        ("Seeing the daughter", "impossible; remnant + T(M) test the core class", 1)]
cm = matplotlib.colors.ListedColormap(["#f2a07e", "#f6dc8a", "#9fd8bd"])
ax[5].imshow(np.array([p[2] for p in prem])[:, None], cmap=cm, vmin=0, vmax=2, aspect="auto")
ax[5].set_yticks(range(len(prem))); ax[5].set_yticklabels([p[0] for p in prem], fontsize=7)
ax[5].set_xticks([]); ax[5].grid(False)
for i, p in enumerate(prem):
    ax[5].text(0, i, p[1], ha="center", va="center", fontsize=7, color=INK)
ax[5].set_title("(f) Round 6l status: green computed and satisfied, yellow conditional or open,\nred not achieved", fontsize=8)
fig.tight_layout()
fig.savefig(f"{FIG}/round6l.png", dpi=150, facecolor="#fcfcfb")
print("done")
