"""Round 8: model U against current data (LHC monojet and black-hole searches, the M_D convention, neutron-star safety) and future
machines.  Writes DATA/round8.json and FIGURES/round8.png."""
import json
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import data8 as E
from babyuniverse import collider6l as C
from babyuniverse import unified6l as U

warnings.simplefilter("ignore")

if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round8.json"))
else:
    L = C.load_lumi(f"{DATA}/round6l_lumi.json")
    r7 = json.load(open(f"{DATA}/round7.json"))
    out = {}
    out["convention"] = E.convention()
    out["monojet"] = E.monojet_bounds()
    Ms = np.linspace(3, 14, 45)
    out["lhc_sigma"] = dict(M=Ms.tolist(), cons=[C.sigma_pp(L, 13.0, 6, 8.0, m, "conservative") * 1e3 for m in Ms],
                            opt=[C.sigma_pp(L, 13.0, 6, 8.0, m, "optimistic") * 1e3 for m in Ms])
    out["lhc_lower_edge"] = [E.lhc_lower_edge(L, n, md, k, 0.5, tr) for n in (4, 6, 7) for md in (3.0, 8.0) for k in ("hayward",)
                             for tr in ("optimistic", "conservative")]
    out["multiobject_acceptance"] = {k: E.multiobject_search_acceptance(v["raw"]) for k, v in r7["benchmarks"].items() if k != "separation"}
    out["remnant_safety"] = [E.remnant_safety(*a) for a in ((6, 8.0, 0.6, "hayward"), (4, 6.0, 0.6, "hayward"),
                                                            (6, 8.0, 0.7, "dymnikova"), (7, 3.0, 0.7, "hayward"))]
    # luminosity-vs-excess curve for the safety figure
    me = U.extremal_mu(10)
    from babyuniverse.pipeline6o import a_n
    l = 0.6 / 8000.0
    rc_fac = (4.5) ** (1 / 7) * np.sqrt(9 / 7)
    curve = []
    for d in np.geomspace(1e-11, 1e-1, 60):
        tp = U.temperature(me * (1 + d), 10)
        if tp["r_plus"] is None:
            continue
        T = tp["T"] / l
        curve.append((float(d), float((np.pi**2 / 120) * 10 * np.pi * (rc_fac * tp["r_plus"] * l) ** 2 * T**4 / 6.582119569e-25)))
    out["safety_curve"] = curve
    out["future_reach"] = E.future_reach(L)
    out["allowed_region"] = {k: E.allowed_region(L, 6, k) for k in ("hayward", "dymnikova")}
    dump("round8.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
ax = ax.ravel()
mj = out["monojet"]
x = np.arange(len(mj))
ax[0].bar(x - 0.2, [m["MD_PDG_TeV"] for m in mj], 0.4, color=MUTED, label="CMS monojet bound, PDG convention")
ax[0].bar(x + 0.2, [m["MD_DL_TeV"] for m in mj], 0.4, color=C1, label="same bound, Dimopoulos-Landsberg convention (this work)")
ax[0].set_xticks(x); ax[0].set_xticklabels([f"n = {m['n']}" for m in mj]); ax[0].set_ylabel("lower bound on M_D (TeV)")
ax[0].set_title("(a) the M_D convention: experiments quote PDG, our formulas use DL;\nM_D(DL) = 0.89-0.38 M_D(PDG) for n = 2-6 (correction F64)", fontsize=8)
ax[0].legend(frameon=False, fontsize=7)
ls = out["lhc_sigma"]
ax[1].semilogy(ls["M"], np.maximum(ls["opt"], 1e-9), color=C2, label="13 TeV, full trapping")
ax[1].semilogy(ls["M"], np.maximum(ls["cons"], 1e-9), color=C1, label="13 TeV, conservative trapping (y = 0.59)")
ax[1].axhline(E.ATLAS_MONOJET_LIMIT_FB / 0.5, color=INK, ls="--", lw=1, label="ATLAS monojet limit 0.3 fb / efficiency 0.5")
ax[1].set_ylim(1e-6, 1e4)
ax[1].set_xlabel("model-U threshold M_ext (TeV)"); ax[1].set_ylabel("sigma(pp -> model-U hole) (fb)")
ax[1].set_title("(b) LHC: model-U holes sit in the monojet high-recoil region; current data exclude\nM_ext below ~8-9 TeV (full trapping) or ~4.7-5.5 TeV (conservative)", fontsize=8)
ax[1].legend(frameon=False, fontsize=7)
acc = out["multiobject_acceptance"]
names = list(acc.keys())
ax[2].bar(range(len(names)), [acc[k] for k in names], color=[C1, C1, C3, C2, C2])
ax[2].set_xticks(range(len(names))); ax[2].set_xticklabels(names, rotation=20, fontsize=7)
ax[2].set_ylabel("fraction with S_T > 5 TeV and >= 3 objects")
ax[2].set_title("(c) the multi-object black-hole searches see GR holes (~90 %) but almost none of the\nnear-threshold model-U holes: their limits do not apply to model U", fontsize=8)
sc = np.array(out["safety_curve"])
rs = out["remnant_safety"][0]
ax[3].loglog(sc[:, 0], sc[:, 1], color=C4, label="Hawking re-emission of the excess above extremality")
ax[3].axhline(rs["accretion_GeV_per_s"], color=C2, ls="--", label="geometric accretion in a neutron star")
ax[3].axvline(rs["equilibrium_excess_fraction"], color=INK, ls=":", lw=1)
ax[3].set_xlabel("excess above extremality (M - M_ext)/M_ext"); ax[3].set_ylabel("power (GeV/s)")
ax[3].set_title(f"(d) safety: a captured remnant re-radiates what it accretes, sitting {rs['equilibrium_excess_fraction']:.0e} above\nextremality at T = {rs['equilibrium_T_GeV']:.2f} GeV: no growth", fontsize=8)
ax[3].legend(frameon=False, fontsize=7)
fr = [r for r in out["future_reach"] if r["N_required"] == 10.0]
ax[4].barh(range(len(fr)), [r["M_ext_reach_TeV"] for r in fr], color=[MUTED, MUTED, C4, C1, C2])
ax[4].axvspan(6.5, 9.3, color="#f2a07e", alpha=0.4, label="current LHC lower edge (trapping-dependent)")
ax[4].set_yticks(range(len(fr))); ax[4].set_yticklabels([r["machine"] for r in fr])
ax[4].set_xlabel("largest M_ext with >= 10 model-U holes (conservative trapping), TeV")
ax[4].set_title("(e) future: the HL-LHC cannot reach the window; a 27 TeV machine reaches ~13 TeV;\nFCC-hh 85 / 100 TeV cover M_ext up to ~41 / 48 TeV", fontsize=8)
ax[4].legend(frameon=False, fontsize=7, loc="lower right")
for k, col in (("hayward", C1), ("dymnikova", C3)):
    ar = out["allowed_region"][k]
    md = [a["M_D_DL_TeV"] for a in ar]
    ax[5].fill_between(md, [a["c_l_min_lhc"] for a in ar], [a["c_l_max_100"] for a in ar], color=col, alpha=0.25, label=f"{k}: testable at 100 TeV, allowed by LHC")
    ax[5].plot(md, [a["c_l_max_85"] for a in ar], color=col, ls="--", lw=1)
mono6 = [m for m in out["monojet"] if m["n"] == 6][0]["MD_DL_TeV"]
ax[5].axvspan(0, mono6, color="#f2a07e", alpha=0.5, label="excluded by monojet (n = 6)")
ax[5].set_xlim(1.5, 15)
ax[5].set_xlabel("M_D (Dimopoulos-Landsberg convention, TeV)"); ax[5].set_ylabel("limiting length c_l = l M_D")
ax[5].set_title("(f) the surviving, testable region (n = 6): allowed by all current data and covered\nby a 100 TeV collider (dashed: 85 TeV upper edge)", fontsize=8)
ax[5].legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(f"{FIG}/round8.png", dpi=150, facecolor="#fcfcfb")
print("done")
