"""Round 6j: what the bridge invariant actually bounds (geometry, stored energy, supplied energy), and fifteen loopholes in four groups
(classical gravity, environment-prepared states, quantum spacetime, different gravity).  Writes DATA/round6j.json and FIGURES/round6j.png;
`--figure-only` redraws the figure."""
import json
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import accessibility6j as A, loopholes6j as L

warnings.simplefilter("ignore")

if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round6j.json"))
else:
    out = {}
    # ---------------------------------------------------------------- Part I
    out["identity"] = [A.compactness_identity(E, 1.0) for E in (1e20, 1e40, A.C_SI**4 / (2 * A.G_SI))]
    out["misner_sharp"] = A.misner_sharp_identity()
    out["which_density"] = A.which_density()
    out["three_bounds"] = [A.three_bounds(R) for R in (1e-18, 1e-3, 1.0)]
    out["precursor"] = [A.precursor_window(), A.precursor_window(M_kg=A.M_SUN), A.precursor_window(M_kg=9e8)]
    out["qei"] = [A.qei_negative_energy(R) for R in (1e-30, 1e-18, 1e-3)]
    # ---------------------------------------------------------------- Group A
    Fs = np.geomspace(1e-60, 0.99, 120)
    out["strong_gravity"] = [A.strong_gravity_core(F) for F in Fs]
    out["strong_gravity_threshold"] = A.strong_gravity_threshold()
    out["dirichlet_check"] = A.dirichlet_check()
    out["pocket"] = [A.pocket_expansion(v, R) for v, R in ((500, 1e-18), (1e4, 1e-19), (1e10, 1e-25), (1.5e16, 1e-31))]
    # ---------------------------------------------------------------- Groups B-D
    out["accretion"] = [L.accretion_requirement(m) for m in (1.4, 10.0, 1e6)]
    out["qcd_bag_partner"] = L.qcd_bag_partner()
    out["bag_star"] = L.bag_star()
    out["oy"] = {str(kw): [L.oy_channel(f, 1.0, kw) for f in (0.02, 0.05, 0.08, 0.12, 0.15, 0.17, 0.185)] for kw in (0.5, 0.1)}
    out["oy_expected_gut"] = L.oy_expected_daughters(4e-4, 0.05, np.sqrt(5.2e10))
    out["oy_expected_planckian"] = [dict(L_over_lP=x, N=L.oy_expected_daughters(4e-4, 0.05, x)) for x in (10, 100, 300, 600, 1000)]
    out["wineglass"] = {str(Lam): L.wineglass_family(Lam, np.geomspace(1e-5, 0.5, 30)) for Lam in (1.0, 0.1, 0.01)}
    out["fmm"] = [L.fmm_junction(M) for M in (0.4, 0.6, 1.0, 10.0, 1e3, 1e38)]
    out["add"] = [L.add_black_holes(n, md, mbh) for n, md, mbh in ((2, 11.0, 40.0), (4, 7.5, 35.0), (6, 6.0, 30.0), (6, 8.0, 30.0),
                                                                   (7, 5.5, 25.0), (6, 15.0, 60.0), (6, 20.0, 80.0))]
    out["add_semiclassical"] = [dict(n=n, M_over_MD_for_RS_3=L.add_semiclassical_mass(n, 3.0), RS_MD_at_5MD=L.add_rs_over_fundamental(5.0, 1.0, n))
                                for n in (2, 3, 4, 5, 6, 7)]
    out["ledger"] = L.ledger6j()
    dump("round6j.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
ax = ax.ravel()
# (a) which density: peak over mean and positive energy needed, by profile
wd = out["which_density"]
x = np.arange(len(wd))
ax[0].bar(x - 0.2, [r["peak_over_mean"] for r in wd], 0.4, color=C1, label="peak density / areal mean")
ax[0].bar(x + 0.2, [r["positive_energy_over_R_half"] for r in wd], 0.4, color=C2, label="positive energy / (R c^4/2G)")
ax[0].axhline(1, color=INK, lw=0.8)
ax[0].set_xticks(x); ax[0].set_xticklabels([r["profile"] for r in wd], rotation=25, ha="right", fontsize=7)
ax[0].set_title("(a) at the first marginal sphere: non-uniformity raises the peak,\nnegative components raise the positive energy needed (Prop. 6j-1)", fontsize=8)
ax[0].legend(frameon=False, fontsize=7)
# (b) the strong-gravity core: what it costs against what it buys
sg = out["strong_gravity"]
F = [r["F"] for r in sg]
ax[1].loglog(F, [r["core_compactness"] for r in sg], color=C2, label="core's own compactness (3/4) ln^2 F")
ax[1].loglog(F, [r["invariant_gain"] for r in sg], color=C1, label="gain in the invariant F^{3/4}")
ax[1].axhline(1, color=INK, lw=0.8)
ax[1].axvline(out["strong_gravity_threshold"], color=MUTED, ls=":")
ax[1].axvline(10 ** (-39.5 / 0.75), color=C3, ls="--", lw=1)
ax[1].text(10 ** (-52.0), 3e-30, "gain needed\n(10^-39.5)", color=C3, fontsize=7)
ax[1].set_xlabel("F = (effective Planck mass / Planck mass)^2 in the core")
ax[1].set_title("(b) field-dependent G: the core is horizon-class for F < 0.32,\nlong before it buys anything (Prop. 6j-4)", fontsize=8)
ax[1].legend(frameon=False, fontsize=7, loc="center right")
# (c) Oshita-Yokoyama: total exponent shape factor versus hole mass
for kw, col in (("0.5", C1), ("0.1", C2)):
    rows = [r for r in out["oy"][kw] if r.get("exists")]
    ax[2].semilogy([r["M_in"] for r in rows], [max(r["Sigma"], 1e-6) for r in rows], "o-", color=col, label=f"wall k_w = {kw}")
ax[2].axvline(1 / (3 * np.sqrt(3)), color=MUTED, ls=":")
ax[2].text(0.15, 2e-1, "Nariai\nbound", fontsize=7, color=MUTED)
ax[2].set_xlabel("hole mass x chi"); ax[2].set_ylabel("Sigma_hat = Sigma / (l_f/l_P)^2")
ax[2].set_title("(c) evaporating-hole catalysis, both steps (thermal + tunnelling):\nGUT vacuum (l_f/l_P)^2 ~ 5e10 gives Sigma >~ 2e7", fontsize=8)
ax[2].legend(frameon=False, fontsize=7)
# (d) wineglass wormholes: action excess against throat area
for Lam, col in (("1.0", C1), ("0.1", C2), ("0.01", C3)):
    rows = out["wineglass"][Lam]["rows"]
    ax[3].loglog([r["a0"] ** 2 for r in rows], [r["dS"] for r in rows], "o", ms=3, color=col, label=f"Lambda = {Lam}")
a2 = np.geomspace(1e-5, 1, 50)
ax[3].loglog(a2, 3 * np.pi**3 * a2, color=INK, lw=0.8, label="3 pi^3 a_0^2")
ax[3].set_xlabel("throat radius^2 a_0^2 (8 pi G = 1)"); ax[3].set_ylabel("Delta S = S(q) - S_no-boundary")
ax[3].set_title("(d) Euclidean wormholes: staying connected costs action ~ throat area;\nthe throat-free limit is the disconnected no-boundary universe", fontsize=8)
ax[3].legend(frameon=False, fontsize=7)
# (e) accretion and precursors: densities against the Nariai partner
acc = out["accretion"]
labels = [f"{a['M_sun']:g} M_sun: Eddington flow" for a in acc[:2]] + [f"{a['M_sun']:g} M_sun: merger debris" for a in acc[:2]]
vals = [np.log10(a["rho_eddington_kg_m3"] / a["rho_f_kg_m3"]) for a in acc[:2]] + [np.log10(a["merger_over_required"]) for a in acc[:2]]
ax[4].barh(range(4), vals, color=[C1, C1, C3, C3])
ax[4].axvline(0, color=INK, lw=0.8)
ax[4].set_yticks(range(4)); ax[4].set_yticklabels(labels, fontsize=7)
ax[4].set_xlabel("log10 (available density / Nariai partner density)")
bs = out["bag_star"]
ax[4].set_title(f"(e) environment-prepared state: accretion falls short of the partner density;\nmost compact bag star: 2GM/c^2R = {bs['compactness']:.2f}, "
                f"R = {bs['R_over_l_f']:.2f} l_f (sub-Hubble)", fontsize=8)
# (f) the fifteen candidates
led = out["ledger"]
score = []
for name, _, verdict in led:
    score.append(2 if "CONDITIONAL" in verdict else (1 if any(w in verdict for w in ("open", "natural", "early", "toy", "Planck")) else 0))
cm = matplotlib.colors.ListedColormap(["#f2a07e", "#f6dc8a", "#9fd8bd"])
ax[5].imshow(np.array(score)[:, None], cmap=cm, vmin=0, vmax=2, aspect="auto")
ax[5].set_yticks(range(len(led))); ax[5].set_yticklabels([l[0] for l in led], fontsize=7)
ax[5].set_xticks([])
ax[5].grid(False)
for i, (_, _, verdict) in enumerate(led):
    ax[5].text(0, i, verdict, ha="center", va="center", fontsize=6.3, color=INK)
ax[5].set_title("(f) fifteen candidates: green = conditional door (new physics, testable),\nyellow = natural / open / Planck-scale only, red = closed", fontsize=8)
fig.tight_layout()
fig.savefig(f"{FIG}/round6j.png", dpi=150, facecolor="#fcfcfb")
print("done")
