"""Round 6g: (A) the false-vacuum region made by the universe itself (a supercooled patch left behind by a first-order transition in a
radiation-dominated universe, evolved with gravity); (B) whether any machine humanity could build can make it.  Writes DATA/round6g.json
and FIGURES/round6g.png; `--figure-only` redraws the figure."""
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from _common import FIG, DATA, dump
from babyuniverse import formation as F, throat_evolution as TE, gsl_solution as GS, accessibility as AC

warnings.simplefilter("ignore")
RHO_F = 3 / 900 / (8 * np.pi)
V = 0.0015


def frw_check(a):
    gauge, h = a
    lam = 4 * RHO_F / (4 * V**4)
    pot = TE.Potential(lam, V, RHO_F)
    pl = F.Plasma(RHO_F, V, lam, f_H=1e-12, T_floor_over_Tc=1e-12)
    eta0 = 20.0
    rho0 = 3 / (8 * np.pi * eta0**2)
    Tout = brentq(lambda T: pl.eos(T, V)[0] - rho0, 1e-6, 5)
    Z0 = F.frw_cone_data(h, 40.0 if gauge == "conformal" else 30.0, pl, pot, eta0, T_out0=Tout, gauge=gauge)
    res = F.evolve_centre(Z0, h, pot, pl, n_rays=int((30 if gauge == "conformal" else 20) / h), clip=0)
    R, S = res["r"], res["s"]
    if gauge == "conformal":
        i, j = np.indices(R.shape)
        eta = eta0 + 0.5 * (i + j) * h
        chi = 0.5 * (j - i) * h
        m = np.isfinite(R) & (chi > 0.5)
        return dict(gauge=gauge, h=h, r_err=float(np.nanmax(np.abs(R[m] / ((eta[m] / eta0) * chi[m]) - 1))),
                    s_err=float(np.nanmax(np.abs(S[m] - np.log(eta[m] / eta0)))))
    ru, rv = np.gradient(R, h, axis=0), np.gradient(R, h, axis=1)
    msh = 0.5 * R * (1 + 4 * ru * rv * np.exp(-2 * S))
    rho = pl.eos(Tout, V)[0] * (np.exp(res["lnTf"]) / Tout) ** 4
    ok = np.isfinite(R) & (R > 1.0)
    ok[:, :3] = False
    ok[:3, :] = False
    rel = np.abs(msh[ok] / (4 * np.pi / 3 * rho[ok] * R[ok] ** 3) - 1)
    return dict(gauge=gauge, h=h, misner_sharp_median=float(np.median(rel)), misner_sharp_p95=float(np.percentile(rel, 95)))


def patch_job(a):
    f, x, h, lam_f, eps = a
    o = F.cosmic_patch(f, x, h=h, lam_f=lam_f, eps=eps)
    cr = F.constraint_residual(o["res"], o["plasma"], r_lo=1.0, r_hi=1e3, s_max=8.0)
    out = {k: o[k] for k in ("eta0", "chi_p", "chi_p_over_lH", "outcome", "r_fv_anti_max", "trapped_cells", "anti_fv_cells",
                             "rows_with_daughter_behind_trapped_true_vacuum", "first_hidden_u")}
    out.update(f=f, chi_p_over_eta0=x, h=h, lam_f=lam_f, eps=eps, Evv_rel=cr["Evv_rel"], Euu_rel=cr["Euu_rel"])
    if (f, x, h, lam_f, eps) == (0.5, 0.8, 0.05, 0.0, 4.0):              # the hidden daughter: GSL and a map for the figure
        g = dict(res=o["res"], cls=o["cls"], r_N=o["r_N"], pot=o["pot"])
        eh = GS.parent_event_horizon(g)
        dh = GS.daughter_horizons(g, n_points=6)
        out.update(u_EH=eh.get("u_EH"), parent_deficit=eh.get("deficit"), daughter_points=dh.get("points"), daughter_max_deficit=dh.get("max_deficit"))
        R, P = o["res"]["r"], o["res"]["phi"]
        out["map"] = dict(cls=np.nan_to_num(o["cls"][::4, ::4], nan=9).astype(int).tolist(),
                          fv=(np.abs(P[::4, ::4] / V + 1) < 0.5).astype(int).tolist(), h=4 * h)
    # false-vacuum radius against conformal time for the figure
    R, P = o["res"]["r"], o["res"]["phi"]
    ok = np.isfinite(R)
    fvm = ok & (np.abs(P / V + 1) < 0.5)
    i, j = np.indices(R.shape)
    eta = o["eta0"] + 0.5 * (i + j) * h
    hist = []
    for e in np.linspace(o["eta0"], o["eta0"] * 5, 33):
        m = ok & (np.abs(eta - e) < h)
        if m.any():
            hist.append((float(e / o["eta0"]), float(R[m & fvm].max()) if (m & fvm).any() else 0.0))
    out["hist"] = hist
    return out


if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round6g.json"))
else:
    out = {}
    with ProcessPoolExecutor(4) as ex:
        out["frw_verification"] = list(ex.map(frw_check, [("conformal", h) for h in (0.2, 0.1, 0.05)] + [("affine", h) for h in (0.1, 0.05)]))
        jobs = [(0.9, x, 0.05, 0.0, 4.0) for x in (0.5, 0.58, 0.6, 0.62, 0.64, 0.7, 1.0, 2.0)] + \
               [(0.9, x, 0.1, 0.0, 4.0) for x in (0.58, 0.6, 0.62, 0.64)] + \
               [(0.5, x, 0.05, 0.0, 4.0) for x in (0.6, 0.8, 1.0)] + [(1.0, x, 0.05, 0.0, 4.0) for x in (0.55, 0.65)] + \
               [(0.9, 0.6, 0.05, 0.0, 2.0), (0.9, 0.6, 0.05, 1e-6, 2.0), (0.9, 0.56, 0.05, 1e-6, 2.0), (0.9, 0.6, 0.05, 1e-6, 4.0)]
        out["patches"] = list(ex.map(patch_job, jobs))
    crit = {}
    for f in (0.9, 0.5, 1.0):
        for h in (0.05, 0.1):
            rows = sorted([p for p in out["patches"] if p["f"] == f and p["h"] == h and p["lam_f"] == 0 and p["eps"] == 4.0],
                          key=lambda p: p["chi_p_over_lH"])
            lo = max([p["chi_p_over_lH"] for p in rows if p["outcome"] == "collapses"], default=None)
            hi = min([p["chi_p_over_lH"] for p in rows if p["outcome"] == "inflates"], default=None)
            if rows:
                crit[f"f={f}, h={h}"] = (lo, hi)
    out["critical_radius_over_lH"] = crit
    out["route_ball_radius_over_lH"] = 0.595
    out["machines"] = AC.machine_audit()
    out["tev_gravity"] = AC.tev_gravity_conditional()
    dump("round6g.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(1, 4, figsize=(19, 4.6))
P_ = out["patches"]
# (a) false-vacuum radius against time: collapsing and inflating patches
for p, col in ((next(q for q in P_ if (q["f"], q["chi_p_over_eta0"], q["h"], q["lam_f"]) == (0.9, 0.6, 0.05, 0.0) and q["eps"] == 4.0), C1),
               (next(q for q in P_ if (q["f"], q["chi_p_over_eta0"], q["h"], q["lam_f"]) == (0.9, 0.62, 0.05, 0.0) and q["eps"] == 4.0), C2),
               (next(q for q in P_ if (q["f"], q["chi_p_over_eta0"], q["h"], q["lam_f"]) == (0.9, 1.0, 0.05, 0.0) and q["eps"] == 4.0), C3)):
    hh = [x for x in p["hist"]]
    if p["outcome"] == "inflates":                      # the computed domain ends where the inflating interior becomes stiff (F45)
        hh = hh[: int(np.argmax([x[1] for x in hh])) + 1]
    ax[0].plot([x[0] for x in hh], [max(x[1], 1e-1) for x in hh], "o-", color=col, ms=4,
               label=f"patch {p['chi_p_over_lH']:.2f} l_H: {p['outcome']}")
ax[0].set_yscale("log")
ax[0].set_xlabel("conformal time / initial Hubble time"); ax[0].set_ylabel("largest radius still in the false vacuum")
ax[0].set_title("(a) a patch left behind by the transition: collapse or inflation\n(inflating curves end where the computed domain ends)", fontsize=8)
ax[0].legend(frameon=False, fontsize=7)
# (b) map of the hidden daughter (f = 0.5, patch 0.57 l_H)
hp = next(q for q in P_ if "map" in q)
cls = np.array(hp["map"]["cls"], float)
fvm = np.array(hp["map"]["fv"], float)
img = np.where(cls == 9, np.nan, cls)
cmap = matplotlib.colors.ListedColormap(["#5e5e5a", "#f4f3ef", "#f2a07e"])
hm = hp["map"]["h"]
ax[1].imshow(img, extent=[0, img.shape[1] * hm, img.shape[0] * hm, 0], aspect="auto", cmap=cmap, vmin=-1, vmax=1)
ax[1].contour(np.arange(img.shape[1]) * hm, np.arange(img.shape[0]) * hm, fvm, levels=[0.5], colors=C1, linewidths=1.2)
ax[1].set_xlabel("v"); ax[1].set_ylabel("u")
ax[1].set_title("(b) hidden daughter (f = 0.5, patch 0.57 l_H): false vacuum inside the blue line,\nanti-trapped (orange), trapped black-hole region (dark); parent horizon at u = %.1f" % hp["u_EH"], fontsize=8)
# (c) outcome against patch radius
for f, col, mk in ((0.9, C1, "o"), (0.5, C2, "s"), (1.0, C3, "^")):
    rows = [p for p in P_ if p["f"] == f and p["h"] == 0.05 and p["lam_f"] == 0 and p["eps"] == 4.0]
    ax[2].plot([p["chi_p_over_lH"] for p in rows], [1 if p["outcome"] == "inflates" else 0 for p in rows], mk, color=col, ms=8,
               label=f"vacuum fraction f = {f}")
ax[2].axvline(0.595, color=INK, lw=0.8, ls="--")
ax[2].text(0.605, 0.82, "the route's ball\n(0.595 l_H)", fontsize=7, color=INK)
ax[2].set_yticks([0, 1]); ax[2].set_yticklabels(["collapses\n(black hole)", "inflates\n(daughter)"])
ax[2].set_xscale("log")
ax[2].set_xlabel("patch radius / false-vacuum Hubble length")
ax[2].set_title("(c) critical size: near-critical patches are the route's ball", fontsize=9)
ax[2].legend(frameon=False, fontsize=7, loc="center left", bbox_to_anchor=(0.45, 0.4))
# (d) the power ladder
asm = out["machines"]["assembly"]
names = ["humanity (all primary power)", "peak of GW150914 (black-hole merger)", "max luminosity c^5/4G (conjecture)", "needed to build the ball"]
vals = [AC.WORLD_POWER, asm["GW150914_peak_W"], asm["max_luminosity_W"], asm["power_W"]]
cols = [C1, C3, MUTED, C2]
ax[3].barh(range(len(names)), np.log10(vals), color=cols)
for k, val in enumerate(vals):
    ax[3].text(np.log10(val) + 0.5, k, f"{val:.1e} W", va="center", fontsize=7, color=INK)
ax[3].set_yticks(range(len(names))); ax[3].set_yticklabels(names, fontsize=7)
ax[3].set_xlim(0, 60)
ax[3].set_xlabel("log10 power (W)")
ax[3].set_title(f"(d) building the ball needs 10^{asm['log10_power_gap']:.0f} x humanity's power\n(Theorem 6g-2; any machine, any vacuum scale)", fontsize=8)
fig.tight_layout()
fig.savefig(f"{FIG}/round6g.png", dpi=150, facecolor="#fcfcfb")
print("done")
