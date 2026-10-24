"""Round 6c (P-R6f): full semiclassical Einstein-Maxwell-scalar-CFT evolution of the magnetic throat, and the obstruction to growth."""
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse import throat_evolution as TE, obstructions as O

warnings.simplefilter("ignore")
out = {}
Q = 10.0

# ------------------------------------------------------------------ 1. Reissner-Nordstrom validation
M, r0 = 10.1, 13.0
f0 = 1 - 2 * M / r0 + Q**2 / r0**2
val = {}
Rs = {}
for h in (0.1, 0.05):
    a, b, _ = TE.characteristic_data(Q, 0.0, r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, h, int(40 / h), int(40 / h))
    res = TE.evolve(a, b, Q, 0.0, TE.Potential(), h)
    Mm = TE.misner_sharp(res)
    val[f"h={h}"] = dict(M_dev=float(np.nanmax(np.abs(Mm - M))), **{k: float(v) for k, v in TE.constraint_residual(res).items()})
    Rs[h] = res["r"]
val["convergence_r_diff_h_vs_h/2"] = float(np.nanmax(np.abs(Rs[0.1] - Rs[0.05][::2, ::2])))
out["rn_validation"] = val

# ------------------------------------------------------------------ 2. QNEC flip, fully backreacted
Mx = 10.0001
rp, rm = Mx + np.sqrt(Mx**2 - Q**2), Mx - np.sqrt(Mx**2 - Q**2)
fq = 1e-4
from scipy.optimize import brentq
rq = brentq(lambda r: (r - rp) * (r - rm) / r**2 - fq, rp * (1 + 1e-12), 2 * Q)
c, Y, up, wp = 50.0, 0.1, 5.0, 0.5
yp = lambda u: 0.5 * Y * (1 + np.tanh((u - up) / wp))
ypp = lambda u: 0.5 * Y / wp / np.cosh((u - up) / wp) ** 2
tu = lambda u: c / (24 * np.pi) * (0.5 * yp(u) ** 2 - ypp(u))
flip = {}
for h in (0.1, 0.05):
    a, b, info = TE.characteristic_data(Q, c, rq, -np.sqrt(fq) / 2, np.sqrt(fq) / 2, h, int(60 / h), int(60 / h), t_u=tu)
    res = TE.evolve(a, b, Q, c, TE.Potential(), h)
    cl = TE.classify(res)
    u = np.arange(int(60 / h) + 1) * h
    flip[f"h={h}"] = dict(n_anti=cl["n_anti"], r_min=cl["r_min"], r_max=cl["r_max"],
                          antitrapped_u=[float(u[info["ru_b"] > 0].min()), float(u[info["ru_b"] > 0].max())],
                          **{k: float(v) for k, v in TE.constraint_residual(res, t_u=tu).items()})
    if h == 0.05:
        flip_map = (cl["cls"], res["r"])
flip["kick_over_need"] = float(c * Y / (24 * np.pi * rq) / (np.sqrt(fq) / 2))
out["qnec_flip"] = flip

# ------------------------------------------------------------------ 3-4. false vacuum from normal data: normality bound and growth scan
x, vphi, eps = 0.5, 0.005, 2.0
rho_f = x / (8 * np.pi * Q**2)
lam = 4 * rho_f / (eps * vphi**4)
w = np.sqrt(2) / (np.sqrt(lam) * vphi)
pot = TE.Potential(lam, vphi, rho_f)
lH = np.sqrt(3 * Q**2 / x)


def fv_data(r0, ru0, rv0, u2, v2, h, Lu, Lv, Nflux=0.0, upos=1.0, wpp=0.3):
    phi_v = lambda vv: vphi * (np.tanh((np.asarray(vv, float) - v2) / w))
    phi_u = lambda uu: vphi * (np.tanh((np.asarray(uu, float) - u2) / w))
    Nu = (lambda uu: -Nflux * np.exp(-((uu - upos) / wpp) ** 2)) if Nflux else None
    a, b, _ = TE.characteristic_data(Q, 0.0, r0, ru0, rv0, h, int(Lv / h), int(Lu / h), phi_v=phi_v, phi_u=phi_u, N_u=Nu)
    return TE.evolve(a, b, Q, 0.0, pot, h, r_inf=1e4)


def max_normal_v2(r0, f0, h=0.05):
    lo, hi = 0.5, 80.0
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        cl = TE.classify(fv_data(r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, 30.0, mid, h, 8, mid + 4))
        lo, hi = (mid, hi) if cl["n_anti"] == 0 else (lo, mid)
    return lo


normal = []
scan = []
for fa in (0.01, 0.05, 0.1, 0.3, 0.6):
    r0 = Q / (1 - np.sqrt(fa))
    C = r0 * (1 + Q**2 / r0**2 - x / Q**2 * r0**2 / 3 - fa)
    rc = brentq(lambda r: O.rnds_f(r, C, Q, x / Q**2), r0, 10 * r0)
    v_c = (rc - r0) / (np.sqrt(fa) / 2)
    v2n = max_normal_v2(r0, fa)
    normal.append(dict(f_a=fa, r0=float(r0), v_c_analytic=float(v_c), v2_numeric=float(v2n)))
    need = np.sqrt(fa) / 2 * r0 / (np.sqrt(np.pi) * 0.3)
    for k in (1.05, 1.5, 3.0, 10.0):
        res = fv_data(r0, -np.sqrt(fa) / 2, np.sqrt(fa) / 2, 30.0, 0.97 * v2n, 0.05, 70, 0.97 * v2n + 60, Nflux=k * need)
        P, R = res["phi"], res["r"]
        scan.append(dict(f_a=fa, flip_factor=k, reached_scri=bool(np.any(res["end"] == 2)),
                         growth=float(np.nanmax(np.where(P < 0, R, np.nan)) / r0),
                         late_false_fraction=float(np.nanmean(P[-1] < 0))))
out["normality_bound"] = normal
out["growth_from_normal_data"] = scan
out["segment_ratio_max"] = {f"x={xx}": O.growth_ratio_max(Q, xx, 1.01 * Q) for xx in (0.3, 0.5, 1.0)}

# ------------------------------------------------------------------ 5-6. super-Hubble anti-trapped false vacuum, and stability
sup = []
hist = {}
for r0s, fneg in ((15.0, -0.05), (15.0, -0.5), (30.0, -0.05), (30.0, -0.5)):
    for L in (10.0, 25.0, 50.0):
        rv0 = 0.3
        res = fv_data(r0s, -fneg / (4 * rv0), rv0, L, L, 0.05, L + 40, L + 40)
        P, R = res["phi"], res["r"]
        reached = bool(np.any(res["end"] == 2))
        sup.append(dict(r0=r0s, f0=fneg, L=L, L_over_lH=float(L / lH), reached_scri=reached,
                        r_max_false=float(np.nanmax(np.where(P < 0, R, np.nan)))))
        if r0s == 15.0 and fneg == -0.05:
            ends = np.where(res["end"] == 2)[0]
            hist[L] = (np.nanmax(np.where(P < 0, R, np.nan), axis=1), int(ends[0]) if len(ends) else None)
out["superhubble"] = sup
stab = []
for dL, df, dw in ((1.0, 1.0, 1.0), (0.9, 1.0, 1.0), (1.1, 1.0, 1.0), (1.0, 0.9, 1.0), (1.0, 1.1, 1.0), (1.0, 1.0, 1.5)):
    w_save = w
    w = w_save * dw
    for L in (25.0, 50.0):
        res = fv_data(15.0, 0.05 * df / 1.2, 0.3, L * dL, L * dL, 0.05, L * dL + 40, L * dL + 40)
        stab.append(dict(L=L, dL=dL, df=df, dw=dw, reached_scri=bool(np.any(res["end"] == 2))))
    w = w_save
out["stability"] = stab
dump("round6c_throat.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED,
                     "ytick.color": MUTED})
fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.2))
a0 = ax[0, 0]
cls, Rq = flip_map
cmap = matplotlib.colors.ListedColormap(["#9cc3ee", "#f4f3ef", "#f2a07e"])
a0.imshow(cls, extent=[0, 60, 60, 0], aspect="auto", cmap=cmap, vmin=-1, vmax=1)
a0.set_xlabel("v (outgoing affine time, throat units)")
a0.set_ylabel("u (ingoing)")
a0.text(30, 7.4, "anti-trapped after the QNEC pulse", color=INK, fontsize=8, ha="center")
a0.text(30, 35, "normal (throat persists, r in [9.7, 10.4])", color=INK, fontsize=8, ha="center")
a0.set_title("(a) Full backreacted flip: c = 50, Q = 10, f = 1e-4 (constraints ~1e-7)", fontsize=9)

a1 = ax[0, 1]
fs = np.geomspace(1e-4, 0.95, 50)
for xx, col in ((0.3, C1), (0.5, C2), (1.0, C3)):
    lHx = np.sqrt(3 * Q**2 / xx)
    a1.plot(fs, [O.segment_Lmax(Q, xx, 1.01 * Q, fa)["L_max"] / (2 * lHx) for fa in fs], color=col, label=f"x = Lambda Q^2 = {xx}")
a1.axhline(1, color=INK, lw=1)
a1.text(1.2e-4, 1.03, "super-Hubble (needed for growth)", fontsize=8, color=INK)
a1.set_xscale("log")
a1.set_ylim(0, 1.15)
a1.set_xlabel("f at the inner end of the segment (depth in the throat)")
a1.set_ylabel("max normal length / (2 l_H)")
a1.set_title("(b) Theorem 6c-1: normal false vacuum is always sub-Hubble", fontsize=9)
a1.legend(frameon=False, fontsize=8, loc="center right")

a2 = ax[1, 0]
for L, col in ((10.0, C1), (25.0, C2), (50.0, C3)):
    hh, iend = hist[L]
    n = len(hh) if iend is None else iend
    uu = np.arange(n) * 0.05
    a2.semilogy(uu, hh[:n], color=col, label=f"null extent {L:g} = {L / lH:.2f} l_H")
    if iend is not None:
        a2.plot([uu[-1]], [hh[n - 1]], "o", color=col, ms=8, mec="#fcfcfb", mew=2)
        a2.text(uu[-1] + 1, hh[n - 1], "r > 10^4 (reaches I+)", fontsize=8, color=INK, va="center")
a2.set_xlabel("u")
a2.set_ylabel("largest r in the false vacuum on ray u")
a2.set_title("(c) Anti-trapped false vacuum: super-Hubble regions inflate to I+", fontsize=9)
a2.legend(frameon=False, fontsize=8, loc="upper left")

a3 = ax[1, 1]
for fa, col in zip((0.01, 0.05, 0.1, 0.3, 0.6), (C1, C2, C3, "#6250d6", "#e87ba4")):
    rows = [s for s in scan if s["f_a"] == fa]
    a3.plot([s["flip_factor"] for s in rows], [s["growth"] for s in rows], "o-", color=col, ms=6, mec="#fcfcfb", label=f"f_a = {fa}")
a3.set_xscale("log")
a3.set_yscale("log")
a3.set_xlabel("flip strength / G1 budget at the corner")
a3.set_ylabel("growth of r in the false vacuum (no run reaches I+)")
a3.set_title("(d) From normal data the flipped false vacuum is eaten or swallowed", fontsize=9)
a3.legend(frameon=False, fontsize=8, loc="upper left", ncol=2)
fig.tight_layout()
fig.savefig(f"{FIG}/round6c_throat.png", dpi=150, facecolor="#fcfcfb")
print("wrote FIGURES/round6c_throat.png")
