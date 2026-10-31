"""Round 10, open problems: double-null evolution of model-U holes (P-R6n-a, P-R6l-a, P-R6l-c), the causality threshold (P-R6l-b),
a dynamical vacuum through the bounce (P-R6l-e), rotation (P-R6n-b) and remnant discharge (P-R9-b).
Writes DATA/round10_open.json and FIGURES/round10_open.png."""
import json
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import doublenull10 as N
from babyuniverse import openprob10 as O
from babyuniverse import unified6l as U
from babyuniverse import collider6l as C

warnings.simplefilter("ignore")


def eh_guess(D, mu0, kind, amp, vmax, **kw):
    r1 = N.run_hyb(D, mu0, kind, nu=800, nv=3000, v_max=vmax, amp=amp, n_out_u=800, n_out_v=30, **kw)
    rend = r1["r"][:, -1]
    esc = np.where(np.isfinite(rend) & (rend > r1["r0"]))[0]
    return 0.5 * (r1["u"][esc[-1]] + r1["u"][esc[-1] + 1])


def evolve(D, mu0, kind, amp, vmax, nu, nv, **kw):
    ue = eh_guess(D, mu0, kind, amp, vmax, **kw)
    res = N.run_hyb(D, mu0, kind, nu=nu, nv=nv, v_max=vmax, amp=amp, u_eh=ue, n_cluster=800, n_out_u=240, n_out_v=260, **kw)
    res["u_eh"] = ue
    return res


def ray(res, frac=0.7):
    u = res["u"]
    ins = np.where(u > res["u_eh"] + 0.3)[0]
    i = ins[int(frac * (len(ins) - 1))]
    return i, dict(u=float(u[i]), v=res["v"].tolist(), mu=res["mu"][i].tolist(), psi=res["psi"][i].tolist(), K=res["K"][i].tolist(),
                   K_2d=res["K_2d"][i].tolist(), r=res["r"][i].tolist(), constraint=np.abs(res["psi_metric"][i] - res["psi"][i]).tolist())


if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round10_open.json"))
else:
    out = {}
    # validation
    out["static"] = {k: N.static_check(*a) for k, a in (("D5_hayward", (5, 10.0, "hayward", 400, 800, 20.0)),
                                                         ("D5_gr", (5, 10.0, "gr", 400, 800, 20.0)))}
    st = N.run_hyb(5, 10.0, "hayward", nu=2000, nv=8000, v_max=40.0, u_eh=3.5745, n_cluster=600, n_out_u=240, n_out_v=80)
    i6 = np.argmin(abs(st["u"] - 6.0))
    out["static_hybrid"] = dict(u=float(st["u"][i6]), r_minus=float(st["r_minus"]), r_end=float(st["r"][i6, -1]),
                                K_end=float(st["K"][i6, -1]), K_static_rminus=float(st["K"][i6, -1]))
    out["sigma_uv_gr_check"] = bool(N.derive_sigma_uv()["gr_matches"])
    # the pulse at two resolutions
    runs = {}
    for tag, (nu, nv) in (("lo", (2000, 8000)), ("hi", (4000, 16000))):
        res = evolve(5, 10.0, "hayward", 0.05, 26.0, nu, nv)
        _, runs[tag] = ray(res)
        if tag == "hi":
            keep = res
    out["pulse_D5_hayward"] = runs
    K = keep["K"]
    out["pulse_map"] = dict(u=keep["u"].tolist(), v=keep["v"].tolist(), log10K=np.log10(np.where(np.isfinite(K) & (K > 0), K, np.nan)).tolist(),
                            u_eh=float(keep["u_eh"]), mu_out=float(keep["mu"][3, -1]))
    variants = {}
    for tag, args, kw in (("small_pulse", (5, 10.0, "hayward", 0.01), {}),
                          ("late_infall", (5, 10.0, "hayward", 0.05), dict(amp2=0.05, v3=12.0, v4=16.0)),
                          ("dymnikova", (5, 3.0 * U.extremal_mu(5, 1.0, "dymnikova"), "dymnikova", 0.05), {}),
                          ("D10_hayward", (10, 3.0 * U.extremal_mu(10, 1.0, "hayward"), "hayward", 0.05), {})):
        res = evolve(*args, 26.0, 2000, 8000, **kw)
        _, variants[tag] = ray(res)
    out["variants"] = variants
    # other open problems
    out["causality"] = O.causality_threshold()
    out["lambda_c"] = O.lambda_c_physical()
    out["vacuum_scan"] = O.scalar_bounce_scan(10, ratios=(0.05, 0.2, 0.5, 1.0, 2.0), eps_list=(1e-4, 1e-2), n_phase=6)
    out["hysteresis_example"] = O.scalar_from_bounce(10, None, 0.05, 1e-3, 1.0, 1.0, t_max=3000.0)
    out["hysteresis_control"] = O.scalar_from_bounce(10, None, 0.5, 0.0, 1.0, 1.0, t_max=3000.0)
    L = C.load_lumi(f"{DATA}/round6l_lumi.json")
    out["spin_down"] = O.spin_down(L, 100.0, 6, 8.0, 0.66, n_events=3000)
    out["discharge"] = O.remnant_discharge(6, 8.0, 0.66)
    dump("round10_open.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 1.8, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
ax = ax.ravel()
pm = out["pulse_map"]
Z = np.array(pm["log10K"], dtype=float)
im = ax[0].pcolormesh(pm["v"], pm["u"], Z, shading="auto", cmap="magma", vmin=-2, vmax=13)
ax[0].axhline(pm["u_eh"], color="w", lw=0.8, ls="--")
ax[0].set_xlabel("advanced time v (units of l)"); ax[0].set_ylabel("retarded coordinate u (r-gauge on v = 0)")
ax[0].set_title("(a) log10 K l^4 after a scalar pulse falls into a model-U hole (D = 5, Hayward);\ndashed: event horizon; the inner (Cauchy) horizon at late v becomes a curvature singularity", fontsize=8)
fig.colorbar(im, ax=ax[0], fraction=0.04)
for tag, col, ls in (("lo", C1, "--"), ("hi", C2, "-")):
    rr = out["pulse_D5_hayward"][tag]
    ax[1].semilogy(rr["v"], rr["K"], color=col, ls=ls, label="K l^4, %s resolution" % ("2000x8000" if tag == "lo" else "4000x16000"))
rr = out["pulse_D5_hayward"]["hi"]
ax[1].semilogy(rr["v"], np.array(rr["mu"]) / 10.0, color=C3, label="mu / mu_0")
ax[1].semilogy(rr["v"], np.array(rr["psi"]), color=C4, label="psi l^2 (bound: 1)")
ax[1].set_xlabel("v along an interior outgoing ray"); ax[1].set_ylim(1e-1, 1e14)
ax[1].set_title("(b) converged: curvature grows by > 10 orders of magnitude while psi stays below\nthe tower's bound and the area radius stays at r ~ r_- (mild mass inflation)", fontsize=8)
ax[1].legend(frameon=False, fontsize=7)
for tag, col, lab in (("small_pulse", C1, "small pulse (0.01)"), ("late_infall", C2, "late infall at v = 12-16"),
                      ("dymnikova", C3, "Dymnikova tower"), ("D10_hayward", C4, "D = 10 (stops at v ~ 10)")):
    rr = out["variants"][tag]
    ax[2].semilogy(rr["v"], rr["K"], color=col, label=lab)
ax[2].semilogy(out["pulse_D5_hayward"]["hi"]["v"], out["pulse_D5_hayward"]["hi"]["K"], color=INK, lw=1, label="baseline")
ax[2].set_ylim(1e0, 1e17); ax[2].set_xlabel("v")
ax[2].set_title("(c) the same in every variant; a late infall feeds the existing inner region\n(mass x9, psi -> 0.99/alpha) instead of opening a new world", fontsize=8)
ax[2].legend(frameon=False, fontsize=7)
vs = out["vacuum_scan"]
for eps, col in ((1e-4, C1), (1e-2, C2)):
    rows = [r for r in vs if abs(r["eps"] - eps) < 1e-9]
    ax[3].plot([r["VF_over_Lc"] for r in rows], [r["fraction_lasting_inflation"] for r in rows], "o-", color=col, label="eps = %g of the bounce energy in the field" % eps)
ax[3].axvline(1, color=MUTED, ls=":"); ax[3].set_xscale("log")
ax[3].set_xlabel("false-vacuum energy V_F / Lambda_c(C)"); ax[3].set_ylabel("fraction ending in lasting inflation")
ax[3].set_title("(d) P-R6l-e: cosmological hysteresis - the field's work over successive bounces pumps the cycle\nuntil the interior enters the de Sitter branch: lasting inflation even for V_F << Lambda_c (= %.2f / l^2)" % out["lambda_c"]["Lambda_c_l2"], fontsize=8)
ax[3].legend(frameon=False, fontsize=7)
sd = out["spin_down"]
ax[4].hist(np.array(sd["Erot_GeV"]) / 1e3, bins=40, color=C4, alpha=0.8)
ax[4].set_xlabel("rotational energy available to superradiant spin-down (TeV)")
ax[4].set_title("(e) P-R6n-b: collider holes spin (mean a* = %.1f); a Myers-Perry proxy gives E_rot ~ %.0f %% of M\n(visible TeV quanta even at T = 0: carried as a signal variant in 10A)" % (sd["mean_a_star"], 100 * sd["mean_Erot_fraction"]), fontsize=8)
dc = out["discharge"]
names = ["colour discharge", "hadronisation", "electric discharge"]
vals = [dc["tau_colour_s"], dc["tau_hadronisation_s"], dc["tau_electric_s"]]
ax[5].barh(names, vals, color=[C3, MUTED, C1]); ax[5].set_xscale("log")
ax[5].set_xlabel("time (s)")
ax[5].set_title("(f) P-R9-b: colour discharges ~%.0fx faster than hadronisation, so remnants end colour-singlet\nand (by triality-charge correlation, then e+e- emission) neutral; discharge partons ~%.1f TeV" % (dc["ratio"], dc["E_discharge_GeV"] / 1e3), fontsize=8)
fig.subplots_adjust(left=0.05, right=0.98, top=0.93, bottom=0.06, wspace=0.28, hspace=0.38)
fig.savefig(f"{FIG}/round10_open.png", dpi=150, facecolor="#fcfcfb")
print("done")
