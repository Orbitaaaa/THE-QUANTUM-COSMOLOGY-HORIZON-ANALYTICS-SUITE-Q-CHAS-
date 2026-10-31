"""Round 11 physics (G11-01 .. G11-05, G11-08): singularity strength, D = 10, the inner-extremal tower, horizon curvature of collider holes,
nonlinear hysteresis.  Writes DATA/round11_physics.json and FIGURES/round11_physics.png."""
import json
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import strength11 as S
from babyuniverse import stable11 as ST
from babyuniverse import hysteresis11 as HY
from babyuniverse import unified6l as U

warnings.simplefilter("ignore")


def pack(ri, keys=("v", "lam", "I_K", "I_T", "sigma", "r")):
    return {k: [float(x) for x in ri[k]] for k in keys}


if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round11_physics.json"))
else:
    out = {}
    # G11-01: strength along an interior outgoing ray, D = 5 Hayward, pulse only and pulse + power-law tail, two resolutions each
    for tag, kw in (("pulse", {}), ("pulse_tail_p3", dict(tail_amp=0.01, tail_p=3.0))):
        out["strength_D5_" + tag] = {}
        for nu, nv in ((2000, 8000), (4000, 16000)):
            res = S.evolve(5, 10.0, "hayward", 0.05, 26.0, nu, nv, **kw)
            out["strength_D5_" + tag]["%dx%d" % (nu, nv)] = pack(S.ray_integrals(res, 5, 0.7))
        print("strength", tag, flush=True)
    # G11-03: D = 10 (mu = 3 mu_ext), two resolutions, converged range
    D = 10; mu0 = 3 * U.extremal_mu(10)
    out["D10"] = {}
    for nu, nv in ((4000, 16000), (8000, 32000)):
        res = S.refined_run(D, mu0, amp=0.02, vmax=22.0, nu=nu, nv=nv)
        u = res["u"]; ins = np.where(u > res["u_eh"] + 0.15)[0]; i = ins[int(0.3 * (len(ins) - 1))]
        ri = S.ray_integrals(res, D, u_value=u[i])
        K = res["K"][i]
        out["D10"]["%dx%d" % (nu, nv)] = dict(u=float(u[i]), v=res["v"].tolist(), K=np.where(np.isfinite(K), K, np.nan).tolist(),
                                              mu=res["mu"][i].tolist(), psi=res["psi"][i].tolist(), **pack(ri, ("I_K", "I_T")), v_ray=ri["v"])
    print("D10", flush=True)
    # classification map
    out["classification"] = []
    for D_, fac in ((5, 10.0 / U.extremal_mu(5)), (5, 25.0 / U.extremal_mu(5)), (5, 100.0 / U.extremal_mu(5)), (10, 1.3), (10, 3.0), (10, 10.0)):
        try:
            out["classification"].append(S.classify(D_, fac * U.extremal_mu(D_), amp=0.05 if D_ == 5 else 0.02, vmax=24.0 if D_ == 5 else 16.0))
        except Exception as e:
            out["classification"].append(dict(D=D_, error=str(e)))
    print("classification", flush=True)
    # G11-02: inner-extremal tower
    out["ie_checks"] = ST.ie_checks()
    r1 = ST.run_ie(amp=0.02, v_max=36.0, nu=800, nv=3000, n_out_u=800, n_out_v=30)
    rend = r1["r"][:, -1]; k = int(np.where(~(np.isfinite(rend) & (rend > r1["r0"])))[0][0]); ue = 0.5 * (r1["u"][k - 1] + r1["u"][k])
    out["ie"] = {}
    for nu, nv in ((2000, 8000), (4000, 16000)):
        res = ST.run_ie(amp=0.02, v_max=36.0, nu=nu, nv=nv, u_eh=ue)
        u = res["u"]; ins = np.where(u > ue + 0.3)[0]; i = ins[int(0.8 * (len(ins) - 1))]
        out["ie"]["%dx%d" % (nu, nv)] = dict(u=float(u[i]), v=res["v"].tolist(), K=res["K"][i].tolist(), mu=res["mu"][i].tolist(), r=res["r"][i].tolist(),
                                             r_minus=float(res["r_minus"]))
    print("ie", flush=True)
    # G11-05: horizon curvature of collider holes (D = 10)
    rows = []
    for kind in ("hayward", "dymnikova"):
        me = U.extremal_mu(10, 1.0, kind)
        for fac in (1.01, 1.1, 1.5, 3.0, 10.0):
            rm, rp = U.static_horizons(fac * me, 10, 1.0, kind)
            rows.append(dict(kind=kind, mu_over_ext=fac, psi_plus_alpha=float(1 / rp**2), psi_minus_alpha=float(1 / rm**2), r_plus_over_l=float(rp)))
    out["horizon_curvature"] = rows
    # G11-08: nonlinear hysteresis (test-field inhomogeneities + energy accounting)
    out["hysteresis"] = [HY.run(VF_over_Lc=0.5, eps=1e-3, phi_b=0.3), HY.run(VF_over_Lc=2.0, eps=1e-2, phi_b=0.0), HY.run(VF_over_Lc=0.05, eps=1e-3, phi_b=1.0)]
    print("hysteresis", flush=True)
    dump("round11_physics.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 1.8, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
ax = ax.ravel()
for tag, col in (("pulse", C1), ("pulse_tail_p3", C2)):
    d = out["strength_D5_" + tag]["4000x16000"]
    ax[0].semilogy(d["v"], np.maximum(d["I_K"], 1e-6), color=col, label="Krolak integral, %s" % tag.replace("_", " "))
    ax[0].semilogy(d["v"], np.maximum(d["I_T"], 1e-6), color=col, ls="--", label="Tipler integral, %s" % tag.replace("_", " "))
ax[0].set_xlim(0, 24); ax[0].set_xlabel("v along an interior outgoing ray (D = 5)")
ax[0].set_title("(a) G11-01: Tipler integral converges (deformationally weak), Krolak integral grows\nbut decelerates; area radius finite (Raychaudhuri: r'' = -r R_kk/(D-2))", fontsize=8)
ax[0].legend(frameon=False, fontsize=7)
for tag, col in (("4000x16000", C1), ("8000x32000", C2)):
    d = out["D10"][tag]
    v = np.array(d["v"]); K = np.array(d["K"], dtype=float)
    ax[1].semilogy(v, K, color=col, label="K l^4, D = 10, %s" % tag)
    ax[1].semilogy(v, np.array(d["mu"], dtype=float) / 32.5, color=col, ls=":", label="mu / mu_0, %s" % tag)
ax[1].axvline(14, color=MUTED, ls="--", lw=1)
ax[1].set_xlim(0, 21); ax[1].set_ylim(0.5, 1e7); ax[1].set_xlabel("v (D = 10, mu = 3 mu_ext)")
ax[1].set_title("(b) G11-03: in D = 10 the mass inflates x3 and psi -> 0.99/alpha with bounded curvature\nthrough v ~ 14 (converged); later times unresolved", fontsize=8)
ax[1].legend(frameon=False, fontsize=7)
for tag, col in (("2000x8000", C3), ("4000x16000", C4)):
    d = out["ie"][tag]
    ax[2].semilogy(d["v"], np.maximum(np.array(d["K"], dtype=float), 1e-3), color=col, label="inner-extremal tower, %s" % tag)
d5 = out["strength_D5_pulse"]["4000x16000"]
ax[2].set_ylim(1, 1e14); ax[2].set_xlabel("v (D = 5)")
cls = [c for c in out["classification"] if c.get("D") == 5]
ax[2].axhline(cls[0]["K_max"] if cls else 1e12, color=C2, ls=":", lw=1, label="Hayward D = 5: K reaches %.0e" % (cls[0]["K_max"] if cls else 1e12))
ax[2].set_title("(c) G11-02: an inner-extremal tower (kappa_- = 0, tuned mass) has no mass inflation and\nbounded curvature; rays cross r_- at finite v (converged)", fontsize=8)
ax[2].legend(frameon=False, fontsize=7)
ax[3].axis("off")
txt = ["G11-01/03 classification of a scalar pulse into a model-U hole (Hayward)", ""]
for c in out["classification"]:
    if "error" in c:
        continue
    txt.append("D=%d  psi_-alpha=%.3f  K growth x%.1e  inner AH crossed: %s  -> %s" % (c["D"], c["psi_minus"], c["growth"], "yes" if c["crossed_inner_AH"] else "no", c["verdict"]))
txt += ["", "G11-05 horizon curvature of collider holes (D = 10):"]
for r in out["horizon_curvature"]:
    txt.append("  %-9s mu/mu_ext=%5.2f  psi_+ alpha=%.2f  psi_- alpha=%.3f" % (r["kind"], r["mu_over_ext"], r["psi_plus_alpha"], r["psi_minus_alpha"]))
ax[3].text(-0.05, 1.0, "\n".join(txt), va="top", family="monospace", fontsize=6.2, transform=ax[3].transAxes)
ax[3].set_title("(d) outcomes and horizon curvature (EFT breaks down at psi alpha ~ 1)", fontsize=8)
hy = out["hysteresis"]
for k, (h, col) in enumerate(zip(hy, (C1, C2, C3))):
    ax[4].semilogy(h["ells"], h["max_amplification"], "o-", color=col, label="V_F/Lambda_c = %g: amplification" % h["VF_over_Lc"])
    ax[4].semilogy(h["ells"], np.array(h["max_energy_fraction"]) * 1e8, "s:", color=col, label="V_F/Lambda_c = %g: E_l/E_hom x 1e8" % h["VF_over_Lc"])
ax[4].set_xscale("log"); ax[4].set_xlabel("angular mode l on S^9")
ax[4].set_title("(e) G11-08: inhomogeneities through the hysteresis cycles (test field, initial 1e-6):\namplified <~ 1e4, energy fraction <~ 2e-8 before escape -> survives for delta phi <~ 1e-3", fontsize=8)
ax[4].legend(frameon=False, fontsize=6.5)
for k, (h, col) in enumerate(zip(hy, (C1, C2, C3))):
    w = np.array(h["work_per_cycle"])
    ax[5].plot(range(1, len(w) + 1), np.sign(w) * np.log10(1 + np.abs(w)), "o-", color=col, label="V_F/Lambda_c = %g (constraint %.0e)" % (h["VF_over_Lc"], h["constraint_max_abs"]))
ax[5].set_xlabel("bounce cycle"); ax[5].set_ylabel("sign(W) log10(1 + |W|)")
ax[5].set_title("(f) G11-08 energy accounting: net positive work per cycle (the pumping), Friedmann\nconstraint held to <= 1e-7", fontsize=8)
ax[5].legend(frameon=False, fontsize=7)
fig.subplots_adjust(left=0.05, right=0.98, top=0.93, bottom=0.06, wspace=0.3, hspace=0.4)
fig.savefig(f"{FIG}/round11_physics.png", dpi=150, facecolor="#fcfcfb")
print("done")
