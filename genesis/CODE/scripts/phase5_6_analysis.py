"""Phases V and VI: reduced quantum system (domain of validity, channels, sectors, information,
topology, semiclassical validity) and observability (laboratory bounds, characteristic masses,
primordial-black-hole windows).  Writes DATA/phase5_quantum_system.json,
DATA/phase6_observability.json, FIGURES/phase5_quantum_system.png, FIGURES/phase6_observability.png."""
import json
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import DATA, FIG, dump
from babyuniverse.model import ThinWallModel
from babyuniverse import quantum_system as Q, observability as O, causal

warnings.simplefilter("ignore")

# ===================================================================== Phase V
p5 = dict(channels=[], topology=[], wkb=[], information=[], critical_mass_limit=[])
for M in (0.02, 0.08, 0.1545, 0.25, 0.3, 0.36):
    m = ThinWallModel(1.0, 1.0, M)
    sd = Q.static_domain(m)
    tp = Q.asymptotic_time(m, "parent")
    td = Q.asymptotic_time(m, "daughter")
    p5["channels"].append(dict(M=M, **sd, parent_time=tp, daughter_time=td,
                               parent_log_slope=float(np.mean(tp["increments"][1:]) / np.log(10)),
                               two_GM=2 * M,
                               P_parent=Q.momentum_near_edge(m, "parent"), P_daughter=Q.momentum_near_edge(m, "daughter")))
    for lab, Rt in (("parent", sd["R1"]), ("daughter", sd["R2"])):
        s = causal.time_symmetric_slice(m, Rt)
        p5["topology"].append(dict(M=M, branch=lab, R=Rt, eps_in=s["eps_in"], eps_out=s["eps_out"],
                                   extremal=[(e.kind, e.r, e.where) for e in s["extremal"]], closed=s["closed"],
                                   spatial_topology="S^3" if s["closed"] else "R^3"))
    print(f"M={M}: parent t grows {p5['channels'][-1]['parent_log_slope']:.4f} per e-fold (2GM = {2*M}); "
          f"daughter t -> {td['t'][-1]:.5f} (finite); edges {sd['parent_edge']:.3f}, {sd['daughter_edge']:.3f}")

# M_cr -> 1/(2 G L) for kappa/chi -> 0 (Phase II formulas: g_min ~ 4 kappa^2 / L^2)
for k in (1.0, 0.3, 0.1, 0.03, 0.01, 0.0036):
    m = ThinWallModel(1.0, k, 1e-9)
    L = np.sqrt(1 + k**2)
    p5["critical_mass_limit"].append(dict(kappa_over_chi=k, M_cr=float(m.M_critical), one_over_2GL=float(1 / (2 * L)),
                                          ratio=float(m.M_critical * 2 * L)))

# physical (GUT-like) samples from Stage One: WKB validity, information, evaporation
phys = json.load(open(os.path.join(DATA, "physical_scales.json")))
for row in phys:
    chi_pl = np.sqrt(8 * np.pi * (row["vacuum_scale_GeV"] / 1.22e19) ** 4 / 3)
    k = row["kappa_over_chi"]
    for smp in row["samples"]:
        mm = ThinWallModel(1.0, k, smp["M_over_Mcr"] * ThinWallModel(1.0, k, 1e-9).M_critical)
        w = Q.wkb_parameter(mm, chi_pl)
        info = Q.information_budget(smp["B"], smp["M_kg"])
        ev = Q.evaporation_competition(smp["B"], smp["M_kg"])
        p5["wkb"].append(dict(E_vac=row["vacuum_scale_GeV"], E_wall=row["wall_scale_GeV"], M_over_Mcr=smp["M_over_Mcr"],
                              M_kg=smp["M_kg"], B=smp["B"], **w))
        p5["information"].append(dict(E_vac=row["vacuum_scale_GeV"], E_wall=row["wall_scale_GeV"],
                                      M_over_Mcr=smp["M_over_Mcr"], M_kg=smp["M_kg"], **info, **ev))
print("WKB eps (physical) max over GUT samples:", max(x["eps_physical"] for x in p5["wkb"]))
print("log10 total daughter probability bound (max, all samples):", max(x["log10_total_probability_bound"] for x in p5["information"]))
print("   restricted to semiclassical black holes (M > 10 m_Pl):",
      max(x["log10_total_probability_bound"] for x in p5["information"] if x["semiclassical_black_hole"]))
for k in p5["critical_mass_limit"]:
    print("   M_cr * 2L at kappa/chi =", k["kappa_over_chi"], ":", round(k["ratio"], 6))
dump("phase5_quantum_system.json", p5)

# ===================================================================== Phase VI
p6 = dict(grid=[], masses=[], gut_reference=None)
Ev = np.logspace(3, 18.5, 32)
Ew = np.logspace(3, 18.5, 32)
B0 = np.full((len(Ev), len(Ew)), np.nan)
for a, ev in enumerate(Ev):
    for b, ew in enumerate(Ew):
        try:
            B0[a, b] = O.B_zero(ev, ew)
        except Exception:
            pass
p6["grid"] = dict(E_vac_GeV=Ev.tolist(), E_wall_GeV=Ew.tolist(), log10_B_lab_lower=np.log10(0.56 * B0).tolist())
p6["gut_reference"] = dict(B_zero=O.B_zero(1e16, 1e16), B_zero_flat=O.B_zero_flat(1e16, 1e16), **O.masses(1e16, 1e16))
# where could B_lab be O(1)?
mask = 0.56 * B0 < 10
p6["B_lab_below_10"] = [dict(E_vac=float(Ev[a]), E_wall=float(Ew[b]), B0=float(B0[a, b]))
                        for a, b in zip(*np.nonzero(mask))][:40]
for ev in (1e3, 1e5, 1e7, 1e9, 1e11, 1e13, 1e16, 1e18):
    ms = O.masses(ev, ev)
    p6["masses"].append(dict(E_vac=ev, E_wall=ev, **ms, lifetime_s=O.hawking_lifetime_s(ms["M_cr_g"]),
                             window=O.window_of(ms["M_cr_g"])))
    print(f"E = {ev:.0e} GeV: M_cr = {ms['M_cr_g']:.3g} g, M_H = {ms['M_H_g']:.3g} g, "
          f"lifetime {O.hawking_lifetime_s(ms['M_cr_g']):.3g} s, window: {O.window_of(ms['M_cr_g'])}")
# the energy scales mapped to each PBH window (sigma^(1/3) = rho^(1/4) = E):  M_cr ~ M_GUT (1e16/E)^2
M_ref = p6["gut_reference"]["M_cr_g"]
p6["window_energy_ranges"] = []
for name, lo, hi in O.PBH_WINDOWS:
    E_hi = 1e16 * np.sqrt(M_ref / lo) if lo > 0 else np.inf
    E_lo = 1e16 * np.sqrt(M_ref / hi)
    p6["window_energy_ranges"].append(dict(window=name, M_g=[lo, hi], E_vac_GeV=[float(E_lo), float(E_hi)]))
    print(f"   {name}: E_vac from {E_lo:.2g} to {E_hi:.2g} GeV")
# loophole search: the most favourable laboratory case within the validity of the calculation
p6["best_laboratory"] = {}
for Mp in (10.0, 100.0):
    res = O.best_laboratory_case(np.logspace(10, 18.7, 36), np.linspace(0.5, 1.6, 45), M_planck=Mp)
    adm = [x for x in res["rows"] if x.get("valid") and x["metastable"] and x["wkb_ok"]]
    prof = []
    for ev in sorted(set(x["E_vac"] for x in adm)):
        rr = [x for x in adm if x["E_vac"] == ev]
        bb = min(rr, key=lambda x: x["B_best"])
        prof.append(dict(E_vac=ev, B_min=bb["B_best"], ratio=bb["E_wall"] / ev, B_CdL=bb["B_CdL"]))
    p6["best_laboratory"][f"M={Mp:g} m_Pl"] = dict(best=res["best"], n_scanned=res["n_scanned"],
                                                    n_admissible=res["n_admissible"], profile=prof)
    print(f"best laboratory case, M = {Mp:g} m_Pl: B = {res['best']['B_best']:.4g} "
          f"(log10 P = {res['best']['log10_P_best']:.4g}) at E_vac = {res['best']['E_vac']:.3g} GeV, "
          f"ratio {res['best']['E_wall'] / res['best']['E_vac']:.3f}")
dump("phase6_observability.json", p6)

# ===================================================================== figures
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for c in p5["channels"]:
    ax[0].semilogx(c["parent_time"]["cutoffs"], c["parent_time"]["t"], "o-", label=f"parent, M = {c['M']}")
    ax[0].semilogx(c["daughter_time"]["cutoffs"], c["daughter_time"]["t"], "s--", label=f"daughter, M = {c['M']}")
ax[0].invert_xaxis()
ax[0].set_xlabel("distance to the edge of the static domain (fraction)")
ax[0].set_ylabel("asymptotic time t_inf from the turning point (units 1/chi)")
ax[0].set_title("parent reaches its edge at t = inf (channel); daughter in finite t (boundary)", fontsize=8)
ax[0].legend(fontsize=6, ncol=2)
for (evac, ewall), mk in (((1e16, 1e16), "o"), ((1e16, 1e17), "^"), ((1e18, 1e18), "s")):
    w = [x for x in p5["wkb"] if x["E_vac"] == evac and x["E_wall"] == ewall]
    lab = f"vac {evac:.0e}, wall {ewall:.0e} GeV"
    ln = ax[1].loglog([x["M_kg"] for x in w], [x["eps_physical"] for x in w], mk + "-", label="WKB parameter, " + lab)
    ax[1].loglog([x["M_kg"] for x in w], [1 / x["B"] for x in w], mk + "--", color=ln[0].get_color(), label="1/B, " + lab)
ax[1].axvline(2.176e-8, color="k", lw=0.7, ls=":", label="Planck mass")
ax[1].set_xlabel("bubble mass (kg)")
ax[1].set_title("semiclassical validity in the GUT examples", fontsize=9)
ax[1].legend(fontsize=6)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "phase5_quantum_system.png"), dpi=150)
plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(19, 5))
for key, mk in (("M=10 m_Pl", "o-"), ("M=100 m_Pl", "s-")):
    prof = p6["best_laboratory"][key]["profile"]
    ax[2].loglog([x["E_vac"] for x in prof], [x["B_min"] for x in prof], mk, label=f"min B over wall scale, {key}")
ax[2].set_xlabel("vacuum energy scale (GeV)")
ax[2].set_ylabel("smallest admissible exponent B")
ax[2].set_title("loophole search: metastable vacuum, semiclassical BH, WKB valid", fontsize=9)
ax[2].legend(fontsize=7)
X, Y = np.meshgrid(np.log10(Ew), np.log10(Ev))
lb = np.clip(np.log10(0.56 * B0), -30, 79)
cs = ax[0].contourf(X, Y, lb, levels=np.arange(-30, 81, 5), cmap="viridis")
ax[0].contour(X, Y, lb, levels=[0, 1, 2], colors="w", linewidths=1)
ax[0].contour(X, Y, np.log10(B0), levels=[0.0], colors="r", linewidths=1.5, linestyles="--")
ax[0].text(4.0, 16.5, "B_CdL < 1:\nfalse vacuum not metastable\n(no laboratory)", color="r", fontsize=7)
fig.colorbar(cs, ax=ax[0], label="log10 (lower bound on B_lab)")
ax[0].set_xlabel("log10 wall scale sigma^(1/3) / GeV")
ax[0].set_ylabel("log10 vacuum scale rho^(1/4) / GeV")
ax[0].set_title("laboratory: P <= exp(-B_lab); white contours B = 1, 10, 100", fontsize=9)
E = np.logspace(0, 19, 200)
Mc = [O.masses(e, e)["M_cr_g"] for e in E]
ax[1].loglog(E, Mc, "k-", lw=1.5, label="M_cr (= 1/(2GL)), sigma^(1/3) = rho^(1/4) = E")
cols = ["#eeeeee", "#f6c1c1", "#f3d9a4", "#c9e9c0", "#c6d8f0", "#ddc9ef"]
for (name, lo, hi), cc in zip(O.PBH_WINDOWS, cols):
    ax[1].axhspan(max(lo, 1e-6), hi, color=cc, alpha=0.8, label=name)
ax[1].set_xlabel("vacuum energy scale E (GeV)")
ax[1].set_ylabel("characteristic daughter-hosting black-hole mass (g)")
ax[1].set_ylim(1e-6, 1e40)
ax[1].legend(fontsize=6, loc="upper right")
ax[1].set_title("early universe: characteristic mass vs primordial-black-hole windows", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "phase6_observability.png"), dpi=150)
plt.close(fig)
