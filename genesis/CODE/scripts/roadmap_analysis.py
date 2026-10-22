"""Roadmap Stages 1, 2 and 4: data and figure for ROADMAP_STAGE1_2_4.md.

Part A  analytic suppression bounds (bounds.py)
Part B  exact Stokes multiplier (stokes_exact.py)
Part C/D quantum sectors, daughter edge and clock (daughter.py)
Part E  alternative mechanisms (alternatives.py)
Part F  theoretical requirements ("five numbers")
"""
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse import tunnelling as T, contours as Cc
from babyuniverse import bounds as Bd, stokes_exact as SE, daughter as D, alternatives as Al

warnings.simplefilter("ignore")
out = {}

# ---------------------------------------------------------------- Part A
grid = []
for gam in np.linspace(-0.95, 0.95, 9):
    for fac in (1.01, 1.1, 1.5, 3.0, 10.0, 100.0):
        m = model_from_dimensionless(gam, (2 - 2 * gam) * fac)
        Ban, Bsq = T.B_canonical(m), Cc.B_sector_blind(m)
        grid.append(dict(gamma=gam, E_over_EAB=fac, B_an=Ban, B_sq=Bsq, B_mem=Bd.sine_bound(m),
                         B_elem=Bd.elementary_bound(m), flip=Cc.flip_class(m), predicted=Bd.flip_class_predicted(m)))
out["bounds_grid"] = grid
out["bounds_grid_min_ratio_Bsq_over_Belem"] = min(g["B_sq"] / g["B_elem"] for g in grid)
out["flip_class_theorem_agreement"] = all(g["flip"] == g["predicted"] for g in grid)

thin = []
for eps in (0.1, 0.03, 0.01, 0.003):
    Mcr = ThinWallModel(1.0, eps, 1e-9).M_critical
    for f in (1e-4, 1e-3, 1e-2, 0.1, 0.5, 0.9):
        m = ThinWallModel(1.0, eps, f * Mcr)
        B0 = T.B_canonical(ThinWallModel(1.0, eps, 0.0))
        Rf = Bd.R_flip(m)
        thin.append(dict(kappa_over_chi=eps, M_over_Mcr=f, E=m.E, x=Rf / (2 * eps / m.L2), B_over_B0=T.B_canonical(m) / B0,
                         Phi=Bd.Phi(m.E), elementary_over_B=Bd.elementary_bound(m) / T.B_canonical(m)))
out["thin_barrier"] = thin
out["thin_barrier_max_rel_err_eps0.003"] = max(abs(r["B_over_B0"] / r["Phi"] - 1) for r in thin if r["kappa_over_chi"] == 0.003)

thr = []
for gam in (-0.99, -0.5, 0.0, 0.2, 0.4, 0.6, 0.8, 0.99):
    lo, hi = 1.00001, 4.0
    for _ in range(50):
        mid = np.sqrt(lo * hi)
        if Bd.omega_condition(gam, (2 - 2 * gam) * mid):
            hi = mid
        else:
            lo = mid
    thr.append(dict(gamma=gam, omega_threshold_E_over_EAB=hi, E_sd_over_EAB=Bd.E_sd(gam) / (2 - 2 * gam) if gam > 0 else None))
out["delta_tE_criterion"] = thr

# ---------------------------------------------------------------- Part B
hs = (0.01, 0.005, 0.0025, 0.00125, 0.000625)
conv = SE.convergence_table(1.0, 1.5, hs)
out["stokes_convergence"] = conv
rq = np.sqrt(1.5**4 - 1)
out["stokes_controls"] = dict(nu0=SE.SquaredConstraintODE(1.0, 1.5, 0.5, beta=rq).multiplier(0.005),
                              nu_m1=SE.SquaredConstraintODE(1.0, 1.5, 0.5, beta=-rq).multiplier(0.005),
                              nu_m1_pred=SE.SquaredConstraintODE(1.0, 1.5, 0.5, beta=-rq).predicted(0.005))
out["stokes_other_parameters"] = {f"{c},{k}": SE.SquaredConstraintODE(c, k, 0.5).multiplier(0.0025) / np.sqrt(2)
                                  for c, k in ((1.0, 1.2), (1.0, 2.5), (0.5, 1.0))}

# ---------------------------------------------------------------- Parts C/D
mods = dict(single=ThinWallModel(1.0, 1.0, 0.1545), double=model_from_dimensionless(0.6, 2.4))
out["edge"] = {k: D.edge_report(m) for k, m in mods.items()}
out["clock"] = {k: D.clock_orientation(m) for k, m in mods.items()}
out["charges"] = [{kk: (vv if not isinstance(vv, np.ndarray) else None)
                   for kk, vv in D.charge_of_wave(0.7, n, 0.05).items()} for n in (-1, 0, 1, 2)]

# ---------------------------------------------------------------- Part E
M = 0.3 * ThinWallModel(1.0, 1.0, 1e-9).M_critical
B0 = T.B_canonical(ThinWallModel(1.0, 1.0, M))
out["seeded"] = [dict(m_over_M=f, B_over_B=Al.SeededModel(1.0, 1.0, M, f * M).B() / B0,
                      equator=Al.SeededModel(1.0, 1.0, M, f * M).initial_data_contains_equator(),
                      eps_out_R1=float(np.sign(Al.SeededModel(1.0, 1.0, M, f * M).A_out(
                          Al.SeededModel(1.0, 1.0, M, f * M).turning_points()[0]))))
                 for f in (0.0, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9)]
out["qei"] = {f"{E:.0e}": Al.species_required(E) for E in (1e15, 1e16, 1e17, 1e18)}
out["ds_parent"] = Al.ds_parent_upward_exponent()
out["fmm"] = [Al.fmm_junction(m) for m in (1.0, 10.0, 1e3)]
out["ledger"] = Al.mechanism_ledger()

# ---------------------------------------------------------------- Part F: requirements
def E_req(log10P, B_survive, M_planck):
    X = -log10P * np.log(10) / (4 * np.pi * np.sqrt(2 * B_survive / 27) * (3 * M_planck / (4 * np.pi)) ** (2 / 3))
    return Al.M_PL_GEV * X ** -1.5


req = []
for P in (-10, -100, -1000):
    for Bs in (10.0, 100.0):
        for Mp in (1.0, 10.0):
            E = E_req(P, Bs, Mp)
            rho = (E / Al.M_PL_GEV) ** 4
            RM = (3 * Mp / (4 * np.pi * rho)) ** (1 / 3)
            req.append(dict(log10_P_target=P, B_survive=Bs, M_planck=Mp, E_vac_min_GeV=E, R_M_lP=RM,
                            R_M_m=RM * 1.616e-35, tau_s=RM * 5.391e-44,
                            kappa_over_chi_needed=float(np.sqrt(6 * np.pi * np.sqrt(2 * Bs * rho**3 / (27 * np.pi**2)) / rho))))
out["requirements"] = req
for r in req:
    r["thin_barrier_valid"] = bool(r["kappa_over_chi_needed"] < 0.3 and r["E_vac_min_GeV"] < 0.1 * Al.M_PL_GEV and r["R_M_lP"] > 10)
out["frontier"] = [dict(B_survive=Bs, M_min=Mm, R_min=Rm, **Bd.frontier_optimum(Bs, Mm, Rm))
                   for Bs in (10.0, 100.0) for Mm in (1.0, 10.0) for Rm in (100.0, 1000.0)]
dump("roadmap_stage1_2_4.json", out)

# ---------------------------------------------------------------- figure
C1, C2, C3, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
                     "lines.linewidth": 2.0, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(2, 2, figsize=(10.5, 8.0))

a = ax[0, 0]
x = np.geomspace(0.05, 300, 200)
a.loglog(x, [Bd.Phi(1 / xx**2) for xx in x], color=C1, label=r"$\Phi(R_{\rm CdL}^2/R_M^2)$ (universal)")
a.loglog(x, 8 / 3 * x**2 + 2 / 3, color=C2, lw=1.2, ls="--", label=r"$\frac{8}{3}x^2+\frac{2}{3}$")
a.loglog(x, 1 + 16 / np.pi * x**3, color=C3, lw=1.2, ls=":", label=r"$1+\frac{16}{\pi}x^3$")
for eps, mk in ((0.03, "o"), (0.003, "s")):
    pts = [r for r in thin if r["kappa_over_chi"] == eps]
    a.loglog([r["x"] for r in pts], [r["B_over_B0"] for r in pts], mk, ms=8, mfc="none", mec=INK,
             label=rf"exact $B/B_{{\rm CdL}}$, $\kappa/\chi={eps}$")
a.set_ylim(0.8, 3e5)
a.set_xlabel(r"$x=R_M/R_{\rm CdL}$ (bubble radius / CdL radius)")
a.set_ylabel(r"$B_{\rm lab}/B_{\rm CdL}$")
a.set_title("A. Laboratory exponent: thin-barrier law", loc="left", fontsize=10)
a.legend(fontsize=7, frameon=False)

a = ax[0, 1]
for lam, col, lab in ((1.0, C1, "g left"), (0.5, C2, "Weyl"), (0.0, C3, "g right")):
    e = np.abs(conv[str(lam)]["rel_err"])
    a.loglog(hs, e, "-o", color=col, ms=8, label=f"ordering: {lab}")
a.set_xlabel(r"$h=\hbar G$ (units $\chi=1$)")
a.set_ylabel(r"$|s/\sqrt{2}-1|$")
a.set_title(r"B. Exact Stokes multiplier $\to\sqrt{2}$ (P-Q2)", loc="left", fontsize=10)
a.legend(fontsize=8, frameon=False)

a = ax[1, 0]
Ev = np.geomspace(1e14, 1e19, 80)
for Bs, Mp, col in ((10.0, 1.0, C1), (100.0, 1.0, C2), (10.0, 10.0, C3)):
    rho = (Ev / Al.M_PL_GEV) ** 4
    a.loglog(Ev, [Bd.minimum_exponent(Bs, Mp, r) for r in rho], color=col,
             label=rf"$B_{{\rm CdL}}\geq{Bs:.0f}$, $M\geq{Mp:.0f}\,m_P$")
a.axhline(np.log(1e100), color=MUTED, lw=0.8)
a.text(1.2e14, np.log(1e100) * 1.15, r"$P=10^{-100}$", color=MUTED, fontsize=8)
a.set_xlabel(r"false-vacuum scale $E_{\rm vac}$ [GeV]")
a.set_ylabel(r"minimum $B_{\rm lab}$")
a.set_title("F. Minimum laboratory exponent (thin barrier)", loc="left", fontsize=10)
a.legend(fontsize=8, frameon=False)

a = ax[1, 1]
m = mods["single"]
R2d = m.turning_points()[-1]
Rs = np.concatenate([np.linspace(R2d * (1 + 1e-4), 0.99999 / m.chi, 80), np.linspace(1.00001 / m.chi, 1.1 / m.chi, 60)])
P = [D.lorentzian_momentum(m, R) for R in Rs]
Pt = [D.regularized_momentum(m, R) for R in Rs]
a.plot(Rs * m.chi, P, color=C1, label=r"static-frame $P$ (Phase V 'edge')")
a.plot(Rs * m.chi, Pt, color=C2, label=r"$\tilde P=P-(R/2G)\ln|f_{\rm in}|$")
a.axvline(1.0, color=MUTED, lw=0.8)
a.text(1.002, -4.5, r"interior horizon $f_{\rm in}=0$", color=MUTED, fontsize=8)
a.set_xlabel(r"daughter wall radius $\chi R$")
a.set_ylabel("canonical momentum")
a.set_title("D. The daughter edge is a coordinate singularity", loc="left", fontsize=10)
a.legend(fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "roadmap_stage1_2_4.png"), dpi=150)
print("thin-barrier max rel err (eps=0.003):", out["thin_barrier_max_rel_err_eps0.003"])
print("flip-class theorem agreement:", out["flip_class_theorem_agreement"])
print("stokes:", {k: v["rel_err"][-1] for k, v in conv.items()})
print("requirements:", [(r["log10_P_target"], r["B_survive"], r["M_planck"], f"{r['E_vac_min_GeV']:.2e}") for r in req])
