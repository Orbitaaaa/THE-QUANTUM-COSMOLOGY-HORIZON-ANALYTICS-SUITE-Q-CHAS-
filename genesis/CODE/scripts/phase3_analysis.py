"""Phase III: monodromy theorem, sheet theorem, Lorentzian complex-time saddle,
covariant action rule for the laboratory pseudomanifold."""
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, dump
from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse import contours as C, tunnelling as T

warnings.simplefilter("ignore")
out = dict(monodromy=[], sheets=[], saddle=[], covariant=[], branch_locus=[], sector_sheets=[],
           zero_net=[], prescriptions=[])

models = [("gamma=-0.5 A", model_from_dimensionless(-0.5, 1.6 * 3.0)),
          ("gamma=0 A", model_from_dimensionless(0.0, 3.0 * 2.0)),
          ("gamma=0.6 A", model_from_dimensionless(0.6, 1.5 * 0.8)),
          ("SdS", ThinWallModel(1.0, 1.0, 0.1, 0.4))]

# 1. monodromy theorem: windings around every horizon point
for name, m in models:
    hz = list(m.horizons_out()) + ([1 / m.chi] if m.chi > 0 else [])
    for rh in hz:
        for t in (-1, 1, 2):
            r = C.monodromy_shift(m, rh, t)
            r.update(model=name, ratio=r["delta_B"] / r["prediction_per_turn"])
            out["monodromy"].append(r)
worst = max(abs(abs(r["ratio"]) - abs(r["turns"])) for r in out["monodromy"])
print("monodromy: max | |dB/pred| - |n| | =", worst, " max |dRe| =",
      max(abs(r["delta_Re"]) for r in out["monodromy"]))

# 2. sheet theorem: analytic branch <-> sign-free constraint, principal <-> sign-dressed
for name, m in models[:3]:
    s = C.principal_branch_constraint_violation(m)
    s["model"] = name
    out["sheets"].append(s)
    print("sheets", name, s)

# 2b. branch locus of the true Hamiltonian, sector sheets, zero-net periods
for name, m in models:
    b = C.branch_locus_check(m)
    b.update(model=name, flip_class=C.flip_class(m))
    out["branch_locus"].append(b)
    sh = C.sector_hamiltonians(m)
    sh["model"] = name
    out["sector_sheets"].append(sh)
    z = C.zero_net_shift(m, m.horizons_out()[0], 1 / m.chi)
    z["model"] = name
    out["zero_net"].append(z)
    print(f"branch locus {name}: class={b['flip_class']} radicand(R_f)={b['relative_radicand']:.1e} "
          f"|H_-(R_f)-M|={b['H_minus_minus_M']:.1e} |H_+(R_f)-M|={b['H_plus_minus_M']:.1e}; "
          f"zero-net dB={z['delta_B']:.8f} period={z['period']:.8f} endpoint dP={z['endpoint_dP']:.1e}")
    for side in ("R1", "R2"):
        print(f"   {side} eps={sh[side]['eps']} artanh |H_-+ - M|={sh[side]['artanh']}  sorkin={sh[side]['sorkin']}")

# 2c. the three connection rules on a grid: analytic (sign-free constraint),
#     principal (splice at the flip), sector switch at C = 0 (squared constraint)
for gam in np.linspace(-0.9, 0.9, 9):
    for fac in (1.05, 1.5, 3.0, 8.0):
        m = model_from_dimensionless(gam, fac * (2 - 2 * gam))
        out["prescriptions"].append(dict(gamma=float(gam), E_over_EAB=fac, flip_class=C.flip_class(m),
                                         B_analytic=T.B_canonical(m), B_principal=T.B_canonical(m, "principal"),
                                         B_switch=C.B_sector_switch(m),
                                         A2_over_4G=float(np.pi * m.turning_points()[-1] ** 2 / m.G)))
pr = out["prescriptions"]
print("prescriptions: single-flip models", sum(p["flip_class"] == "single" for p in pr),
      " double-flip", sum(p["flip_class"] == "double" for p in pr),
      " B_principal<0 in", sum(p["B_principal"] < 0 for p in pr), "of", len(pr),
      " 0<B_switch<B_analytic in all single:",
      all(0 < p["B_switch"] < p["B_analytic"] for p in pr if p["flip_class"] == "single"))

# 3. Lorentzian complex-time saddle: T* = dS/dM, S = int P dR along the barrier.
#    Under the barrier P is purely imaginary, so T*_barrier = i dImS/dM and
#    2 Im S = -B~ - A/4G (orientation fixed by B > 0): Im T* = -(1/2) dB/dM.
for name, m in models[:3]:
    h = 1e-5 * m.M
    def ImS(Mv):
        mm = m.with_mass(Mv)
        I, _ = C.continue_P(mm, C.real_barrier_path(mm))
        return abs(np.imag(I))
    dImS = (ImS(m.M + h) - ImS(m.M - h)) / (2 * h)
    pcs = T.euclidean_pieces(m)
    R2 = m.turning_points()[-1]
    dend = (T.endpoint_area_term(m.with_mass(m.M + h)) - T.endpoint_area_term(m.with_mass(m.M - h))) / (2 * h)
    out["saddle"].append(dict(model=name, Im_Tstar=dImS, minus_half_dtW_plus_half_dA=(-pcs["dt_W"] + dend) / 2,
                              beta_half=4 * np.pi * m.G * m.M,
                              shift_per_BH_winding=(np.pi / m.G) * 2 * R2 * 0 - 4 * np.pi * m.G * m.M))
    print(f"saddle {name}: Im T* = {dImS:.8f}  vs  (-dt_W + dA/dM/4G)/2 = {(-pcs['dt_W'] + dend) / 2:.8f};"
          f"  beta/2 = {4*np.pi*m.G*m.M:.5f}")

# 4. covariant action rule for the laboratory (regime A) pseudomanifold:
#    B = [S_in + S_wall] + [(M/2) dt_inf - M dt_inf] + A_BH/4G + (folds and dS bolts)/4G,
#    dt_inf = dt_W + 8 pi G M (one half-turn of the exterior at the tip per half-bounce)
for name, m in models[:3]:
    s = T.smarr_check(m)
    p = s["pieces"]
    dt_inf = p["dt_W"] + 8 * np.pi * m.G * m.M
    A_BH = 4 * np.pi * (2 * m.G * m.M) ** 2
    R2 = m.turning_points()[-1]
    fold = 4 * np.pi * R2**2 / (4 * m.G) if (m.A_in(R2) > 0 > m.A_out(R2)) else 0.0
    dsbolt = np.pi / (m.G * m.chi**2) if m.A_in(R2) < 0 else 0.0
    rhs = p["S_in"] + p["S_wall"] + (m.M / 2) * dt_inf - m.M * dt_inf + A_BH / (4 * m.G) + fold + dsbolt
    out["covariant"].append(dict(model=name, B=s["B"], S_in=p["S_in"], S_wall=p["S_wall"],
                                 GHY_infinity=(m.M / 2) * dt_inf, background=-m.M * dt_inf,
                                 BH_bolt=A_BH / (4 * m.G), wall_fold=fold, dS_bolt=dsbolt,
                                 sum=rhs, residual=s["B"] - rhs))
    print(f"covariant rule {name}: B={s['B']:.10f} sum={rhs:.10f} residual={s['B']-rhs:.1e}")

dump("phase3_quantum_transition.json", out)

# figure: complex R-plane structure and exponents of the winding sectors
m = models[0][1]
R1, R2 = m.turning_points()
rh, rc = 2 * m.G * m.M, 1 / m.chi
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for path, lab, col in ((C.real_barrier_path(m, 400), "real barrier contour (n=0)", "k"),
                       (C.horizon_loop_path(m, rh, 1, 4000), "loop around 2GM (n_BH=+1)", "C0"),
                       (C.horizon_loop_path(m, rc, 1, 4000), "loop around 1/chi (n_dS=+1)", "C3")):
    ax[0].plot(path.real, path.imag, color=col, lw=1.2, label=lab)
ax[0].plot([R1, R2], [0, 0], "o", color="#c0392b", label="turning points (sqrt)")
ax[0].plot([rh, rc], [0, 0], "x", color="#2c3e50", ms=9, mew=2, label="horizons (log): P -> P +- i pi R/G")
ax[0].plot([m.flip_radius_out()], [0], "s", color="#8e44ad", label="eps_out flip (regular point)")
ax[0].axhline(0, color="k", lw=0.4)
ax[0].set_xlabel("Re R")
ax[0].set_ylabel("Im R")
ax[0].legend(fontsize=7)
ax[0].set_title("Riemann-surface structure of the reduced momentum P(R)", fontsize=9)
for fac, ls in ((1.5, "-"), (8.0, "--")):
    rows = [p for p in pr if p["E_over_EAB"] == fac]
    g = [p["gamma"] for p in rows]
    ax[1].plot(g, [p["B_analytic"] for p in rows], "k" + ls, marker="o", ms=3,
               label=f"analytic (sign-free constraint), E/E_AB={fac}")
    ax[1].plot(g, [p["B_principal"] for p in rows], color="#888888", ls=ls, marker="s", ms=3,
               label=f"principal (splice at flip), E/E_AB={fac}")
    ax[1].plot(g, [p["B_switch"] for p in rows], color="#c0392b", ls=ls, marker="^", ms=3,
               label=f"sector switch at C=0 (single flip), E/E_AB={fac}")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xlabel("gamma")
ax[1].set_ylabel("B (units G = 1, kappa = 1)")
ax[1].legend(fontsize=6)
ax[1].set_title("exponents of the three connection rules", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "phase3_monodromy.png"), dpi=150)
plt.close(fig)
