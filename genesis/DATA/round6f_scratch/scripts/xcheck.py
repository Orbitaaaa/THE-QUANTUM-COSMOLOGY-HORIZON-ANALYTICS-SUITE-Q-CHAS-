"""Item 4 cross-check: the route (ball + shell + wall, no black hole on the data) evolved by the independent formation solver (fluid +
scalar + condensate + clock, Preissmann fluid box, implicit scalar) with a negligible plasma, against condensate.evolve."""
import numpy as np, warnings, json, sys, time
warnings.simplefilter("ignore")
from scipy.optimize import brentq
from babyuniverse import condensate as C, formation as F, throat_evolution as TE
lam_f = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-6
h = float(sys.argv[2]) if len(sys.argv) > 2 else 0.04
E, l_H, r0, eps, vphi, Lu, Lv, v_shell, w_shell, wall_frac, delta_target = 7.258540773391722, 30.0, 15.0, 4.0, 0.0015, 30.0, 600.0, 4.5, 0.25, 0.6, 4e-3
Lam = 3 / l_H**2; rho_f = Lam / (8 * np.pi); lam = 4 * rho_f / (eps * vphi**4); w = np.sqrt(2) / (np.sqrt(lam) * vphi)
pot = TE.Potential(lam, vphi, rho_f)
rN, MN = 1 / np.sqrt(Lam), 1 / (3 * np.sqrt(Lam))
rv0 = np.sqrt(1 - Lam * r0**2 / 3) / 2
M = MN * (1 - delta_target)
rc = brentq(lambda r: 1 - 2 * M / r - Lam * r**2 / 3, rN, 3 * rN)
shell = lambda v: E / (np.sqrt(np.pi) * w_shell) * np.exp(-((np.asarray(v, float) - v_shell) / w_shell) ** 2)
Nv, Nu = int(Lv / h), int(Lu / h)
a0, _, _ = TE.characteristic_data(0.0, 0.0, r0, -rv0, rv0, h, Nv, Nu, phi_v=lambda v: -vphi + 0 * np.asarray(v, float),
                                  phi_u=lambda u: -vphi + 0 * np.asarray(u, float), N_v=shell)
j = np.where(a0[0] >= rN + wall_frac * (rc - rN))[0]
v_wall = j[0] * h
phi_v = lambda v: vphi * np.tanh((np.asarray(v, float) - v_wall) / w)
phi_u = lambda u: -vphi + 0 * np.asarray(u, float)
z = lambda x: 0 * np.asarray(x, float)
t0 = time.time()
A, B, info = C.condensate_data(0.0, 0.0, r0, -rv0, rv0, h, Nv, Nu, z, z, phi_v=phi_v, phi_u=phi_u, pot=pot, frame="corner", extra_Nv=shell)
ra = C.evolve(A, B, 0.0, 0.0, pot, h, r_inf=1e4, phi_max=10 * vphi, lam_f=lam_f)
t1 = time.time()
pl = F.Plasma(rho_f, vphi, lam, g_star=1e-12)
A2, B2, info2 = F.corner_data(r0, -rv0, rv0, h, Nv, Nu, pl, phi_v=phi_v, phi_u=phi_u, extra_Nv=shell)
rb = F.evolve_corner(A2, B2, h, pot, pl, lam_f=lam_f, clip=0, r_inf=1e4)
t2 = time.time()
ga = TE.inflated(ra, TE.classify(ra), rN); gb = TE.inflated(rb, TE.classify(rb), rN)
both = np.isfinite(ra["r"]) & np.isfinite(rb["r"])
reg = both & (ra["r"] > 0.5 * rN) & (ra["r"] < 3 * rN)
out = dict(lam_f=lam_f, h=h, data_r_maxdiff=float(np.max(np.abs(A[0] - A2[0]))), grew_condensate=bool(ga["grew"]), grew_formation=bool(gb["grew"]),
           n_big=(int(ga["n_big"]), int(gb["n_big"])), finite_cells=(int(np.isfinite(ra["r"]).sum()), int(np.isfinite(rb["r"]).sum())),
           r_rel_med=float(np.median(np.abs(rb["r"][reg] / ra["r"][reg] - 1))), r_rel_p95=float(np.percentile(np.abs(rb["r"][reg] / ra["r"][reg] - 1), 95)),
           phi_med_over_v=float(np.median(np.abs(rb["phi"][reg] - ra["phi"][reg])) / vphi), phi_p95_over_v=float(np.percentile(np.abs(rb["phi"][reg] - ra["phi"][reg]), 95) / vphi),
           D_rel_p95=float(np.percentile(np.abs(rb["D"][reg] - ra["D"][reg]), 95) / max(np.nanmax(np.abs(ra["D"][reg])), 1e-300)),
           secs=(round(t1 - t0, 1), round(t2 - t1, 1)))
print(json.dumps(out), flush=True)
