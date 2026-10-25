"""Item 4: one continuous evolution from the route's ball (false-vacuum ball + imploding shell + wall, no black hole, zero condensate
charge, cold plasma remnant) through the self-generated NEC phase, tipping, a hidden daughter, its graceful exit and recollapse."""
import numpy as np, warnings, json, sys, time, pickle
warnings.simplefilter("ignore")
from scipy.optimize import brentq
from babyuniverse import condensate as C, formation as F, throat_evolution as TE, gsl_solution as GS
lam_f = float(sys.argv[1]); Lv = float(sys.argv[2]); h = float(sys.argv[3]); solver = sys.argv[4]; dphi0 = float(sys.argv[5]) if len(sys.argv) > 5 else 1.5e-9
Lu = float(sys.argv[6]) if len(sys.argv) > 6 else 30.0
clip = int(sys.argv[7]) if len(sys.argv) > 7 else 1
eps_arg = float(sys.argv[8]) if len(sys.argv) > 8 else 4.0
gstar = float(sys.argv[9]) if len(sys.argv) > 9 else 1e-12
rinf = float(sys.argv[10]) if len(sys.argv) > 10 else 1e5
E, l_H, r0, eps, vphi, v_shell, w_shell, wall_frac, delta_target = 7.258540773391722, 30.0, 15.0, eps_arg, 0.0015, 4.5, 0.25, 0.6, 4e-3
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
# false vacuum displaced by dphi0 inside the wall (the clipped potential holds phi = -v exactly; F42)
phi_v = lambda v: vphi * np.tanh((np.asarray(v, float) - v_wall) / w) + dphi0 * 0.5 * (1 - np.tanh((np.asarray(v, float) - v_wall) / w))
phi_u = lambda u: -vphi + dphi0 + 0 * np.asarray(u, float)
z = lambda x: 0 * np.asarray(x, float)
t0 = time.time()
if solver == "formation":
    pl = F.Plasma(rho_f, vphi, lam, g_star=gstar)                     # the thermal sector of item 1, at its floor temperature (cold remnant)
    A, B, info = F.corner_data(r0, -rv0, rv0, h, Nv, Nu, pl, phi_v=phi_v, phi_u=phi_u, extra_Nv=shell)
    res = F.evolve_corner(A, B, h, pot, pl, lam_f=lam_f, clip=clip, r_inf=rinf, phi_max=10 * vphi)
    cr = F.constraint_residual(res, pl, extra_Nv=shell, r_lo=0.5 * rN, r_hi=3 * rN, s_max=5.0)
else:
    A, B, info = C.condensate_data(0.0, 0.0, r0, -rv0, rv0, h, Nv, Nu, z, z, phi_v=phi_v, phi_u=phi_u, pot=pot, frame="corner", extra_Nv=shell)
    res = C.evolve(A, B, 0.0, 0.0, pot, h, r_inf=rinf, phi_max=10 * vphi, lam_f=lam_f, clip=clip)
    cr = C.constraint_residual(res, r_lo=0.5 * rN, r_hi=3 * rN, s_max=5.0)
secs = time.time() - t0
cl = TE.classify(res)
g = TE.inflated(res, cl, rN)
R, P, T = res["r"], res["phi"], res["T"]
ok = np.isfinite(R)
fv = ok & (np.abs(P / vphi + 1) < 0.5)
tv = ok & (np.abs(P / vphi - 1) < 0.5)
big = ok & (R > 3 * rN)
o = dict(res=res, cls=cl["cls"], r_N=rN, M_N=MN, pot=pot)
eh = GS.parent_event_horizon(o)
# exit: true vacuum reached at large r (inside the daughter); recollapse: daughter rows whose r peaks and then ends at r -> 0
exit_pts = big & tv
r_max_daughter = float(np.nanmax(np.where(big, R, np.nan))) if big.any() else None
rows = np.where(big.any(axis=1))[0]
recollapse = 0
for i in rows:
    js = np.where(ok[i])[0]
    rr = R[i, js]
    if rr.max() > 3 * rN and rr[-1] < 0.1 * rr.max():
        recollapse += 1
out = dict(eps=eps, clip=clip, solver=solver, lam_f=lam_f, Lv=Lv, Lu=Lu, h=h, dphi0=dphi0, secs=round(secs, 1), grew=bool(g["grew"]), n_big=int(g["n_big"]),
           r_max_daughter=r_max_daughter, T_max_false_vacuum=float(np.nanmax(np.where(fv & big, T, np.nan))) if (fv & big).any() else None,
           T_exit_predicted=float(rho_f * (16 / (3 * eps) - 1) / lam_f),
           exit_cells=int(exit_pts.sum()), first_exit_T=float(np.nanmin(np.where(exit_pts, T, np.nan))) if exit_pts.any() else None,
           daughter_rows=int(len(rows)), recollapsing_rows=int(recollapse), u_EH=eh.get("u_EH"), parent_deficit=eh.get("deficit"),
           Evv_rel=cr.get("Evv_rel"), Euu_rel=cr.get("Euu_rel"), normal_on_data=bool(np.all(info["ru_b"] < 0) and np.all(info["rv_a"] > 0)))
print(json.dumps(out), flush=True)
pickle.dump(dict(out=out, R=R[::4, ::4], P=P[::4, ::4], T=T[::4, ::4], h=h * 4), open(f"/tmp/claude-0/-home-user-THE-QUANTUM-COSMOLOGY-HORIZON-ANALYTICS-SUITE-Q-CHAS-/5d1f4f78-e2f6-5e44-a419-ab5b29ff8104/scratchpad/r7/cont_{solver}_{lam_f:g}_{Lv:g}_c{clip}_d{dphi0:g}_e{eps:g}_ri{rinf:g}.pkl", "wb"))
