import numpy as np, warnings
warnings.simplefilter("ignore")
from scipy.integrate import solve_ivp
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; eps = 4.0
lam = 4*rho_f/(eps*v**4)
pot = TE.Potential(lam, v, rho_f)
pl = F.Plasma(rho_f, v, lam, alpha_pt=3.0)
# 1. constraint convergence with gravity: hot lump (both species) in Schwarzschild data, phi at +v
M, r0 = 1.0, 8.0; f0 = 1 - 2*M/r0
for h in (0.1, 0.05, 0.025):
    N = int(6/h)
    lumpT = lambda x: pl.T_floor + 1.4*pl.Tc*np.exp(-((np.asarray(x, float)-2.0)/1.0)**2)
    A, B, _ = F.corner_data(r0, -np.sqrt(f0)/2, np.sqrt(f0)/2, h, N, N, pl, T_v=lumpT, T_u=lumpT, phi_v=lambda x: v+0*x, phi_u=lambda x: v+0*x)
    res = F.evolve_corner(A, B, h, pot, pl)
    cr = F.constraint_residual(res, pl)
    print("h", h, "Evv %.2e Euu %.2e rel %.1e %.1e" % (cr["Evv"], cr["Euu"], cr["Evv_rel"], cr["Euu_rel"]), "min phi/v", np.nanmin(res["phi"])/v)
# 2. restoration in a uniform hot region, G -> 0: compare with the ODE phi'' + Gamma phi' + U' + F_phi = 0 (T fixed to leading order)
h = 0.02; N = 500
T0 = 1.4*pl.Tc
A, B, _ = F.corner_data(60.0, -0.5, 0.5, h, N, N, pl, T_v=lambda x: T0+0*x, T_u=lambda x: T0+0*x, phi_v=lambda x: v+0*x, phi_u=lambda x: v+0*x, G=1e-10)
res = F.evolve_corner(A, B, h, pot, pl, G=1e-10)
# at fixed v - u (same r), the field is a function of t = (u+v)/2 (data are homogeneous along both rays only at t = 0 ... compare far from edges)
ts, ph, Ts = [], [], []
for k in range(0, N + 1, 10):
    i = j = k
    ts.append(k*h); ph.append(res["phi"][i, j]/v); Ts.append(np.exp(res["lnTf"][i, j]))
def rhs(t, y):
    T = y[2]
    rho, p, Fp, Fpp = pl.eos(T, y[0])
    dU = lam*y[0]*(y[0]**2 - v**2) - rho_f*(0.75 - 0.75*(y[0]/v)**2)/v
    # energy: total rho_fluid + phi'^2/2 + U conserved (flat, homogeneous, no expansion)
    return [y[1], -pl.Gamma*y[1] - dU - Fp, 0.0]
# with T held fixed the ODE is only the leading-order comparison; report the restoration time and final value from both
sol = solve_ivp(rhs, (0, ts[-1]), [v, 0.0, T0], max_step=0.005, rtol=1e-9)
print("PDE: phi/v at t =", [(round(ts[q],1), round(ph[q],3)) for q in (0, 5, 10, 20, 30, 40, 50)])
print("ODE: phi/v at t =", [(round(t,1), round(np.interp(t, sol.t, sol.y[0])/v,3)) for t in (0, 1, 2, 4, 6, 8, 10)])
print("plasma T/T0 change along diagonal:", Ts[-1]/T0)
