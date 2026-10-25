import numpy as np, warnings, time
warnings.simplefilter("ignore")
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
pot = TE.Potential(lam, v, rho_f)
def U(p): x = p/v; return lam/4*(p*p - v*v)**2 + rho_f*(0.5 - 0.75*x + 0.25*x**3)
def dU(p): return lam*p*(p*p - v*v) - rho_f*(0.75 - 0.75*(p/v)**2)/v
pl = F.Plasma(rho_f, v, lam)
T0 = 1.1*pl.spinodal()*pl.Tc
E = pl.eos(T0, v)[0] + U(v)
Tof = lambda p, pd: brentq(lambda T: pl.eos(T, p)[0] - (E - 0.5*pd*pd - U(p)), 1e-6, 10)
sol = solve_ivp(lambda t, y: [y[1], -pl.Gamma*y[1] - dU(y[0]) - pl.eos(Tof(y[0], y[1]), y[0])[2]], (0, 12), [v, 0.0],
                max_step=0.0005, rtol=1e-11, atol=1e-16, dense_output=True)
ph = lambda t: sol.sol(np.asarray(t, float))[0]
Tt = lambda t: Tof(*sol.sol(float(t)))
# 1. homogeneous restoration through the centre (G -> 0): on the cone, t = r (alpha = 1 gauge, r = v/2)
for h in (0.02, 0.01):
    Z0 = F.centre_cone_data(h, 16.0, pl, phi_r=lambda r: ph(r), T_r=lambda r: Tt(r), G=1e-10)
    t0 = time.time()
    res = F.evolve_centre(Z0, h, pot, pl, G=1e-10, clip=0)
    # axis: u = v = t;  off-axis point (i, i+20): t = (u + v)/2
    errs = []
    for t in (1.0, 2.0, 4.0, 7.0):
        i = int(t/h); k = i + 40
        tt = 0.5*(i + k)*h
        errs.append((res["phi"][i, k] - ph(tt))/v)
    print(f"h={h}: phi error off-axis (in v) at t~1,2,4,7: {np.round(errs, 5)}; axis T/T_ode at t=4: {np.exp(res['lnTf'][int(4/h), int(4/h)])/Tt(4):.5f}; {time.time()-t0:.0f}s")
    out = []
    for off in (0, 1, 2, 5, 20, 80):
        i = int(4.0/h); k = i + off
        out.append((off, round(float(np.exp(res["lnTf"][i, k]) / Tt(0.5*(i + k)*h) - 1), 5)))
    print("   T rel. error at t~4 vs distance from axis (grid points):", out)
