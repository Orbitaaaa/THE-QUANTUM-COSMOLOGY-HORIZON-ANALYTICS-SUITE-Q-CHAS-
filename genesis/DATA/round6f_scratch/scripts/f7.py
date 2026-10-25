import numpy as np, warnings
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
Tof = lambda p, pd: brentq(lambda T: pl.eos(T, p)[0] - (E - 0.5*pd*pd - U(p)), 1e-5, 10)
sol = solve_ivp(lambda t, y: [y[1], -pl.Gamma*y[1] - dU(y[0]) - pl.eos(Tof(y[0], y[1]), y[0])[2]], (0, 12), [v, 0.0],
                max_step=0.001, rtol=1e-11, atol=1e-16, dense_output=True)
phi_t = lambda t: sol.sol(np.asarray(t, float))[0]
T_t = np.vectorize(lambda t: Tof(*sol.sol(t)))
for h in (0.02, 0.01, 0.005):
    N = int(12/h)
    A, B, _ = F.corner_data(60.0, -0.5, 0.5, h, N, N, pl, T_v=lambda x: T_t(np.asarray(x)/2), T_u=lambda x: T_t(np.asarray(x)/2),
                            phi_v=lambda x: phi_t(np.asarray(x)/2), phi_u=lambda x: phi_t(np.asarray(x)/2), G=1e-10)
    res = F.evolve_corner(A, B, h, pot, pl, G=1e-10, clip=0)
    ts = np.array([1.0, 3.0, 6.0, 9.0, 12.0])
    dphi = [res["phi"][int(t/h), int(t/h)]/v - phi_t(t)/v for t in ts]
    dT = [np.exp(res["lnTf"][int(t/h), int(t/h)])/T_t(t) - 1 for t in ts]
    print(f"h={h}: phi error (in v) {np.round(dphi, 5)}   T rel error {np.round(dT, 5)}")
print("ODE phi/v at t=1,3,6,9,12:", np.round(phi_t(ts)/v, 4), " T/Tc:", np.round(T_t(ts)/pl.Tc, 4))
