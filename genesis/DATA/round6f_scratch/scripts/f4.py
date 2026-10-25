import numpy as np, warnings
warnings.simplefilter("ignore")
from scipy.integrate import solve_ivp
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
pot = TE.Potential(lam, v, rho_f)
pl = F.Plasma(rho_f, v, lam)
T1 = pl.spinodal(); T0 = 1.1*T1*pl.Tc
rho0 = pl.eos(T0, v)[0]
print("T1/Tc", T1, "T0/Tc", T0/pl.Tc, "rho(T0,+v)/rho_f", rho0/rho_f)
h = 0.02; N = 600
A, B, _ = F.corner_data(60.0, -0.5, 0.5, h, N, N, pl, T_v=lambda x: T0+0*x, T_u=lambda x: T0+0*x, phi_v=lambda x: v+0*x, phi_u=lambda x: v+0*x, G=1e-10)
res = F.evolve_corner(A, B, h, pot, pl, G=1e-10)
ph = [res["phi"][k, k]/v for k in range(0, N+1, 50)]
Tt = [np.exp(res["lnTf"][k, k])/pl.Tc for k in range(0, N+1, 50)]
print("PDE diag t:", [round(k*h, 1) for k in range(0, N+1, 50)])
print(" phi/v:", np.round(ph, 3)); print(" T/Tc :", np.round(Tt, 3))
# ODE with energy conservation: rho_plasma(T, phi) + phi'^2/2 + U(phi) = const, friction heats the plasma
def U(p): x = p/v; return lam/4*(p*p - v*v)**2 + rho_f*(0.5 - 0.75*x + 0.25*x**3)
def dU(p): return lam*p*(p*p - v*v) - rho_f*(0.75 - 0.75*(p/v)**2)/v
Etot = rho0 + U(v)
def rhs(t, y):
    p, pd = y
    from scipy.optimize import brentq
    target = Etot - 0.5*pd*pd - U(p)
    T = brentq(lambda T: pl.eos(T, p)[0] - target, 1e-4, 10)
    return [pd, -pl.Gamma*pd - dU(p) - pl.eos(T, p)[2]]
sol = solve_ivp(rhs, (0, N*h), [v, 0.0], max_step=0.01, rtol=1e-8)
print("ODE phi/v:", np.round([np.interp(k*h, sol.t, sol.y[0])/v for k in range(0, N+1, 50)], 3))
