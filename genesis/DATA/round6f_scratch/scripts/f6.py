import numpy as np, warnings, copy
warnings.simplefilter("ignore")
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
pot = TE.Potential(lam, v, rho_f)
def U(p): x = p/v; return lam/4*(p*p - v*v)**2 + rho_f*(0.5 - 0.75*x + 0.25*x**3)
def dU(p): return lam*p*(p*p - v*v) - rho_f*(0.75 - 0.75*(p/v)**2)/v
def compare(tag, pl, phi0, T0, h=0.01, tmax=3.0):
    N = int(tmax/h)
    A, B, _ = F.corner_data(60.0, -0.5, 0.5, h, N, N, pl, T_v=lambda x: T0+0*x, T_u=lambda x: T0+0*x, phi_v=lambda x: phi0+0*x, phi_u=lambda x: phi0+0*x, G=1e-10)
    res = F.evolve_corner(A, B, h, pot, pl, G=1e-10)
    E = pl.eos(T0, phi0)[0] + U(phi0)
    def rhs(t, y):
        p, pd = y
        T = brentq(lambda T: pl.eos(T, p)[0] - (E - 0.5*pd*pd - U(p)), 1e-5, 10)
        return [pd, -pl.Gamma*pd - dU(p) - pl.eos(T, p)[2]]
    sol = solve_ivp(rhs, (0, tmax), [phi0, 0.0], max_step=0.002, rtol=1e-10)
    ts = (0.25, 0.5, 1.0, 2.0, 3.0)
    print(tag, "PDE", np.round([res["phi"][int(t/h), int(t/h)]/v for t in ts], 4), " ODE", np.round([np.interp(t, sol.t, sol.y[0])/v for t in ts], 4))
pl = F.Plasma(rho_f, v, lam)
T0 = 1.1*pl.spinodal()*pl.Tc
p0 = copy.copy(pl); p0.aH = 0.0; p0.Gamma = 0.0
compare("no plasma force, no friction, phi0=0.5v", p0, 0.5*v, T0)
p1 = copy.copy(pl); p1.aH = 0.0
compare("no plasma force, friction,    phi0=0.5v", p1, 0.5*v, T0)
p2 = copy.copy(pl); p2.Gamma = 0.0
compare("plasma force, no friction,    phi0=+v  ", p2, v, T0)
compare("full,                         phi0=+v  ", pl, v, T0)
