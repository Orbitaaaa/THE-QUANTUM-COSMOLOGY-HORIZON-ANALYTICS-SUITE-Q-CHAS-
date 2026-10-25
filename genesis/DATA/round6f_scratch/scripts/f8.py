import numpy as np, warnings, copy
warnings.simplefilter("ignore")
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
pot = TE.Potential(lam, v, rho_f)
def U(p): x = p/v; return lam/4*(p*p - v*v)**2 + rho_f*(0.5 - 0.75*x + 0.25*x**3)
def dU(p): return lam*p*(p*p - v*v) - rho_f*(0.75 - 0.75*(p/v)**2)/v
def test(tag, pl, phi0, T0, h=0.01, tmax=4.0):
    E = pl.eos(T0, phi0)[0] + U(phi0)
    Tof = lambda p, pd: brentq(lambda T: pl.eos(T, p)[0] - (E - 0.5*pd*pd - U(p)), 1e-6, 10)
    sol = solve_ivp(lambda t, y: [y[1], -pl.Gamma*y[1] - dU(y[0]) - pl.eos(Tof(y[0], y[1]), y[0])[2]], (0, tmax), [phi0, 0.0],
                    max_step=0.0005, rtol=1e-11, atol=1e-16, dense_output=True)
    ph = lambda t: sol.sol(np.asarray(t, float))[0]
    Tt = np.vectorize(lambda t: Tof(*sol.sol(t)))
    N = int(tmax/h)
    A, B, _ = F.corner_data(60.0, -0.5, 0.5, h, N, N, pl, T_v=lambda x: Tt(np.asarray(x)/2), T_u=lambda x: Tt(np.asarray(x)/2),
                            phi_v=lambda x: ph(np.asarray(x)/2), phi_u=lambda x: ph(np.asarray(x)/2), G=1e-10)
    res = F.evolve_corner(A, B, h, pot, pl, G=1e-10, clip=0)
    ts = np.array([0.1, 0.25, 0.5, 1.0, 2.0, 4.0])
    print(f"{tag}: max|dphi|/v over diag {np.nanmax(np.abs([res['phi'][k,k]-ph(k*h) for k in range(N+1)]))/v:.2e}; at t={ts}: ",
          np.round([(res['phi'][int(t/h), int(t/h)] - ph(t))/v for t in ts], 4))
pl = F.Plasma(rho_f, v, lam)
T0 = 1.1*pl.spinodal()*pl.Tc
p0 = copy.copy(pl); p0.aH = 1e-30; p0.Gamma = 0.0
test("scalar only (phi0=0.5v)", p0, 0.5*v, T0)
p1 = copy.copy(pl); p1.aH = 1e-30
test("scalar+friction (phi0=0.5v)", p1, 0.5*v, T0)
p2 = copy.copy(pl); p2.Gamma = 0.0
test("scalar+plasma force (phi0=v)", p2, v, T0)
