import numpy as np, warnings, json
warnings.simplefilter("ignore")
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from babyuniverse import formation as F
from concurrent.futures import ProcessPoolExecutor
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
def U(p): x = p/v; return lam/4*(p*p - v*v)**2 + rho_f*(0.5 - 0.75*x + 0.25*x**3)
def dU(p): return lam*p*(p*p - v*v) - rho_f*(0.75 - 0.75*(p/v)**2)/v
def restores(pl, E):
    Tof = lambda p, pd: brentq(lambda T: pl.eos(T, p)[0] - (E - 0.5*pd*pd - U(p)), 1e-8, 20)
    try:
        sol = solve_ivp(lambda t, y: [y[1], -pl.Gamma*y[1] - dU(y[0]) - pl.eos(Tof(y[0], y[1]), y[0])[2]], (0, 80), [v, 0.0], max_step=0.02, rtol=1e-8)
    except Exception:
        return False
    return sol.y[0][-1] < -0.9*v
def scan(a):
    fH, xc = a
    pl = F.Plasma(rho_f, v, lam, f_H=fH, x_c=xc)
    lo, hi = 1.0, 6.0
    if not restores(pl, hi*rho_f):
        return dict(f_H=fH, x_c=xc, E_min=None)
    for _ in range(12):
        mid = 0.5*(lo + hi)
        if restores(pl, mid*rho_f): hi = mid
        else: lo = mid
    return dict(f_H=fH, x_c=xc, E_min_over_rho_f=hi, thermo_bound=1/(1 - (1 - fH)**(1/3)))
jobs = [(fH, xc) for fH in (0.95, 0.99, 0.999) for xc in (3.0, 8.0, 16.0, 30.0)]
with ProcessPoolExecutor(4) as ex:
    for r in ex.map(scan, jobs):
        print(json.dumps(r), flush=True)
print("normality at R = 0.62 l_H requires E < (l_H/R)^2 rho_f = %.2f rho_f (spacelike estimate; the cone is stricter)" % (1/0.62**2))
