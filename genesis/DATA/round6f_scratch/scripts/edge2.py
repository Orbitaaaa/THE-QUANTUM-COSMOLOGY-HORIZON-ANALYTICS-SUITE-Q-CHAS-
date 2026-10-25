import numpy as np, warnings, time
warnings.simplefilter("ignore")
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; eps = 5.2
lam = 4*rho_f/(eps*v**4)
pot = TE.Potential(lam, v, rho_f)
pl = F.Plasma(rho_f, v, lam, f_H=0.999, x_c=30.0, Gamma_over_m=0.03)
U = lambda p: lam/4*(p*p - v*v)**2 + rho_f*(0.5 - 0.75*(p/v) + 0.25*(p/v)**3)
dU = lambda p: lam*p*(p*p - v*v) - rho_f*(0.75 - 0.75*(p/v)**2)/v
def restores(E):
    Tof = lambda p, pd: brentq(lambda T: pl.eos(T, p)[0] - (E - 0.5*pd*pd - U(p)), 1e-8, 20)
    s = solve_ivp(lambda t, y: [y[1], -pl.Gamma*y[1] - dU(y[0]) - pl.eos(Tof(y[0], y[1]), y[0])[2]], (0, 150), [v, 0.0], max_step=0.02, rtol=1e-8)
    k = np.where(s.y[0] < -0.9*v)[0]
    return (s.y[0][-1] < -0.9*v), (s.t[k[0]] if len(k) else None)
lo, hi = 1.0, 4.0
for _ in range(11):
    mid = 0.5*(lo + hi)
    if restores(mid*rho_f)[0]: hi = mid
    else: lo = mid
Emin = hi
print(f"eps={eps}, weak friction: E_min = {Emin:.3f} rho_f; restoration time at 1.03 E_min: {restores(1.03*Emin*rho_f)[1]}, at 1.1 E_min: {restores(1.1*Emin*rho_f)[1]}")
h = 0.04
for fac in (1.03, 1.1):
    Ef = fac*Emin
    Th = brentq(lambda T: pl.eos(T, v)[0] - Ef*rho_f, pl.Tc*0.3, pl.Tc*6)
    Rmax = None
    for R in np.arange(15.0, 22.0, 0.25):
        ball = lambda r, R=R: pl.T_floor + (Th - pl.T_floor)*0.5*(1 - np.tanh((np.asarray(r, float) - R)/0.5))
        if np.min(np.gradient(F.centre_cone_data(0.05, 3*R + 20, pl, phi_r=lambda r: v + 0*np.asarray(r, float), T_r=ball)[0], 0.05)) > 0:
            Rmax = R
    R = Rmax
    ball = lambda r: pl.T_floor + (Th - pl.T_floor)*0.5*(1 - np.tanh((np.asarray(r, float) - R)/0.5))
    Z0 = F.centre_cone_data(h, 140.0, pl, phi_r=lambda r: v + 0*np.asarray(r, float), T_r=ball)
    res = F.evolve_centre(Z0, h, pot, pl, n_rays=int(70/h))
    Rr, P = res["r"], res["phi"]
    ru = np.gradient(Rr, h, axis=0); rv = np.gradient(Rr, h, axis=1)
    hist = []
    for t in (10, 20, 30, 40, 50, 60, 70):
        pts = [(i, int(round(2*t/h)) - i) for i in range(Rr.shape[0]) if i <= int(round(2*t/h)) - i < Rr.shape[1]]
        ii, jj = np.array(pts).T
        r_, p_ = Rr[ii, jj], P[ii, jj]/v
        ok = np.isfinite(r_)
        fv = ok & (p_ < -0.5)
        hist.append((t, round(float(r_[fv].max()), 1) if fv.any() else 0.0, bool(np.any((ru[ii, jj] < 0) & (rv[ii, jj] < 0) & ok))))
    print(f"E = {fac} E_min = {Ef:.2f} rho_f, largest normal ball R = {R}: (t, false-vacuum radius, trapped):", hist)
