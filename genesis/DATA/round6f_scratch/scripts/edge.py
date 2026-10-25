import numpy as np, warnings, time, json
warnings.simplefilter("ignore")
from scipy.optimize import brentq
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
pot = TE.Potential(lam, v, rho_f)
pl = F.Plasma(rho_f, v, lam, f_H=0.999, x_c=30.0)
h = 0.04
for Ef, R, w in ((2.05, 17.5, 0.5), (2.05, 17.0, 0.5), (2.3, 16.0, 0.5)):
    Th = brentq(lambda T: pl.eos(T, v)[0] - Ef*rho_f, pl.Tc*0.5, pl.Tc*5)
    ball = lambda r: pl.T_floor + (Th - pl.T_floor)*0.5*(1 - np.tanh((np.asarray(r, float) - R)/w))
    Z0 = F.centre_cone_data(h, 120.0, pl, phi_r=lambda r: v + 0*np.asarray(r, float), T_r=ball)
    rvmin = np.min(np.gradient(Z0[0], h))
    t0 = time.time()
    res = F.evolve_centre(Z0, h, pot, pl, n_rays=int(60/h))
    Rr, P, Tf = res["r"], res["phi"], np.exp(res["lnTf"])
    ru = np.gradient(Rr, h, axis=0); rv = np.gradient(Rr, h, axis=1)
    print(f"E={Ef} R={R}: cone min r_v {rvmin:+.3f}  ({time.time()-t0:.0f}s)")
    for t in (5, 10, 20, 30, 40, 50):
        # on the slice t = (u+v)/2: radius of the false-vacuum region and its plasma; trapped spheres
        pts = [(i, int(round(2*t/h)) - i) for i in range(Rr.shape[0]) if 0 <= int(round(2*t/h)) - i < Rr.shape[1] and int(round(2*t/h)) - i >= i]
        if not pts: continue
        ii, jj = np.array(pts).T
        r_, p_, T_ = Rr[ii, jj], P[ii, jj]/v, Tf[ii, jj]/pl.Tc
        ok = np.isfinite(r_)
        fv = ok & (p_ < -0.5)
        trapped = np.any((ru[ii, jj] < 0) & (rv[ii, jj] < 0) & ok)
        print(f"   t={t}: false-vacuum radius {r_[fv].max() if fv.any() else 0:.2f}, frac FV inside r<10: {np.mean(p_[ok & (r_ < 10)] < -0.5) if (ok & (r_<10)).any() else float('nan'):.2f}, "
              f"T/Tc at centre {T_[ok][0] if ok.any() else float('nan'):.2f}, max r on slice {np.nanmax(r_):.1f}, trapped {trapped}, computed pts {ok.sum()}/{len(ok)}")
