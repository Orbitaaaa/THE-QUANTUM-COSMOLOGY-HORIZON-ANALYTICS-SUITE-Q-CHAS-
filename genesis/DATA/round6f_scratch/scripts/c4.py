import numpy as np, warnings, time
warnings.simplefilter("ignore")
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
pot = TE.Potential(lam, v, rho_f)
pl = F.Plasma(rho_f, v, lam)
T1 = pl.spinodal()
Th = 1.05*T1*pl.Tc
ball = lambda r, R=18.5, w=1.0: pl.T_floor + (Th - pl.T_floor)*0.5*(1 - np.tanh((np.asarray(r, float) - R)/w))
for h in (0.08, 0.04, 0.02):
    t0 = time.time()
    Z0 = F.centre_cone_data(h, 60.0, pl, phi_r=lambda r: v + 0*np.asarray(r, float), T_r=ball)
    res = F.evolve_centre(Z0, h, pot, pl, n_rays=int(25/h))
    cr = F.constraint_residual(res, pl, r_lo=2.0, s_max=3.0)
    R = res["r"]
    m = 0.5*R*(1 + 4*np.gradient(R, h, axis=0)*np.gradient(R, h, axis=1)*np.exp(-2*res["s"]))
    print(f"h={h}: Evv {cr['Evv']:.2e} Euu {cr['Euu']:.2e} rel {cr['Evv_rel']:.1e} {cr['Euu_rel']:.1e}; max 2m/r {np.nanmax(2*m/R):.3f}; "
          f"phi at centre region u=20: {np.nanmin(res['phi'][int(20/h), int(20/h)+5:int(20/h)+50])/v:.3f}; {time.time()-t0:.0f}s")
