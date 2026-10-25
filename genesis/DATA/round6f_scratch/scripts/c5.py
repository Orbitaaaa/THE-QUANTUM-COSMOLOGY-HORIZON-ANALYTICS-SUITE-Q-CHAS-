import numpy as np, warnings, time
warnings.simplefilter("ignore")
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
pot = TE.Potential(lam, v, rho_f)
pl = F.Plasma(rho_f, v, lam, f_H=0.98, x_c=8.0)
Th = 1.02*pl.spinodal()*pl.Tc
ball = lambda r, R=19.0, w=1.0: pl.T_floor + (Th - pl.T_floor)*0.5*(1 - np.tanh((np.asarray(r, float) - R)/w))
for h in (0.04, 0.02):
    t0 = time.time()
    Z0 = F.centre_cone_data(h, 80.0, pl, phi_r=lambda r: v + 0*np.asarray(r, float), T_r=ball)
    resA = F.evolve_centre(Z0, h, pot, pl, n_rays=int(30/h), lam_f=1e-6)
    i1, j1 = int(2.0/h), int(32.0/h)
    A, B = F.handoff(resA, i1, j1, Lu=26.0, Lv=40.0)
    resB = F.evolve_corner(A, B, h, pot, pl, lam_f=1e-6)
    ov = F.overlap_check(resA, resB, i1, j1)
    print(f"h={h}: cone min r_v {np.min(np.gradient(Z0[0], h)):.3f}; corner r={A[0,0]:.3f} phi/v={A[2,0]/v:.3f}; overlap:",
          {k: (f"{val:.2e}" if isinstance(val, float) else val) for k, val in ov.items()}, f"{time.time()-t0:.0f}s")
    print("   phi/v in ball core at u=20 (min,max over v<32):", np.nanmin(resA["phi"][int(20/h), int(21/h):j1])/v, np.nanmax(resA["phi"][int(20/h), int(21/h):j1])/v)
