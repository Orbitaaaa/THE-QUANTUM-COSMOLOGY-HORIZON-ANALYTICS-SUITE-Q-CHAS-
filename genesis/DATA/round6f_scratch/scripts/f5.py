import numpy as np, warnings
warnings.simplefilter("ignore")
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
pot = TE.Potential(lam, v, rho_f)
pl = F.Plasma(rho_f, v, lam)
T0 = 1.1*pl.spinodal()*pl.Tc
for h, th in ((0.02, 0.6), (0.02, 0.5), (0.01, 0.6), (0.01, 0.5), (0.005, 0.5)):
    N = int(12/h)
    A, B, _ = F.corner_data(60.0, -0.5, 0.5, h, N, N, pl, T_v=lambda x: T0+0*x, T_u=lambda x: T0+0*x, phi_v=lambda x: v+0*x, phi_u=lambda x: v+0*x, G=1e-10)
    res = F.evolve_corner(A, B, h, pot, pl, G=1e-10, theta=th)
    idx = [int(t/h) for t in (0.5, 1.0, 1.5, 2.0, 6.0, 11.0, 12.0)]
    print(f"h={h} th={th}: phi/v", np.round([res["phi"][k, k]/v for k in idx], 3), " T/Tc", np.round([np.exp(res["lnTf"][k, k])/pl.Tc for k in idx], 3))
