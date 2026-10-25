import numpy as np, warnings, time
warnings.simplefilter("ignore")
from babyuniverse import formation as F, throat_evolution as TE
pot0 = TE.Potential(0.0, 0.0015, 0.0)
pl = F.Plasma(rho_f=3/900/(8*np.pi), vphi=0.0015, lam=4*(3/900/(8*np.pi))/(4*0.0015**4), alpha_pt=3.0)
pl.aH = 0.0                                                     # light species only for the pure-fluid tests
# 1. uniform static plasma in (near-)flat space, G tiny: should stay uniform and at rest
T0 = pl.Tc
for G in (1e-8,):
    h = 0.05; N = 200
    A, B, _ = F.corner_data(20.0, -0.5, 0.5, h, N, N, pl, T_v=lambda x: T0 + 0*x, T_u=lambda x: T0 + 0*x, phi_v=lambda x: 0.0015+0*x, phi_u=lambda x: 0.0015+0*x, G=G)
    t = time.time(); res = F.evolve_corner(A, B, h, pot0, pl, G=G); 
    Tf = np.exp(res["lnTf"])
    print("uniform: max|T/T0-1|", np.nanmax(np.abs(Tf/T0-1)), "max|eta_f|", np.nanmax(np.abs(res["etaf"])), "time", round(time.time()-t, 1))
# 2. sound pulse: bump in T on the outgoing ray; track the peak speed in flat space (expect c_s = 1/sqrt(3))
h = 0.025; N = 800
bump = lambda x: T0 * (1 + 0.01 * np.exp(-((np.asarray(x, float) - 5.0) / 0.5) ** 2))
A, B, _ = F.corner_data(40.0, -0.5, 0.5, h, N, N, pl, T_v=bump, T_u=lambda x: T0 + 0*x, phi_v=lambda x: 0.0015+0*x, phi_u=lambda x: 0.0015+0*x, G=1e-8)
res = F.evolve_corner(A, B, h, pot0, pl, G=1e-8)
Tf = np.exp(res["lnTf"]) / T0 - 1
# along the line t = (u+v)/2 = const, find the two pulses in r = r0 + (v - u)/2
for tt in (6.0, 10.0):
    pts = []
    for i in range(N + 1):
        j = int(round(2 * tt / h)) - i
        if 0 <= j <= N:
            pts.append((40 + (j - i) * h / 2, Tf[i, j]))
    pts = np.array(pts)
    k = np.argsort(pts[:, 1])[-1]
    print("t", tt, "peak at r", pts[k, 0], "amp", pts[k, 1])
