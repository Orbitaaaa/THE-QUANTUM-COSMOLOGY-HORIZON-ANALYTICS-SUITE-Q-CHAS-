import numpy as np, warnings
warnings.simplefilter("ignore")
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
pot = TE.Potential(lam, v, rho_f)
pl = F.Plasma(rho_f, v, lam, alpha_pt=3.0)
M, r0 = 1.0, 8.0; f0 = 1 - 2*M/r0
h = 0.05; N = int(6/h)
lumpT = lambda x: pl.T_floor + 1.4*pl.Tc*np.exp(-((np.asarray(x, float)-2.0)/1.0)**2)
A, B, _ = F.corner_data(r0, -np.sqrt(f0)/2, np.sqrt(f0)/2, h, N, N, pl, T_v=lumpT, T_u=lumpT, phi_v=lambda x: v+0*x, phi_u=lambda x: v+0*x)
res = F.evolve_corner(A, B, h, pot, pl)
P = res["phi"]/v; T = np.exp(res["lnTf"])/pl.Tc; ef = res["etaf"]
bad = np.argwhere(np.abs(P - 1) > 0.5)
print("first bad (i,j):", bad[np.argsort(bad.sum(axis=1))[:3]] if len(bad) else None)
if len(bad):
    i, j = bad[np.argsort(bad.sum(axis=1))[0]]
    sl = (slice(max(i-3,0), i+2), slice(max(j-3,0), j+2))
    print("phi/v\n", np.round(P[sl], 3)); print("T/Tc\n", np.round(T[sl], 4)); print("eta_f\n", np.round(ef[sl], 3))
print("T range", np.nanmin(T), np.nanmax(T), "eta_f range", np.nanmin(ef), np.nanmax(ef))
