import numpy as np, warnings
warnings.simplefilter("ignore")
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
pot = TE.Potential(lam, v, rho_f)
h = 0.02
for tag, kw, fl in (("default", {}, 0.02), ("floor 1e-6", {}, 1e-6), ("k_axis=1", dict(k_axis=1), 0.02), ("k_axis=4", dict(k_axis=4), 0.02),
                    ("theta=0.75", dict(theta=0.75), 0.02), ("n_iter=5", dict(n_iter=5), 0.02)):
    pl = F.Plasma(rho_f, v, lam, T_floor_over_Tc=fl)
    Z0 = F.centre_cone_data(h, 8.0, pl, phi_r=lambda r: v + 0*r, T_r=lambda r: pl.T_floor + 0*r, G=1e-10)
    res = F.evolve_centre(Z0, h, pot, pl, G=1e-10, clip=0, **kw)
    S = res["s"]; n = np.sum(np.isfinite(res["r"]).any(axis=1))
    ax = [abs(S[i, i]) for i in range(1, min(n, 390))]
    print(f"{tag}: rows {n}, max|s_axis| {max(ax):.2e}, first i with |s_axis|>1e-6: {next((i+1 for i, a in enumerate(ax) if a > 1e-6), None)}, max|eta_f| {np.nanmax(np.abs(res['etaf'])):.2e}")
