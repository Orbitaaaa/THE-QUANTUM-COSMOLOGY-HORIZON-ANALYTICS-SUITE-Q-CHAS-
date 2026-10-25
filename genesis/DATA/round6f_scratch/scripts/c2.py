import numpy as np, warnings
warnings.simplefilter("ignore")
from babyuniverse import formation as F, throat_evolution as TE
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
pot = TE.Potential(lam, v, rho_f)
pl = F.Plasma(rho_f, v, lam)
h = 0.02
# simplest: pure vacuum (field at +v, plasma at floor), G tiny: must stay trivial
Z0 = F.centre_cone_data(h, 8.0, pl, phi_r=lambda r: v + 0*r, T_r=lambda r: pl.T_floor + 0*r, G=1e-10)
res = F.evolve_centre(Z0, h, pot, pl, G=1e-10, clip=0)
R = res["r"]; fin = np.isfinite(R)
last = [np.where(fin[i])[0].max() if fin[i].any() else -1 for i in range(R.shape[0])]
print("vacuum test: rows computed", sum(1 for l in last if l >= 0), "of", R.shape[0], "; row 50 extends to j", last[50] if len(last) > 50 else None)
print("max |phi/v - 1|", np.nanmax(np.abs(res["phi"]/v - 1)), "s range", np.nanmin(res["s"]), np.nanmax(res["s"]), "max|r - (v-u)/2|",
      np.nanmax(np.abs(R - (np.arange(R.shape[1])[None, :] - np.arange(R.shape[0])[:, None])*h/2)))
i = 5
print("row 5 near axis: r", np.round(R[i, i:i+5], 4), "s", np.round(res["s"][i, i:i+5], 5), "lnT", np.round(res["lnTf"][i, i:i+5], 3), "etaf", np.round(res["etaf"][i, i:i+5], 3))
S = res["s"]
print("axis s for i=1..12:", np.round([S[i, i] for i in range(1, 13)], 5))
print("axis s at i=50,100,110,115:", [round(S[i, i], 4) for i in (50, 100, 110, 115) if i < S.shape[0]])
dev = np.abs(R - (np.arange(R.shape[1])[None, :] - np.arange(R.shape[0])[:, None])*h/2)
k = np.unravel_index(np.nanargmax(dev), dev.shape)
print("largest r deviation at (i,j)", k, "value", dev[k], "r", R[k], "s", S[k])
print("Z0 cone r first values", np.round(Z0[0, :6], 5), "lnT", np.round(Z0[6, :4], 3))
