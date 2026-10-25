import numpy as np, warnings
warnings.simplefilter("ignore")
from babyuniverse import condensate as C, throat_evolution as TE
for n0, lf in ((0.0, 0.0), (0.0, 4e-7)):
    h = 0.04
    o = C.nariai_condensate(n0, frame="corner", h=h, lam_f=lf, Lu=40.0, Lv=40.0)
    res = o["res"]; R, S, P = res["r"], res["s"], res["phi"]
    T = C.stress(res)
    rvv = (R[:, 2:] - 2 * R[:, 1:-1] + R[:, :-2]) / h**2
    rv = (R[:, 2:] - R[:, :-2]) / (2 * h); sv = (S[:, 2:] - S[:, :-2]) / (2 * h); pv = (P[:, 2:] - P[:, :-2]) / (2 * h)
    r = R[:, 1:-1]
    E = rvv - 2 * sv * rv + 4 * np.pi * r * pv**2 + T["T2vv"][:, 1:-1] / r
    E[~((r > 0.5 * o["r_N"]) & (r < 3 * o["r_N"]))] = np.nan
    i, j = np.unravel_index(np.nanargmax(np.abs(E)), E.shape)
    print("n0", n0, "lam_f", lf, "max|E|", np.nanmax(np.abs(E)), "at u", i * h, "v", (j + 1) * h, "r", r[i, j], "phi", P[i, j + 1], "s", S[i, j + 1],
          "v_wall", o["v_wall"], "p99", np.nanpercentile(np.abs(E), 99))
