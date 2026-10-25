import numpy as np, warnings
warnings.simplefilter("ignore")
from babyuniverse import condensate as C
for n0, lf in ((-6e-6, 0.0), (0.0, 4e-7)):
    for h in (0.08, 0.04, 0.02):
        o = C.nariai_condensate(n0, frame="corner", h=h, lam_f=lf, Lu=40.0, Lv=40.0)
        res = o["res"]; R, S, P = res["r"], res["s"], res["phi"]
        T = C.stress(res); out = []
        for ax in (1, 0):
            sl = lambda a, k: np.take(a, np.arange(1 + k, a.shape[ax] - 1 + k), axis=ax)
            rxx = (sl(R, 1) - 2 * sl(R, 0) + sl(R, -1)) / h**2; rx = (sl(R, 1) - sl(R, -1)) / (2 * h)
            sx = (sl(S, 1) - sl(S, -1)) / (2 * h); px = (sl(P, 1) - sl(P, -1)) / (2 * h); r = sl(R, 0)
            E = rxx - 2 * sx * rx + 4 * np.pi * r * px**2 + sl(T["T2vv" if ax == 1 else "T2uu"], 0) / r
            m = (r > 0.5 * o["r_N"]) & (r < 3 * o["r_N"]) & (np.abs(sl(S, 0)) < 3) & (np.abs(sl(P, 0) + 0.0015) < 1.5e-5)
            out.append(np.nanmax(np.abs(np.where(m, E, np.nan))))
        print(f"n0={n0} lam_f={lf} h={h} false-vacuum region: max|Evv|={out[0]:.3e} max|Euu|={out[1]:.3e}", flush=True)
