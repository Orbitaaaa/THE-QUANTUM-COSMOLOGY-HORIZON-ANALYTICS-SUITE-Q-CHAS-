import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate as C, gsl_solution as GS
from concurrent.futures import ProcessPoolExecutor
def one(a):
    tag, eps, dphi = a
    o = C.nariai_condensate(-8e-6, frame="corner", h=0.08, eps=eps, dphi0=dphi)
    res = o["res"]; R, P, Cc = res["r"], res["phi"], o["cls"]; h = 0.08
    big = np.isfinite(R) & (R > 5 * o["r_N"])
    out = dict(tag=tag, grew=o["grew"], n_big=o["n_big"])
    # late-time state of the large-r region (the daughter): rows u = 100, 150, 199
    for u in (60, 100, 150, 199):
        i = int(u / h); m = big[i]
        out[f"u{u}"] = dict(n_big=int(m.sum()), frac_true_vac=float(np.mean(P[i][m] > 0)) if m.any() else None,
                           frac_anti=float(np.mean(Cc[i][m] == 1)) if m.any() else None, v_max_big=float(np.where(m)[0].max() * h) if m.any() else None,
                           r_max=float(np.nanmax(R[i])) if np.isfinite(R[i]).any() else None)
    # past boundaries of the latest large-r points
    idx = np.argwhere(big)
    late = idx[np.argsort(idx[:, 0])[-5:]]
    out["late_points"] = [(float(i * h), float(j * h), round(GS._deficit(R[: i + 1, j]), 2), round(GS._deficit(R[i, : j + 1]), 2)) for i, j in late]
    return out
jobs = [("eternal", 4.0, 0.0), ("hilltop_mu10", 5.334003623997308, 1.5e-6), ("hilltop_mu30", 5.33534471090311, 1.5e-6)]
with ProcessPoolExecutor(3) as ex:
    for r in ex.map(one, jobs):
        print(json.dumps(r), flush=True)
