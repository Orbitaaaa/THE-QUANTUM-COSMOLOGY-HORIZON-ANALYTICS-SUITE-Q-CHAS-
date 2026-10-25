import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate as C, gsl_solution as GS
from concurrent.futures import ProcessPoolExecutor
def one(seed):
    h, L = 0.08, 520.0
    o = C.nariai_condensate(0.0, frame="corner", h=h, lam_f=2e-7, Lu=L, Lv=L, clip=1, dphi0=seed)
    res, Cc = o["res"], o["cls"]; R, P, T = res["r"], res["phi"], res["T"]
    big = np.isfinite(R) & (R > 5 * o["r_N"])
    rows = {}
    for u in range(80, int(L) + 1, 40):
        i = min(int(u / h), R.shape[0] - 1); fin = np.isfinite(R[i])
        dau = big[i] & (P[i] < 0) & (Cc[i] == 1)
        rows[u] = dict(r_max=round(float(np.nanmax(R[i])), 1) if fin.any() else None, n_daughter=int(dau.sum()),
                       max_phi_v=round(float(np.nanmax(P[i]) / 0.0015), 3) if fin.any() else None,
                       T_max=round(float(np.nanmax(T[i])), 1) if fin.any() else None)
    tv = (P > 0) & np.isfinite(R)
    eh = GS.parent_event_horizon(o)
    return dict(seed=seed, grew=bool(o["grew"]), T_exit_pred=3 / 900 / (8 * np.pi) / (3 * 2e-7), rows=rows, parent=eh,
                TV_untrapped_far=int(np.sum(tv & (Cc >= 0) & (R > 5 * o["r_N"]))))
with ProcessPoolExecutor(2) as ex:
    for r in ex.map(one, [1.5e-9, 0.0]):
        print(json.dumps(r), flush=True)
