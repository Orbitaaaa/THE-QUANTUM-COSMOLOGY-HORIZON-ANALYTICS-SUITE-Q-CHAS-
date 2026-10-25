import numpy as np, warnings, json, sys
warnings.simplefilter("ignore")
from babyuniverse import condensate as C, gsl_solution as GS
from concurrent.futures import ProcessPoolExecutor
S = "/tmp/claude-0/-home-user-THE-QUANTUM-COSMOLOGY-HORIZON-ANALYTICS-SUITE-Q-CHAS-/5d1f4f78-e2f6-5e44-a419-ab5b29ff8104/scratchpad/r7/"
def one(a):
    lf, L = a
    h = 0.08
    o = C.nariai_condensate(0.0, frame="corner", h=h, lam_f=lf, Lu=L, Lv=L, clip=1, dphi0=1.5e-9)
    res = o["res"]; R, P, Cc, T = res["r"], res["phi"], o["cls"], res["T"]
    big = np.isfinite(R) & (R > 5 * o["r_N"])
    rho_f = 3 / 900 / (8 * np.pi)
    out = dict(lam_f=lf, L=L, grew=o["grew"], T_exit_pred=rho_f / (3 * lf), rows={})
    for u in np.arange(40, L + 1, 40):
        i = min(int(u / h), R.shape[0] - 1); m = big[i]
        fin = np.isfinite(R[i])
        out["rows"][int(u)] = dict(n_big=int(m.sum()), frac_TV_big=float(np.mean(P[i][m] > 0)) if m.any() else None,
                                   v_max_finite=float(np.where(fin)[0].max() * h) if fin.any() else None,
                                   v_max_big=float(np.where(m)[0].max() * h) if m.any() else None,
                                   T_at_vmin_big=float(T[i][m][0]) if m.any() else None, r_max=float(np.nanmax(R[i])) if fin.any() else None)
    idx = np.argwhere(big)
    if len(idx):
        late = idx[np.argsort(idx[:, 0])[-4:]]
        out["late_points"] = [(float(i * h), float(j * h), round(GS._deficit(R[: i + 1, j]), 2), round(GS._deficit(R[i, : j + 1]), 2),
                               float(P[i, j])) for i, j in late]
    np.savez_compressed(S + f"exit_{lf:.0e}.npz", R=R[::4, ::4], P=P[::4, ::4], C=Cc[::4, ::4])
    return out
with ProcessPoolExecutor(2) as ex:
    for r in ex.map(one, [(4e-7, 320.0), (1e-6, 240.0)]):
        print(json.dumps(r), flush=True)
