import numpy as np, warnings, json, sys
warnings.simplefilter("ignore")
from babyuniverse import collapse as CO
from concurrent.futures import ProcessPoolExecutor
pot, info = CO.potential(eps=2.0)
lam, vphi, rho_f = pot
mt = np.sqrt(2 * lam * vphi**2 + 1.5 * rho_f / vphi**2)
G = float(sys.argv[1]); n = int(sys.argv[2])
def one(a):
    A, R0, w, kf, seed = a
    h = 0.02
    phi = lambda r: vphi * (1 - A * np.exp(-((r - R0) / w) ** 2) * np.cos(kf * mt * (r - R0)))
    try:
        res = CO.evolve(phi, lambda r: 0 * r, pot, h, v_max=2 * R0 + 6 * w + 60, G=G, gauge="alpha1")
    except Exception as e:
        return dict(A=A, err=str(e))
    R, P = res["r"], res["phi"]
    u_pass = 2 * R0 + 4 * w + 5
    late = []
    for i in range(int(u_pass / h), R.shape[0]):
        fb = np.isfinite(P[i]) & (P[i] < -0.375 * vphi)
        if fb.any():
            late.append((i * h, float(R[i][fb].max())))
    u_last_row = float(np.max(np.where(np.isfinite(R).any(axis=1))[0]) * h)
    return dict(A=round(A, 3), R0=round(R0, 2), w=round(w, 2), kf=round(kf, 3), n_late_rows=len(late), u_pass=u_pass,
                last_u_false=late[-1][0] if late else None, r_false_max_late=max(x[1] for x in late) if late else 0.0, u_last_row=u_last_row)
rng = np.random.default_rng(int(G * 1000) + n)
jobs = [(rng.uniform(1.0, 1.8), rng.uniform(5, 20), rng.uniform(1, 5), rng.uniform(0.15, 1.2), i) for i in range(n)]
with ProcessPoolExecutor(4) as ex:
    res = list(ex.map(one, jobs))
res = [r for r in res if "err" not in r]
res.sort(key=lambda r: (r["n_late_rows"], r["r_false_max_late"]), reverse=True)
for r in res[:10]:
    print(json.dumps(r))
print("runs", len(res), "with any late false-basin row:", sum(r["n_late_rows"] > 0 for r in res))
