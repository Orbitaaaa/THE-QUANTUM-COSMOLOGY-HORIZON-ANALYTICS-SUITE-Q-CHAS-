import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import collapse as CO
from concurrent.futures import ProcessPoolExecutor
pot, info = CO.potential(eps=2.0)
lam, vphi, rho_f = pot
mt = np.sqrt(2 * lam * vphi**2 + 1.5 * rho_f / vphi**2)
def one(a):
    A, R0, w, kf, G, extra = a
    h = 0.02
    phi = lambda r: vphi * (1 - A * np.exp(-((r - R0) / w) ** 2) * np.cos(kf * mt * (r - R0)))
    res = CO.evolve(phi, lambda r: 0 * r, pot, h, v_max=2 * R0 + 6 * w + 60 + extra, G=G, gauge="alpha1")
    R, P = res["r"], res["phi"]
    m, ru, rv = CO.misner_sharp(res)
    rows = []
    for u in np.arange(40, R.shape[0] * h, 8):
        i = int(u / h)
        ok = np.isfinite(P[i]) & np.isfinite(R[i])
        if not ok.any():
            break
        x = P[i][ok] / vphi; rr = R[i][ok]
        inner = rr < 12
        rows.append(dict(u=round(float(u), 1), r_range=(round(float(rr.min()), 2), round(float(rr.max()), 2)),
                         frac_false_inner=round(float(np.mean(x[inner] < -0.375)), 3) if inner.any() else None,
                         x_inner_mean=round(float(np.mean(x[inner])), 3) if inner.any() else None,
                         x_inner_std=round(float(np.std(x[inner])), 3) if inner.any() else None,
                         r_false_max=round(float(rr[x < -0.375].max()), 2) if (x < -0.375).any() else 0.0,
                         max_2m_r=round(float(np.nanmax((2 * m[i] / R[i])[ok])), 3) if np.isfinite(m[i]).any() else None,
                         trapped=bool(np.nanmin((ru[i] * rv[i])[ok]) < 0) if np.isfinite(ru[i]).any() else None))
    return dict(A=A, R0=R0, G=G, rows=rows)
jobs = [(1.07, 16.13, 3.89, 0.32, 0.001, 60), (1.07, 16.13, 3.89, 0.32, 1.0, 60), (1.058, 11.31, 1.18, 0.381, 1.0, 60)]
with ProcessPoolExecutor(3) as ex:
    for r in ex.map(one, jobs):
        print(json.dumps(r), flush=True)
