import numpy as np, warnings, json, time
warnings.simplefilter("ignore")
from babyuniverse import collapse as CO
from concurrent.futures import ProcessPoolExecutor
pot, info = CO.potential(eps=2.0)
lam, vphi, rho_f = pot
mt = np.sqrt(2 * lam * vphi**2 + 1.5 * rho_f / vphi**2)
def one(a):
    A, R0, w, kf = a
    h = 0.02
    phi = lambda r: vphi * (1 - A * np.exp(-((r - R0) / w) ** 2) * np.cos(kf * mt * (r - R0)))
    chi = lambda r: 0 * r
    t = time.time()
    res = CO.evolve(phi, chi, pot, h, v_max=2 * R0 + 6 * w + 40, gauge="alpha1")
    R, P = res["r"], res["phi"]
    out = dict(A=A, R0=R0, w=w, kf=kf, min_phi_data=float(np.min(phi(np.linspace(0, 3 * R0, 4000))) / vphi), t=round(time.time() - t))
    # along rays u = const (ingoing index i): largest areal radius of the false-basin region (phi < -0.375 v) at late times
    sizes = []
    for i in range(0, R.shape[0], R.shape[0] // 16):
        row = P[i]; rr = R[i]
        fb = np.isfinite(row) & (row < -0.375 * vphi)
        sizes.append((round(i * h, 1), float(rr[fb].max()) if fb.any() else 0.0, int(fb.sum())))
    out["false_basin_rmax_by_u"] = sizes
    m, ru, rv = CO.misner_sharp(res)
    out["max_2m_over_r"] = float(np.nanmax(2 * m / R))
    return out
jobs = [(1.2, 10.0, 2.0, 1.0), (1.2, 10.0, 2.0, 0.5), (1.3, 15.0, 3.0, 0.7), (1.0, 15.0, 4.0, 0.7)]
with ProcessPoolExecutor(4) as ex:
    for r in ex.map(one, jobs):
        print(json.dumps(r), flush=True)
