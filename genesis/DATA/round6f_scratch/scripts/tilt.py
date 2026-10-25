import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate as C
from concurrent.futures import ProcessPoolExecutor
def one(lf):
    o = C.nariai_condensate(0.0, frame="corner", h=0.08, lam_f=lf)
    res = o["res"]; n = C.stress(res)["n"]; R = res["r"]
    band = np.isfinite(R) & (np.abs(R - o["r_N"]) < 0.5 * (o["r_c"] - o["r_b"]))
    return dict(lam_f=lf, grew=o["grew"], n_big=o["n_big"], hidden=o["true_vacuum_untrapped_large_r"] == 0, r_max=o["r_max"],
                n_corner_u10=float(n[int(10 / 0.08), int(10 / 0.08)]), n_band_min=float(np.nanmin(np.where(band, n, np.nan))))
with ProcessPoolExecutor(2) as ex:
    for r in ex.map(one, [0.0, 1e-7, 2e-7, 4e-7, 8e-7, 1.6e-6]):
        print(json.dumps(r), flush=True)
