import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate as C
from concurrent.futures import ProcessPoolExecutor
def thr(h):
    lo, hi = 1e-8, 4e-7
    for _ in range(7):
        mid = np.sqrt(lo * hi)
        if C.nariai_condensate(0.0, frame="corner", h=h, lam_f=mid)["grew"]:
            hi = mid
        else:
            lo = mid
    return dict(h=h, lam_th=(lo, hi))
with ProcessPoolExecutor(2) as ex:
    for r in ex.map(thr, [0.08, 0.04]):
        print(json.dumps(r), flush=True)
