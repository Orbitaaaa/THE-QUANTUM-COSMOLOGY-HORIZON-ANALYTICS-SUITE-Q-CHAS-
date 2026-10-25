import numpy as np, warnings, sys, json
warnings.simplefilter("ignore")
from babyuniverse import condensate as C
from concurrent.futures import ProcessPoolExecutor

def thr(kw):
    lo, hi = kw.pop("lo"), kw.pop("hi")
    for _ in range(kw.pop("n", 6)):
        mid = -np.sqrt(lo * hi)
        if C.nariai_condensate(mid, frame="corner", **kw)["grew"]:
            hi = -mid
        else:
            lo = -mid
    return dict(kw=kw, n_th=(lo, hi))

jobs = [dict(h=0.08, lo=5e-6, hi=1e-5), dict(h=0.04, lo=5e-6, hi=1e-5),
        dict(h=0.08, profile="lump", lump_width=4.0, lo=1e-6, hi=1e-3), dict(h=0.08, profile="lump", lump_width=8.0, lo=1e-6, hi=1e-3),
        dict(h=0.08, u_end_frac=0.3, lo=1e-6, hi=1e-3), dict(h=0.08, v_end=20.0, lo=1e-6, hi=1e-3),
        dict(h=0.08, profile="lump", lump_width=16.0, lo=1e-6, hi=1e-3), dict(h=0.08, u_end_frac=0.9, lo=1e-6, hi=1e-3)]
with ProcessPoolExecutor(4) as ex:
    for r in ex.map(thr, jobs):
        print(json.dumps(r), flush=True)
