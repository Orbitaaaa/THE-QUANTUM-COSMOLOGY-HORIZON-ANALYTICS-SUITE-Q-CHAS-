import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate as C
from concurrent.futures import ProcessPoolExecutor
def one(a):
    h, L, n = a
    o = C.nariai_condensate(n, frame="corner", h=h, Lu=L, Lv=L)
    return dict(h=h, L=L, n=n, grew=o["grew"], n_big=o["n_big"], hidden=o["true_vacuum_untrapped_large_r"] == 0)
jobs = [(0.04, 120.0, -5.45e-6), (0.04, 120.0, -5.55e-6), (0.02, 120.0, -4.6e-6), (0.02, 120.0, -5.0e-6)]
with ProcessPoolExecutor(2) as ex:
    for r in ex.map(one, jobs):
        print(json.dumps(r), flush=True)
