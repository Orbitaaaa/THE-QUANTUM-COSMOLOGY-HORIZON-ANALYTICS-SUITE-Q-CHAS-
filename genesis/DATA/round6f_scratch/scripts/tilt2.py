import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate as C
from concurrent.futures import ProcessPoolExecutor
def one(a):
    tag, lf, ls, h, L = a
    o = C.nariai_condensate(0.0, frame="corner", h=h, lam_f=lf, lam_s=ls, Lu=L, Lv=L)
    return dict(tag=tag, lam_f=lf, lam_s=ls, h=h, grew=o["grew"], n_big=o["n_big"], hidden=o["true_vacuum_untrapped_large_r"] == 0,
                constraint=o["constraint"])
jobs = [("cons_h08", 4e-7, 4e-7, 0.08, 40.0), ("cons_h04", 4e-7, 4e-7, 0.04, 40.0),
        ("drift_only", 4e-7, 0.0, 0.08, 200.0), ("source_only", 0.0, 4e-7, 0.08, 200.0),
        ("drift_only", 1.6e-6, 0.0, 0.08, 200.0), ("source_only", 0.0, 1.6e-6, 0.08, 200.0)]
with ProcessPoolExecutor(2) as ex:
    for r in ex.map(one, jobs):
        print(json.dumps(r), flush=True)
