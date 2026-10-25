import numpy as np, warnings, json, gc
warnings.simplefilter("ignore")
from babyuniverse import condensate as C
for n in (-4.7e-6, -5.1e-6):
    o = C.nariai_condensate(n, frame="corner", h=0.02)
    print(json.dumps(dict(h=0.02, n=n, grew=o["grew"], n_big=o["n_big"], hidden=o["true_vacuum_untrapped_large_r"] == 0)), flush=True)
    del o; gc.collect()
