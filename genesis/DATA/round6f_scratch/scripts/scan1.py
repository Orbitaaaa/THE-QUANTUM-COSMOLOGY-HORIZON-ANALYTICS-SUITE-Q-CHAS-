import numpy as np, warnings, time, sys
warnings.simplefilter("ignore")
from babyuniverse import condensate as C
h = float(sys.argv[1]) if len(sys.argv) > 1 else 0.08
for n0 in (0.0, -2e-6, -5e-6, -1e-5, -2e-5, -4e-5):
    t = time.time()
    o = C.nariai_condensate(n0, h=h, frame="corner")
    print(f"n0={n0:.1e} grew={o['grew']} n_big={o['n_big']} r_max={o['r_max']} hidden_ok={o['true_vacuum_untrapped_large_r']==0} "
          f"W_u={o['W_u']:.4f} W_v={o['W_v']:.4f} t={time.time()-t:.0f}s", flush=True)
