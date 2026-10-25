import numpy as np, warnings
warnings.simplefilter("ignore")
from babyuniverse import condensate as C
for lf, ls in ((4e-7, 4e-7), (4e-7, 0.0)):
    for h in (0.08, 0.04, 0.02):
        o = C.nariai_condensate(0.0, frame="corner", h=h, lam_f=lf, lam_s=ls, Lu=40.0, Lv=40.0)
        rN = o["r_N"]
        cr = C.constraint_residual(o["res"], r_lo=0.5 * rN, r_hi=3 * rN, s_max=5.0)
        print(f"lam_f={lf} lam_s={ls} h={h} Evv={cr['Evv']:.3e} Euu={cr['Euu']:.3e} rel {cr['Evv_rel']:.2e} {cr['Euu_rel']:.2e}", flush=True)
