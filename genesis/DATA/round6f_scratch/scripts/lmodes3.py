import numpy as np, warnings, json, sys
warnings.simplefilter("ignore")
from babyuniverse import condensate as C, condensate_eft as E
tag = sys.argv[1]
if tag == "band_self_tipping":
    o = C.nariai_condensate(-1e-5, frame="corner", h=0.08)
elif tag == "control_no_tip":
    o = C.nariai_condensate(-4e-6, frame="corner", h=0.08)
elif tag == "route_self_made_push":
    o = C.route_condensate(lam_f=1e-6)
elif tag == "route_control":
    o = C.route_condensate(lam_f=0.0)
else:
    o = C.nariai_condensate(0.0, frame="corner", h=0.08, lam_f=2e-7, Lu=400.0, Lv=400.0, clip=1, dphi0=1.5e-9)
out = dict(tag=tag, grew=bool(o["grew"]), r_N=o["r_N"])
for kappa in (0.0, 1.0):
    ells = (0, 1, 2, 3, 5, 10, 20) if kappa == 0 else (2, 3, 5, 10, 20)
    out[f"k{int(kappa)}_all"] = E.lmode_amplification(o["res"], ells=ells, kappa=kappa, pot=o["pot"])
    out[f"k{int(kappa)}_regular"] = E.lmode_amplification(o["res"], ells=ells, kappa=kappa, pot=o["pot"], r_cut=0.5 * o["r_N"])
print(json.dumps(out), flush=True)
