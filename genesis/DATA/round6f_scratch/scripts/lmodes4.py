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
rN = o["r_N"]
out = dict(tag=tag, grew=bool(o["grew"]))
out["k0_regular"] = E.lmode_amplification(o["res"], ells=(0, 1, 2, 3, 5, 10, 20), kappa=0.0, pot=o["pot"], r_cut=0.5 * rN)
out["k1_weyl_regular"] = E.lmode_amplification(o["res"], ells=(2, 3, 5, 10, 20), kappa=1.0, pot=o["pot"], r_cut=0.5 * rN, power=3)
out["k1_weyl_parent"] = E.lmode_amplification(o["res"], ells=(2, 3, 5, 10, 20), kappa=1.0, pot=o["pot"], r_cut=0.5 * rN, power=3)
R = o["res"]["r"]
print(json.dumps(out), flush=True)
