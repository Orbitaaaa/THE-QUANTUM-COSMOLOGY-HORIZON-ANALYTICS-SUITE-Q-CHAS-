import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate as C, condensate_eft as E
from concurrent.futures import ProcessPoolExecutor
def case(tag):
    if tag == "band_self_tipping":
        o = C.nariai_condensate(-1e-5, frame="corner", h=0.08)
    elif tag == "route_self_made_push":
        o = C.route_condensate(lam_f=1e-6)
    else:
        o = C.nariai_condensate(0.0, frame="corner", h=0.08, lam_f=2e-7, Lu=400.0, Lv=400.0, clip=1, dphi0=1.5e-9)
    out = {}
    for kappa in (0.0, 1.0):
        out[f"kappa={kappa}"] = E.lmode_amplification(o["res"], ells=(1, 2, 3, 5, 10, 20), kappa=kappa)
    return tag, out
with ProcessPoolExecutor(2) as ex:
    for tag, out in ex.map(case, ["band_self_tipping", "route_self_made_push", "transient_daughter"]):
        print(tag, json.dumps({k: {l: round(v, 3) for l, v in d.items()} for k, d in out.items()}), flush=True)
