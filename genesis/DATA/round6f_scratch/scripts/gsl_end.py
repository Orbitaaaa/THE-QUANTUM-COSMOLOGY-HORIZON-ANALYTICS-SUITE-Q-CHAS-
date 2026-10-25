import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate as C, gsl_solution as GS
from concurrent.futures import ProcessPoolExecutor
def case(tag):
    if tag == "plain_black_hole":                       # the route data with no push: a black hole, no daughter
        o = C.route_condensate(lam_f=0.0)
    elif tag == "transient_daughter":                   # eps = 4, lam_f = 2e-7, seeded exit: inflates, exits, recollapses
        o = C.nariai_condensate(0.0, frame="corner", h=0.08, lam_f=2e-7, Lu=520.0, Lv=520.0, clip=1, dphi0=1.5e-9)
    else:                                                # Schwarzschild-de Sitter near-Nariai band with no push: collapses
        o = C.nariai_condensate(0.0, frame="corner", h=0.08)
    R, P, Cc = o["res"]["r"], o["res"]["phi"], o["cls"]
    # observers whose rows end next to r -> 0 (singularity / crunch), anywhere in the domain
    region = np.isfinite(R) & (R < 0.1 * o["r_N"])
    d = GS.ending_observer_deficits(o, region, n_points=10)
    return tag, dict(points=d.get("points"), max_deficit=d.get("max_deficit"), S_N=float(np.pi * o["r_N"]**2),
                     rows=[(round(x["u"], 1), round(x["v"], 1), round(x["deficit_ingoing"], 1), round(x["deficit_outgoing"], 1)) for x in d.get("rows", [])][:5])
with ProcessPoolExecutor(3) as ex:
    for tag, r in ex.map(case, ["plain_black_hole", "near_nariai_collapse", "transient_daughter"]):
        print(tag, json.dumps(r), flush=True)
