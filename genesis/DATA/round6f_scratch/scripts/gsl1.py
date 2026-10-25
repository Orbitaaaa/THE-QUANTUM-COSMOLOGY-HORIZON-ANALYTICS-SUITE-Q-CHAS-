import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate as C, gsl_solution as GS, preparation as PR
from concurrent.futures import ProcessPoolExecutor
def summ(tag, o):
    dh = GS.daughter_horizons(o, n_points=12)
    eh = GS.parent_event_horizon(o)
    return dict(tag=tag, grew=bool(o.get("grew")), parent_deficit=eh.get("deficit"), max_daughter_deficit=dh.get("max_deficit"), S_N=dh.get("S_N"),
                rows=[(r["u"], r["v"], round(r["deficit_ingoing_piece"], 3), round(r["deficit_outgoing_piece"], 3)) for r in dh.get("rows", [])][:6])
def one(job):
    tag = job[0]
    if tag == "fill":
        o = C.nariai_condensate(job[1], frame="corner", h=0.04)
    elif tag == "tilt":
        o = C.nariai_condensate(0.0, frame="corner", h=0.04, lam_f=job[1])
    elif tag == "route_tilt":
        o = C.route_condensate(lam_f=job[1])
    else:
        res, cl, g, mt = PR.run(7.258540773391722, N_amp=10.0, u_p=1.0, w_p=0.1, v_shell=4.5, w_shell=0.25, Lv=600.0, Lu=30.0)
        o = dict(res=res, cls=cl["cls"], r_N=mt["r_N"], grew=g["grew"])
    return summ(f"{tag}:{job[1]}", o)
jobs = [("fill", -6e-6), ("fill", -2e-5), ("tilt", 1.3e-7), ("tilt", 4e-7), ("route_tilt", 1e-6), ("route_tilt", 3e-6), ("route_push", 10.0)]
with ProcessPoolExecutor(3) as ex:
    for r in ex.map(one, jobs):
        print(json.dumps(r), flush=True)
