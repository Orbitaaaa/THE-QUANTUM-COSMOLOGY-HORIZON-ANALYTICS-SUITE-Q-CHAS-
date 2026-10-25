import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate as C, gsl_solution as GS
from concurrent.futures import ProcessPoolExecutor
def one(lf):
    o = C.route_condensate(lam_f=lf)
    eh = GS.parent_event_horizon(dict(res=o["res"], cls=o["cls"], r_N=o["r_N"]))
    return dict(lam_f=lf, grew=o["grew"], n_big=o["n_big"], r_max=o["r_max"], normal=o["normal_on_data"], rows=o["r_at_end_of_rows_u2_u3p2"],
                parent_EH_deficit=eh.get("deficit"), u_EH=eh.get("u_EH"))
with ProcessPoolExecutor(3) as ex:
    for r in ex.map(one, [0.0, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5]):
        print(json.dumps(r), flush=True)
