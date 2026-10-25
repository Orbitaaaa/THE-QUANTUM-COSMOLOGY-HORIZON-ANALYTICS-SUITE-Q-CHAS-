import numpy as np, warnings, json, sys
warnings.simplefilter("ignore")
from babyuniverse import condensate as C
from concurrent.futures import ProcessPoolExecutor
def run(a):
    eps, dt, lam, wf = a
    o = C.route_condensate(lam_f=lam, eps=eps, clip=1, dphi0=1.5e-9, Lv=600.0, r_inf=1e40, delta_target=dt, wall_frac=wf)
    hr = C.hiding_report(o)
    return dict(eps=eps, delta=dt, wall_frac=wf, lam_f=lam, normal=o["normal_on_data"], grew=bool(o["grew"]), n_big=int(o["n_big"]),
                u_EH=hr["u_EH"], hidden=hr["hidden"], edge_rows=hr["daughter_rows_reaching_outer_edge"], first_dau=hr["first_daughter_u"])
jobs = [(2.0, 6e-3, 1e-6, 0.6), (2.0, 1e-2, 1e-6, 0.6), (2.0, 4e-3, 1e-6, 0.5), (2.0, 2e-2, 1e-6, 0.6)]
with ProcessPoolExecutor(4) as ex:
    for r in ex.map(run, jobs):
        print(json.dumps(r), flush=True)
