"""Recompute the route and transient entries of DATA/round6e.json with the hiding report and the in-band exit metric (same code as the
updated scripts/round6e.py)."""
import json, warnings, gc, sys
import numpy as np
warnings.simplefilter("ignore")
sys.path.insert(0, ".")
from babyuniverse import condensate as C, gsl_solution as GS
from concurrent.futures import ProcessPoolExecutor
RHO_F = 3 / 900 / (8 * np.pi)
src = open("round6e.py").read()
ns = {}
exec(src[src.index("def transient(seed):"):src.index("out = {}")], dict(np=np, C=C, GS=GS, RHO_F=RHO_F), ns)
def tr(seed):
    return ns["transient"](seed)


def route(lf):
    o = C.route_condensate(lam_f=lf)
    return dict(lam_f=lf, grew=bool(o["grew"]), n_big=o["n_big"], normal_on_data=o["normal_on_data"], **C.hiding_report(o))
p = "../../DATA/round6e.json"
d = json.load(open(p))
with ProcessPoolExecutor(2) as ex:
    d["problem1"]["route_runs"] = list(ex.map(route, [0.0, 3e-7, 1e-6, 3e-6]))
    d["problem2"]["transient"] = list(ex.map(tr, [1.5e-9, 0.0]))
json.dump(d, open(p, "w"), indent=1, default=float)
for r in d["problem1"]["route_runs"]: print(r)
for t in d["problem2"]["transient"]: print({k: t[k] for k in ("seed", "grew", "T_exit_pred", "T_at_first_exit", "parent_deficit")}, t["hiding"], {k: (v["r_max"], v["max_phi_over_v"]) for k, v in t["windows"].items()})
