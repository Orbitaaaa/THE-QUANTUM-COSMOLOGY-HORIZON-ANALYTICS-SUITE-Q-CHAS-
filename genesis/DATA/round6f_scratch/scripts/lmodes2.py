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
res = o["res"]; R, S, h = res["r"], res["s"], res["h"]; rN = o["r_N"]
v = np.arange(R.shape[1]) * h
pa = np.exp(-((v - 2.0) / 2.0) ** 2) * R[0]; pb = np.zeros(R.shape[0]); pb[0] = pa[0]
# trapped flags from the background
ru = np.gradient(R, h, axis=0); rv = np.gradient(R, h, axis=1)
out = dict(tag=tag, r_N=rN, grew=bool(o["grew"]))
for kappa in (0.0, 1.0):
    for ell in (0, 1, 2, 3, 5, 10, 20):
        P = E._lmode(np.ascontiguousarray(R), np.ascontiguousarray(S), h, float(ell), float(kappa), pa, pb)
        F = np.abs(P / R)
        good = np.isfinite(F)
        big = good & (R > 0.5 * rN)
        untr = good & (ru < 0) & (rv > 0)
        i, j = np.unravel_index(np.nanargmax(np.where(good, F, -1)), F.shape)
        out[f"k{int(kappa)}_l{ell}"] = dict(all=float(F[good].max()), r_gt_half_rN=float(F[big].max()) if big.any() else None,
                                           untrapped=float(F[untr].max()) if untr.any() else None,
                                           r_at_max=float(R[i, j]), trapped_at_max=bool(not (ru[i, j] < 0 and rv[i, j] > 0)))
print(json.dumps(out), flush=True)
