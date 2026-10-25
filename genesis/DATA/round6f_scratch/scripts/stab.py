import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate as C, ghost_condensate as GC
out = {}
for tag, n0, lf in (("initial_charge", -6e-6, 0.0), ("self_generated", 0.0, 4e-7)):
    o = C.nariai_condensate(n0, frame="corner", h=0.04, lam_f=lf)
    res = o["res"]; rN = o["r_N"]
    n = C.stress(res)["n"]
    k = np.arange(min(n.shape))
    T, nd, rd = res["T"][k, k], n[k, k], res["r"][k, k]
    ok = np.isfinite(T) & np.isfinite(nd)
    x, nu, rr = T[ok] / rN, nd[ok] * rN**2, rd[ok] / rN
    sel = x <= 9.6
    I = float(np.trapezoid(np.abs(nu[sel]), x[sel]))
    out[tag] = dict(x_max=float(x.max()), nu_min=float(nu.min()), I_3tau=I, r_end_over_rN=float(rr[sel][-1]),
                    window_lo=I / (2 * np.sqrt(2 * 1.0 * 2.0)), window_hi=GC.window(1e6)["M3rN_max"])
    np.savez(f"/tmp/claude-0/-home-user-THE-QUANTUM-COSMOLOGY-HORIZON-ANALYTICS-SUITE-Q-CHAS-/5d1f4f78-e2f6-5e44-a419-ab5b29ff8104/scratchpad/r7/hist_{tag}.npz", x=x, nu=nu, rr=rr)
    print(tag, json.dumps(out[tag]), flush=True)
