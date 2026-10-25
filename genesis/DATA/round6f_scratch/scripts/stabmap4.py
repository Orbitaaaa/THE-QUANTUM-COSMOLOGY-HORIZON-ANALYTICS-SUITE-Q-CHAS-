import numpy as np, warnings, json, pickle
warnings.simplefilter("ignore")
from babyuniverse import condensate as C, condensate_eft as E
from concurrent.futures import ProcessPoolExecutor
def history(tag):
    n0, lf = (-6e-6, 0.0) if tag == "initial_charge" else (0.0, 4e-7)
    o = C.nariai_condensate(n0, frame="corner", h=0.04, lam_f=lf)
    res, rN = o["res"], o["r_N"]
    n = C.stress(res)["n"]
    k = np.arange(min(n.shape))
    T, nd, rd = res["T"][k, k], n[k, k], res["r"][k, k]
    ok = np.isfinite(T) & np.isfinite(nd)
    return tag, T[ok] / rN, nd[ok] * rN**2, rd[ok] / rN
with ProcessPoolExecutor(2) as ex:
    H = {t: (x, nu, rr) for t, x, nu, rr in ex.map(history, ["initial_charge", "self_generated"])}
pickle.dump(H, open("/tmp/claude-0/-home-user-THE-QUANTUM-COSMOLOGY-HORIZON-ANALYTICS-SUITE-Q-CHAS-/5d1f4f78-e2f6-5e44-a419-ab5b29ff8104/scratchpad/r7/histories.pkl", "wb"))
for t, (x, nu, rr) in H.items():
    print(t, "x_max", round(float(x.max()), 2), "int|nu| dx", float(np.trapezoid(np.abs(nu), x)), "nu range", float(nu.min()), float(nu.max()), flush=True)
grid = np.geomspace(1e-5, 1.0, 61)
def scan(a):
    tag, rN, al = a
    x, nu, rr = H[tag]
    G = [E.growth_exponent_history(m, rN, al, x, nu, rr, decay_after=(tag == "initial_charge"))["G"] for m in grid]
    ok = [float(m) for m, g in zip(grid, G) if g < 1]
    return dict(tag=tag, r_N=rN, alpha=al, lo=min(ok) if ok else None, hi=max(ok) if ok else None, G=[float(g) for g in G])
jobs = [(t, rN, al) for t in H for rN in (1e6, 1e20, 1e40) for al in (0.1, 1.0, 10.0)]
with ProcessPoolExecutor(4) as ex:
    out = list(ex.map(scan, jobs))
for r in out:
    print(f"{r['tag']:15s} r_N={r['r_N']:.0e} alpha={r['alpha']:5.1f}: G<1 for M^3 r_N in " + (f"[{r['lo']:.2e}, {r['hi']:.2e}] width {r['hi']/r['lo']:.1f}" if r['lo'] else "none"), flush=True)
json.dump(dict(grid=list(map(float, grid)), scans=out), open("/tmp/claude-0/-home-user-THE-QUANTUM-COSMOLOGY-HORIZON-ANALYTICS-SUITE-Q-CHAS-/5d1f4f78-e2f6-5e44-a419-ab5b29ff8104/scratchpad/r7/stabmap4.json", "w"))
