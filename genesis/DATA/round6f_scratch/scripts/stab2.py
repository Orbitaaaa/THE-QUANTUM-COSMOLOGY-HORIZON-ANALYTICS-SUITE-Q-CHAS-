import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import ghost_condensate as GC
from concurrent.futures import ProcessPoolExecutor
S = "/tmp/claude-0/-home-user-THE-QUANTUM-COSMOLOGY-HORIZON-ANALYTICS-SUITE-Q-CHAS-/5d1f4f78-e2f6-5e44-a419-ab5b29ff8104/scratchpad/r7/"
rN = 1e6
def hist(tag):
    d = np.load(S + f"hist_{tag}.npz"); x, nu, rr = d["x"], d["nu"], d["rr"]
    xm = x.max()
    if tag == "initial_charge":
        nu_f = lambda X: np.where(X <= xm, np.interp(X, x, nu), nu[-1] * np.exp(-(X - xm)))
        r_f = lambda X: rN * np.where(X <= xm, np.interp(X, x, rr), rr[-1])
    else:
        nu_f = lambda X: np.where(X <= xm, np.interp(X, x, nu), nu[-1])
        r_f = lambda X: rN * np.where(X <= xm, np.interp(X, x, rr), rr[-1] * np.exp((X - xm) / np.sqrt(3)))
    xs = np.linspace(0, 9.6, 2001)
    return nu_f, r_f, float(np.trapezoid(np.abs(nu_f(xs)), xs))
def one(a):
    tag, f = a
    nu_f, r_f, I = hist(tag)
    lo = I / (2 * np.sqrt(4.0))
    M3rN = lo * f
    M = (M3rN / rN) ** (1 / 3)
    eps_t = lambda t: nu_f(t / rN) / rN**2 / (2 * M**4)
    r_t = lambda t: float(r_f(t / rN))
    kmax = np.sqrt(abs(nu_f(np.linspace(0, 9.6, 200))).max() / rN**2 / (2 * M**4) / (1 / (2 * M * M))) * 3
    g = GC.mode_growth(M, 0.0, rN, t_end_over_rN=9.6, eps_t=eps_t, r_t=r_t, ks=np.geomspace(1e-2 / rN, kmax, 30),
                       ls=np.unique(np.concatenate([np.arange(6), np.geomspace(6, max(7, kmax * rN), 15).astype(int)])))
    return dict(tag=tag, I_3tau=I, M3rN=M3rN, factor=f, window_hi=GC.window(rN)["M3rN_max"], **g)
jobs = [(t, f) for t in ("initial_charge", "self_generated") for f in (0.5, 1.5, 5.0)]
with ProcessPoolExecutor(4) as ex:
    for r in ex.map(one, jobs):
        print(json.dumps(r), flush=True)
