import numpy as np, warnings, json, time
warnings.simplefilter("ignore")
from babyuniverse import condensate_eft as E
from concurrent.futures import ProcessPoolExecutor
K = 7.5e-4            # |eps0| M^4 r_N^2 needed to tip (Round 6e: n_th r_N^2 ~ 1.5e-3, n = 2 eps M^4)
def point(a):
    rN, M3rN, alpha, mode = a
    mu = (M3rN * rN**2) ** (1/3)
    eps0 = -K * rN**2 / mu**4
    J = 8 * np.pi * alpha * (mu / rN) ** 2
    sigma = abs(eps0) if mode == "self" else 0.0
    t_end = 9.6 if mode == "dilute" else 20.0
    lam_unst = np.sqrt(2 * mu**2 * (abs(eps0) + J) / alpha)       # largest unstable comoving wavenumber at t = 0
    L = 30.0 * 2 * np.pi / lam_unst
    li = E.linear_rms_1p1(mu, alpha, eps0, J=J, sigma=sigma, L=L, N=512, t_end=t_end, seed_amp=1e-8)
    return dict(r_N=rN, M3rN=M3rN, alpha=alpha, mode=mode, mu=mu, eps0=eps0, J=J, amplification=li["amplification"])
grid = np.geomspace(1e-6, 1.0, 31)
jobs = [(rN, m, al, mode) for rN in (1e6, 1e20, 1e40) for al in (0.1, 1.0, 10.0) for mode in ("dilute", "self") for m in grid]
t0 = time.time()
with ProcessPoolExecutor(4) as ex:
    res = list(ex.map(point, jobs, chunksize=4))
json.dump(res, open("/tmp/claude-0/-home-user-THE-QUANTUM-COSMOLOGY-HORIZON-ANALYTICS-SUITE-Q-CHAS-/5d1f4f78-e2f6-5e44-a419-ab5b29ff8104/scratchpad/r7/stabmap.json", "w"))
print(f"{len(res)} points in {time.time()-t0:.0f}s")
for rN in (1e6, 1e20, 1e40):
    for al in (0.1, 1.0, 10.0):
        for mode in ("dilute", "self"):
            rows = [r for r in res if r["r_N"] == rN and r["alpha"] == al and r["mode"] == mode]
            ok = [r["M3rN"] for r in rows if r["amplification"] < np.e]
            s = f"[{min(ok):.2e}, {max(ok):.2e}]" if ok else "none"
            print(f"r_N={rN:.0e} alpha={al:5.1f} {mode:6s}: amp<e for M^3 r_N in {s} | min amp {min(r['amplification'] for r in rows):.3f}")
