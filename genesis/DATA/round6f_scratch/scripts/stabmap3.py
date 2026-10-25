import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate_eft as E
K = 7.5e-4
grid = np.geomspace(1e-7, 10.0, 81)
res = []
for rN in (1e6, 1e20, 1e40, 1e60):
    for al in (0.01, 0.1, 1.0, 10.0, 100.0):
        for mode in ("dilute", "self"):
            row = []
            for m in grid:
                mu = (m * rN**2) ** (1/3)
                eps0 = -K * rN**2 / mu**4
                J = 8 * np.pi * al * (mu / rN) ** 2
                g = E.growth_exponent(mu, al, eps0, J=J, sigma=abs(eps0) if mode == "self" else 0.0, t_end=9.6 if mode == "dilute" else 20.0)
                row.append(dict(M3rN=float(m), G=g["G"], eps0=float(eps0), valid=g["valid"] and abs(eps0) < 1))
            ok = [r["M3rN"] for r in row if r["G"] < 1 and r["valid"]]
            res.append(dict(r_N=rN, alpha=al, mode=mode, lo=min(ok) if ok else None, hi=max(ok) if ok else None, row=row))
            print(f"r_N={rN:.0e} alpha={al:6.2f} {mode:6s}: G<1 for M^3 r_N in " + (f"[{min(ok):.2e}, {max(ok):.2e}]  width {max(ok)/min(ok):.0f}" if ok else "none"), flush=True)
json.dump(res, open("/tmp/claude-0/-home-user-THE-QUANTUM-COSMOLOGY-HORIZON-ANALYTICS-SUITE-Q-CHAS-/5d1f4f78-e2f6-5e44-a419-ab5b29ff8104/scratchpad/r7/stabmap3.json", "w"))
