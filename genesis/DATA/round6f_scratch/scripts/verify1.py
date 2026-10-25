import numpy as np, warnings, time
warnings.simplefilter("ignore")
from babyuniverse import throat_evolution as TE, condensate as C
# 1. zero condensate reproduces the diamond solver (Schwarzschild data, r0 = 6, M = 1)
M, r0 = 1.0, 6.0
f0 = 1 - 2 * M / r0
for h in (0.1, 0.05, 0.025):
    Nv = Nu = int(10 / h)
    z = lambda x: 0 * np.asarray(x, float)
    A, B, info = C.condensate_data(0.0, 0.0, r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, h, Nv, Nu, z, z)
    res = C.evolve(A, B, 0.0, 0.0, TE.Potential(), h)
    a, b, _ = TE.characteristic_data(0.0, 0.0, r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, h, Nv, Nu)
    ref = TE.evolve(a, b, 0.0, 0.0, TE.Potential(), h)
    # 2. positive dust lump (proper density 1e-3) centred on the corner, static frame
    lump = lambda x: 2e-3 * np.exp(-((np.asarray(x, float) - 2.0) / 0.8) ** 2)
    A2, B2, info2 = C.condensate_data(0.0, 0.0, r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, h, Nv, Nu, lump, lump)
    t = time.time()
    r2 = C.evolve(A2, B2, 0.0, 0.0, TE.Potential(), h)
    cr = C.constraint_residual(r2)
    print(h, "D=0 vs TE max|dr| =", np.nanmax(np.abs(res["r"] - ref["r"])), " dust: Evv", cr["Evv"], "Euu", cr["Euu"],
          "max|eta|", np.nanmax(np.abs(r2["eta"])), "t", round(time.time() - t, 2))
