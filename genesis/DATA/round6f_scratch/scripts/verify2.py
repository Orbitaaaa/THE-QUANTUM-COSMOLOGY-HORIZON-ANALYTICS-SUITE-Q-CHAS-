import numpy as np, warnings
warnings.simplefilter("ignore")
from babyuniverse import throat_evolution as TE, condensate as C
M, r0 = 1.0, 6.0
f0 = 1 - 2 * M / r0
lump = lambda x: 2e-3 * np.exp(-((np.asarray(x, float) - 2.0) / 0.8) ** 2)
for frame in ("corner", "static"):
    for h in (0.1, 0.05, 0.025):
        Nv = Nu = int(6 / h)
        A2, B2, info2 = C.condensate_data(0.0, 0.0, r0, -np.sqrt(f0) / 2, np.sqrt(f0) / 2, h, Nv, Nu, lump, lump, frame=frame)
        r2 = C.evolve(A2, B2, 0.0, 0.0, TE.Potential(), h)
        cr = C.constraint_residual(r2)
        R = r2["r"]; T = C.stress(r2)
        # Misner-Sharp mass at the far corner vs. M + dust energy
        d = TE.derivatives(r2); m = TE.misner_sharp(dict(r=R, s=r2["s"], Q=0.0, h=h), d)
        print(frame, h, "Evv", f"{cr['Evv']:.2e}", "Euu", f"{cr['Euu']:.2e}", "rel", f"{cr['Evv_rel']:.1e} {cr['Euu_rel']:.1e}",
              "eta range", f"{np.nanmin(r2['eta']):.2f} {np.nanmax(r2['eta']):.2f}", "m(far)", f"{m[-2, -2]:.8f}", "m(0,far)", f"{m[1, -2]:.8f}")
