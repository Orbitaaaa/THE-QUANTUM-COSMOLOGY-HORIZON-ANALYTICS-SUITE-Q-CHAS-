import numpy as np, warnings
warnings.simplefilter("ignore")
from babyuniverse import throat_evolution as TE, condensate as C
pot = TE.Potential(0.0, 1.0, 0.0)
z = lambda x: 0 * np.asarray(x, float)
m1 = lambda x: -1 + 0 * np.asarray(x, float)
for h in (0.1, 0.05):
    N = int(4 / h)
    A, B, info = C.condensate_data(0.0, 0.0, 5.0, -0.5, 0.5, h, N, N, z, z, phi_v=m1, phi_u=m1, pot=pot, frame="corner")
    res = C.evolve(A, B, 0.0, 0.0, pot, h, lam_f=1e-6)
    T = res["T"]; n = C.stress(res)["n"]
    d = TE.derivatives(res)
    Tu = np.full_like(T, np.nan); Tv = np.full_like(T, np.nan)
    Tu[1:-1] = (T[2:] - T[:-2]) / (2 * h); Tv[:, 1:-1] = (T[:, 2:] - T[:, :-2]) / (2 * h)
    eik = -4 * np.exp(-2 * res["s"]) * Tu * Tv
    print(h, "clock eikonal |g(dT,dT)+1| max", np.nanmax(np.abs(eik + 1)), " n/(-lam T) at far corner", n[-1, -1] / (-1e-6 * T[-1, -1]),
          "T far", T[-1, -1], "(u+v)/2", 4.0)
