import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse import condensate_eft as E
from concurrent.futures import ProcessPoolExecutor
K, rN, al = 7.5e-4, 1e6, 1.0
def run(a):
    M3rN, mode, seed_amp = a
    mu = (M3rN * rN**2) ** (1/3); eps0 = -K * rN**2 / mu**4; J = 8 * np.pi * al * (mu / rN) ** 2
    qb = np.sqrt(2 * (abs(eps0) + J) * mu**2 / al); L = 30 * 2 * np.pi / qb
    nl = E.evolve_1p1(mu, al, eps0, J=J, sigma=abs(eps0) if mode == "self" else 0.0, L=L, N=512, t_end=9.6, seed_amp=seed_amp, k_seed=2 * qb)
    # size of the perturbation of X on the data: |2 pi_t - pi_x^2| ~ pi_x^2 for a static seed
    return dict(M3rN=M3rN, mode=mode, seed=seed_amp, dX_over_eps=float((seed_amp * qb) ** 2 / abs(eps0)), nonlinear=nl["amplification"], broke=bool(nl["broke"]), t_break=float(nl["t"][-1]))
jobs = [(m, mode, s) for m in (3.8e-4, 3.2e-2) for mode in ("dilute", "self") for s in (1e-6, 1e-5, 1e-4)]
with ProcessPoolExecutor(4) as ex:
    for r in ex.map(run, jobs):
        print(json.dumps(r), flush=True)
