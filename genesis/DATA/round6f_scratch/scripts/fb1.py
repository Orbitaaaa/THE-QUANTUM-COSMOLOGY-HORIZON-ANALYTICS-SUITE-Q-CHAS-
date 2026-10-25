import numpy as np, warnings, json, time
warnings.simplefilter("ignore")
from babyuniverse.fireball import Fireball
for alpha, Gm in ((1.0, 0.3), (3.0, 0.3), (10.0, 0.3)):
    fb = Fireball(alpha_pt=alpha, Gamma_over_m=Gm)
    for Ti in (1.3, 2.0):
        T_i = Ti * fb.Tc
        rho_i = fb.thermo(fb.v, 0.0, T_i)[0]
        for f in (0.5, 0.9):                              # (8 pi/3) rho_i a_i^2 = 1 + adot_i^2 ; choose adot_i = f * 4/3 (sin chi_b = 0.6 edge stays normal at start)
            ad_i = f * 4 / 3
            a_i = np.sqrt((1 + ad_i**2) / (8 * np.pi / 3 * rho_i))
            t = time.time()
            try:
                sol = fb.run(T_i, a_i, t_max=600.0)
                sm = fb.summary(sol, sin_chi_b=0.6)
            except Exception as e:
                sm = dict(err=str(e))
            print(json.dumps(dict(alpha_pt=alpha, Ti_over_Tc=Ti, a_i=round(a_i, 2), adot_i=round(ad_i, 3), t=round(time.time() - t, 1),
                                  **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in sm.items()})), flush=True)
