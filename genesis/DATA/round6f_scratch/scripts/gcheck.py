import numpy as np, warnings
warnings.simplefilter("ignore")
from babyuniverse import condensate_eft as E
for (mu, al, e0, J, sg, t_end) in [(30.0, 1.0, -0.05, 0.0, 0.0, 9.6), (30.0, 1.0, -0.2, 0.0, 0.0, 9.6), (30.0, 1.0, -0.05, 0.0, 0.05, 20.0),
                                   (60.0, 1.0, -0.1, 0.0, 0.1, 20.0), (30.0, 1.0, -1e-3, 0.05, 0.0, 9.6), (20.0, 0.3, -0.3, 0.0, 0.0, 9.6)]:
    g = E.growth_exponent(mu, al, e0, J=J, sigma=sg, t_end=t_end)
    qb = np.sqrt(2 * (max(abs(e0), sg) + J) * mu**2 / al)
    L = 30 * 2 * np.pi / qb
    a = E.linear_rms_1p1(mu, al, e0, J=J, sigma=sg, L=L, N=256, t_end=t_end)["amplification"]
    print(f"mu={mu} eps0={e0} J={J} sigma={sg}:  G={g['G']:.3f} e^G={np.exp(g['G']):.3f}   ODE rms amplification={a:.3f}")
