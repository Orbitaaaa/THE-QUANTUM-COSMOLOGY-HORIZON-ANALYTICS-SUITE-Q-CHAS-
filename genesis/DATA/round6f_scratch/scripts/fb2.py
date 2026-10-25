import numpy as np, warnings, json
warnings.simplefilter("ignore")
from babyuniverse.fireball import Fireball
from concurrent.futures import ProcessPoolExecutor
def spinodal(fb):
    """T_1: lowest T at which phi = +v side has no local minimum of F (field must leave the true vacuum)."""
    xs = np.linspace(-0.37, 1.3, 4000) * fb.v
    for r in np.linspace(1.0, 3.0, 401):
        T = r * fb.Tc
        Fp = fb.cT * T * T * (xs + fb.v) + fb.dU(xs)
        if not np.any((Fp[:-1] < 0) & (Fp[1:] >= 0)):
            return r
    return None
def one(a):
    alpha, Ti, sb, frac = a
    fb = Fireball(alpha_pt=alpha, Gamma_over_m=0.3)
    T_i = Ti * fb.Tc
    rho_i = fb.thermo(fb.v, 0.0, T_i)[0]
    tanb = sb / np.sqrt(1 - sb * sb)
    ad_i = frac / tanb
    a_i = np.sqrt((1 + ad_i**2) / (8 * np.pi / 3 * rho_i))
    try:
        sm = fb.summary(fb.run(T_i, a_i, t_max=400.0), sin_chi_b=sb)
    except Exception as e:
        return dict(alpha=alpha, Ti=Ti, err=str(e))
    return dict(alpha=alpha, Ti=Ti, sin_chi_b=sb, a_i=round(a_i, 2), T1_over_Tc=spinodal(fb), **{k: (round(v, 3) if isinstance(v, float) else v) for k, v in sm.items()
                 if k in ("turnaround", "x_end", "T_end_over_Tc", "R_ball_turn_over_lH", "rad_over_rho_f", "normal_throughout", "max_adot_tanchib", "a_max_over_lH")})
jobs = [(al, Ti, sb, 0.95) for al in (0.5, 1.0, 2.0) for Ti in (1.4, 1.6, 1.8) for sb in (0.6, 0.75)]
with ProcessPoolExecutor(4) as ex:
    for r in ex.map(one, jobs):
        print(json.dumps(r), flush=True)
