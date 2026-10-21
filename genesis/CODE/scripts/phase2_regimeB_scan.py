"""Phase II: regime-B Euclidean angle theta_D(gamma, E) on a grid.
Checks (i) |theta_D| > pi everywhere, (ii) monotonic increase of |theta_D| with E
from the critical end (closed form, proved > pi) to the A/B boundary."""
import warnings

import numpy as np

from _common import dump
from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse import tunnelling as T

warnings.simplefilter("ignore")
rows, mono_ok, min_ratio = [], True, np.inf
for gam in np.linspace(-0.95, 0.95, 39):
    gmin, g1 = ThinWallModel.g_min(gam), 2 - 2 * gam
    fr = np.array([1e-4, 0.02, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95, 0.999])
    th = []
    for f in fr:
        m = model_from_dimensionless(gam, gmin + f * (g1 - gmin))
        Ra, Rb = T.barrier(m)
        D = T._quad_s(m, lambda r: m.A_out(r) / m.f_out(r), Ra, Rb)   # half-bounce Killing time
        th.append(abs(D / (4 * m.G * m.M)) / np.pi)
    th = np.array(th)
    mono = bool(np.all(np.diff(th) > -1e-9))
    mono_ok &= mono
    min_ratio = min(min_ratio, th.min())
    rows.append(dict(gamma=gam, fractions=fr, theta_over_pi=th, monotone=mono,
                     critical_closed_form=T.theta_D_critical(gam)))
print("regime B: min |theta_D|/pi =", min_ratio, " monotone in E for all gamma:", mono_ok)
dump("phase2_regimeB_scan.json", dict(rows=rows, min_ratio=min_ratio, monotone_all=mono_ok))
