"""Reduced Hamiltonian checks and inequivalent quantizations."""
import warnings

import numpy as np

from _common import dump
from babyuniverse.model import ThinWallModel, model_from_dimensionless
from babyuniverse import quantum as Q, tunnelling as T

warnings.simplefilter("ignore")
checks = []
for gam, Ef in ((-0.5, 1.5), (0.6, 1.5), (0.0, 3.0)):
    m = model_from_dimensionless(gam, Ef * (2 - 2 * gam))
    R1, R2 = m.turning_points()
    cands = list(np.linspace(0.3 * R1, R1 * 0.999, 6)) + list(np.linspace(R2 * 1.001, 3 * R2, 8))
    for R in cands:
        if m.f_in(R) > 0 and m.f_out(R) > 0 and m.V(R) < 0:
            Rd = np.sqrt(-m.V(R))
            if abs(Rd) < abs(m.A_in(R)) and abs(Rd) < abs(m.A_out(R)):
                for sg in (1, -1):
                    c = Q.hamilton_equation_check(m, R, sign=sg)
                    c.update(gamma=gam, rel_err=abs(c["dHdP"] / c["dRdt"] - 1))
                    checks.append(c)
print("Hamilton equation: max rel err", max(c["rel_err"] for c in checks), "over", len(checks), "points;",
      "sectors:", sorted(set(c["sector"] for c in checks)))
comp = []
for gam in (-0.6, 0.0, 0.6):
    for Ef in (1.2, 3.0, 10.0):
        m = model_from_dimensionless(gam, Ef * (2 - 2 * gam))
        comp.append(dict(gamma=gam, E=Ef * (2 - 2 * gam), M=m.M, B_canonical=T.B_canonical(m),
                         B_naive_membrane=Q.B_naive(m, "membrane"), B_naive_unit=Q.B_naive(m, "unit"),
                         eikonal_residual=Q.eikonal_check(m)))
for c in comp:
    print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in c.items()})
dump("quantum_checks.json", dict(hamilton=checks, quantization_comparison=comp))
