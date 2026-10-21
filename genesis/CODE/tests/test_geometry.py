"""Independent check of the explicit trajectories: finite differences of the
integrated null coordinates reproduce the normalisation u.u = -1 and the
unsquared junction value n^r = A(R)."""
import warnings

import numpy as np
import pytest

from babyuniverse.model import model_from_dimensionless
from babyuniverse.branches import classify
from babyuniverse.geometry import branch_trajectories

warnings.simplefilter("ignore")


@pytest.mark.parametrize("gam,E", [(-0.5, 4.5), (0.6, 0.791)])
def test_kruskal_normalisation_and_junction(gam, E):
    m = model_from_dimensionless(gam, E)
    GM = m.G * m.M
    for b in classify(m):
        tr = branch_trajectories(m, b)["out"]
        t, R, V, U = tr.tau, tr.R, tr.p, tr.q
        sel = (R > 0.05 * GM) & (R < 20 * GM) & (np.abs(U) < 50) & (np.abs(V) < 50)
        Vd_all, Ud_all = np.gradient(V, t), np.gradient(U, t)
        # skip the joint at tau = 0 and the end points (one-sided differences)
        sel &= np.abs(t) > 3 * np.max(np.diff(t))
        idx = np.where(sel)[0][5:-5]
        Vd, Ud = Vd_all[idx], Ud_all[idx]
        r = R[idx]
        k = 4 * GM**2 / r * np.exp(-r / (2 * GM))
        F = 8 * GM * k
        assert np.median(np.abs(F * Ud * Vd - 1)) < 1e-3
        nr = (-U[idx] * k) * Vd - (-V[idx] * k) * Ud
        assert np.median(np.abs(nr - m.A_out(r))) < 1e-3


def test_desitter_interior_normalisation():
    m = model_from_dimensionless(0.6, 1.2)
    for b in classify(m):
        tr = branch_trajectories(m, b)["in"]
        t, p, q, R = tr.tau, tr.p, tr.q, tr.R
        T = 0.5 * (p + q)
        good = (np.cos(T) > 0.05) & (np.abs(t) > 3 * np.max(np.diff(t)))
        idx = np.where(good)[0][5:-5]
        pd, qd = np.gradient(p, t)[idx], np.gradient(q, t)[idx]
        Om2 = 1 / (m.chi**2 * np.cos(T[idx]) ** 2)
        assert np.median(np.abs(Om2 * pd * qd - 1)) < 1e-3
        rp = np.cos(q[idx]) / (2 * m.chi * np.cos(T[idx]) ** 2)
        rq = -np.cos(p[idx]) / (2 * m.chi * np.cos(T[idx]) ** 2)
        assert np.median(np.abs(rp * pd - rq * qd - m.A_in(R[idx]))) < 1e-3
