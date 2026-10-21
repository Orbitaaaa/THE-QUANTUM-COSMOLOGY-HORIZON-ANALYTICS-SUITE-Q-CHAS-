"""Phase IV: double-null Einstein-scalar evolution."""
import warnings

import numpy as np
import pytest

from babyuniverse import doublenull as D, stability as S
from babyuniverse.thickwall import QuarticPotential
from babyuniverse.model import ThinWallModel

warnings.simplefilter("ignore")
RHO_F = 3 / (8 * np.pi)


def bubble_potential(w, kappa=1.0):
    q = QuarticPotential.from_width_tension(w, kappa / (4 * np.pi), RHO_F)
    return q, D.Potential.from_quartic(q)


def test_symbolic_field_equations():
    ratios = D.symbolic_check()
    for k, v in ratios.items():
        assert "Derivative" not in v and v not in ("0", "nan", "zoo"), (k, v)


def test_flat_space_scalar_is_exact():
    pot0 = D.Potential(lam=0.0, v=1.0, rho_f=0.0)
    prof = lambda r: 0.01 * np.exp(-((r - 0.5) / 0.1) ** 2)
    res = D.evolve(prof, pot0, 0.01, v_max=2.5, G=1e-12)
    U, V = np.meshgrid(res["u"], res["v"], indexing="ij")
    F = lambda x: (x / 2) * prof(x / 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        ex = (F(V) - F(U)) / ((V - U) / 2)
    assert np.nanmax(np.abs(res["phi"] - ex)[(V - U) > 0.05]) < 1e-12


def test_de_sitter_convergence():
    pot = D.Potential(lam=50.0, v=0.1, rho_f=RHO_F)
    errs = []
    for h in (0.02, 0.01, 0.005):
        k = int(round(0.04 / h))
        res = D.evolve(lambda r: -pot.v * np.ones_like(r), pot, h, v_max=1.6, ks=k, kv=k)
        U, V = np.meshgrid(res["u"], res["v"], indexing="ij")
        rex, _ = D.de_sitter_exact(1.0, U, V)
        errs.append(np.nanmax(np.abs(res["r"] - rex)))
    assert errs[0] / errs[1] > 2.4 and errs[1] / errs[2] > 2.7
    assert errs[-1] < 1e-4


def test_self_gravitating_pulse_second_order():
    pot0 = D.Potential(lam=0.0, v=1.0, rho_f=0.0)
    sols = {}
    for h in (0.004, 0.002, 0.001):
        k = int(round(0.016 / h))
        sols[h] = D.evolve(lambda r: 0.15 * np.exp(-((r - 0.5) / 0.1) ** 2), pot0, h, v_max=2.4, ks=k, kv=k)
    mask = sols[0.004]["r"] > 0.05
    for hh in sols:
        mask &= np.isfinite(sols[hh]["r"]) & (np.gradient(sols[hh]["r"], 0.016, axis=1) > 0.1)
    d1 = np.nanmax(np.abs(sols[0.004]["r"] - sols[0.002]["r"])[mask])
    d2 = np.nanmax(np.abs(sols[0.002]["r"] - sols[0.001]["r"])[mask])
    assert d1 / d2 > 3.3


def test_collapse_threshold():
    pot0 = D.Potential(lam=0.0, v=1.0, rho_f=0.0)
    lo = D.evolve(lambda r: 0.03 * np.exp(-((r - 0.5) / 0.1) ** 2), pot0, 0.004, v_max=3.0, ks=10, kv=10)
    hi = D.evolve(lambda r: 0.15 * np.exp(-((r - 0.5) / 0.1) ** 2), pot0, 0.004, v_max=3.0, ks=10, kv=10)
    assert np.nanmax(lo["n_trapped"]) == 0
    assert np.nanmax(hi["n_trapped"]) > 0
    M = np.nanmedian(hi["m_boundary"][1:20])
    assert abs(np.nanmax(hi["r_ah"]) / (2 * M) - 1) < 0.03


def test_laboratory_bubble_collapses_like_thin_wall():
    w = 0.02
    q, pot = bubble_potential(w)
    res = D.evolve(lambda r: q.v * np.tanh((r - 0.4) / w), pot, w / 16, v_max=2.0, ks=8, kv=8)
    tr = D.wall_trajectory(res)
    M = np.nanmedian(res["m_boundary"][1:50])
    R1 = ThinWallModel(1.0, 1.0, M).turning_points()[0]
    assert abs(np.nanmax(tr["r"]) / R1 - 1) < 0.015
    r_ah = res["r_ah"][np.isfinite(res["r_ah"])][0]
    assert abs(r_ah / (2 * M) - 1) < 0.03
    assert np.nanmax(res["n_antitrapped"]) == 0           # Farhi-Guth in field theory
    early = tr["r"] > 0.95 * np.nanmax(tr["r"])
    ratio = (tr["A_in"] - tr["A_out"]) / tr["r"]           # Israel: A_in - A_out = kappa R
    assert abs(np.nanmedian(ratio[early]) - 1) < 0.03


def test_bounce_branch_daughter_inflates_beyond_bridge():
    w = 0.02
    q, pot = bubble_potential(w)
    R2 = ThinWallModel(1.0, 1.0, 0.1545).turning_points()[-1]
    wv = 2 * w * np.sqrt(1 - R2**2)
    h = wv / 12
    r0, s0, p0 = D.initial_cone_alpha1(lambda v: q.v * np.tanh((v - 2 * R2) / wv), h, v_max=2 * R2 + 0.25)
    res = D.evolve_cone(r0, s0, p0, pot, h, ks=16, kv=16)
    tr = D.wall_trajectory(res)
    g = np.isfinite(tr["A_out"]) & tr["resolved"] & (tr["r"] < 1.8)
    assert tr["r"][g].max() > 1.5 * tr["r"][g].min()
    assert np.all(tr["A_out"][g] < 0) and np.all(tr["A_in"][g] > 0)
    assert np.nanmax(res["n_antitrapped"]) > 0
    M = np.nanmedian(tr["m_out"][g][:20])
    tm = ThinWallModel(1.0, 1.0, M)
    sel = g & (tr["r"] > 1.1)
    dev = (np.abs(tr["Rdot"][sel]) - np.sqrt(-tm.V(tr["r"][sel]))) / np.abs(tr["Rdot"][sel])
    assert abs(np.median(dev)) < 0.04


def test_probe_jacobi_zero_modes():
    tr = S.cmc_trajectory(1 / (4 * np.pi), RHO_F, 0.41, 0.05)
    z = S.zero_mode_check(tr)
    assert z["l1_translation"] < 1e-4 and z["l0_time_translation"] < 1e-4
