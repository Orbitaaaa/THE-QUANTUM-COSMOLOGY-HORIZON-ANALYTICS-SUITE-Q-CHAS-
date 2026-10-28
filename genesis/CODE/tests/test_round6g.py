"""Round 6g: the false-vacuum region made by the universe (a supercooled patch in a radiation-dominated universe, with gravity), and why
no machine can make it."""
import warnings

import numpy as np
from scipy.optimize import brentq

from babyuniverse import formation as F, throat_evolution as TE, gsl_solution as GS, accessibility as AC

warnings.simplefilter("ignore")
RHO_F = 3 / 900 / (8 * np.pi)
V = 0.0015


def _frw_setup():
    lam = 4 * RHO_F / (4 * V**4)
    pot = TE.Potential(lam, V, RHO_F)
    pl = F.Plasma(RHO_F, V, lam, f_H=1e-12, T_floor_over_Tc=1e-12)
    eta0 = 20.0
    Tout = brentq(lambda T: pl.eos(T, V)[0] - 3 / (8 * np.pi * eta0**2), 1e-6, 5)
    return pot, pl, eta0, Tout


def test_frw_is_reproduced_in_both_gauges():
    pot, pl, eta0, Tout = _frw_setup()
    h = 0.1
    Z0 = F.frw_cone_data(h, 40.0, pl, pot, eta0, T_out0=Tout)                  # conformal null coordinates: s = ln a, r = a chi
    res = F.evolve_centre(Z0, h, pot, pl, n_rays=300, clip=0)
    i, j = np.indices(res["r"].shape)
    eta, chi = eta0 + 0.5 * (i + j) * h, 0.5 * (j - i) * h
    m = np.isfinite(res["r"]) & (chi > 0.5)
    assert np.nanmax(np.abs(res["r"][m] / ((eta[m] / eta0) * chi[m]) - 1)) < 1e-2
    assert np.nanmax(np.abs(res["s"][m] - np.log(eta[m] / eta0))) < 2e-3
    Z0 = F.frw_cone_data(h, 30.0, pl, pot, eta0, T_out0=Tout, gauge="affine")   # s = 0 on the cone, comoving fluid at eta_f = ln a
    res = F.evolve_centre(Z0, h, pot, pl, n_rays=200, clip=0)
    R, S = res["r"], res["s"]
    ru, rv = np.gradient(R, h, axis=0), np.gradient(R, h, axis=1)
    msh = 0.5 * R * (1 + 4 * ru * rv * np.exp(-2 * S))
    rho = pl.eos(Tout, V)[0] * (np.exp(res["lnTf"]) / Tout) ** 4
    ok = np.isfinite(R) & (R > 1.0)
    ok[:, :3] = False
    ok[:3, :] = False
    assert np.median(np.abs(msh[ok] / (4 * np.pi / 3 * rho[ok] * R[ok] ** 3) - 1)) < 1.5e-2


def test_patches_collapse_or_inflate():
    small = F.cosmic_patch(0.9, 0.5, h=0.1)
    big = F.cosmic_patch(0.9, 1.0, h=0.1)
    assert small["outcome"] == "collapses" and small["anti_fv_cells"] == 0               # a black hole with no daughter
    assert big["outcome"] == "inflates" and big["r_fv_anti_max"] > 100 * big["chi_p"] / 10  # the false vacuum inflates
    lo, hi = F.cosmic_patch(0.9, 0.6, h=0.1), F.cosmic_patch(0.9, 0.62, h=0.1)
    assert lo["outcome"] == "collapses" and hi["outcome"] == "inflates"                   # critical radius 0.57-0.59 l_H


def test_cosmic_daughter_is_hidden_and_gsl_consistent():
    o = F.cosmic_patch(0.5, 0.8, h=0.1)
    assert o["outcome"] == "inflates" and o["rows_with_daughter_behind_trapped_true_vacuum"] > 0
    g = dict(res=o["res"], cls=o["cls"], r_N=o["r_N"], pot=o["pot"])
    eh = GS.parent_event_horizon(g)
    dh = GS.daughter_horizons(g, n_points=4)
    assert eh["u_EH"] is not None and eh["deficit"] == 0.0                              # parent horizon area never decreases
    assert dh["points"] > 0 and dh["max_deficit"] == 0.0                                 # no daughter observer loses horizon area


def test_condensate_turns_a_collapsing_patch_into_a_daughter():
    assert F.cosmic_patch(0.9, 0.6, h=0.1, eps=2.0)["outcome"] == "collapses"
    assert F.cosmic_patch(0.9, 0.6, h=0.1, eps=2.0, lam_f=1e-6)["outcome"] == "inflates"


def test_external_assembly_bounds_are_conditional():
    a = AC.assembly_bounds()
    assert a["log10_power_gap"] > 38                                                     # external delivery within one crossing
    assert AC.assembly_bounds(alpha=1.5)["power_W"] < a["max_luminosity_W"]             # depends on the assembly time (Round 6h)
    assert a["force_N"] < a["max_tension_N"]                                             # below Gibbons' conjectured maximum
    m = AC.machine_audit()
    assert all(r["log10_shortfall"] > 1.5 for r in m["machines"])                        # every audited machine, ideally focused
    assert m["single_collision_quantum_GeV"] > 1e18 and m["linac_100GV_per_m_m"] > 3e16


def test_background_requirements_and_trigger():
    b = AC.background_requirements(1.5e16)                                               # a GUT-scale false vacuum
    assert 1e-31 < b["r_N_m"] < 1e-29 and 1e-4 < b["M_N_kg"] < 1e-2 and b["quantum_GeV"] > 1e13
    assert AC.penrose_farhi_guth()["applies_to_any_assembly_strategy"]
    lo, hi = F.cosmic_patch(0.9, 0.60, h=0.1), F.cosmic_patch(0.9, 0.62, h=0.1)         # a positive-energy trigger of ~10% of the
    assert lo["outcome"] == "collapses" and hi["outcome"] == "inflates"                  # patch's vacuum energy, in a cosmological
    assert (0.62 / 0.60) ** 3 - 1 < 0.11                                                 # (anti-trapped) background
