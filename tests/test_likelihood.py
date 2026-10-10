"""Orientation-scanned torus likelihood and the pre-registered resonance search, [TOP] Secs. 9 and 11.6."""

import json
import time
from types import SimpleNamespace

import numpy as np
import pytest

from qchas.cmb import covariance as C
from qchas.cmb import harmonics as H
from qchas.cmb import likelihood as L
from qchas.cmb import resonance as Rz

SMALL = L.SO3Grid(12, 7, 24)


@pytest.fixture(scope="module")
def sky_lcdm():
    rng = np.random.default_rng(2026)
    return np.sqrt(C.lcdm_covariance()) * rng.standard_normal(C.lcdm_covariance().size)


@pytest.fixture(scope="module")
def eng_small():
    return L.Engine([("FCC", 1.05), ("BCC", 1.0), ("cubic", 1.2)], SMALL, search=True, keep_matrices=True)


# --- (1) Fourier / FFT scan equals brute-force rotation of the data ---------------------------------

def test_fourier_scan_equals_bruteforce_on_grid(eng_small, sky_lcdm):
    rng = np.random.default_rng(1)
    for t in range(len(eng_small.models)):
        v = eng_small.grid_values(sky_lcdm, t)
        for _ in range(8):
            j, ia, k = rng.integers(v.shape[0]), rng.integers(v.shape[1]), rng.integers(v.shape[2])
            rot = H.rotation_from_euler(eng_small.alpha[ia], eng_small.beta[j], eng_small.gamma[k])
            assert abs(v[j, ia, k] - eng_small.dlnl_direct(sky_lcdm, t, rot)) < 1e-9


def test_fourier_series_equals_bruteforce_at_random_rotations(eng_small, sky_lcdm):
    ac = L._to_complex(sky_lcdm, eng_small.lay)
    rng = np.random.default_rng(2)
    for _ in range(6):
        rot = H.random_rotation(rng)
        al, be, ga = H.euler_zyz(rot)
        for t in range(len(eng_small.models)):
            m = eng_small.models[t]
            q = L._q_coeffs(L._split(L._ct_points(ac, eng_small.lay, al, be), eng_small.lay), m, eng_small.lay)
            four = -0.5 * (L._series_at(q, np.array([ga]))[0] + m.logdet_ratio)
            assert abs(four - eng_small.dlnl_direct(sky_lcdm, t, rot)) < 1e-9
        # the series is real and only multiples of the detected period s survive (cube quarter turn)
        assert np.max(np.abs(q[:, 1:][:, np.arange(1, L.DMAX + 1) % eng_small.s != 0])) < 1e-9 * np.max(np.abs(q))


def test_symmetry_reduction_is_exact(sky_lcdm):
    models = [("FCC", 1.05), ("cubic", 1.2)]
    e1 = L.Engine(models, SMALL, search=True)
    e0 = L.Engine(models, SMALL, search=True, symmetry=False)
    assert e1.s == 4 and e1.reflect and e0.s == 1 and not e0.reflect
    r1, r0 = e1.scan(sky_lcdm, refine=False), e0.scan(sky_lcdm, refine=False)
    assert np.allclose(r1["lnB"], r0["lnB"], rtol=0, atol=1e-10)
    assert np.allclose(r1["max_grid"], r0["max_grid"], rtol=0, atol=1e-10)
    assert np.allclose(r1["line_lnB"], r0["line_lnB"], rtol=0, atol=1e-10)
    assert e1.symmetry_info()["scanned_fraction"] < 0.15


# --- (2) Woodbury line statistic equals dense evaluation --------------------------------------------

def test_woodbury_line_statistic_equals_dense(eng_small, sky_lcdm):
    t = [m.shape for m in eng_small.models].index("cubic")
    m = eng_small.models[t]
    v = eng_small.line_grid_values(sky_lcdm, t)
    c0 = C.torus_covariance("cubic", 1.2)
    cl = C.lcdm_covariance()
    rng = np.random.default_rng(3)
    for _ in range(6):
        j, ia, k, i = (rng.integers(n) for n in v.shape)
        md = m.line.models[i]
        rot = H.rotation_from_euler(eng_small.alpha[ia], eng_small.beta[j], eng_small.gamma[k])
        u, w = Rz.line_lowrank("cubic", 1.2, md["shell_m2"], md["f"])
        iso = cl + Rz.line_isotropic_diag("cubic", 1.2, md["shell_m2"], md["f"])
        dense = Rz.gaussian_loglike(sky_lcdm, C.rotate_covariance(c0 + (u * w) @ u.T, rot)) - \
            Rz.gaussian_loglike(sky_lcdm, np.diag(iso))
        assert abs(v[j, ia, k, i] - dense) < 1e-8


# --- (3) injection: a hidden FCC orientation is recovered ----------------------------------------------

def test_injection_recovers_orientation():
    """FCC torus at d = D with a hidden Haar orientation (paper: recovered to 0.17 deg at its grid).

    Over 15 injections (FCC 0.7-1.0 D) our maximum lies 0.26-0.92 deg from the truth, and always at a
    HIGHER likelihood than the true orientation (by 0.3-4.4 ~ chi^2_3 / 2): the offset is the statistical
    scatter of the maximum-likelihood orientation, not a limit of the search (zoom resolution 0.05 deg).
    """
    eng = L.Engine([("FCC", 1.0)], L.SO3Grid(*L.PAPER_GRID), keep_matrices=True)
    chol = np.linalg.cholesky(C.torus_covariance("FCC", 1.0))
    rng = np.random.default_rng([42, 0])
    r_true = H.random_rotation(rng)
    a = H.rotate_real_vector(chol @ rng.standard_normal(chol.shape[0]), r_true, C.LMIN, C.LMAX)
    out = eng.scan(a, refine=True, lines=False)
    err = L.orientation_error_deg(out["rot"][0], r_true)
    assert err < 1.5, err
    assert out["max"][0] >= eng.dlnl_direct(a, 0, r_true) - 1e-6  # the search found a maximum at least as high
    assert abs(out["max"][0] - eng.dlnl_direct(a, 0, out["rot"][0])) < 1e-8
    assert out["lnB"][0] > 5.0  # strong evidence when the torus is real and d <= D


def test_orientation_error_modulo_cubic_group():
    g = L.cubic_rotation_group()
    assert g.shape == (24, 3, 3) and np.allclose(np.linalg.det(g), 1.0)
    rng = np.random.default_rng(7)
    r = H.random_rotation(rng)
    for gm in g[::5]:
        assert L.orientation_error_deg(r @ gm, r) < 1e-6
    assert [len(L.lattice_axes(s)) for s in ("cubic", "BCC", "FCC")] == [3, 4, 6]
    axes = L.axes_galactic("FCC", r)
    assert len(axes) == 6 and all(b >= 0 for _, b in axes)


# --- (4) null skies ------------------------------------------------------------------------------------

def test_null_skies_lnB_centred_at_or_below_zero():
    """E_null[B] = 1 exactly (each torus law is normalized), hence E[ln B] <= 0 by Jensen."""
    eng = L.Engine([("FCC", 1.05), ("BCC", 1.1), ("cubic", 1.2)], L.SO3Grid(*L.PAPER_SEARCH_GRID))
    rng = np.random.default_rng(11)
    lnb = np.array([eng.scan(np.sqrt(eng.c_l) * rng.standard_normal(eng.c_l.size), refine=False, lines=False)["lnB"]
                    for _ in range(16)])
    assert np.all(lnb.mean(axis=0) < 0.05)
    assert np.median(lnb[:, 0]) < 0.0
    assert abs(np.exp(lnb[:, 2]).mean() - 1.0) < 0.1  # normalization check at cubic 1.2 D (sd of mean ~ 0.03)
    assert np.isclose(eng.w_rows.sum(), 1.0) and np.isclose(SMALL.w_beta.sum(), 1.0)


# --- (5) bundle contract ---------------------------------------------------------------------------------

def test_export_fields_contract(sky_lcdm):
    from qchas.export import _clean

    sky = SimpleNamespace(alm_real=sky_lcdm, source_name="simulated-LCDM")
    out = L.export_fields(sky, fast=True, null_sims=3, workers=2, true_sims=2, grid=(12, 7, 24))
    text = json.dumps(_clean(out), allow_nan=False)
    back = json.loads(text)
    assert set(back) == {"likelihood", "likelihood_meta", "search"}
    row = back["likelihood"][0]
    for k in ("shape", "d_over_D", "max_dlnL", "lnB", "null_q95", "p", "n_null", "best_axes_gal", "paper", "computed"):
        assert k in row
    assert row["computed"] is True and row["n_null"] == 3 and "exclusion" in row
    fcc105 = next(r for r in back["likelihood"] if r["shape"] == "FCC" and r["d_over_D"] == 1.05)
    assert fcc105["paper"]["lnB"] == -1.20 and fcc105["paper"]["max_dlnL"] == 7.20 and fcc105["paper"]["verdict"] == "allowed"
    fcc0985 = next(r for r in back["likelihood"] if r["shape"] == "FCC" and r["d_over_D"] == 0.985)
    assert fcc0985["paper"]["lnB_planck"] == -22.14 and fcc0985["paper"]["verdict"] == "excluded (95%)"
    meta = back["likelihood_meta"]
    for k in ("map", "n_null", "grid", "statistics", "runtime_s", "notes"):
        assert k in meta
    assert {"n_alpha", "n_beta", "n_gamma"} <= set(meta["grid"])
    assert any("independent re-analysis" in n.lower() for n in meta["notes"])
    s = back["search"]
    for k in ("map", "models", "S_max", "S_avg", "best", "null", "p_max", "p_avg", "n_null", "paper", "computed"):
        assert k in s
    assert len(s["models"]) == 160 and len(s["null"]["S_max"]) == 3
    assert {"shape", "d_over_D", "shell_m2", "k_chi", "n_vectors", "f", "lnB"} <= set(s["models"][0])
    # the pre-registered grid of [TOP] Sec. 11.6: 80 shells (<= 12 vectors, 2 <= k0 chi_* <= 14) x f in {0.2, 0.4}
    assert {m["f"] for m in s["models"]} == {0.2, 0.4}
    assert all(m["n_vectors"] <= 12 and 2.0 <= m["k_chi"] <= 14.0 for m in s["models"])
    assert {(m["shape"], m["d_over_D"]) for m in s["models"]} == {(x, d) for x in ("FCC", "BCC", "cubic")
                                                                   for d in (1.05, 1.15, 1.2, 1.3)}
    assert s["paper"]["S_max"] == 1.01 and s["paper"]["p_max"] == 0.78
    assert s["S_max"] == max(m["lnB"] for m in s["models"])


@pytest.mark.slow
def test_export_fields_fast_mode_runtime():
    """fast=True (CI smoke build) must finish in under ~2 minutes."""
    from qchas.cmb import maps

    try:
        sky = maps.prepare_sky()
    except Exception:  # pragma: no cover - raw data not on this host
        rng = np.random.default_rng(5)
        sky = SimpleNamespace(alm_real=np.sqrt(C.lcdm_covariance()) * rng.standard_normal(957), source_name="simulated")
    t0 = time.perf_counter()
    out = L.export_fields(sky, fast=True)
    assert time.perf_counter() - t0 < 120.0
    assert out["likelihood_meta"]["n_null"] == 12
