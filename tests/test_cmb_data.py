import json

import healpy as hp
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import eval_legendre

from qchas.cmb import harmonics as H
from qchas.cmb import maps as M
from qchas.cmb import stats as S
from qchas.data.fetch import available

needs_wmap = pytest.mark.skipif(
    available("wmap_ilc") is None or available("wmap_kq75") is None or available("wmap_tt") is None,
    reason="WMAP 9-yr files not in data/raw",
)
needs_planck = pytest.mark.skipif(
    available("planck_commander") is None or available("planck_tt") is None,
    reason="Planck PR3 files not in data/raw",
)

#: healpy-data pixel_window_functions/pixel_window_n0064.fits, temperature column
HEALPY_PIXWIN_64 = {2: 0.99993174405345, 10: 0.998749256132102, 30: 0.9894664229888386, 64: 0.9535939863969998}


def _toy_cl(lmax: int) -> np.ndarray:
    """Sachs-Wolfe-like plateau D_l = 1000 uK^2 (test input only)."""
    ell = np.arange(lmax + 1, dtype=float)
    cl = np.zeros(lmax + 1)
    cl[2:] = 2 * np.pi * 1000.0 / (ell[2:] * (ell[2:] + 1))
    return cl


@pytest.fixture(scope="module")
def sky():
    if available("wmap_ilc") is None:
        pytest.skip("WMAP ILC map not in data/raw")
    return M.prepare_sky("wmap")


@pytest.fixture(scope="module")
def wmap_tt():
    if available("wmap_tt") is None:
        pytest.skip("WMAP TT spectrum not in data/raw")
    return S.load_wmap_tt()


# --- maps -----------------------------------------------------------------------------------

def test_unit_factor_and_beam():
    assert M._unit_factor("mK, thermodynamic") == 1e3
    assert M._unit_factor("K_CMB") == 1e6
    assert M._unit_factor("uK") == 1.0
    for bad in ("K_RJ", "mK, antenna", "MJy/sr", "counts"):
        with pytest.raises(ValueError):
            M._unit_factor(bad)
    assert np.allclose(M.gaussian_beam(60.0, 200), hp.gauss_beam(np.radians(1.0), lmax=200), rtol=1e-12)


def test_degrade_is_nested_pixel_average():
    rng = np.random.default_rng(0)
    m = rng.standard_normal(hp.nside2npix(128))
    d = M.degrade_nested(m, 32)
    assert np.allclose(d, hp.ud_grade(m, 32, order_in="NESTED", order_out="RING"), atol=1e-12)
    m[:16] = hp.UNSEEN
    m[16:20] = np.nan
    d = M.degrade_nested(m, 32)
    first = hp.nest2ring(32, 0)
    assert d[first] == hp.UNSEEN
    assert np.isclose(d[hp.nest2ring(32, 1)], m[20:32].mean())
    assert np.isclose(d[hp.nest2ring(32, 2)], m[32:48].mean())


def test_exact_pixel_window_matches_healpy_file():
    w = M.exact_pixel_window(64, 64)
    for l, ref in HEALPY_PIXWIN_64.items():
        assert abs(w[l] - ref) < 2e-5 * (1 if l <= 30 else 5)
    pw, method = M.pixel_window(64, 30)
    assert method in ("healpy", "exact") and abs(pw[30] - HEALPY_PIXWIN_64[30]) < 2e-5


def test_pipeline_recovers_simulated_alm_fast():
    r = M.simulate_pipeline(256, M.WMAP_ILC_BEAM_FWHM_ARCMIN, seed=1, boosted_lmax=0)
    assert r["relative_residual"] < 0.003
    assert np.all(np.abs(np.array(r["gain_per_l"]) - 1) < 3e-3)
    assert r["alm_in"].shape == (957,)


def test_pipeline_toy_spectrum_planck_beam():
    r = M.simulate_pipeline(256, M.PLANCK_COMMANDER_BEAM_FWHM_ARCMIN, seed=2, cl=_toy_cl(767))
    assert r["relative_residual"] < 0.01


@needs_wmap
def test_wmap_ilc_header_units_and_order():
    full = M.load_wmap_ilc()
    assert full.meta["ordering_in"] == "NESTED" and full.nside == 512
    assert full.meta["unit_factor_to_uK"] == 1e3 and full.beam_fwhm_arcmin == 60.0
    assert 50.0 < np.std(full.map_uK) < 150.0


def test_planck_absent_is_none(tmp_path):
    assert M.load_planck_commander(tmp_path) is None
    assert S.load_planck_tt(tmp_path) is None
    with pytest.raises(FileNotFoundError):
        M.prepare_sky("planck", raw_dir=tmp_path)


def _write_commander_like(path, inp_K: np.ndarray, nside: int) -> None:
    """PR3 Commander column layout (NESTED, K_CMB); I_STOKES deliberately differs from I_STOKES_INP."""
    zero = np.zeros(12 * nside**2, dtype=np.float32)
    names = ["I_STOKES", "Q_STOKES", "U_STOKES", "TMASK", "PMASK", "I_STOKES_INP", "Q_STOKES_INP", "U_STOKES_INP"]
    cols = [zero + 1.0, zero, zero, zero + 1, zero + 1, inp_K.astype(np.float32), zero, zero]
    units = ["K_CMB"] * 3 + ["NONE"] * 2 + ["K_CMB"] * 3
    hp.write_map(str(path), cols, nest=True, coord="G", column_names=names, column_units=units,
                 fits_IDL=False, dtype=np.float32, overwrite=True)


def test_planck_branch_on_synthetic_files(tmp_path, monkeypatch):
    """The Planck code path end to end on small files in the PR3 formats (no KQ75 mask on disk)."""
    from qchas.data.registry import DATASETS

    def on_disk(key, raw_dir=None):  # the registry's size floor (1 GB) would hide the small test files
        p = tmp_path / DATASETS[key].filename
        return p if p.exists() else None

    monkeypatch.setattr(M, "available", on_disk)
    monkeypatch.setattr(S, "available", on_disk)
    nside_in = 128
    rng = np.random.default_rng(11)
    alm = M.gaussian_alm(_toy_cl(3 * nside_in - 1), 3 * nside_in - 1, rng)
    sky_ring = hp.alm2map(hp.almxfl(alm, M.gaussian_beam(5.0, 3 * nside_in - 1)), nside_in, lmax=3 * nside_in - 1)
    inp_K = hp.reorder(sky_ring, r2n=True) * 1e-6
    _write_commander_like(tmp_path / DATASETS["planck_commander"].filename, inp_K, nside_in)
    ell = np.arange(2, 2509)
    rows = np.column_stack([ell, 1000.0 + ell, np.full(ell.size, 100.0), np.full(ell.size, 160.0)])
    np.savetxt(tmp_path / DATASETS["planck_tt"].filename, rows, fmt=["%6d", "%15.8E", "%15.8E", "%15.8E"],
               header="   L    TT             -dTT            +dTT")

    full = M.load_planck_commander(tmp_path)
    assert full.meta["column"] == "I_STOKES_INP" and full.meta["unit_factor_to_uK"] == 1e6
    assert full.beam_fwhm_arcmin == 5.0 and full.nside == nside_in
    sky_p = M.prepare_sky("auto", raw_dir=tmp_path)
    assert sky_p.source_name == "Planck-PR3-Commander"
    a_in = H.healpy_to_real(alm, 2, 30, 3 * nside_in - 1)
    assert np.linalg.norm(sky_p.alm_real - a_in) / np.linalg.norm(a_in) < 0.01
    assert sky_p.mask_fraction is None and any("KQ75" in w for w in sky_p.meta["warnings"])
    meta, arrays = M.sky_map_export(sky_p)
    assert meta["mask_file"] is None and "sky_mask.bin" not in arrays and meta["source"] == "Planck-PR3-Commander"

    planck = S.load_planck_tt(tmp_path)
    assert np.allclose(planck.select(2, 4).err_minus, 100.0) and np.allclose(planck.select(2, 4).err_plus, 160.0)
    model = S.cl_from_dl(np.arange(2, 31), np.full(29, 1000.0))
    f = S.cmb_fields(sky_p, None, planck, model_cl_band=model, n_realizations=5_000, nside_search=8)
    json.dumps(f)
    assert f["source"] == "Planck-PR3-Commander" and f["stats"]["primary"] == "planck"
    assert np.isclose(f["stats"]["D2_obs"], 1002.0) and f["observed"]["wmap"] is None
    assert f["observed"]["planck"]["err_plus"][0] == 160.0
    assert f["stats"]["extraction_check"]["published_D"][0] == 1002.0
    assert len(f["stats"]["extraction_check"]["paper_D_extracted_rel_diff"]) == 5
    assert f["stats"]["paper_check"]["D2_planck_spectrum"] == pytest.approx(1002.0)
    assert not any("I_STOKES_INP" in w for w in sky_p.meta["warnings"])

    inp_K[:4] = hp.UNSEEN  # the four NESTED children of one Nside-64 pixel
    _write_commander_like(tmp_path / DATASETS["planck_commander"].filename, inp_K, nside_in)
    sky_u = M.prepare_sky("planck", raw_dir=tmp_path)
    assert sky_u.meta["extraction"]["unseen_pixels_zeroed"] == 1
    assert any("no valid input" in w for w in sky_u.meta["warnings"])


@needs_wmap
def test_prepare_sky_wmap(sky):
    npix = hp.nside2npix(64)
    assert sky.source_name == "WMAP9-ILC" and sky.nside == 64
    assert sky.map_uK.shape == sky.lowl_map_uK.shape == sky.mask.shape == (npix,)
    assert sky.alm_real.shape == (957,) and list(sky.ell) == list(range(2, 31))
    assert abs(sky.map_uK.mean()) < 1e-6
    assert np.allclose(hp.map2alm(sky.map_uK, lmax=1, iter=4)[:3], 0, atol=1e-2)
    th, ph = hp.pix2ang(64, np.arange(0, npix, 97))
    assert np.allclose(H.real_ylm(2, 30, th, ph) @ sky.alm_real, sky.lowl_map_uK[::97], atol=1e-6)
    assert np.allclose(sky.D_ell, M.dl_from_cl(sky.ell, M.cl_from_real(sky.alm_real, 2, 30)))
    assert 0.6 < sky.fsky < 0.75 and set(np.unique(sky.mask)) <= {0, 1}
    assert sky.meta["extraction"]["iterations"] == 4 and sky.meta["extraction"]["lmax_alm"] == 64


@needs_wmap
def test_wmap_map_matches_direct_nside512_analysis(sky):
    full = M.load_wmap_ilc()
    m = hp.remove_dipole(hp.reorder(full.map_uK, n2r=True))
    alm = hp.almxfl(hp.map2alm(m, lmax=1023, iter=3), 1 / M.gaussian_beam(60.0, 1023))
    direct = H.healpy_to_real(alm, 2, 30, 1023)
    assert np.linalg.norm(direct - sky.alm_real) / np.linalg.norm(direct) < 0.005


@needs_wmap
def test_wmap_map_D_consistent_with_published(sky, wmap_tt):
    w = wmap_tt.select(2, 30)
    assert np.all(np.abs(sky.D_ell - w.D) < w.err)


@needs_wmap
def test_sky_map_export(sky):
    meta, arrays = M.sky_map_export(sky)
    json.dumps(meta)
    assert meta["nside"] == 64 and meta["ordering"] == "RING" and meta["source"] == "WMAP9-ILC"
    assert [m["file"] for m in meta["maps"]] == ["sky_map_full.bin", "sky_map_lowl.bin"]
    for name in ("sky_map_full.bin", "sky_map_lowl.bin"):
        a = arrays[name]
        assert a.dtype == np.dtype("<f4") and a.size == 12 * 64**2
    assert arrays["sky_mask.bin"].dtype == np.uint8
    assert np.array_equal(S.sky_map_export(sky)[1]["sky_mask.bin"], arrays["sky_mask.bin"])


# --- stats ----------------------------------------------------------------------------------

def test_load_wmap_tt(wmap_tt):
    assert wmap_tt.ell[0] == 2 and np.isclose(wmap_tt.D[0], 150.6398)
    b = wmap_tt.select(2, 30)
    assert b.D.size == 29 and np.allclose(b.extra["err_cosmic_variance"][0], 708.2630)


def test_s12_quadrature():
    rng = np.random.default_rng(4)
    cl = rng.uniform(10.0, 1000.0, 29) / np.arange(2, 31) ** 2
    s = float(S.s_half(cl))
    assert np.isclose(s, S.s_half_direct(cl), rtol=1e-12)

    def c_theta(x):
        return sum((2 * l + 1) / (4 * np.pi) * c * eval_legendre(l, x) for l, c in zip(range(2, 31), cl))

    ref, _ = quad(lambda x: c_theta(x) ** 2, -1.0, 0.5, limit=400, epsabs=0, epsrel=1e-12)
    assert np.isclose(s, ref, rtol=1e-9)
    assert np.allclose(S.s_half(np.stack([cl, 2 * cl])), [s, 4 * s])


def test_realizations_and_pvalues():
    model = S.cl_from_dl(np.arange(2, 31), np.full(29, 1000.0))
    sims = S.lcdm_realizations(model, n=40_000, seed=1)
    assert np.allclose(sims.mean(axis=0), model, rtol=0.03)
    obs = model.copy()
    obs[0] *= 0.25
    p = S.anomaly_pvalues(obs, model, sims=sims)
    assert abs(p["p_C2"] - p["p_C2_analytic"]) < 0.005
    assert np.isclose(p["D2"], 250.0)


def test_axis_quadratic_form_and_known_axis():
    rng = np.random.default_rng(5)
    for l in (2, 3):
        a = rng.standard_normal(2 * l + 1)
        v = rng.standard_normal((10, 3))
        v /= np.linalg.norm(v, axis=1)[:, None]
        M3 = S.axis_matrix(a, l)
        assert np.allclose(S.dispersion_by_rotation(a, l, v), np.einsum("ni,ij,nj->n", v, M3, v), atol=1e-10)
        assert np.isclose(np.trace(M3), l * (l + 1) * a @ a)
        sectoral = np.zeros(2 * l + 1)
        sectoral[0], sectoral[-1] = 1.0, -0.3
        rot = H.random_rotation(rng)
        ax = S.multipole_axis(H.real_rotation_block(l, rot) @ sectoral, l, nside_search=16)
        assert S.axis_angle_deg(ax["vector"], rot[:, 2]) < 1e-6
        assert ax["grid_offset_deg"] < 4.0 and np.isclose(ax["dispersion"], 1.0)


@needs_wmap
def test_wmap_axes(sky):
    al = S.quadrupole_octopole_alignment(sky.alm_real, nside_search=16, source=sky.source_name)
    for key in ("quadrupole", "octopole"):
        assert 230 < al[key]["l"] < 250 and 55 < al[key]["b"] < 66
    assert al["angle_deg"] < 6.0
    assert np.isclose(al["p_isotropic"], 1 - np.cos(np.radians(al["angle_deg"])))


@needs_wmap
def test_cmb_fields_structure(sky, wmap_tt):
    model = S.cl_from_dl(np.arange(2, 31), np.full(29, 1000.0))
    f = S.cmb_fields(sky, wmap_tt, None, model_cl_band=model, n_realizations=20_000, nside_search=8)
    json.dumps(f)
    assert f["source"] == "WMAP9-ILC" and f["ell"] == list(range(2, 31))
    assert f["observed"]["planck"] is None and len(f["observed"]["wmap"]["D"]) == 29
    assert len(f["map_D"]) == 29 and np.allclose(f["lcdm_D"], 1000.0)
    st = f["stats"]
    assert st["primary"] == "wmap" and np.isclose(st["D2_obs"], 150.6398)
    assert st["p_C2"]["planck"] is None and 0 < st["p_C2"]["wmap"] < st["p_C2"]["map"] < 1
    assert st["axes"]["source"] == "WMAP9-ILC"
    assert st["extraction_check"]["published_D"][0] == pytest.approx(150.6398)
    for key in ("D2_obs", "S12_obs", "S12_lcdm_median", "p_C2", "p_S12"):  # docs/ARCHITECTURE.md cmb.json stats
        assert key in st
    # keys read by qchas.validation._group_cmb_stats
    pc = st["paper_check"]
    assert 0 < pc["p_S12_planck_value"] < 1 and 0 < pc["p_C2_planck_value"] < 1 and pc["S12_lcdm_median"] > 0


@needs_wmap
def test_statistics_on_a_sub_band(sky, wmap_tt):
    model = S.cl_from_dl(np.arange(3, 21), np.full(18, 1000.0))
    f = S.cmb_fields(sky, wmap_tt, None, lmin=3, lmax=20, model_cl_band=model, n_realizations=2_000, nside_search=4)
    json.dumps(f)
    st = f["stats"]
    assert f["map_D"] == pytest.approx(sky.D_ell[1:19].tolist())
    assert st["D2_obs"] is None and st["p_C2"]["wmap"] is None and st["p_S12"]["map"] is not None
    assert "axes" not in st and st["paper_check"]["p_S12_planck_value"] is None


# --- slow: paper numbers and full-resolution simulations ------------------------------------

@pytest.mark.slow
def test_lcdm_reproduces_paper_statistics():
    model = S.lcdm_cl_band()
    sims = S.lcdm_realizations(model)
    s = S.s_half(sims)
    assert abs(S.dl_from_cl(2, model[0]) / S.PAPER_STATS["D2_lcdm"] - 1) < 0.005
    assert abs(np.median(s) / S.PAPER_STATS["S12_lcdm_median"] - 1) < 0.015
    assert abs(np.mean(s <= S.PAPER_STATS["S12_planck"]) - S.PAPER_STATS["p_S12_planck"]) < 0.005
    assert abs(np.mean(sims[:, 0] <= S.cl_from_dl(2, S.PAPER_STATS["D2_planck"])) - S.PAPER_STATS["p_C2_planck"]) < 0.003


@pytest.mark.slow
def test_pipeline_recovers_simulated_alm_wmap_resolution():
    r = M.simulate_pipeline(512, M.WMAP_ILC_BEAM_FWHM_ARCMIN, seed=0)
    assert r["relative_residual"] < 0.003


@pytest.mark.slow
def test_pipeline_recovers_simulated_alm_planck_like():
    r = M.simulate_pipeline(1024, M.PLANCK_COMMANDER_BEAM_FWHM_ARCMIN, seed=0)
    assert r["relative_residual"] < 0.01


@pytest.mark.slow
@needs_planck
def test_planck_pipeline_against_paper():
    sky_p = M.prepare_sky("auto")
    assert sky_p.source_name == "Planck-PR3-Commander"
    d = np.asarray(sky_p.D_ell[:5])
    # same map and procedure: only the paper's integer rounding and its simulated pixel-window
    # correction separate the two (a 1 deg beam applied by mistake would already move D_6 by 3 %)
    assert np.allclose(d, M.PAPER_PIPELINE["D_extracted"], rtol=0.02)
    planck = S.load_planck_tt()
    assert np.allclose(planck.select(2, 6).D, M.PAPER_PIPELINE["D_published"], rtol=0.01)
    assert abs(float(S.s_half(planck.cl_band())) / S.PAPER_STATS["S12_planck"] - 1) < 0.01
