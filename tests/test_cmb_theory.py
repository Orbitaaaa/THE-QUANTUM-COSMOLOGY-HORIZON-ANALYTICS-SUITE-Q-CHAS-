"""Transfer functions, [TOP] Eq. (2) torus covariance, KL / Pinsker, and topological resonance."""

import numpy as np
import pytest
from scipy.interpolate import CubicSpline
from scipy.special import sph_harm_y, spherical_jn

from qchas.cmb import covariance as C
from qchas.cmb import harmonics as H
from qchas.cmb import resonance as Rz
from qchas.cmb.transfer import load_transfer


@pytest.fixture(scope="module")
def tf():
    return load_transfer()


# --- transfer functions ------------------------------------------------------------------

def test_transfer_reproduces_camb_cl(tf):
    cl = tf.cl_isotropic(2, 30)
    assert np.max(np.abs(cl / tf.cl_camb[2:31] - 1.0)) < 3e-3
    assert np.max(np.abs(cl / tf.cl_camb[2:31] - 1.0)) < 1e-3  # ~2e-4; l = 18 carries the +0.08 % repair
    assert abs(tf.unit_factor_measured / tf.t_cmb_uk**2 - 1.0) < 1e-3
    assert abs(tf.chi_star - 13870.0) < 2.0  # [TOP] Sec. 7
    assert np.array_equal(tf.ells, np.arange(2, 41))


@pytest.mark.slow
def test_transfer_matches_camb_model_lcdm_cl(tf):
    from qchas.cmb.camb_model import lcdm_cl

    cl = lcdm_cl(40)
    assert np.max(np.abs(tf.cl_isotropic(2, 30) / cl[2:31] - 1.0)) < 3e-3
    d2 = 3.0 * cl[2] / np.pi
    assert abs(d2 - 1016.0) < 1.0  # [TOP] Sec. 8


def test_transfer_spline_accuracy(tf):
    """Leave-one-out: a spline through every other CAMB node predicts the others (spacing doubled)."""
    x = tf.q * tf.chi_star
    sel = np.flatnonzero((x > 0.5) & (x < 160.0))
    keep, test = sel[::2], sel[1::2]
    spl = CubicSpline(tf.q[keep], tf.delta_q[:, keep].T, axis=0)
    err = np.abs(spl(tf.q[test]).T - tf.delta_q[:, test])
    scale = np.max(np.abs(tf.delta_q[:, sel]), axis=1, keepdims=True)
    assert np.max(err / scale) < 2e-2  # h^4 scaling: <~1e-3 at CAMB's own spacing
    full = tf.delta(np.arange(2, 41)[:, None], tf.q[None, sel])
    assert np.allclose(full, tf.delta_q[:, sel], rtol=0, atol=1e-9 * np.max(np.abs(tf.delta_q)))


def test_camb_low_k_truncation_repair(tf):
    """CAMB 2.0.4 zeroes Delta_l below k chi ~ 14.5 for 18 <= l <= 27; the repair is small in C_l."""
    raw = load_transfer(repair=False)
    assert set(tf.repaired_ells.tolist()) <= set(range(18, 28)) and 18 in tf.repaired_ells
    i = 18 - 2
    first = np.flatnonzero(tf.delta_q_camb[i] != 0)[0]
    assert abs(tf.delta_q_camb[i, first]) > 0.05 * np.max(np.abs(tf.delta_q_camb[i]))
    x = tf.q * tf.chi_star
    ratio = tf.delta_q[i, first - 3:first + 3] / spherical_jn(18, x[first - 3:first + 3])
    assert np.all(tf.delta_q[i, first - 3:first] != 0)
    assert np.max(np.abs(np.diff(ratio))) < 0.02 * np.abs(ratio).max()  # Delta_l / j_l continuous
    dcl = tf.cl_isotropic(2, 30) / raw.cl_isotropic(2, 30) - 1
    assert 0 < dcl[16] < 1.5e-3 and np.max(np.abs(np.delete(dcl, [16, 17]))) < 1e-4


@pytest.mark.slow
def test_transfer_against_denser_camb_sampling(tf):
    """Independent check of the spline: CAMB rerun with 4x denser k sampling (IntkAccuracyBoost = 4)."""
    import camb

    from qchas.constants import PLANCK2018

    pars = camb.set_params(H0=PLANCK2018["H0"], ombh2=PLANCK2018["ombh2"], omch2=PLANCK2018["omch2"],
                           tau=PLANCK2018["tau"], As=PLANCK2018["As"], ns=PLANCK2018["ns"],
                           mnu=PLANCK2018["mnu"], lmax=200)
    pars.set_accuracy(AccuracyBoost=2.0, lSampleBoost=50.0, lAccuracyBoost=2.0)
    pars.Accuracy.IntkAccuracyBoost = 4.0
    td = camb.get_results(pars).get_cmb_transfer_data(tp="scalar")
    q = np.asarray(td.q)
    dense = np.asarray(td.delta_p_l_k[0][:39]) * tf.t_cmb_uk
    sel = (q * tf.chi_star > 0.5) & (q * tf.chi_star < 160.0)
    ours = tf.delta_table(q[sel]).T
    scale = np.max(np.abs(dense[:, sel]), axis=1, keepdims=True)
    valid = dense[:, sel] != 0  # the dense run has its own low-k truncation (18 <= l <= 26)
    assert np.max(np.where(valid, np.abs(ours - dense[:, sel]), 0.0) / scale) < 1e-3


def test_delta_vectorization_and_limits(tf):
    k = np.array([1e-7, 0.01, 0.1, 1.0, 5.0, 150.0]) / tf.chi_star
    tab = tf.delta_table(k, 2, 40)
    assert tab.shape == (6, 39)
    assert np.allclose(tf.delta(np.arange(2, 41)[None, :], k[:, None]), tab)
    assert np.all(np.abs(tab[0]) < 1e-6 * np.abs(tab[3]).max())
    assert tf.k_fine[-1] * tf.chi_star >= 150.0
    with pytest.raises(ValueError):
        tf.delta(41, 0.001)


def test_continuum_split_is_additive(tf):
    kc = 70.0 / tf.chi_star
    lo = tf.cl_isotropic(2, 30, kmax=kc)
    hi = tf.cl_isotropic(2, 30, kmin=kc)
    assert np.allclose(lo + hi, tf.cl_isotropic(2, 30), rtol=1e-6)


# --- lattices and harmonics ---------------------------------------------------------------

@pytest.mark.parametrize("shape,vol", [("cubic", 1.0), ("BCC", 4 / (3 * np.sqrt(3))), ("FCC", 1 / np.sqrt(2))])
def test_lattice_geometry(shape, vol):
    b = C.lattice_basis(shape)
    assert abs(C.shortest_vector_length(b) - 1.0) < 1e-12
    assert abs(abs(np.linalg.det(b)) - vol) < 1e-12
    k = C.dual_basis(b)
    assert np.allclose(b @ k.T, 2 * np.pi * np.eye(3))
    full = C.dual_vectors(b, 60.0)
    half = C.dual_vectors(b, 60.0, half=True)
    assert full.shape[0] == 2 * half.shape[0]
    expected = 4 / 3 * np.pi * 60.0**3 * vol / (2 * np.pi) ** 3
    assert abs(full.shape[0] / expected - 1) < 0.1


def test_sec8_volumes():
    for r in C.PAPER_SEC8["rows"]:
        assert abs(C.volume_over_D3(r["shape"], r["d_over_D"]) - r["V_over_D3"]) < 6e-4


def test_fast_real_ylm_matches_harmonics():
    rng = np.random.default_rng(7)
    v = rng.standard_normal((300, 3))
    v[0] = (0, 0, 1)
    v[1] = (0, 0, -1)
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    assert np.allclose(C.real_ylm_unit(2, 30, v), H.real_ylm_vec(2, 30, v), atol=1e-11)


def _complex_eq2(kvecs, vol, tf, lmin, lmax):
    """Literal [TOP] Eq. (2) with complex Y_lm, converted to the real basis by harmonics' unitary."""
    k = np.linalg.norm(kvecs, axis=1)
    th = np.arccos(kvecs[:, 2] / k)
    ph = np.arctan2(kvecs[:, 1], kvecs[:, 0])
    pk = 2 * np.pi**2 * C.primordial_power(k) / k**3
    dl = tf.delta_table(k, lmin, lmax)
    cols = []
    for l in range(lmin, lmax + 1):
        y = np.stack([sph_harm_y(l, m, th, ph) for m in range(-l, l + 1)], axis=1)
        cols.append(np.sqrt(16 * np.pi**2 * pk / vol)[:, None] * (1j**l) * dl[:, l - lmin, None] * np.conj(y))
    g = np.concatenate(cols, axis=1)  # a_lm = sum_kappa g[kappa, lm] zeta_kappa / sqrt(P/V)
    cc = g.T @ g.conj()  # <a_lm a*_l'm'>
    u = np.zeros_like(cc)
    for l, sl in H.block_slices(lmin, lmax).items():
        u[sl, sl] = H._c2r(l)
    return u @ cc @ u.conj().T


def test_eq2_real_assembly_matches_complex_formula(tf):
    lmin, lmax = 2, 7
    rng = np.random.default_rng(3)
    rot = H.random_rotation(rng)
    b = C.scaled_basis("BCC", 1.1, tf.chi_star, rot)
    vol = abs(np.linalg.det(b))
    kc = 12.0 / tf.chi_star
    ref = _complex_eq2(C.dual_vectors(b, kc), vol, tf, lmin, lmax)
    assert np.max(np.abs(ref.imag)) < 1e-10 * np.max(np.abs(ref.real))
    ours = C._assemble(C.dual_vectors(b, kc, half=True), vol, tf, lmin, lmax)
    assert np.allclose(ours, ref.real, rtol=0, atol=1e-11 * np.max(np.abs(ours)))
    ells, _ = H.lm_index(lmin, lmax)
    odd = (ells[:, None] + ells[None, :]) % 2 == 1
    assert np.max(np.abs(ref.real[odd])) < 1e-11 * np.max(np.abs(ours))


def test_rotated_lattice_equals_wigner_rotation():
    rng = np.random.default_rng(11)
    rot = H.random_rotation(rng)
    cw = C.torus_covariance("FCC", 1.05, rot, method="wigner")
    cl = C.torus_covariance("FCC", 1.05, rot, method="lattice")
    assert np.max(np.abs(cw - cl)) / np.max(np.abs(cl)) < 1e-10  # [TOP] Sec. 9: 3e-15


def test_orientation_average_is_torus_cl():
    c = C.torus_covariance("FCC", 1.05)
    diag = np.diag(c)
    avg = np.array([diag[sl].mean() for sl in H.block_slices(2, 30).values()])
    assert np.allclose(avg, C.torus_cl("FCC", 1.05), rtol=1e-10)


def test_large_torus_approaches_lcdm():
    d2 = C.torus_dl("cubic", 2.0, 2, 2)[0]
    lc = load_transfer().dl_isotropic(2, 2)[0]
    assert abs(d2 - 1009.0) < 5.0  # [TOP] Sec. 14: D_2 = 1009 against 1016
    assert abs(lc - 1016.0) < 1.0
    d_all = C.torus_dl("cubic", 2.0)
    assert np.max(np.abs(d_all[2:] / load_transfer().dl_isotropic(4, 30) - 1)) < 0.01


# --- KL and Pinsker -------------------------------------------------------------------------

def test_kl_gaussian_matches_eigenvalue_formula():
    rng = np.random.default_rng(0)
    a = rng.standard_normal((20, 20))
    c1 = a @ a.T + 20 * np.eye(20)
    c0 = np.diag(rng.uniform(5, 30, 20))
    lam = np.linalg.eigvals(np.linalg.solve(c0, c1)).real
    ref = 0.5 * np.sum(lam - 1 - np.log(lam))
    assert abs(C.kl_gaussian(c1, c0) - ref) < 1e-10 * ref
    assert abs(C.kl_gaussian(c1, np.diag(c0)) - ref) < 1e-10 * ref
    assert C.kl_gaussian(c0, c0) == pytest.approx(0.0, abs=1e-12)


def test_pinsker_bound_paper_row():
    bound = C.pinsker_bound(C.PAPER_KL["cubic"])
    assert np.allclose(bound, C.PAPER_KL["pinsker"]["cubic"], atol=0.006)


@pytest.mark.parametrize("shape,d,paper", [("cubic", 1.3, 0.12), ("FCC", 1.05, 10.15), ("BCC", 1.1, 2.14)])
def test_kl_spot_values(shape, d, paper):
    assert C.torus_kl(shape, d) == pytest.approx(paper, rel=0.02, abs=0.005)


@pytest.mark.slow
def test_kl_theorem_p_table():
    """[TOP] Theorem P table: every entry within 1.5 % or the rounding of the printed digits.

    Largest deviation: FCC 1.2 D, 0.672 against 0.68 (raw CAMB transfer gives 0.675).
    """
    for shape in C.SHAPES:
        for d, paper in zip(C.PAPER_KL["d_over_D"], C.PAPER_KL[shape]):
            assert C.torus_kl(shape, d) == pytest.approx(paper, rel=0.015, abs=0.005), (shape, d)


@pytest.mark.slow
def test_sec8_table():
    tab = C.sec8_table()
    for row in tab["rows"]:
        assert row["D2"] == pytest.approx(row["paper"]["D2"], abs=1.5), row["shape"]
        assert row["kl"] == pytest.approx(row["paper"]["kl"], rel=0.01, abs=0.01), row["shape"]
    assert tab["D2_at_2D"]["cubic"] == pytest.approx(1009.0, abs=1.0)


# --- Woodbury / determinant lemma ------------------------------------------------------------

def _random_problem(n=957, r=12, seed=1):
    rng = np.random.default_rng(seed)
    c0 = C.torus_covariance("cubic", 1.3)[:n, :n]
    u = rng.standard_normal((n, r)) * 30.0
    w = rng.uniform(0.0, 2.0, r)
    return rng, c0, u, w


def test_woodbury_matches_dense():
    rng, c0, _, _ = _random_problem()
    base = Rz.Factorized.of(c0)
    c0_inv = base.inverse()
    ref_diag = np.diag(c0) * 1.1
    for rank in (6, 8, 12):
        u = rng.standard_normal((c0.shape[0], rank)) * 30.0
        w = rng.uniform(0.0, 2.0, rank)
        c = c0 + (u * w) @ u.T
        inv, logdet = Rz.woodbury_inverse_logdet(c0_inv, base.logdet, u, w)
        sign, ref_logdet = np.linalg.slogdet(c)
        assert sign > 0 and abs(logdet - ref_logdet) < 1e-10 * abs(ref_logdet)
        if rank == 12:
            assert np.max(np.abs(inv @ c - np.eye(c.shape[0]))) < 1e-10
        x = np.linalg.cholesky(c) @ rng.standard_normal(c.shape[0])
        ll = Rz.woodbury_loglike(x, base, u, w)
        assert abs(ll - Rz.gaussian_loglike(x, c)) < 1e-10 * abs(ll)
        kl = Rz.kl_lowrank(base, u, w, ref_diag, c0=c0)
        assert abs(kl - C.kl_gaussian(c, ref_diag)) < 1e-10 * kl


def test_woodbury_singular_and_matrix_weights():
    rng, c0, u, w = _random_problem(n=200, r=6)
    w[2] = 0.0
    base = Rz.Factorized.of(c0)
    wm = np.diag(w)
    wm[0, 1] = wm[1, 0] = 0.3 * np.sqrt(w[0] * w[1])
    for weights in (w, wm):
        c = c0 + u @ (np.diag(weights) if np.ndim(weights) == 1 else weights) @ u.T
        inv, logdet = Rz.woodbury_inverse_logdet(base.inverse(), base.logdet, u, weights)
        assert abs(logdet - np.linalg.slogdet(c)[1]) < 1e-10 * abs(logdet)
        assert np.allclose(inv @ c, np.eye(200), atol=1e-9)


# --- resonance ---------------------------------------------------------------------------------

def test_shell_counts_and_preregistered_grid():
    assert Rz.shell_vectors("cubic", 1.2, 9).shape[0] == 30
    assert Rz.shell_vectors("cubic", 1.2, 16).shape[0] == 6
    grid = Rz.preregistered_grid()
    assert len(grid) == Rz.PAPER_SEARCH["grid"]["n_shells"] == 80
    kchi = {(g["shape"], g["d_over_D"], g["n_vectors"], round(g["k_chi"], 2)) for g in grid}
    assert ("BCC", 1.05, 6, 5.18) in kchi  # [TOP] Sec. 11.6 best model
    for row in Rz.PAPER_SINGLE_FIELD["rows"]:
        assert (row["torus"], row["d_over_D"], row["n_vectors"], row["k0_chi"]) in kchi


def test_line_addition_theorem_and_rank(tf):
    u, w = Rz.line_lowrank("cubic", 1.2, 9, 0.5)
    dc = (u * w) @ u.T
    assert np.linalg.matrix_rank(dc, tol=1e-9 * np.max(np.abs(dc))) <= 30
    iso = Rz.line_isotropic_diag("cubic", 1.2, 9, 0.5)
    diag = np.diag(dc)
    for l, sl in H.block_slices(2, 30).items():
        assert diag[sl].mean() == pytest.approx(iso[sl][0], rel=1e-10)
    k0 = 2 * np.pi * 3 / (1.2 * 2 * tf.chi_star)
    bump = Rz.isotropic_line_cl(k0, Rz.line_amplitude(k0, 0.5)) / C.lcdm_cl_iso()
    assert np.max(bump) == pytest.approx(0.5, rel=1e-12)
    ells, _ = H.lm_index(2, 30)
    odd = (ells[:, None] + ells[None, :]) % 2 == 1
    assert np.max(np.abs(dc[odd])) < 1e-12 * np.max(np.abs(dc))


def test_line_rotates_with_lattice():
    rng = np.random.default_rng(4)
    rot = H.random_rotation(rng)
    u0, w = Rz.line_lowrank("cubic", 1.3, 16, 0.5)
    u1, _ = Rz.line_lowrank("cubic", 1.3, 16, 0.5, R=rot)
    d = C.block_diag_rotation(rot)
    assert np.allclose((u1 * w) @ u1.T, d @ ((u0 * w) @ u0.T) @ d.T, atol=1e-10 * np.max(np.abs(u0)) ** 2 * w[0])


def test_resonance_kl_matches_dense():
    u, w = Rz.line_lowrank("cubic", 1.3, 16, 0.5)
    c1 = C.torus_covariance("cubic", 1.3) + (u * w) @ u.T
    c0 = C.lcdm_covariance() + Rz.line_isotropic_diag("cubic", 1.3, 16, 0.5)
    assert Rz.resonance_kl("cubic", 1.3, 16, 0.5) == pytest.approx(C.kl_gaussian(c1, c0), rel=1e-10)


def test_resonance_row_1p3():
    row = Rz.PAPER_RESONANCE["rows"][1.3]
    got = [Rz.resonance_kl("cubic", 1.3, None, 0.0), Rz.resonance_kl("cubic", 1.3, 9, 0.25),
           Rz.resonance_kl("cubic", 1.3, 9, 0.5), Rz.resonance_kl("cubic", 1.3, 16, 0.25),
           Rz.resonance_kl("cubic", 1.3, 16, 0.5)]
    assert np.allclose(got, row[:5], rtol=0.01, atol=0.005)
    assert got[4] / got[0] == pytest.approx(row[5], rel=0.03)  # paper prints the gain floored: 64x


@pytest.mark.slow
def test_resonance_table():
    for r in Rz.resonance_table():
        assert r["kl"] == pytest.approx(r["paper"]["kl"], rel=0.01, abs=0.005)
        assert r["no_line"] == pytest.approx(r["paper"]["no_line"], abs=0.005)
        assert abs(r["gain"] - r["paper"]["gain"]) < 1.0


def test_woodbury_benchmark_small():
    out = Rz.woodbury_benchmark(ranks=(6, 12), n_orient=2)
    assert out["n"] == 957
    assert max(out["max_rel_diff"]) < 1e-10
