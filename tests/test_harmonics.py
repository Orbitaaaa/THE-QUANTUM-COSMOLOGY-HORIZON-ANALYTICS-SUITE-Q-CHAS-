import healpy as hp
import numpy as np

from qchas.cmb import harmonics as H


def test_real_ylm_orthonormal():
    nside = 32
    th, ph = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    y = H.real_ylm(0, 8, th, ph)
    gram = y.T @ y * (4 * np.pi / th.size)
    assert np.allclose(gram, np.eye(gram.shape[0]), atol=2e-3)


def test_healpy_roundtrip_and_synthesis():
    rng = np.random.default_rng(1)
    lmin, lmax, nside = 2, 12, 32
    vec = rng.standard_normal(H.n_modes(lmin, lmax))
    alm = H.real_to_healpy(vec, lmin, lmax)
    assert np.allclose(H.healpy_to_real(alm, lmin, lmax), vec)
    m = hp.alm2map(alm, nside, lmax=lmax)
    th, ph = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    direct = H.real_ylm(lmin, lmax, th, ph) @ vec
    assert np.allclose(m, direct, atol=1e-10)


def test_rotation_law_and_group():
    rng = np.random.default_rng(2)
    lmin, lmax = 2, 10
    rot = H.random_rotation(rng)
    n = rng.standard_normal((50, 3))
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    y = H.real_ylm_vec(lmin, lmax, n)
    y_rot = H.real_ylm_vec(lmin, lmax, n @ rot)  # rows R^{-1} n = R^T n
    for l, sl in H.block_slices(lmin, lmax).items():
        r = H.real_rotation_block(l, rot)
        assert np.allclose(r @ r.T, np.eye(2 * l + 1), atol=1e-12)
        assert np.allclose(y_rot[:, sl], y[:, sl] @ r, atol=1e-10)
    r2 = H.random_rotation(rng)
    for l in (3, 7):
        assert np.allclose(H.real_rotation_block(l, rot @ r2),
                           H.real_rotation_block(l, rot) @ H.real_rotation_block(l, r2), atol=1e-10)


def test_euler_roundtrip():
    rng = np.random.default_rng(3)
    for _ in range(20):
        r = H.random_rotation(rng)
        assert np.allclose(H.rotation_from_euler(*H.euler_zyz(r)), r, atol=1e-10)
