"""Real spherical harmonics, healpy <-> real a_lm conversion, and real Wigner rotations.

Real basis (Condon-Shortley complex Y_lm underneath)::

    Y^r_{l0}  = Y_{l0}
    Y^r_{lm}  = sqrt(2) (-1)^m Re Y_{lm}      (m > 0)
    Y^r_{l,-m} = sqrt(2) (-1)^m Im Y_{lm}     (m > 0)

For a real field T = sum a_lm Y_lm (healpy convention a_{l,-m} = (-1)^m conj(a_lm)) the
real coefficients are a^r_{l0} = a_{l0}, a^r_{lm} = sqrt(2)(-1)^m Re a_lm and
a^r_{l,-m} = -sqrt(2)(-1)^m Im a_lm. Within each l the real vector is ordered
m = -l..l. An analysis vector stacks l = lmin..lmax (957 entries for 2..30).

Rotations act actively: f'(n) = f(R^{-1} n) maps a^r_l -> R_l a^r_l with R_l real
orthogonal, so that Y^r_l(R^{-1} n) = R_l^T Y^r_l(n).
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.special import sph_harm_y


def lm_index(lmin: int, lmax: int) -> tuple[np.ndarray, np.ndarray]:
    """(ell, m) arrays for the stacked real vector over lmin..lmax."""
    ells = np.concatenate([np.full(2 * l + 1, l) for l in range(lmin, lmax + 1)])
    ms = np.concatenate([np.arange(-l, l + 1) for l in range(lmin, lmax + 1)])
    return ells, ms


def n_modes(lmin: int, lmax: int) -> int:
    return (lmax + 1) ** 2 - lmin**2


def block_slices(lmin: int, lmax: int) -> dict[int, slice]:
    out, off = {}, 0
    for l in range(lmin, lmax + 1):
        out[l] = slice(off, off + 2 * l + 1)
        off += 2 * l + 1
    return out


def real_ylm(lmin: int, lmax: int, theta, phi) -> np.ndarray:
    """Real spherical harmonics, shape (n_dirs, n_modes(lmin, lmax)).

    ``theta`` is the polar angle (0 at the north pole), ``phi`` the azimuth.
    """
    theta = np.atleast_1d(np.asarray(theta, dtype=float))
    phi = np.atleast_1d(np.asarray(phi, dtype=float))
    out = np.empty((theta.size, n_modes(lmin, lmax)))
    col = 0
    for l in range(lmin, lmax + 1):
        y0 = sph_harm_y(l, 0, theta, phi).real
        block = np.empty((theta.size, 2 * l + 1))
        block[:, l] = y0
        for m in range(1, l + 1):
            y = sph_harm_y(l, m, theta, phi)
            s = np.sqrt(2.0) * (-1.0) ** m
            block[:, l + m] = s * y.real
            block[:, l - m] = s * y.imag
        out[:, col:col + 2 * l + 1] = block
        col += 2 * l + 1
    return out


def real_ylm_vec(lmin: int, lmax: int, vecs: np.ndarray) -> np.ndarray:
    """Real Y evaluated at unit vectors (n, 3)."""
    v = np.asarray(vecs, dtype=float)
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    theta = np.arccos(np.clip(v[..., 2], -1.0, 1.0))
    phi = np.arctan2(v[..., 1], v[..., 0])
    return real_ylm(lmin, lmax, theta, phi)


def healpy_to_real(alm: np.ndarray, lmin: int, lmax: int, alm_lmax: int | None = None) -> np.ndarray:
    """Convert a healpy complex alm array to the stacked real vector over lmin..lmax."""
    import healpy as hp

    alm_lmax = alm_lmax if alm_lmax is not None else hp.Alm.getlmax(alm.size)
    out = np.empty(n_modes(lmin, lmax))
    col = 0
    for l in range(lmin, lmax + 1):
        out[col + l] = alm[hp.Alm.getidx(alm_lmax, l, 0)].real
        for m in range(1, l + 1):
            a = alm[hp.Alm.getidx(alm_lmax, l, m)]
            s = np.sqrt(2.0) * (-1.0) ** m
            out[col + l + m] = s * a.real
            out[col + l - m] = -s * a.imag
        col += 2 * l + 1
    return out


def real_to_healpy(vec: np.ndarray, lmin: int, lmax: int, alm_lmax: int | None = None) -> np.ndarray:
    """Inverse of :func:`healpy_to_real` (multipoles outside lmin..lmax set to zero)."""
    import healpy as hp

    alm_lmax = alm_lmax if alm_lmax is not None else lmax
    alm = np.zeros(hp.Alm.getsize(alm_lmax), dtype=complex)
    col = 0
    for l in range(lmin, lmax + 1):
        alm[hp.Alm.getidx(alm_lmax, l, 0)] = vec[col + l]
        for m in range(1, l + 1):
            s = np.sqrt(2.0) * (-1.0) ** m
            alm[hp.Alm.getidx(alm_lmax, l, m)] = (vec[col + l + m] - 1j * vec[col + l - m]) / s
        col += 2 * l + 1
    return alm


# --- Wigner rotations ----------------------------------------------------------------

@lru_cache(maxsize=None)
def _jy_eig(l: int) -> tuple[np.ndarray, np.ndarray]:
    """Eigen-decomposition of J_y in the |l m> basis (m = -l..l)."""
    m = np.arange(-l, l)  # J+ |l m> = sqrt((l-m)(l+m+1)) |l m+1>
    jp = np.diag(np.sqrt((l - m) * (l + m + 1.0)), -1).astype(complex)
    jy = (jp - jp.conj().T) / 2j
    w, v = np.linalg.eigh(jy)
    return w, v


@lru_cache(maxsize=None)
def _c2r(l: int) -> np.ndarray:
    """Unitary U with a^r = U a (complex a over m = -l..l, real convention above)."""
    u = np.zeros((2 * l + 1, 2 * l + 1), dtype=complex)
    u[l, l] = 1.0
    for m in range(1, l + 1):
        # with a_-m = (-1)^m conj(a_m): Re a_m = (a_m + (-1)^m a_-m)/2, Im a_m = (a_m - (-1)^m a_-m)/2i
        s = (-1.0) ** m / np.sqrt(2.0)
        u[l + m, l + m] = s
        u[l + m, l - m] = s * (-1.0) ** m
        u[l - m, l + m] = -s / 1j
        u[l - m, l - m] = s * (-1.0) ** m / 1j
    return u


def wigner_d_small(l: int, beta: float) -> np.ndarray:
    """d^l_{m m'}(beta) = <l m| exp(-i beta J_y) |l m'>."""
    w, v = _jy_eig(l)
    return (v * np.exp(-1j * beta * w)) @ v.conj().T


def wigner_D(l: int, alpha: float, beta: float, gamma: float) -> np.ndarray:
    """D^l(alpha, beta, gamma) = e^{-i alpha J_z} e^{-i beta J_y} e^{-i gamma J_z} (zyz, active)."""
    m = np.arange(-l, l + 1)
    return np.exp(-1j * m * alpha)[:, None] * wigner_d_small(l, beta) * np.exp(-1j * m * gamma)[None, :]


def euler_zyz(rot: np.ndarray) -> tuple[float, float, float]:
    """Euler angles with rot = Rz(alpha) Ry(beta) Rz(gamma)."""
    r = np.asarray(rot, dtype=float)
    beta = float(np.arccos(np.clip(r[2, 2], -1.0, 1.0)))
    if abs(np.sin(beta)) > 1e-12:
        alpha = float(np.arctan2(r[1, 2], r[0, 2]))
        gamma = float(np.arctan2(r[2, 1], -r[2, 0]))
    else:  # gimbal lock: only alpha + gamma (or alpha - gamma) is defined
        alpha = float(np.arctan2(r[1, 0], r[0, 0])) if r[2, 2] > 0 else float(np.arctan2(-r[1, 0], -r[0, 0]))
        gamma = 0.0
    return alpha, beta, gamma


def rotation_from_euler(alpha: float, beta: float, gamma: float) -> np.ndarray:
    def rz(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])

    c, s = np.cos(beta), np.sin(beta)
    ry = np.array([[c, 0, s], [0, 1.0, 0], [-s, 0, c]])
    return rz(alpha) @ ry @ rz(gamma)


def real_rotation_block(l: int, rot: np.ndarray) -> np.ndarray:
    """Real orthogonal R_l with a^r -> R_l a^r for the active rotation ``rot``."""
    d = wigner_D(l, *euler_zyz(rot))
    u = _c2r(l)
    r = u @ d @ u.conj().T
    return r.real


def rotate_real_vector(vec: np.ndarray, rot: np.ndarray, lmin: int, lmax: int) -> np.ndarray:
    """Apply the active rotation to a stacked real a_lm vector."""
    out = np.empty_like(vec)
    for l, sl in block_slices(lmin, lmax).items():
        out[sl] = real_rotation_block(l, rot) @ vec[sl]
    return out


def random_rotation(rng: np.random.Generator) -> np.ndarray:
    """Haar-random rotation matrix (uniform on SO(3))."""
    q, r = np.linalg.qr(rng.standard_normal((3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def gal_to_vec(l_deg, b_deg) -> np.ndarray:
    l, b = np.radians(l_deg), np.radians(b_deg)
    return np.stack([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], axis=-1)


def vec_to_gal(v) -> tuple[np.ndarray, np.ndarray]:
    v = np.asarray(v, dtype=float)
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360.0, np.degrees(np.arcsin(np.clip(v[..., 2], -1, 1)))
