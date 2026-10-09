"""Torus covariance of the low-l CMB multipoles, KL divergence and Pinsker bound, [TOP] Secs. 8 and 10.

[TOP] Eq. (2): on the 3-torus R^3/Lambda the allowed wavevectors are the dual lattice Lambda*, and

    <a_lm a*_l'm'> = (16 pi^2 / V) i^{l-l'} sum_{kappa in Lambda* \\ 0} P(kappa) Delta_l(kappa) Delta_l'(kappa)
                     Y*_lm(kappa_hat) Y_l'm'(kappa_hat),        P(k) = 2 pi^2 P_R(k) / k^3.

In the real basis of :mod:`qchas.cmb.harmonics` this is real, and pairing kappa with -kappa
kills every block with l + l' odd. With i^l = r_l (l even) or i r_l (l odd), r_l = +-1,

    C_even = G_e^T G_e,   C_odd = G_o^T G_o,   G[kappa, lm] = sqrt(2 * 16 pi^2 P(kappa) / V) r_l Delta_l Y^r_lm,

summed over one representative of each +-kappa pair, in chunks with BLAS. As in [TOP] Sec. 14 the
sum is exact for |kappa| chi_* <= 70, capped at 24 and at least 12 fundamental shells
(|kappa| <= 24 k_1, >= 12 k_1 with k_1 the shortest dual vector), and the isotropic continuum
4 pi int_{k > k_cut} dln k P_R Delta_l^2 is added on the diagonal. Lattice sizes are d/D with d
the shortest lattice vector and D = 2 chi_*.

The KL divergence of two zero-mean Gaussians, [TOP] Theorem P, is
KL(P_T || P_inf) = 1/2 [tr(C_inf^{-1} C_T) - n - ln det(C_inf^{-1} C_T)], and every test has
power minus size at most min(1, sqrt(KL/2)).
"""

from __future__ import annotations

import hashlib
import itertools
import os
from functools import lru_cache

import numpy as np
from scipy.linalg import cho_factor, solve_triangular
from scipy.optimize import brentq

from . import harmonics as H
from .camb_model import primordial_power
from .transfer import CACHE_DIR, TransferFunctions, cache_key as transfer_key, load_npz, load_transfer, save_npz_atomic

KCUT_CHI = 70.0
SHELLS_MIN = 12.0
SHELLS_MAX = 24.0
LMIN, LMAX = 2, 30
SHAPES = ("cubic", "BCC", "FCC")
COV_CACHE_VERSION = 2

#: [TOP] Theorem P table (page 28): KL(P_T || P_inf) for the l <= 30 temperature sky
PAPER_KL = {
    "d_over_D": [1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.5],
    "cubic": [19.09, 5.05, 1.48, 0.64, 0.33, 0.12, 0.03],
    "BCC": [26.01, 7.03, 2.14, 0.95, 0.50, 0.17, 0.03],
    "FCC": [38.15, 10.15, 3.01, 1.30, 0.68, 0.23, 0.04],
    "pinsker": {"cubic": [1.00, 1.00, 0.86, 0.56, 0.41, 0.25, 0.12],
                "BCC": [1.00, 1.00, 1.00, 0.69, 0.50, 0.29, 0.13],
                "FCC": [1.00, 1.00, 1.00, 0.81, 0.58, 0.34, 0.15]},
    "source": "[TOP] Theorem P, Sec. 10",
}
#: [TOP] Sec. 8 table (page 21): D_2 (uK^2) and KL (l <= 30); LCDM D_2 = 1016
PAPER_SEC8 = {
    "lcdm_D2": 1016.0,
    "rows": [
        {"shape": "cubic", "d_over_D": 0.985, "V_over_D3": 0.956, "D2": 870, "p_C2": 0.078, "p_S12": 0.105, "S12_median": 26749, "kl": 19.3},
        {"shape": "BCC", "d_over_D": 0.985, "V_over_D3": 0.736, "D2": 813, "p_C2": 0.078, "p_S12": 0.111, "S12_median": 26310, "kl": 26.4},
        {"shape": "FCC", "d_over_D": 0.985, "V_over_D3": 0.676, "D2": 778, "p_C2": 0.084, "p_S12": 0.115, "S12_median": 25329, "kl": 38.5},
        {"shape": "FCC", "d_over_D": 1.10, "V_over_D3": 0.941, "D2": 888, "p_C2": 0.063, "p_S12": 0.099, "S12_median": 28679, "kl": 3.0},
        {"shape": "FCC", "d_over_D": 1.25, "V_over_D3": 1.381, "D2": 946, "p_C2": 0.055, "p_S12": 0.088, "S12_median": 30910, "kl": 0.38},
        {"shape": "FCC", "d_over_D": 0.70, "V_over_D3": 0.243, "D2": 407, "p_C2": 0.271, "p_S12": 0.376, "S12_median": 8897, "kl": 63.8},
        {"shape": "cubic", "d_over_D": 0.50, "V_over_D3": 0.125, "D2": 510, "p_C2": 0.204, "p_S12": 0.335, "S12_median": 10311, "kl": 120.8},
    ],
    "lcdm": {"D2": 1016, "p_C2": 0.047, "p_S12": 0.076, "S12_median": 34576},
    "d2_at_2D": {"torus": 1009, "lcdm": 1016, "source": "[TOP] Sec. 14"},
    "kl_above_1_up_to": {"cubic": 1.13, "BCC": 1.16, "FCC": 1.18, "source": "[TOP] Sec. 8"},
    "source": "[TOP] Sec. 8",
}


# --- lattices ------------------------------------------------------------------------

def lattice_basis(shape) -> np.ndarray:
    """Basis rows of a named lattice (cubic, BCC, FCC) or a given 3x3 basis, shortest vector 1.

    Named lattices have their cube axes along x, y, z (the "lattice frame" of [TOP] Sec. 7):
    FCC volume 1/sqrt(2), BCC volume 4/(3 sqrt(3)) in units of d^3.
    """
    if isinstance(shape, str):
        key = shape.upper()
        if key in ("CUBIC", "SC"):
            b = np.eye(3)
        elif key == "BCC":
            b = np.array([[-1.0, 1, 1], [1, -1, 1], [1, 1, -1]]) / np.sqrt(3.0)
        elif key == "FCC":
            b = np.array([[0.0, 1, 1], [1, 0, 1], [1, 1, 0]]) / np.sqrt(2.0)
        else:
            raise ValueError(f"unknown lattice shape {shape!r}")
        return b
    b = np.asarray(shape, dtype=float).reshape(3, 3)
    return b / shortest_vector_length(b)


def _check_basis(basis) -> np.ndarray:
    b = np.asarray(basis, dtype=float).reshape(3, 3)
    norms = np.linalg.norm(b, axis=1)
    if not np.all(np.isfinite(b)) or abs(np.linalg.det(b)) <= 1e-10 * np.prod(norms):
        raise ValueError("degenerate lattice basis")
    return b


def shortest_vector_length(basis: np.ndarray, max_points: int = 10_000_000) -> float:
    """Length of the shortest non-zero lattice vector, exact for any (also unreduced) basis.

    A vector v = n B shorter than the shortest basis row r0 has |n_i| <= r0 |B^{-1} e_i|.
    """
    b = _check_basis(basis)
    r0 = float(np.min(np.linalg.norm(b, axis=1)))
    nb = np.ceil(r0 * np.linalg.norm(np.linalg.inv(b), axis=0) + 1e-9).astype(int)
    if np.prod(2 * nb + 1) > max_points:
        raise ValueError("lattice basis too skewed; reduce it (e.g. LLL) first")
    n = np.array(list(itertools.product(*(np.arange(-m, m + 1) for m in nb))), dtype=float)
    n = n[np.any(n != 0, axis=1)]
    return float(np.min(np.linalg.norm(n @ b, axis=1)))


def shape_name(shape) -> str:
    if isinstance(shape, str):
        return {"SC": "cubic", "CUBIC": "cubic"}.get(shape.upper(), shape.upper())
    return "basis-" + hashlib.sha256(np.round(lattice_basis(shape), 12).tobytes()).hexdigest()[:10]


def volume_over_D3(shape, d_over_D: float) -> float:
    """V / D^3 of the torus (paper Sec. 8 column V/D^3)."""
    return float(abs(np.linalg.det(lattice_basis(shape))) * d_over_D**3)


def scaled_basis(shape, d_over_D: float, chi_star: float | None = None, R: np.ndarray | None = None) -> np.ndarray:
    """Basis rows in Mpc with shortest vector d = (d/D) 2 chi_*, rotated by R (rows b -> R b)."""
    if not (np.isfinite(d_over_D) and d_over_D > 0):
        raise ValueError(f"d/D must be positive, got {d_over_D}")
    chi = load_transfer().chi_star if chi_star is None else chi_star
    b = lattice_basis(shape) * d_over_D * 2.0 * chi
    return b if R is None else b @ np.asarray(R, dtype=float).T


def dual_basis(basis: np.ndarray) -> np.ndarray:
    """Rows of the dual lattice Lambda* = 2 pi inv(B)^T (so B K^T = 2 pi I)."""
    return 2.0 * np.pi * np.linalg.inv(np.asarray(basis, dtype=float)).T


def dual_vectors(basis: np.ndarray, kmax: float, half: bool = False) -> np.ndarray:
    """All non-zero dual-lattice vectors with |kappa| <= kmax (one of each +-pair if ``half``)."""
    K = dual_basis(basis)
    kinv = np.linalg.inv(K)
    nb = np.ceil(kmax * np.linalg.norm(kinv, axis=0) + 1e-9).astype(int)
    out = []
    for n1 in range(0 if half else -nb[0], nb[0] + 1):
        g2, g3 = np.meshgrid(np.arange(-nb[1], nb[1] + 1), np.arange(-nb[2], nb[2] + 1), indexing="ij")
        n = np.stack([np.full(g2.size, n1), g2.ravel(), g3.ravel()], axis=1)
        if half and n1 == 0:
            n = n[(n[:, 1] > 0) | ((n[:, 1] == 0) & (n[:, 2] > 0))]
        v = n @ K
        r2 = np.einsum("ij,ij->i", v, v)
        keep = (r2 <= kmax * kmax * (1 + 1e-12)) & np.any(n != 0, axis=1)
        out.append(v[keep])
    return np.concatenate(out) if out else np.zeros((0, 3))


def fundamental_k(basis: np.ndarray) -> float:
    """|kappa| of the shortest dual vector, Mpc^-1."""
    return shortest_vector_length(dual_basis(basis))


def k_cut(basis: np.ndarray, chi_star: float, kcut_chi: float = KCUT_CHI,
          shells: tuple[float, float] = (SHELLS_MIN, SHELLS_MAX)) -> float:
    """Cutoff of the exact sum: kcut_chi / chi_*, clipped to [shells[0], shells[1]] x k_1 ([TOP] Sec. 14)."""
    k1 = fundamental_k(basis)
    return float(np.clip(kcut_chi / chi_star, shells[0] * k1, shells[1] * k1))


# --- fast real harmonics (same convention as qchas.cmb.harmonics) ----------------------

def real_ylm_unit(lmin: int, lmax: int, n: np.ndarray) -> np.ndarray:
    """Real Y^r_lm at unit vectors n (rows), shape (len(n), n_modes); harmonics.py convention.

    Normalized associated Legendre recursion without the Condon-Shortley phase, since
    Y^r_{l,+-m} = sqrt(2) N_lm Pbar_l^m(cos theta) (cos m phi, sin m phi).
    """
    n = np.asarray(n, dtype=float)
    npt = n.shape[0]
    z = np.clip(n[:, 2], -1.0, 1.0)
    s = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    phi = np.arctan2(n[:, 1], n[:, 0])
    out = np.empty((npt, H.n_modes(lmin, lmax)))
    offs = {l: sl.start for l, sl in H.block_slices(lmin, lmax).items()}
    pmm = np.full(npt, np.sqrt(1.0 / (4.0 * np.pi)))
    for m in range(0, lmax + 1):
        if m > 0:
            pmm = pmm * s * np.sqrt((2.0 * m + 1.0) / (2.0 * m))
        cm = np.sqrt(2.0) * np.cos(m * phi) if m else None
        sm = np.sqrt(2.0) * np.sin(m * phi) if m else None
        p_prev, p_cur = None, pmm
        for l in range(m, lmax + 1):
            if l == m + 1:
                p_prev, p_cur = p_cur, np.sqrt(2.0 * m + 3.0) * z * p_cur
            elif l > m + 1:
                a = np.sqrt((4.0 * l * l - 1.0) / (l * l - m * m))
                b = np.sqrt(((l - 1.0) ** 2 - m * m) / (4.0 * (l - 1.0) ** 2 - 1.0))
                p_prev, p_cur = p_cur, a * (z * p_cur - b * p_prev)
            if l >= lmin:
                o = offs[l] + l
                if m == 0:
                    out[:, o] = p_cur
                else:
                    out[:, o + m] = p_cur * cm
                    out[:, o - m] = p_cur * sm
    return out


def ell_phase(lmin: int, lmax: int) -> np.ndarray:
    """r_l per mode: i^l = r_l (l even) or i r_l (l odd), r_l = +-1."""
    ells, _ = H.lm_index(lmin, lmax)
    return np.where(ells % 2 == 0, (-1.0) ** (ells // 2), (-1.0) ** ((ells - 1) // 2))


def mode_vectors(kvecs: np.ndarray, tf: TransferFunctions, lmin: int = LMIN, lmax: int = LMAX,
                 k_norm: float | None = None) -> np.ndarray:
    """Rows v_kappa[lm] = r_l Delta_l(|kappa|) Y^r_lm(kappa_hat) (uK); ``k_norm`` overrides |kappa|."""
    kv = np.asarray(kvecs, dtype=float)
    norms = np.linalg.norm(kv, axis=1)
    y = real_ylm_unit(lmin, lmax, kv / norms[:, None])
    ells, _ = H.lm_index(lmin, lmax)
    if k_norm is None:
        dt = tf.delta_table(norms, lmin, lmax)
    else:
        dt = np.broadcast_to(tf.delta_table(np.array([float(k_norm)]), lmin, lmax), (norms.size, lmax - lmin + 1))
    return y * dt[:, ells - lmin] * ell_phase(lmin, lmax)[None, :]


def block_diag_rotation(R: np.ndarray, lmin: int = LMIN, lmax: int = LMAX) -> np.ndarray:
    """Dense real Wigner matrix D(R) = blockdiag(R_l) with Y^r(R n) = D Y^r(n)."""
    n = H.n_modes(lmin, lmax)
    d = np.zeros((n, n))
    for l, sl in H.block_slices(lmin, lmax).items():
        d[sl, sl] = H.real_rotation_block(l, R)
    return d


def rotate_covariance(C: np.ndarray, R: np.ndarray, lmin: int = LMIN, lmax: int = LMAX) -> np.ndarray:
    """C(R) = D(R) C D(R)^T: covariance of the lattice rotated actively by R."""
    sl = H.block_slices(lmin, lmax)
    blocks = {l: H.real_rotation_block(l, R) for l in sl}
    out = np.empty_like(C)
    for l, si in sl.items():
        left = blocks[l] @ C[si, :]
        for lp, sj in sl.items():
            out[si, sj] = left[:, sj] @ blocks[lp].T
    return out


# --- torus covariance -----------------------------------------------------------------

def _parity_split(lmin: int, lmax: int) -> tuple[np.ndarray, np.ndarray]:
    ells, _ = H.lm_index(lmin, lmax)
    return np.flatnonzero(ells % 2 == 0), np.flatnonzero(ells % 2 == 1)


def _assemble(kvecs_half: np.ndarray, volume: float, tf: TransferFunctions, lmin: int, lmax: int,
              power=None, chunk: int = 4096) -> np.ndarray:
    """Exact Eq. (2) over the given half-lattice (each +-pair once), real basis."""
    pr_fn = primordial_power if power is None else power
    ie, io = _parity_split(lmin, lmax)
    ce = np.zeros((ie.size, ie.size))
    co = np.zeros((io.size, io.size))
    for s in range(0, kvecs_half.shape[0], chunk):
        kv = kvecs_half[s:s + chunk]
        k = np.linalg.norm(kv, axis=1)
        pk = 2.0 * np.pi**2 * pr_fn(k) / k**3
        amp = np.sqrt(2.0 * 16.0 * np.pi**2 * pk / volume)
        g = mode_vectors(kv, tf, lmin, lmax) * amp[:, None]
        ge, go = g[:, ie], g[:, io]
        ce += ge.T @ ge
        co += go.T @ go
    n = H.n_modes(lmin, lmax)
    c = np.zeros((n, n))
    c[np.ix_(ie, ie)] = ce
    c[np.ix_(io, io)] = co
    return c


def continuum_cl(kcut: float, lmin: int = LMIN, lmax: int = LMAX, tf: TransferFunctions | None = None) -> np.ndarray:
    """Isotropic continuum 4 pi int_{k > kcut} dln k P_R Delta_l^2 (uK^2) added beyond the exact sum."""
    tf = load_transfer() if tf is None else tf
    return tf.cl_isotropic(lmin, lmax, kmin=kcut)


def _cov_cache_path(shape, d_over_D: float, lmin: int, lmax: int, kcut_chi: float, shells) -> os.PathLike:
    tag = (f"{shape_name(shape)}_{d_over_D:.6f}_{lmin}_{lmax}_{kcut_chi:g}_{shells[0]:g}-{shells[1]:g}"
           f"_v{COV_CACHE_VERSION}_{transfer_key()}")
    return CACHE_DIR / "covariance" / f"cov_{tag}.npz"


def _lattice_frame_covariance(shape, d_over_D: float, lmin: int, lmax: int, kcut_chi: float,
                              shells: tuple[float, float], cache: bool) -> np.ndarray:
    path = _cov_cache_path(shape, d_over_D, lmin, lmax, kcut_chi, shells)
    ie, io = _parity_split(lmin, lmax)
    z = load_npz(path) if cache else None
    if z is not None and z.get("even", np.zeros(0)).shape == (ie.size,) * 2 and z.get("odd", np.zeros(0)).shape == (io.size,) * 2:
        c = np.zeros((H.n_modes(lmin, lmax),) * 2)
        c[np.ix_(ie, ie)] = z["even"]
        c[np.ix_(io, io)] = z["odd"]
        return c
    c = _direct_covariance(shape, d_over_D, None, lmin, lmax, kcut_chi, shells)
    if cache:
        save_npz_atomic(path, even=c[np.ix_(ie, ie)], odd=c[np.ix_(io, io)])
    return c


def _direct_covariance(shape, d_over_D: float, R, lmin: int, lmax: int, kcut_chi: float,
                       shells: tuple[float, float]) -> np.ndarray:
    tf = load_transfer()
    b = scaled_basis(shape, d_over_D, tf.chi_star, R)
    vol = abs(np.linalg.det(b))
    kc = k_cut(b, tf.chi_star, kcut_chi, shells)
    c = _assemble(dual_vectors(b, kc, half=True), vol, tf, lmin, lmax)
    ells, _ = H.lm_index(lmin, lmax)
    c[np.diag_indices_from(c)] += continuum_cl(kc, lmin, lmax, tf)[ells - lmin]
    return c


def torus_covariance(shape, d_over_D: float, R: np.ndarray | None = None, lmin: int = LMIN, lmax: int = LMAX, *,
                     kcut_chi: float = KCUT_CHI, shells: tuple[float, float] = (SHELLS_MIN, SHELLS_MAX),
                     method: str = "wigner", cache: bool = True) -> np.ndarray:
    """[TOP] Eq. (2) covariance of the real a_lm vector (uK^2), lattice rotated actively by R.

    ``method="wigner"`` rotates the (cached) lattice-frame covariance with real Wigner blocks;
    ``method="lattice"`` sums directly over the rotated dual lattice.
    """
    if method not in ("wigner", "lattice"):
        raise ValueError(f"method must be 'wigner' or 'lattice', got {method!r}")
    shells = (float(shells[0]), float(shells[1]))
    if R is not None and method == "lattice":
        return _direct_covariance(shape, d_over_D, R, lmin, lmax, kcut_chi, shells)
    c0 = _lattice_frame_covariance(shape, d_over_D, lmin, lmax, kcut_chi, shells, cache)
    return c0 if R is None else rotate_covariance(c0, R, lmin, lmax)


def torus_cl(shape, d_over_D: float, lmin: int = LMIN, lmax: int = LMAX, *, kcut_chi: float = KCUT_CHI,
             shells: tuple[float, float] = (SHELLS_MIN, SHELLS_MAX)) -> np.ndarray:
    """Orientation-averaged C_l = (4 pi / V) sum_kappa P(kappa) Delta_l(kappa)^2 + continuum (uK^2)."""
    tf = load_transfer()
    b = scaled_basis(shape, d_over_D, tf.chi_star)
    vol = abs(np.linalg.det(b))
    kc = k_cut(b, tf.chi_star, kcut_chi, shells)
    k = np.linalg.norm(dual_vectors(b, kc, half=True), axis=1)
    pk = 2.0 * np.pi**2 * primordial_power(k) / k**3
    d = tf.delta_table(k, lmin, lmax)
    return 2.0 * 4.0 * np.pi / vol * np.einsum("k,kl->l", pk, d**2) + continuum_cl(kc, lmin, lmax, tf)


def torus_dl(shape, d_over_D: float, lmin: int = LMIN, lmax: int = LMAX, **kw) -> np.ndarray:
    ell = np.arange(lmin, lmax + 1)
    return ell * (ell + 1) * torus_cl(shape, d_over_D, lmin, lmax, **kw) / (2.0 * np.pi)


def lcdm_cl_iso(lmin: int = LMIN, lmax: int = LMAX) -> np.ndarray:
    """Infinite-space C_l from the same transfer functions (isotropic integral), uK^2."""
    return load_transfer().cl_isotropic(lmin, lmax)


def lcdm_covariance(lmin: int = LMIN, lmax: int = LMAX) -> np.ndarray:
    """Diagonal of the infinite-space LCDM covariance of the real a_lm vector (uK^2)."""
    ells, _ = H.lm_index(lmin, lmax)
    return lcdm_cl_iso(lmin, lmax)[ells - lmin]


# --- divergences ------------------------------------------------------------------------

def _whiten(c1: np.ndarray, c0: np.ndarray) -> np.ndarray:
    """M = L0^{-1} C1 L0^{-T} with C0 = L0 L0^T (C0 may be a diagonal given as 1-d array)."""
    if c0.ndim == 1:
        s = 1.0 / np.sqrt(c0)
        return c1 * s[:, None] * s[None, :]
    l0 = np.linalg.cholesky(c0)
    tmp = solve_triangular(l0, c1, lower=True)
    return solve_triangular(l0, tmp.T, lower=True).T


def kl_gaussian(c1: np.ndarray, c0: np.ndarray) -> float:
    """KL(N(0, C1) || N(0, C0)) = 1/2 [tr(C0^{-1} C1) - n - ln det(C0^{-1} C1)] via Cholesky."""
    c1 = np.asarray(c1, dtype=float)
    c0 = np.asarray(c0, dtype=float)
    m = _whiten(c1, c0)
    m = 0.5 * (m + m.T)
    lm, _ = cho_factor(m, lower=True)
    logdet = 2.0 * np.sum(np.log(np.diag(lm)))
    return float(0.5 * (np.trace(m) - m.shape[0] - logdet))


def pinsker_bound(kl) -> np.ndarray:
    """[TOP] Theorem P: power minus size <= min(1, sqrt(KL/2))."""
    return np.minimum(1.0, np.sqrt(np.asarray(kl, dtype=float) / 2.0))


@lru_cache(maxsize=256)
def torus_kl(shape: str, d_over_D: float, lmin: int = LMIN, lmax: int = LMAX, kcut_chi: float = KCUT_CHI,
             shells: tuple[float, float] = (SHELLS_MIN, SHELLS_MAX)) -> float:
    """KL(torus || infinite LCDM) for the l <= lmax temperature sky (any orientation)."""
    c = torus_covariance(shape, d_over_D, None, lmin, lmax, kcut_chi=kcut_chi, shells=shells)
    return kl_gaussian(c, lcdm_covariance(lmin, lmax))


def kl_table(shapes=SHAPES, d_over_D=tuple(PAPER_KL["d_over_D"]), lmax: int = LMAX) -> dict:
    """KL and Pinsker bound per shape over a size grid, with the paper's Theorem P values (``cmb.json -> kl``)."""
    out = {"d_over_D": list(d_over_D), "paper": {}, "pinsker_bound": {}, "computed": True, "source": PAPER_KL["source"]}
    for s in shapes:
        kl = [torus_kl(s, float(x), LMIN, lmax) for x in d_over_D]
        out[s] = kl
        out["pinsker_bound"][s] = pinsker_bound(kl).tolist()
        if s in PAPER_KL:
            out["paper"][s] = {str(float(x)): v for x, v in zip(PAPER_KL["d_over_D"], PAPER_KL[s])}
            out["pinsker_bound_paper"] = PAPER_KL["pinsker"]
    return out


def sec8_table(lmax: int = LMAX) -> dict:
    """[TOP] Sec. 8 models: orientation-averaged D_l, D_2 and KL, computed next to the paper's values."""
    ell = np.arange(LMIN, lmax + 1)
    lcdm_d = ell * (ell + 1) * lcdm_cl_iso(LMIN, lmax) / (2.0 * np.pi)
    rows = []
    for r in PAPER_SEC8["rows"]:
        dl = torus_dl(r["shape"], r["d_over_D"], LMIN, lmax)
        rows.append({"shape": r["shape"], "d_over_D": r["d_over_D"], "V_over_D3": volume_over_D3(r["shape"], r["d_over_D"]),
                     "D2": float(dl[0]), "kl": torus_kl(r["shape"], r["d_over_D"], LMIN, lmax), "D": dl.tolist(),
                     "computed": True, "paper": {k: r[k] for k in ("V_over_D3", "D2", "kl")}})
    d2_2d = {s: float(torus_dl(s, 2.0, LMIN, LMIN)[0]) for s in SHAPES}
    return {"ell": ell.tolist(), "lcdm_D": lcdm_d.tolist(), "lcdm_D2_paper": PAPER_SEC8["lcdm_D2"], "rows": rows,
            "D2_at_2D": d2_2d, "D2_at_2D_paper": PAPER_SEC8["d2_at_2D"], "computed": True, "source": PAPER_SEC8["source"]}


def _contract_size(d_over_D: float) -> str:
    """Size label of the ``cmb.json -> torus_D`` keys: at least two decimals (0.70, 0.985)."""
    s = f"{d_over_D:.3f}".rstrip("0")
    return s if len(s.split(".")[1]) >= 2 else f"{d_over_D:.2f}"


def torus_D_bundle(models=(("FCC", 0.985), ("cubic", 0.985), ("BCC", 0.985), ("FCC", 0.70)),
                   lmax: int = LMAX) -> dict:
    """Orientation-averaged D_l (uK^2, l = 2..lmax) keyed as ``cmb.json -> torus_D`` ("FCC_0.985", ...)."""
    return {f"{shape_name(s)}_{_contract_size(d)}": torus_dl(s, d, LMIN, lmax).tolist() for s, d in models}


def kl_unit_crossing(shape: str, bracket: tuple[float, float] = (1.05, 1.3), lmax: int = LMAX,
                     xtol: float = 5e-4) -> float:
    """Size d/D at which KL(l <= lmax) = 1 (Brent root of ln KL; covariances not cached).

    [TOP] Sec. 8 quotes d ~ 1.13 D (cubic), 1.16 D (BCC), 1.18 D (FCC).
    """
    c_inf = lcdm_covariance(LMIN, lmax)

    def log_kl(d: float) -> float:
        return float(np.log(kl_gaussian(torus_covariance(shape, d, None, LMIN, lmax, cache=False), c_inf)))

    return float(brentq(log_kl, *bracket, xtol=xtol))
