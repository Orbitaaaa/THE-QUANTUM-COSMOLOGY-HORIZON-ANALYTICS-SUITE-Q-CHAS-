"""Orientation-scanned torus likelihood and the pre-registered resonance search, [TOP] Secs. 9 and 11.6.

Torus likelihood ([TOP] Sec. 9). For a lattice of shape X, size d and orientation R the real
multipoles a (2 <= l <= 30) are Gaussian with covariance C(R) = D(R) C0 D(R)^T, D the real
Wigner matrix and C0 the [TOP] Eq. (2) covariance in the lattice frame. Against infinite LCDM
(diagonal C_L)

    Delta lnL(R) = -1/2 b^T K b - 1/2 ln(det C0 / det C_L),   b = D(R)^T a,   K = C0^{-1} - C_L^{-1},

because a^T C(R)^{-1} a = b^T C0^{-1} b (the DATA are rotated into the lattice frame) and C_L is
isotropic. Per model we report max over SO(3) of Delta lnL and the orientation-averaged Bayes
factor ln B = ln < exp(Delta lnL(R)) >_Haar, the paper's optimal statistic.

SO(3) scan. With R = Rz(alpha) Ry(beta) Rz(gamma) the complex-basis coefficients of b are
c_lm = e^{i m gamma} ct_lm(alpha, beta), ct_l = d^l(beta)^T (e^{i m' alpha} a^c_lm'), so

    b^T K b = sum_{Delta = -60}^{60} q_Delta(alpha, beta) e^{-i Delta gamma},
    q_Delta = sum_{m - m' = Delta} conj(ct_lm) K^c_{lm, l'm'} ct_l'm',

an exact finite Fourier series in gamma ("two Euler angles are scanned on a grid; the third is
handled exactly by a Fourier series"). q is computed for every (alpha, beta) of the grid with
BLAS-3 products over m-sorted blocks of K^c; blocks with l + l' odd vanish (Proposition 9.1) and
the symmetries q_{-Delta} = conj(q_Delta), B(m, m') = B(-m', -m) cut the work to a quarter; one
FFT then evaluates the series on the gamma grid (exact at the grid points for any n_gamma).
alpha is uniform, cos(beta) on Gauss-Legendre nodes (Haar weight sin beta dbeta; exact for the
band-limited quadratic form itself), gamma uniform. The maximum is refined from the best grid
cells by a nested zoom in (alpha, beta) with gamma maximized exactly on the trigonometric series.

Resonance search ([TOP] Sec. 11.6, pre-registered). Each of the 80 grid shells (FCC/BCC/cubic at
d/D = 1.05, 1.15, 1.2, 1.3, at most 12 vectors, 2 <= k0 chi_* <= 14) carries a line of C_l bump
f = 0.2 or 0.4: 160 models. Statistic: orientation-averaged ln B of torus + line (C0 + U W U^T,
:func:`qchas.cmb.resonance.line_lowrank`) against infinite space + the same isotropic line. By
the Woodbury identity

    b^T (C0 + U W U^T)^{-1} b = b^T C0^{-1} b - y^T M y,   y = V^T b,  V = C0^{-1} U,
    M = (W^{-1} + U^T V)^{-1},   ln det(C0 + U W U^T) = ln det C0 + ln det(I + W U^T V),

so each torus needs one quadratic orientation scan (shared with the likelihood) plus linear
rotational correlations y_j(R); these are linear in a, hence also short Fourier series in gamma.
S_max = max over the 160 models of ln B, S_avg = ln(mean over models of B).

Calibration: isotropic LCDM skies a ~ N(0, C_L), processed with the identical code and grid;
p = fraction of null skies with statistic >= observed. Torus-true skies (random Haar
orientation) give the exclusion columns of [TOP] Table 9b.
"""

from __future__ import annotations

import ctypes
import math
import multiprocessing as mp
import os
import time
from dataclasses import dataclass, field

import numpy as np

from . import covariance as cov
from . import harmonics as H
from . import resonance as res

LMIN, LMAX = cov.LMIN, cov.LMAX
DMAX = 2 * LMAX  # largest |m - m'| in the gamma series

SHAPES = ("FCC", "BCC", "cubic")
LIKELIHOOD_SIZES = (0.985, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3)
SEARCH_SIZES = tuple(res.PAPER_SEARCH["grid"]["d_over_D"])  # 1.05, 1.15, 1.2, 1.3
SEARCH_F = tuple(res.PAPER_SEARCH["grid"]["f"])  # 0.2, 0.4

#: grids by measured accuracy (see :func:`grid_convergence` and likelihood_meta.grid_check)
FULL_GRID = (60, 31, 120)  # the paper's (alpha, beta, gamma) grid, [TOP] Sec. 14 "Calibration"
FAST_GRID = (36, 19, 72)  # the paper's resonance-search grid, [TOP] Sec. 14 "Resonance search"
FAST_MODELS = (("FCC", 0.985), ("FCC", 1.05), ("BCC", 1.05), ("cubic", 1.05), ("FCC", 1.2))

NULL_SEED, TRUE_SEED = 911, 931  # [TOP] Sec. 14 used seeds 911-914 (null) and 931-942 (resonant)


# --- SO(3) grid ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SO3Grid:
    """Euler grid R = Rz(alpha) Ry(beta) Rz(gamma) with Haar quadrature weights."""

    n_alpha: int
    n_beta: int
    n_gamma: int

    @property
    def alpha(self) -> np.ndarray:
        return 2.0 * np.pi * np.arange(self.n_alpha) / self.n_alpha

    @property
    def gamma(self) -> np.ndarray:
        return 2.0 * np.pi * np.arange(self.n_gamma) / self.n_gamma

    def _gl(self) -> tuple[np.ndarray, np.ndarray]:
        x, w = np.polynomial.legendre.leggauss(self.n_beta)
        return np.arccos(x[::-1]), w[::-1] / 2.0  # beta ascending; weights sum to 1

    @property
    def beta(self) -> np.ndarray:
        return self._gl()[0]

    @property
    def w_beta(self) -> np.ndarray:
        return self._gl()[1]

    @property
    def n_points(self) -> int:
        return self.n_alpha * self.n_beta * self.n_gamma

    def as_dict(self) -> dict:
        return {"n_alpha": self.n_alpha, "n_beta": self.n_beta, "n_gamma": self.n_gamma,
                "alpha": "uniform on [0, 2 pi)", "beta": "Gauss-Legendre nodes in cos(beta) (Haar weight sin beta)",
                "gamma": "uniform on [0, 2 pi), exact Fourier series in gamma evaluated by FFT"}


# --- mode layout --------------------------------------------------------------------------

@dataclass
class _Layout:
    """Mode bookkeeping for an m-period ``s``: K^c couples only modes with m = m' (mod s) and equal parity."""

    n: int
    s: int
    ells: np.ndarray
    ms: np.ndarray
    lslices: dict
    u_blocks: dict  # l -> c2r unitary (a^r = U a^c)
    v_blocks: dict  # l -> eigenvectors of J_y
    w_blocks: dict  # l -> eigenvalues of J_y
    class_idx: list  # per (parity, m mod s) class: mode indices sorted by (m, l)
    steps: list  # per class: list of (m', cols slice, tail start, reduce offsets, deltas, weights)
    mpos_idx: np.ndarray  # modes with m >= 0 sorted by (m, l)
    mpos_bounds: list  # (m, start, stop) within mpos_idx
    u_full: np.ndarray = field(repr=False, default=None)


_LAYOUT_CACHE: dict = {}


def _layout(lmin: int = LMIN, lmax: int = LMAX, s: int = 1) -> _Layout:
    key = (lmin, lmax, s)
    if key in _LAYOUT_CACHE:
        return _LAYOUT_CACHE[key]
    ells, ms = H.lm_index(lmin, lmax)
    lsl = H.block_slices(lmin, lmax)
    u_blocks = {l: H._c2r(l) for l in lsl}
    w_blocks, v_blocks = {}, {}
    for l in lsl:
        w, v = H._jy_eig(l)
        w_blocks[l], v_blocks[l] = w, v
    class_idx, steps = [], []
    for par in (0, 1):
        for r in range(s):
            idx = np.flatnonzero((ells % 2 == par) & (ms % s == r))
            if idx.size == 0:
                continue
            idx = idx[np.lexsort((ells[idx], ms[idx]))]
            msort = ms[idx]
            present = [int(m) for m in np.unique(msort)]
            start = {m: int(np.searchsorted(msort, m, "left")) for m in present}
            stop = {m: int(np.searchsorted(msort, m, "right")) for m in present}
            st = []
            for mp_ in present:
                tail = [m for m in present if m >= abs(mp_)]  # B(m, m') = B(-m', -m): keep m + m' >= 0
                if not tail:
                    continue
                t0 = start[tail[0]]
                st.append((mp_, slice(start[mp_], stop[mp_]), t0,
                           np.array([start[m] - t0 for m in tail], dtype=np.intp),
                           np.array([m - mp_ for m in tail], dtype=np.intp),
                           np.array([1.0 if m == -mp_ else 2.0 for m in tail])))
            class_idx.append(idx)
            steps.append(st)
    sel = np.flatnonzero(ms >= 0)
    sel = sel[np.lexsort((ells[sel], ms[sel]))]
    msel = ms[sel]
    bounds = [(int(m), int(np.searchsorted(msel, m, "left")), int(np.searchsorted(msel, m, "right")))
              for m in np.unique(msel)]
    u_full = np.zeros((ells.size, ells.size), dtype=complex)
    for l, sl in lsl.items():
        u_full[sl, sl] = u_blocks[l]
    lay = _Layout(n=ells.size, s=s, ells=ells, ms=ms, lslices=lsl, u_blocks=u_blocks, v_blocks=v_blocks,
                  w_blocks=w_blocks, class_idx=class_idx, steps=steps, mpos_idx=sel, mpos_bounds=bounds, u_full=u_full)
    _LAYOUT_CACHE[key] = lay
    return lay


def _to_complex(a: np.ndarray, lay: _Layout) -> dict:
    """Real-basis vector -> complex-basis blocks a^c_l = U_l^dagger a^r_l."""
    return {l: lay.u_blocks[l].conj().T @ a[sl] for l, sl in lay.lslices.items()}


def _complex_matrix(k: np.ndarray, lay: _Layout) -> np.ndarray:
    """K^c = U^dagger K U (U block-diagonal in l)."""
    return lay.u_full.conj().T @ k @ lay.u_full


def _m_period(kc: np.ndarray, ms: np.ndarray, rtol: float = 1e-11) -> int:
    """Largest s in (4, 2, 1) with K^c_{lm,l'm'} = 0 unless m = m' (mod s): a 2 pi / s z-rotation symmetry."""
    scale = float(np.max(np.abs(kc)))
    dm = ms[:, None] - ms[None, :]
    for s in (4, 2):
        if float(np.max(np.abs(kc[dm % s != 0]), initial=0.0)) <= rtol * scale:
            return s
    return 1


def _has_beta_reflection(k: np.ndarray, lmin: int, lmax: int, rtol: float = 1e-11) -> bool:
    """True when D(Ry(pi))^T K D(Ry(pi)) = K, so F(alpha, beta, gamma) = F(alpha + pi, pi - beta, pi - gamma)."""
    ry = H.rotation_from_euler(0.0, np.pi, 0.0)
    d = cov.block_diag_rotation(ry, lmin, lmax)
    return bool(np.max(np.abs(d.T @ k @ d - k)) <= rtol * np.max(np.abs(k)))


def _shell_symmetric(shape: str, d_over_D: float, m2: float, s: int, reflect: bool) -> bool:
    """The shell (hence its line term) is invariant under Rz(2 pi / s) and, if ``reflect``, under Ry(pi)."""
    kv = res.shell_vectors(shape, d_over_D, m2)
    ok = True
    rots = [H.rotation_from_euler(2.0 * np.pi / s, 0.0, 0.0)] + ([H.rotation_from_euler(0.0, np.pi, 0.0)] if reflect else [])
    for r in rots:
        rv = kv @ r.T
        dist = np.min(np.linalg.norm(rv[:, None, :] - kv[None, :, :], axis=2), axis=1)
        ok = ok and bool(np.all(dist <= 1e-9 * np.linalg.norm(kv, axis=1)))
    return ok


# --- rotated coefficients ct(alpha, beta) ---------------------------------------------------

def _ct_grid(ac: dict, lay: _Layout, alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """ct for the product grid (beta outer, alpha inner): shape (n_beta * n_alpha, n_modes)."""
    na, nb = alpha.size, beta.size
    out = np.empty((nb, na, lay.n), dtype=complex)
    for l, sl in lay.lslices.items():
        m = np.arange(-l, l + 1)
        xv = (np.exp(1j * np.outer(alpha, m)) * ac[l][None, :]) @ lay.v_blocks[l]  # (na, 2l+1)
        e = np.exp(-1j * np.outer(beta, lay.w_blocks[l]))  # (nb, 2l+1)
        out[:, :, sl] = (xv[None, :, :] * e[:, None, :]) @ lay.v_blocks[l].conj().T
    return out.reshape(nb * na, lay.n)


def _ct_points(ac: dict, lay: _Layout, alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """ct at paired points (alpha_i, beta_i): shape (n_points, n_modes)."""
    alpha, beta = np.atleast_1d(alpha).astype(float), np.atleast_1d(beta).astype(float)
    out = np.empty((alpha.size, lay.n), dtype=complex)
    for l, sl in lay.lslices.items():
        m = np.arange(-l, l + 1)
        x = np.exp(1j * np.outer(alpha, m)) * ac[l][None, :]
        e = np.exp(-1j * np.outer(beta, lay.w_blocks[l]))
        out[:, sl] = ((x @ lay.v_blocks[l]) * e) @ lay.v_blocks[l].conj().T
    return out


def _split(ct: np.ndarray, lay: _Layout) -> list:
    return [np.ascontiguousarray(ct[:, idx]) for idx in lay.class_idx]


# --- quadratic-form Fourier coefficients -------------------------------------------------------

@dataclass
class _TorusModel:
    shape: str
    d_over_D: float
    logdet_ratio: float  # ln det C0 - ln det C_L
    blocks: list  # per class: contiguous (K^c)^T[cols, tail] blocks aligned with the layout steps
    s: int = 1  # detected m-period
    reflect: bool = False  # detected beta reflection symmetry
    k_real: np.ndarray | None = None  # K = C0^{-1} - C_L^{-1} (kept for checks)
    c0_inv: np.ndarray | None = None
    line: "_LineData | None" = None


def _q_coeffs(cp: list, model: _TorusModel, lay: _Layout) -> np.ndarray:
    """q_Delta (Delta = 0..60) for every row of the class-split, m-sorted ct arrays ``cp``."""
    q = np.zeros((cp[0].shape[0], DMAX + 1), dtype=complex)
    for c, steps, blocks in zip(cp, lay.steps, model.blocks):
        cc = c.conj()
        for (_mp, cols, t0, offs, deltas, wts), blk in zip(steps, blocks):
            z = c[:, cols] @ blk
            q[:, deltas] += np.add.reduceat(cc[:, t0:] * z, offs, axis=1) * wts
    return q


def _series_on_grid(q: np.ndarray, n_gamma: int, s: int = 1) -> np.ndarray:
    """F(gamma_k), gamma_k = 2 pi k / n_gamma for k < n_gamma / s (one period 2 pi / s of F).

    F(gamma) = Re q_0 + 2 Re sum_{Delta >= 1} q_Delta e^{-i Delta gamma}; with q_Delta = 0 unless
    s | Delta this is an FFT of length n_gamma / s over j = Delta / s (exact at the grid points).
    """
    ng = n_gamma // s
    coef = q[:, ::s].copy()
    coef[:, 1:] *= 2.0
    h = np.zeros((q.shape[0], ng), dtype=complex)
    np.add.at(h.T, np.arange(coef.shape[1]) % ng, coef.T)
    return np.fft.fft(h, axis=1).real


def _series_at(q: np.ndarray, gamma: np.ndarray, deriv: int = 0) -> np.ndarray:
    """Derivative ``deriv`` of F(gamma) at one gamma per row of q."""
    d = np.arange(DMAX + 1)
    coef = q.copy()
    coef[:, 1:] *= 2.0
    return np.sum(coef * np.exp(-1j * np.outer(gamma, d)) * (-1j * d[None, :]) ** deriv, axis=1).real


def _min_over_gamma(q: np.ndarray, n_fine: int = 720, newton: int = 4) -> tuple[np.ndarray, np.ndarray]:
    """Minimum of F over gamma per row: dense FFT, then guarded Newton steps on the trigonometric series."""
    f = _series_on_grid(q, n_fine)
    k = np.argmin(f, axis=1)
    g = 2.0 * np.pi * k / n_fine
    best = f[np.arange(f.shape[0]), k]
    for _ in range(newton):
        d1, d2 = _series_at(q, g, 1), _series_at(q, g, 2)
        step = np.where(d2 > 0, -d1 / np.where(d2 > 0, d2, 1.0), 0.0)
        g_new = g + np.clip(step, -np.pi / n_fine, np.pi / n_fine)
        f_new = _series_at(q, g_new, 0)
        ok = f_new <= best
        g = np.where(ok, g_new, g)
        best = np.where(ok, f_new, best)
    return g % (2.0 * np.pi), best


# --- line (resonance) data -----------------------------------------------------------------

@dataclass
class _LineData:
    vt_mpos: list  # per m >= 0 group: rows of Vt' = U_c^T (V Q) (complex, group rows x ncols)
    weights: np.ndarray  # (ncols, n_models): 1 / (1/w_f + lambda) on the shell's columns
    const: np.ndarray  # (n_models,) a-independent constant: -1/2 ln det(I + W G) + 1/2 sum ln(C_iso / C_L)
    dinv: np.ndarray  # (n_modes, n_models): 1/C_L - 1/C_iso, for the a-dependent constant
    models: list  # dicts: shape, d_over_D, shell_m2, k_chi, n_vectors, f, grid_index


def _line_data(shape: str, d_over_D: float, c0_inv: np.ndarray, c_l: np.ndarray, lay: _Layout,
               s: int, reflect: bool, fs=SEARCH_F) -> _LineData:
    """Woodbury pieces of the 2 x (shells) line models of one torus ([TOP] Sec. 11.6 grid)."""
    shells = res.list_shells(shape, d_over_D)
    vt_cols, consts, dinvs, models, blocks = [], [], [], [], []
    ncol = 0
    for sh in shells:
        if (s > 1 or reflect) and not _shell_symmetric(shape, d_over_D, sh["m2"], s, reflect):
            raise ValueError(f"shell {sh['m2']} of {shape} {d_over_D} breaks the scan symmetry")
        u, w1 = res.line_lowrank(shape, d_over_D, sh["m2"], 1.0)  # w is linear in f
        v = c0_inv @ u
        g = u.T @ v
        lam, qv = np.linalg.eigh(0.5 * (g + g.T))
        lam = np.clip(lam, 0.0, None)
        vt_cols.append(lay.u_full.T @ (v @ qv))  # y_hat = Q^T V^T b = (U_c^T V Q)^T c
        blocks.append((ncol, ncol + u.shape[1], lam, float(w1[0])))
        ncol += u.shape[1]
    vt = np.concatenate(vt_cols, axis=1)
    weights = []
    for sh, (c0, c1, lam, w1) in zip(shells, blocks):
        for f in fs:
            wf = f * w1
            col = np.zeros(ncol)
            col[c0:c1] = 1.0 / (1.0 / wf + lam)  # M = (W^-1 + U^T V)^-1 = Q diag(1 / (1/w + lambda)) Q^T
            weights.append(col)
            iso = c_l + res.line_isotropic_diag(shape, d_over_D, sh["m2"], f)
            consts.append(-0.5 * float(np.sum(np.log1p(wf * lam))) + 0.5 * float(np.sum(np.log(iso / c_l))))
            dinvs.append(1.0 / c_l - 1.0 / iso)
            models.append({"shape": cov.shape_name(shape), "d_over_D": float(d_over_D), "shell_m2": float(sh["m2"]),
                           "k_chi": float(sh["k_chi"]), "n_vectors": int(sh["n_vectors"]), "f": float(f),
                           "grid_index": int(sh["grid_index"])})
    vt_sorted = vt[lay.mpos_idx]
    return _LineData(vt_mpos=[np.ascontiguousarray(vt_sorted[a:b]) for (_, a, b) in lay.mpos_bounds],
                     weights=np.array(weights).T.copy(), const=np.array(consts), dinv=np.array(dinvs).T.copy(),
                     models=models)


def _gamma_dft(lay: _Layout, gamma: np.ndarray) -> np.ndarray:
    """Real matrix E with y(gamma_k) = [Re h, -Im h] @ E for h_m = Y_m (m = 0) or 2 Y_m (m > 0)."""
    m = np.array([b[0] for b in lay.mpos_bounds])
    ph = np.exp(1j * np.outer(m, gamma))
    return np.concatenate([ph.real, ph.imag], axis=0)


def _line_quad(ct: np.ndarray, line: _LineData, lay: _Layout, edft: np.ndarray) -> np.ndarray:
    """y^T M y on the (rows, gamma) grid for every line model of one torus: (rows, n_gamma, n_models)."""
    rows = ct.shape[0]
    cs = ct[:, lay.mpos_idx]
    ncol = line.weights.shape[0]
    nm = len(lay.mpos_bounds)
    h = np.empty((rows, ncol, 2 * nm))
    for j, ((m, a, b), vt) in enumerate(zip(lay.mpos_bounds, line.vt_mpos)):
        y_m = cs[:, a:b] @ vt
        if m:
            y_m *= 2.0
        h[:, :, j] = y_m.real
        h[:, :, nm + j] = -y_m.imag
    y = h.reshape(rows * ncol, 2 * nm) @ edft  # Re sum_m h_m e^{i m gamma}
    ng = edft.shape[1]
    s2 = np.square(y).reshape(rows, ncol, ng).transpose(0, 2, 1).reshape(rows * ng, ncol)
    return (s2 @ line.weights).reshape(rows, ng, -1)


# --- the engine --------------------------------------------------------------------------------

class _LSE:
    """Running weighted log-sum-exp per model, and the best grid cells (for refinement)."""

    def __init__(self, n_models: int, keep: int = 0):
        self.mx = np.full(n_models, -np.inf)
        self.s = np.zeros(n_models)
        self.keep = keep
        self.top_val = [np.zeros(0) for _ in range(n_models)]
        self.top_pos = [np.zeros((0, 2), dtype=np.intp) for _ in range(n_models)]

    def add(self, i: int, dl: np.ndarray, logw_rows: np.ndarray, row0: int = 0):
        """dl: (rows, n_gamma) for model i, or (rows, n_gamma, k) for models i .. i + k - 1."""
        if dl.ndim == 2:
            dl = dl[:, :, None]
        cm = dl.max(axis=(0, 1))
        cs = np.einsum("r,rgk->k", np.exp(logw_rows), np.exp(dl - cm[None, None, :]))
        for j in range(dl.shape[2]):
            ii = i + j
            m_new = max(self.mx[ii], cm[j])
            prev = self.s[ii] * math.exp(self.mx[ii] - m_new) if np.isfinite(self.mx[ii]) else 0.0
            self.s[ii] = prev + cs[j] * math.exp(cm[j] - m_new)
            self.mx[ii] = m_new
            if self.keep:
                cell = dl[:, :, j]
                kbest = np.argmax(cell, axis=1)
                vals = np.concatenate([self.top_val[ii], cell[np.arange(cell.shape[0]), kbest]])
                pos = np.concatenate([self.top_pos[ii], np.stack([np.arange(cell.shape[0]) + row0, kbest], 1)])
                order = np.argsort(-vals, kind="stable")[: self.keep]
                self.top_val[ii], self.top_pos[ii] = vals[order], pos[order]

    def value(self) -> np.ndarray:
        return self.mx + np.log(self.s)


class Engine:
    """Precomputed torus (and line) models for orientation scans on one SO(3) grid.

    ``symmetry=True`` uses the lattice symmetries detected numerically from K (never assumed):
    a 2 pi / s rotation about the lattice z axis (K^c couples only m = m' mod s, F has period
    2 pi / s in gamma) and the half-turn about y (F(alpha, beta, gamma) = F(alpha + pi, pi - beta,
    pi - gamma), so only beta <= pi / 2 rows are scanned). Both are exact reductions of the same
    quadrature; ``symmetry=False`` scans the full grid (used by the tests to prove equality).
    """

    def __init__(self, models, grid: SO3Grid, *, search: bool = False, search_sizes=SEARCH_SIZES,
                 keep_matrices: bool = False, symmetry: bool = True, lmin: int = LMIN, lmax: int = LMAX,
                 rows_per_chunk: int = 600):
        self.grid = grid
        self.lmin, self.lmax = lmin, lmax
        self.c_l = cov.lcdm_covariance(lmin, lmax)
        self.logdet_l = float(np.sum(np.log(self.c_l)))
        base = _layout(lmin, lmax, 1)
        prepared = []
        for shape, d in models:
            c0 = cov.torus_covariance(shape, float(d), None, lmin, lmax)
            fac = res.Factorized.of(c0)
            k = fac.inv - np.diag(1.0 / self.c_l)
            k = 0.5 * (k + k.T)
            kc = _complex_matrix(k, base)
            s = _m_period(kc, base.ms) if symmetry else 1
            refl = _has_beta_reflection(k, lmin, lmax) if symmetry else False
            prepared.append((shape, float(d), fac, k, kc, s, refl))
        s_all = math.gcd(*[p[5] for p in prepared]) if prepared else 1
        while s_all > 1 and grid.n_gamma % s_all:
            s_all //= 2
        self.s = s_all
        self.reflect = bool(prepared) and all(p[6] for p in prepared) and grid.n_alpha % 2 == 0 and \
            (grid.n_gamma // self.s) % 2 == 0
        self.lay = _layout(lmin, lmax, self.s)
        self.models: list[_TorusModel] = []
        for shape, d, fac, k, kc, s, refl in prepared:
            kct = kc.T  # Z = C @ K^T  <=>  Z_i = sum_j K_ij c_j
            blocks = [[np.ascontiguousarray(kct[np.ix_(idx[cols], idx[t0:])]) for (_, cols, t0, _, _, _) in steps]
                      for idx, steps in zip(self.lay.class_idx, self.lay.steps)]
            m = _TorusModel(shape=cov.shape_name(shape), d_over_D=d, logdet_ratio=fac.logdet - self.logdet_l,
                            blocks=blocks, s=s, reflect=refl, k_real=k if keep_matrices else None,
                            c0_inv=fac.inv if keep_matrices else None)
            if search and any(abs(d - x) < 1e-9 for x in search_sizes):
                m.line = _line_data(shape, d, fac.inv, self.c_l, self.lay, self.s, self.reflect)
            self.models.append(m)
        self.line_models = [lm for t in self.models if t.line is not None for lm in t.line.models]
        self.alpha, self.beta = grid.alpha, grid.beta
        self.n_gamma_eff = grid.n_gamma // self.s
        self.gamma = grid.gamma[: self.n_gamma_eff]
        # beta rows actually scanned and their quadrature weights (reflection pairs row j with n_beta - 1 - j)
        if self.reflect:
            half = (grid.n_beta + 1) // 2
            wb = grid.w_beta[:half].copy()
            wb[: grid.n_beta // 2] *= 2.0
        else:
            half, wb = grid.n_beta, grid.w_beta.copy()
        self.rows_beta, self.w_rows = np.arange(half), wb
        self.edft = _gamma_dft(self.lay, self.gamma)
        self.rows_per_chunk = rows_per_chunk

    def symmetry_info(self) -> dict:
        return {"gamma_period_deg": 360.0 / self.s, "beta_reflection": self.reflect,
                "scanned_fraction": len(self.rows_beta) / self.grid.n_beta / self.s,
                "note": ("cube-aligned lattices: a quarter turn about z gives K^c_{lm,l'm'} = 0 unless m = m' (mod 4); "
                         "the half turn about y maps (alpha, beta, gamma) -> (alpha + pi, pi - beta, pi - gamma). "
                         "Both are detected numerically from K and reduce the same grid quadrature exactly.")}

    # -- full-grid scan ----------------------------------------------------------------------
    def scan(self, a: np.ndarray, *, refine: bool = True, top_cells: int = 3, zoom_levels: int = 7,
             lines: bool = True, only: list[int] | None = None) -> dict:
        """ln B and max Delta lnL per torus model (and ln B per line model) for the sky ``a``."""
        lay, g = self.lay, self.grid
        a = np.asarray(a, dtype=float)
        ac = _to_complex(a, lay)
        sel = list(range(len(self.models))) if only is None else list(only)
        lse = _LSE(len(self.models), keep=top_cells if refine else 0)
        line_off, n_line = {}, 0
        for t in sel:
            if lines and self.models[t].line is not None:
                line_off[t] = n_line
                n_line += len(self.models[t].line.models)
        lse_line = _LSE(n_line) if n_line else None
        a2 = a * a
        nrows = self.rows_beta.size
        nb_chunk = max(1, self.rows_per_chunk // g.n_alpha)
        for j0 in range(0, nrows, nb_chunk):
            j1 = min(nrows, j0 + nb_chunk)
            ct = _ct_grid(ac, lay, self.alpha, self.beta[self.rows_beta[j0:j1]])
            cp = _split(ct, lay)
            logw = np.repeat(np.log(self.w_rows[j0:j1] / (g.n_alpha * self.n_gamma_eff)), g.n_alpha)
            for t in sel:
                m = self.models[t]
                dl = -0.5 * (_series_on_grid(_q_coeffs(cp, m, lay), g.n_gamma, self.s) + m.logdet_ratio)
                lse.add(t, dl, logw, j0 * g.n_alpha)
                if t in line_off:
                    quad = _line_quad(ct, m.line, lay, self.edft)
                    const = m.line.const - 0.5 * (a2 @ m.line.dinv)
                    lse_line.add(line_off[t], dl[:, :, None] + 0.5 * quad + const[None, None, :], logw)
        out = {"lnB": lse.value()[sel], "max_grid": lse.mx[sel], "models": sel}
        if refine:
            best, rots = [], []
            for t in sel:
                v, rot = self._refine(ac, t, lse.top_pos[t], zoom_levels)
                best.append(max(v, lse.mx[t]))
                rots.append(rot)
            out["max"], out["rot"] = np.array(best), rots
        if n_line:
            vals = lse_line.value()
            out["line_lnB"] = np.concatenate([vals[line_off[t]:line_off[t] + len(self.models[t].line.models)]
                                              for t in sel if t in line_off])
        return out

    # -- refinement of the maximum ----------------------------------------------------------------
    def profile(self, ac: dict, t: int, alpha, beta) -> tuple[np.ndarray, np.ndarray]:
        """max over gamma of Delta lnL at points (alpha, beta), and the maximizing gamma."""
        lay, m = self.lay, self.models[t]
        cp = _split(_ct_points(ac, lay, alpha, beta), lay)
        gam, fmin = _min_over_gamma(_q_coeffs(cp, m, lay))
        return -0.5 * (fmin + m.logdet_ratio), gam

    def _refine(self, ac: dict, t: int, top_pos: np.ndarray, levels: int) -> tuple[float, np.ndarray]:
        """Nested 5 x 5 zoom in (alpha, beta) from the best grid cells, gamma maximized exactly."""
        g = self.grid
        ha, hb = 2.0 * np.pi / g.n_alpha, np.pi / g.n_beta
        offs = np.arange(-2, 3)
        best_v, best_r = -np.inf, np.eye(3)
        for row, _k in top_pos:
            j, ia = divmod(int(row), g.n_alpha)
            a0, b0 = self.alpha[ia], self.beta[self.rows_beta[j]]
            v0, gam0 = self.profile(ac, t, np.array([a0]), np.array([b0]))
            cur_v, cur = float(v0[0]), (a0, b0, float(gam0[0]))
            sa, sb = ha, hb
            for _ in range(levels):
                sa, sb = sa / 2.0, sb / 2.0
                aa, bb = np.meshgrid(cur[0] + offs * sa, np.clip(cur[1] + offs * sb, 0.0, np.pi), indexing="ij")
                v, gam = self.profile(ac, t, aa.ravel(), bb.ravel())
                i = int(np.argmax(v))
                if v[i] > cur_v:
                    cur_v, cur = float(v[i]), (float(aa.ravel()[i]), float(bb.ravel()[i]), float(gam[i]))
            if cur_v > best_v:
                best_v, best_r = cur_v, H.rotation_from_euler(*cur)
        return best_v, best_r

    # -- checks ----------------------------------------------------------------------------------
    def grid_values(self, a: np.ndarray, t: int) -> np.ndarray:
        """Delta lnL on all beta rows and the scanned gamma period: (n_beta, n_alpha, n_gamma / s)."""
        ac = _to_complex(np.asarray(a, dtype=float), self.lay)
        cp = _split(_ct_grid(ac, self.lay, self.alpha, self.beta), self.lay)
        f = _series_on_grid(_q_coeffs(cp, self.models[t], self.lay), self.grid.n_gamma, self.s)
        return (-0.5 * (f + self.models[t].logdet_ratio)).reshape(self.grid.n_beta, self.grid.n_alpha, -1)

    def line_grid_values(self, a: np.ndarray, t: int) -> np.ndarray:
        """Delta lnL of every line model of torus t: (n_beta, n_alpha, n_gamma / s, n_models)."""
        m = self.models[t]
        a = np.asarray(a, dtype=float)
        ct = _ct_grid(_to_complex(a, self.lay), self.lay, self.alpha, self.beta)
        dl = self.grid_values(a, t).reshape(ct.shape[0], -1)
        quad = _line_quad(ct, m.line, self.lay, self.edft)
        const = m.line.const - 0.5 * ((a * a) @ m.line.dinv)
        out = dl[:, :, None] + 0.5 * quad + const[None, None, :]
        return out.reshape(self.grid.n_beta, self.grid.n_alpha, self.n_gamma_eff, -1)

    def dlnl_direct(self, a: np.ndarray, t: int, rot: np.ndarray) -> float:
        """Brute force: b = D(R)^T a with real Wigner blocks, Delta lnL = -1/2 (b^T K b + ln det ratio)."""
        m = self.models[t]
        if m.k_real is None:
            raise ValueError("engine built without keep_matrices=True")
        b = H.rotate_real_vector(np.asarray(a, dtype=float), np.asarray(rot).T, self.lmin, self.lmax)
        return float(-0.5 * (b @ m.k_real @ b + m.logdet_ratio))


# --- symmetry helpers -----------------------------------------------------------------------------

def cubic_rotation_group() -> np.ndarray:
    """The 24 proper rotations of the cube (signed permutation matrices with det +1)."""
    import itertools

    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1.0, -1.0), repeat=3):
            m = np.zeros((3, 3))
            for i, (j, s) in enumerate(zip(perm, signs)):
                m[i, j] = s
            if np.linalg.det(m) > 0:
                out.append(m)
    return np.array(out)


def orientation_error_deg(r_est: np.ndarray, r_true: np.ndarray) -> float:
    """Smallest rotation angle between r_est and r_true g, g in the cubic group (O_h acts via +-g)."""
    best = np.inf
    for gm in cubic_rotation_group():
        m = r_true.T @ r_est @ gm
        c = np.clip((np.trace(m) - 1.0) / 2.0, -1.0, 1.0)
        best = min(best, float(np.degrees(np.arccos(c))))
    return best


def lattice_axes(shape) -> np.ndarray:
    """Unit shortest lattice vectors, one per +- pair (cubic 3, BCC 4, FCC 6), lattice frame."""
    b = cov.lattice_basis(shape)
    n = np.stack(np.meshgrid(*(np.arange(-2, 3),) * 3, indexing="ij"), -1).reshape(-1, 3)
    v = n @ b
    r = np.linalg.norm(v, axis=1)
    v = v[np.abs(r - 1.0) < 1e-9]
    keep = []
    for x in v:
        if not any(np.allclose(x, -y) or np.allclose(x, y) for y in keep):
            keep.append(x)
    return np.array(keep)


def axes_galactic(shape, rot: np.ndarray) -> list[list[float]]:
    """Sky-frame lattice axes R v in Galactic (l, b), each pair represented by its b >= 0 member."""
    out = []
    for v in lattice_axes(shape) @ np.asarray(rot).T:
        if v[2] < 0 or (abs(v[2]) < 1e-12 and np.arctan2(v[1], v[0]) < 0):
            v = -v
        l, b = H.vec_to_gal(v)
        out.append([round(float(l), 2), round(float(b), 2)])
    return sorted(out, key=lambda x: -x[1])


# --- multiprocessing (null and torus-true skies) ---------------------------------------------------

_WORKER: dict = {}


def _limit_blas_threads(n: int = 1) -> None:
    """Set OpenBLAS threads in this process (env vars are read only at library load)."""
    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = str(n)
    try:
        with open("/proc/self/maps") as fh:
            libs = {line.split()[-1] for line in fh if "openblas" in line.lower() and line.split()[-1].startswith("/")}
    except OSError:
        return
    for path in libs:
        try:
            lib = ctypes.CDLL(path)
        except OSError:
            continue
        for name in ("openblas_set_num_threads", "openblas_set_num_threads64_", "scipy_openblas_set_num_threads64_",
                     "scipy_openblas_set_num_threads"):
            fn = getattr(lib, name, None)
            if fn is not None:
                fn(ctypes.c_int(n))
                break


def _worker_init():
    _limit_blas_threads(1)


def _null_task(i: int) -> dict:
    eng: Engine = _WORKER["engine"]
    rng = np.random.default_rng([NULL_SEED, i])
    a = np.sqrt(eng.c_l) * rng.standard_normal(eng.c_l.size)
    out = eng.scan(a, refine=_WORKER["refine"])
    return {"i": i, "lnB": out["lnB"], "max": out.get("max", out["max_grid"]), "line_lnB": out.get("line_lnB")}


def _true_task(job: tuple[int, int]) -> dict:
    t, i = job
    eng: Engine = _WORKER["engine"]
    rng = np.random.default_rng([TRUE_SEED, t, i])
    m = eng.models[t]
    chol = _WORKER.setdefault("chol", {})
    if t not in chol:
        chol[t] = np.linalg.cholesky(cov.torus_covariance(m.shape, m.d_over_D, None, eng.lmin, eng.lmax))
    x = chol[t] @ rng.standard_normal(chol[t].shape[0])
    rot = H.random_rotation(rng)
    a = H.rotate_real_vector(x, rot, eng.lmin, eng.lmax)
    out = eng.scan(a, refine=False, lines=False, only=[t])
    return {"t": t, "i": i, "lnB": float(out["lnB"][0])}


def _run_pool(fn, jobs, workers: int):
    if workers <= 1 or len(jobs) <= 1:
        return [fn(j) for j in jobs]
    ctx = mp.get_context("fork")
    with ctx.Pool(workers, initializer=_worker_init) as pool:
        return list(pool.imap_unordered(fn, jobs, chunksize=1))


# --- bundle -----------------------------------------------------------------------------------------

def _paper_row(shape: str, d: float) -> dict:
    from ..topology.manifolds import EVIDENCE_TORUS, EXCLUSION_TORUS

    ev = next((r for r in EVIDENCE_TORUS if r["shape"] == shape and abs(r["d_over_D"] - d) < 1e-9), None)
    ex = next((r for r in EXCLUSION_TORUS if r["shape"] == shape and abs(r["d_over_D"] - d) < 1e-9), None)
    out: dict = {}
    if ev is not None:
        out.update({"max_dlnL": ev["max_dlnL"], "lnB": ev["lnB"], "p": ev["p"], "null95": ev["null95"],
                    "source": ev["source"]})
    if ex is not None:
        out.update({"lnB_planck": ex["lnB_planck"], "verdict": ex["verdict"], "median_if_true": ex["median_if_true"],
                    "q05_if_true": ex["q05_if_true"], "p_le_planck": ex["p_le_planck"], "sims_if_true": ex["sims"]})
        out.setdefault("source", ex["source"])
        if ev is not None:
            out["source"] = "[TOP] Tables 9a, 9b"
    if out:
        out["map"] = "Planck 2018 Commander (PR3)"
    return out


def _pvalue(null: np.ndarray, obs: float) -> float:
    null = np.asarray(null, dtype=float)
    return float(np.mean(null >= obs)) if null.size else float("nan")


def _lse_mean(x: np.ndarray, axis=-1) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    mx = np.max(x, axis=axis, keepdims=True)
    return np.squeeze(mx, axis=axis) + np.log(np.mean(np.exp(x - mx), axis=axis))


def export_fields(sky, fast: bool = False, null_sims: int = 100, *, workers: int | None = None,
                  true_sims: int | None = None) -> dict:
    """``cmb.json -> likelihood, likelihood_meta`` and ``resonance.json -> search`` for the map in ``sky``."""
    t_start = time.perf_counter()
    workers = workers or min(4, os.cpu_count() or 1)
    grid = SO3Grid(*(FAST_GRID if fast else FULL_GRID))
    models = list(FAST_MODELS) if fast else [(s, d) for s in SHAPES for d in LIKELIHOOD_SIZES]
    search_models = [(s, d) for s in SHAPES for d in SEARCH_SIZES]
    for sm in search_models:
        if sm not in models:
            models.append(sm)
    n_null = int(min(null_sims, 12) if fast else null_sims)
    n_true = int(true_sims if true_sims is not None else (0 if fast else 40))
    a = np.asarray(sky.alm_real, dtype=float)
    map_name = getattr(sky, "source_name", "unknown")

    t0 = time.perf_counter()
    eng = Engine(models, grid, search=True)
    t_pre = time.perf_counter() - t0

    t0 = time.perf_counter()
    obs = eng.scan(a, refine=True)
    t_obs = time.perf_counter() - t0

    _WORKER.clear()
    _WORKER.update(engine=eng, refine=True)
    t0 = time.perf_counter()
    null = sorted(_run_pool(_null_task, list(range(n_null)), workers), key=lambda r: r["i"])
    t_null = time.perf_counter() - t0
    lik_idx = [i for i, m in enumerate(models) if fast or m in [(s, d) for s in SHAPES for d in LIKELIHOOD_SIZES]]
    t0 = time.perf_counter()
    true = _run_pool(_true_task, [(t, i) for t in lik_idx for i in range(n_true)], workers) if n_true else []
    t_true = time.perf_counter() - t0
    _WORKER.clear()

    null_lnb = np.array([r["lnB"] for r in null]).reshape(n_null, len(models))
    null_max = np.array([r["max"] for r in null]).reshape(n_null, len(models))

    rows = []
    for t in lik_idx:
        m = eng.models[t]
        lnb, mx = float(obs["lnB"][t]), float(obs["max"][t])
        row = {"shape": m.shape, "d_over_D": m.d_over_D, "max_dlnL": mx, "lnB": lnb,
               "null_q95": float(np.quantile(null_lnb[:, t], 0.95)) if n_null else None,
               "p": _pvalue(null_lnb[:, t], lnb) if n_null else None,
               "null_q95_max": float(np.quantile(null_max[:, t], 0.95)) if n_null else None,
               "p_max": _pvalue(null_max[:, t], mx) if n_null else None,
               "null_median": float(np.median(null_lnb[:, t])) if n_null else None,
               "n_null": n_null, "max_dlnL_grid": float(obs["max_grid"][t]),
               "best_axes_gal": axes_galactic(m.shape, obs["rot"][t]),
               "best_euler_deg": [round(float(np.degrees(x)), 3) for x in H.euler_zyz(obs["rot"][t])],
               "map": map_name, "paper": _paper_row(m.shape, m.d_over_D), "computed": True}
        tl = np.array([r["lnB"] for r in true if r["t"] == t])
        if tl.size:
            p_le = float(np.mean(tl <= lnb))
            row["exclusion"] = {"median_if_true": float(np.median(tl)), "q05_if_true": float(np.quantile(tl, 0.05)),
                                "p_le_obs": p_le, "n_true": int(tl.size),
                                "verdict": "excluded (95%)" if p_le < 0.05 else "allowed"}
        rows.append(row)

    # global statistics over our likelihood models (the paper averaged over its own 27-model grid)
    sel = np.array(lik_idx)
    glob_obs = float(_lse_mean(obs["lnB"][sel]))
    glob_null = _lse_mean(null_lnb[:, sel], axis=1) if n_null else np.zeros(0)
    fcc = np.array([t for t in lik_idx if eng.models[t].shape == "FCC" and eng.models[t].d_over_D >= 1.0])
    fcc_obs = float(_lse_mean(obs["lnB"][fcc])) if fcc.size else None
    fcc_null = _lse_mean(null_lnb[:, fcc], axis=1) if (n_null and fcc.size) else np.zeros(0)
    best_t = int(sel[np.argmax(obs["max"][sel])])

    # resonance search
    line_obs = obs["line_lnB"]
    line_null = np.array([r["line_lnB"] for r in null]) if n_null else np.zeros((0, line_obs.size))
    s_max, s_avg = float(np.max(line_obs)), float(_lse_mean(line_obs))
    null_smax = line_null.max(axis=1) if n_null else np.zeros(0)
    null_savg = _lse_mean(line_null, axis=1) if n_null else np.zeros(0)
    smodels = []
    for i, md in enumerate(eng.line_models):
        smodels.append({**md, "lnB": float(line_obs[i]),
                        "p_local": _pvalue(line_null[:, i], float(line_obs[i])) if n_null else None})
    ib = int(np.argmax(line_obs))
    paper_search = {k: res.PAPER_SEARCH[k] for k in ("S_max", "p_max", "S_avg", "p_avg", "null_median_S_max",
                                                      "best_model", "power", "source")}
    paper_search["map"] = "Planck 2018 SMICA (inpainted)"
    search = {
        "map": map_name, "models": smodels, "n_models": len(smodels), "S_max": s_max, "S_avg": s_avg,
        "best": smodels[ib], "null": {"S_max": null_smax.tolist(), "S_avg": null_savg.tolist()},
        "null_median_S_max": float(np.median(null_smax)) if n_null else None,
        "p_max": _pvalue(null_smax, s_max) if n_null else None, "p_avg": _pvalue(null_savg, s_avg) if n_null else None,
        "n_null": n_null, "grid": grid.as_dict(), "paper": paper_search, "computed": True,
        "statistic": ("ln B per model = ln <exp(Delta lnL(R))>_Haar of torus + line (C0 + U W U^T) against infinite "
                      "space + the same isotropic line; S_max = max over models, S_avg = ln(mean over models of B)"),
        "notes": (f"Independent re-analysis on {map_name}, not a reproduction of the Planck numbers. "
                  "Woodbury identity with linear rotational correlations y_j(R) as exact Fourier series in gamma. "
                  "Resonance gives up to 64x information only at KNOWN orientation ([TOP] Sec. 11.6); averaging over "
                  "SO(3) and the 160-model grid consumes most of it. A null result here neither detects nor excludes "
                  "resonance."),
    }

    runtime = time.perf_counter() - t_start
    meta = {
        "map": map_name, "n_null": n_null, "n_true": n_true, "grid": grid.as_dict(),
        "models": [{"shape": eng.models[t].shape, "d_over_D": eng.models[t].d_over_D} for t in lik_idx],
        "statistics": {
            "dlnL": "Delta lnL(R) = -1/2 a^T (C(R)^-1 - C_L^-1) a - 1/2 ln(det C0 / det C_L), C(R) = D(R) C0 D(R)^T",
            "max_dlnL": "max over SO(3): best grid cells refined by nested zoom in (alpha, beta), gamma maximized exactly",
            "lnB": "ln <exp(Delta lnL(R))>_Haar on the grid (orientation-averaged Bayes factor, the paper's optimal statistic)",
            "p": "fraction of isotropic LCDM null skies with ln B >= observed",
            "p_max": "fraction of null skies with max Delta lnL >= observed",
            "null_q95": "95% quantile of ln B over the null skies",
            "exclusion": "torus-true skies (Haar-random orientation): P(ln B <= observed | torus true); < 0.05 -> excluded (95%)",
        },
        "global": {"lnB_avg": glob_obs, "p": _pvalue(glob_null, glob_obs) if n_null else None,
                   "n_models": int(sel.size), "definition": "ln(mean over the likelihood models of B)",
                   "paper": {"lnB_avg": -0.06, "p": 0.15, "n_models": 27, "source": "[TOP] Sec. 9"}},
        "fcc_family": {"lnB": fcc_obs, "p": _pvalue(fcc_null, fcc_obs) if (n_null and fcc_obs is not None) else None,
                       "sizes": [eng.models[t].d_over_D for t in fcc],
                       "definition": "ln(mean of B over the FCC sizes d/D >= 1 of this table)",
                       "paper": {"lnB": -0.23, "p": 0.25, "source": "[TOP] Sec. 9 (pre-registered FCC family, its own size grid)"}},
        "best_model": {"shape": eng.models[best_t].shape, "d_over_D": eng.models[best_t].d_over_D,
                       "max_dlnL": float(obs["max"][best_t]),
                       "axes_gal": axes_galactic(eng.models[best_t].shape, obs["rot"][best_t]),
                       "note": "largest max Delta lnL of the table; orientations are not significant, not detected structures"},
        "runtime_s": {"total": round(runtime, 1), "precompute": round(t_pre, 1), "observed": round(t_obs, 1),
                      "null": round(t_null, 1), "torus_true": round(t_true, 1), "workers": workers},
        "notes": [
            f"Independent re-analysis on {map_name} (full-sky a_lm, 2 <= l <= 30); the paper's Table 9 used the Planck 2018 "
            "Commander map, so the numbers are not expected to agree. Foreground residuals of a full-sky map near the "
            "Galactic plane enter both maps differently; the paper's masked check moved FCC 1.05 D from -1.20 to -0.76.",
            "No detection is claimed: ln B and max Delta lnL are calibrated on isotropic skies processed with the identical "
            "code and grid. Best-fit orientations are reported for the record only; they are not significant.",
            "Fourier scan in gamma is exact; it matches brute-force Wigner rotation of the data to < 1e-9 (tests).",
        ],
        "computed": True, "source": "[TOP] Sec. 9, Tables 9a/9b, Sec. 14",
    }
    return {"likelihood": rows, "likelihood_meta": meta, "search": search}


# --- grid accuracy check ---------------------------------------------------------------------------

def grid_convergence(a: np.ndarray, models, grids=((24, 13, 48), (36, 19, 72), (60, 31, 120), (90, 46, 180))) -> dict:
    """ln B and grid max per model on several grids (used to choose FULL_GRID / FAST_GRID)."""
    out = {"grids": [list(g) for g in grids], "models": [list(m) for m in models], "lnB": [], "max_grid": [], "seconds": []}
    for g in grids:
        eng = Engine(models, SO3Grid(*g))
        t0 = time.perf_counter()
        r = eng.scan(a, refine=False, lines=False)
        out["seconds"].append(time.perf_counter() - t0)
        out["lnB"].append(r["lnB"].tolist())
        out["max_grid"].append(r["max_grid"].tolist())
    return out
