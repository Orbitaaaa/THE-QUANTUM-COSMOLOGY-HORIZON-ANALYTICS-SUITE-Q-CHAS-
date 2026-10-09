"""Topological resonance: a primordial line on a dual-lattice shell, [TOP] Secs. 11.6-11.8 and Sec. 14.

A narrow line P_line(k) = a delta(ln k - ln k0) in the curvature spectrum P_R adds, in infinite
space, the isotropic diagonal Delta C_l = 4 pi a Delta_l(k0)^2. On the torus the same power sits
only on the N dual vectors of the shell |kappa| = k0, normalized by the addition theorem so that
the orientation-averaged C_l are identical ([TOP] Sec. 14, "Resonance search"):

    Delta C_T = (4 pi a)(4 pi / N) sum_{kappa in shell} Delta_l(k0) Delta_l'(k0) i^{l-l'} Y*_lm Y_l'm'
              = U diag(w) U^T,   U[lm, kappa] = r_l Delta_l(k0) Y^r_lm(kappa_hat),  w = 16 pi^2 a / N,

a real term of rank <= N (the shell is closed under kappa -> -kappa, so blocks with l + l'
odd cancel). The amplitude a is fixed by the fractional bump f = max_l Delta C_l / C_l^LCDM.
Inverses, determinants, KL divergences and Gaussian log-likelihoods of C0 + U W U^T use the
Woodbury identity and the matrix determinant lemma, O(n^2 N) once C0 is factorized.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.linalg import cho_factor, cholesky, solve_triangular

from . import covariance as cov
from . import harmonics as H
from .transfer import load_transfer

LMIN, LMAX = cov.LMIN, cov.LMAX

#: [TOP] Sec. 11.6 table (page 33): KL for the cubic torus at known orientation, l <= 30
PAPER_RESONANCE = {
    "columns": ["no_line", "m2_9_f0.25", "m2_9_f0.5", "m2_16_f0.25", "m2_16_f0.5", "gain"],
    "shells": {9: 30, 16: 6},
    "rows": {
        1.05: [5.05, 6.16, 8.03, 10.32, 16.93, 3],
        1.15: [0.64, 1.21, 2.48, 4.88, 10.49, 16],
        1.2: [0.33, 0.81, 1.93, 4.14, 9.28, 28],
        1.3: [0.12, 0.47, 1.29, 3.33, 7.77, 64],
    },
    "source": "[TOP] Sec. 11.6",
}
#: [TOP] Sec. 11.6 pre-registered grid and Planck result
PAPER_SEARCH = {
    "grid": {"shapes": ["FCC", "BCC", "cubic"], "d_over_D": [1.05, 1.15, 1.2, 1.3], "max_vectors": 12,
             "k0_chi_range": [2.0, 14.0], "n_shells": 80, "f": [0.2, 0.4], "n_models": 160},
    "S_max": 1.01, "p_max": 0.78, "S_avg": -0.03, "p_avg": 0.53, "null_median_S_max": 1.65,
    "best_model": {"shape": "BCC", "d_over_D": 1.05, "k0_chi": 5.18, "n_vectors": 6, "lnB": 1.01, "p_local": 0.11},
    "power": {"cubic_1.2D_f0.4_shell7": 0.35, "FCC_1.2D_f0.4_shell5": 0.30, "level": 0.05, "skies": 20},
    "planck_isotropic_line_cap_f95": [0.3, 0.9],
    "source": "[TOP] Sec. 11.6",
}
#: [TOP] Sec. 11.8 Exit 1: single-field limits on the lowest shells (CMB, l <= 30, known orientation)
PAPER_SINGLE_FIELD = {
    "growth_bound": {"dlnP_dlnk_max": 4, "refined": 5, "a_max": "h (e^{5h} - 1)", "f": "a G(k0)"},
    "shells_with_f_below_0.05": {"count": 32, "of": 80},
    "rows": [
        {"torus": "cubic", "d_over_D": 1.2, "k0_chi": 2.62, "n_vectors": 6, "single_field_max_f": 0.95, "planck_cap": 0.71, "kl_no_line": 0.33, "kl_allowed_line": 1.12, "gain": 3.4},
        {"torus": "cubic", "d_over_D": 1.3, "k0_chi": 2.42, "n_vectors": 6, "single_field_max_f": 0.84, "planck_cap": 0.75, "kl_no_line": 0.12, "kl_allowed_line": 0.81, "gain": 6.7},
        {"torus": "BCC", "d_over_D": 1.15, "k0_chi": 4.73, "n_vectors": 6, "single_field_max_f": 0.35, "planck_cap": 0.85, "kl_no_line": 0.95, "kl_allowed_line": 2.12, "gain": 2.2},
        {"torus": "BCC", "d_over_D": 1.2, "k0_chi": 4.53, "n_vectors": 6, "single_field_max_f": 0.34, "planck_cap": 0.79, "kl_no_line": 0.50, "kl_allowed_line": 1.43, "gain": 2.9},
    ],
    "source": "[TOP] Sec. 11.8",
}
#: [TOP] Sec. 11.8: perfect 3D map of the causal ball (cubic torus, l <= 20, six-vector shells, known orientation)
PAPER_EXITS_3D = {
    "rows": [
        {"spectrum": "smooth (Theorem L applies)", "d1.5": 0.063, "d2": 0.006, "d3": 0.0003},
        {"spectrum": "single-field bump on the fundamental shell, a = 1", "d1.5": 0.50, "d2": 0.12, "d3": 0.009},
        {"spectrum": "single-field narrow line, shell |m|^2 = 16, a = 0.005", "d1.5": 0.089, "d2": 0.012, "d3": 0.0009},
        {"spectrum": "multi-field narrow line, shell |m|^2 = 16, a = 0.05", "d1.5": None, "d2": 0.42, "d3": 0.051},
        {"spectrum": "multi-field narrow line, shell |m|^2 = 16, a = 0.2", "d1.5": 9.07, "d2": 3.18, "d3": 0.54},
        {"spectrum": "multi-field narrow line, shell |m|^2 = 64, a = 0.2", "d1.5": None, "d2": 32.16, "d3": 9.01},
    ],
    "source": "[TOP] Sec. 11.8",
}
#: [TOP] Sec. 11.2 / Theorem L ceiling (cubic torus), quoted
PAPER_CEILING = {
    "d_over_D": [1.05, 1.1, 1.2, 1.3, 1.5, 2.0],
    "causal_ball": [10.48, 3.04, 0.74, 0.31, 0.087, 0.007],
    "boundary_sphere": [8.74, 2.29, 0.55, 0.21, 0.052, 0.005],
    "lss_ball": [4.79, 1.76, 0.50, 0.20, 0.055, 0.006],
    "cmb_l30": [5.05, 1.48, 0.33, 0.12, 0.029, None],
    "source": "[TOP] Sec. 11.2",
}


# --- shells ---------------------------------------------------------------------------

def _dual_by_norm(shape, d_over_D: float, kmax_chi: float):
    tf = load_transfer()
    b = cov.scaled_basis(shape, d_over_D, tf.chi_star)
    kv = cov.dual_vectors(b, kmax_chi / tf.chi_star * (1 + 1e-9))
    x = np.linalg.norm(kv, axis=1) * tf.chi_star
    order = np.argsort(x, kind="stable")
    return kv[order], x[order], tf.chi_star


def _group(x: np.ndarray, rtol: float = 1e-9) -> list[np.ndarray]:
    if x.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(x) > rtol * x[1:]) + 1
    return np.split(np.arange(x.size), breaks)


def shell_label(k_chi: float, d_over_D: float) -> float:
    """m2 = (|kappa| d / 2 pi)^2; for the cubic torus the integer |m|^2 of kappa = 2 pi m / d."""
    return float((k_chi * d_over_D / np.pi) ** 2)


def all_shells(shape, d_over_D: float, kmax_chi: float = 14.0) -> list[dict]:
    """Every dual-lattice shell with k chi_* <= kmax_chi: index (1 = fundamental), m2, k_chi, N."""
    _, x, _ = _dual_by_norm(shape, d_over_D, kmax_chi)
    out = []
    for i, g in enumerate(_group(x), start=1):
        kc = float(x[g].mean())
        out.append({"index": i, "m2": round(shell_label(kc, d_over_D), 9), "k_chi": kc, "n_vectors": int(g.size)})
    return out


def list_shells(shape, d_over_D: float, max_vectors: int = 12, kchi_range: tuple[float, float] = (2.0, 14.0)) -> list[dict]:
    """Shells with at most ``max_vectors`` vectors and k0 chi_* in ``kchi_range`` ([TOP] Sec. 11.6 grid)."""
    lo, hi = kchi_range
    sel = [s for s in all_shells(shape, d_over_D, hi) if s["n_vectors"] <= max_vectors and lo <= s["k_chi"] <= hi]
    for j, s in enumerate(sel, start=1):
        s["grid_index"] = j
    return sel


def preregistered_grid(shapes=("FCC", "BCC", "cubic"), sizes=(1.05, 1.15, 1.2, 1.3)) -> list[dict]:
    """The pre-registered resonance grid of [TOP] Sec. 11.6 (paper: 80 shells)."""
    return [dict(shape=s, d_over_D=d, **sh) for s in shapes for d in sizes for sh in list_shells(s, d)]


def shell_vectors(shape, d_over_D: float, m2: float, R: np.ndarray | None = None) -> np.ndarray:
    """Dual vectors (Mpc^-1) of the shell with label m2 = (|kappa| d / 2 pi)^2, rotated by R."""
    tf = load_transfer()
    k0 = 2.0 * np.pi * np.sqrt(m2) / (d_over_D * 2.0 * tf.chi_star)
    b = cov.scaled_basis(shape, d_over_D, tf.chi_star, R)
    kv = cov.dual_vectors(b, k0 * (1 + 1e-7))
    kv = kv[np.abs(np.linalg.norm(kv, axis=1) / k0 - 1.0) < 1e-7]
    if kv.shape[0] == 0:
        raise ValueError(f"no dual vectors with m2 = {m2} for {shape} at d/D = {d_over_D}")
    return kv


# --- line amplitude and covariance terms -------------------------------------------------

def line_gain(k0: float, lmin: int = LMIN, lmax: int = LMAX) -> np.ndarray:
    """4 pi Delta_l(k0)^2 / C_l^LCDM: fractional bump per unit line strength a (in P_R units)."""
    tf = load_transfer()
    return 4.0 * np.pi * tf.delta_table(np.array([k0]), lmin, lmax)[0] ** 2 / cov.lcdm_cl_iso(lmin, lmax)


def line_amplitude(k0: float, f: float, lmin: int = LMIN, lmax: int = LMAX) -> float:
    """Line strength a (absolute, P_R units) whose isotropic bump peaks at max_l Delta C_l / C_l = f."""
    return float(f / np.max(line_gain(k0, lmin, lmax)))


def isotropic_line_cl(k0: float, a: float, lmin: int = LMIN, lmax: int = LMAX) -> np.ndarray:
    """Infinite-space line: Delta C_l = 4 pi a Delta_l(k0)^2 (uK^2)."""
    tf = load_transfer()
    return 4.0 * np.pi * a * tf.delta_table(np.array([k0]), lmin, lmax)[0] ** 2


def line_lowrank(shape, d_over_D: float, shell: float, f: float, R: np.ndarray | None = None,
                 lmin: int = LMIN, lmax: int = LMAX) -> tuple[np.ndarray, np.ndarray]:
    """Torus line term Delta C_T = U diag(w) U^T for shell label m2 = ``shell`` at bump ``f``.

    Returns U (n_modes, N) in uK and weights w (N,). Rotating by R rotates the shell with the lattice.
    """
    tf = load_transfer()
    kv = shell_vectors(shape, d_over_D, shell, R)
    k0 = float(np.linalg.norm(kv, axis=1).mean())
    a = line_amplitude(k0, f, lmin, lmax)
    u = cov.mode_vectors(kv, tf, lmin, lmax, k_norm=k0).T
    return u, np.full(kv.shape[0], 16.0 * np.pi**2 * a / kv.shape[0])


def line_isotropic_diag(shape, d_over_D: float, shell: float, f: float, lmin: int = LMIN, lmax: int = LMAX) -> np.ndarray:
    """Diagonal (per mode) of the infinite-space isotropic line matching :func:`line_lowrank`."""
    tf = load_transfer()
    k0 = 2.0 * np.pi * np.sqrt(shell) / (d_over_D * 2.0 * tf.chi_star)
    ells, _ = H.lm_index(lmin, lmax)
    return isotropic_line_cl(k0, line_amplitude(k0, f, lmin, lmax), lmin, lmax)[ells - lmin]


# --- Woodbury identity and determinant lemma ----------------------------------------------

@dataclass
class Factorized:
    """Cholesky factor and inverse of a base covariance C0, reused for many low-rank updates."""

    chol: np.ndarray  # lower triangular L, C0 = L L^T
    logdet: float
    inv: np.ndarray  # C0^{-1}, so that each solve is one BLAS-3 product

    @classmethod
    def of(cls, c0: np.ndarray) -> "Factorized":
        l, _ = cho_factor(np.asarray(c0, dtype=float), lower=True)
        l = np.tril(l)
        linv = solve_triangular(l, np.eye(l.shape[0]), lower=True)
        return cls(l, float(2.0 * np.sum(np.log(np.diag(l)))), linv.T @ linv)

    def solve(self, b: np.ndarray) -> np.ndarray:
        return self.inv @ b

    def inverse(self) -> np.ndarray:
        return self.inv


def _capacitance(z: np.ndarray, u: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """W as matrix and I + W U^T C0^{-1} U, with z = C0^{-1} U."""
    wm = np.diag(w) if np.ndim(w) == 1 else np.asarray(w, dtype=float)
    return wm, np.eye(wm.shape[0]) + wm @ (u.T @ z)


def woodbury_inverse_logdet(c0_inv: np.ndarray, logdet_c0: float, U: np.ndarray, w) -> tuple[np.ndarray, float]:
    """(C0 + U W U^T)^{-1} and ln det(C0 + U W U^T) from C0^{-1} in O(n^2 N).

    C^{-1} = C0^{-1} - Z (I + W U^T Z)^{-1} W Z^T with Z = C0^{-1} U; ln det C = ln det C0 + ln det(I + W U^T Z).
    ``w`` is a weight vector (W = diag w) or a symmetric N x N matrix; W may be singular.
    """
    z = c0_inv @ U
    wm, cap = _capacitance(z, U, w)
    sign, ld = np.linalg.slogdet(cap)
    if sign <= 0:
        raise np.linalg.LinAlgError("updated covariance is not positive definite")
    inv = c0_inv - z @ np.linalg.solve(cap, wm @ z.T)
    return 0.5 * (inv + inv.T), float(logdet_c0 + ld)


def woodbury_logdet(base: Factorized, U: np.ndarray, w) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """ln det(C0 + U W U^T) via the determinant lemma; also returns Z = C0^{-1} U, W and the capacitance."""
    z = base.solve(U)
    wm, cap = _capacitance(z, U, w)
    sign, ld = np.linalg.slogdet(cap)
    if sign <= 0:
        raise np.linalg.LinAlgError("updated covariance is not positive definite")
    return float(base.logdet + ld), z, wm, cap


def gaussian_loglike(x: np.ndarray, c: np.ndarray) -> float:
    """ln N(x; 0, C) by one dense Cholesky factorization and one triangular solve, O(n^3 / 3)."""
    l = cholesky(np.asarray(c, dtype=float), lower=True, check_finite=False)
    y = solve_triangular(l, x, lower=True, check_finite=False)
    return float(-0.5 * (y @ y) - np.sum(np.log(np.diag(l))) - 0.5 * x.size * np.log(2.0 * np.pi))


def woodbury_loglike(x: np.ndarray, base: Factorized, U: np.ndarray, w, y0: np.ndarray | None = None) -> float:
    """ln N(x; 0, C0 + U W U^T) with C0 already factorized: O(n^2 N).

    ``y0 = C0^{-1} x`` may be passed when x is fixed and only the line term changes (orientation scans).
    """
    logdet, z, wm, cap = woodbury_logdet(base, U, w)
    y = base.solve(x) if y0 is None else y0
    corr = z.T @ x
    quad = x @ y - corr @ np.linalg.solve(cap, wm @ corr)
    return float(-0.5 * quad - 0.5 * logdet - 0.5 * x.size * np.log(2.0 * np.pi))


def kl_lowrank(base: Factorized, U: np.ndarray, w, c_ref_diag: np.ndarray, base_trace_term: float | None = None,
               c0: np.ndarray | None = None) -> float:
    """KL(N(0, C0 + U W U^T) || N(0, diag(c_ref))) by the determinant lemma.

    Needs tr(diag^{-1} C0) (``base_trace_term``) or C0 itself.
    """
    if base_trace_term is None:
        if c0 is None:
            raise ValueError("kl_lowrank needs base_trace_term or c0")
        base_trace_term = float(np.sum(np.diag(c0) / c_ref_diag))
    logdet1, _, wm, _ = woodbury_logdet(base, U, w)
    tr = base_trace_term + float(np.trace(wm @ (U.T @ (U / c_ref_diag[:, None]))))
    n = c_ref_diag.size
    return 0.5 * (tr - n - logdet1 + float(np.sum(np.log(c_ref_diag))))


# --- resonance divergences -------------------------------------------------------------------

@lru_cache(maxsize=8)
def _torus_factor(shape: str, d_over_D: float, lmin: int, lmax: int) -> tuple[Factorized, np.ndarray]:
    c = cov.torus_covariance(shape, d_over_D, None, lmin, lmax)
    return Factorized.of(c), np.diag(c).copy()


def resonance_kl(shape: str, d_over_D: float, shell: float | None, f: float, lmin: int = LMIN, lmax: int = LMAX) -> float:
    """KL(torus + shell line || infinite space + isotropic line) at known orientation.

    ``shell=None`` or ``f=0`` gives the no-line divergence of [TOP] Theorem P.
    """
    base, cdiag = _torus_factor(shape, float(d_over_D), lmin, lmax)
    ref = cov.lcdm_covariance(lmin, lmax)
    if shell is None or f == 0:
        u, w = np.zeros((ref.size, 0)), np.zeros(0)
    else:
        u, w = line_lowrank(shape, d_over_D, shell, f, None, lmin, lmax)
        ref = ref + line_isotropic_diag(shape, d_over_D, shell, f, lmin, lmax)
    return kl_lowrank(base, u, w, ref, float(np.sum(cdiag / ref)))


def resonance_table(sizes=(1.05, 1.15, 1.2, 1.3), shells=(9, 16), fs=(0.25, 0.5), shape: str = "cubic") -> list[dict]:
    """[TOP] Sec. 11.6 table: KL with and without a line, gain = max KL / no-line KL."""
    rows = []
    for d in sizes:
        k_none = resonance_kl(shape, d, None, 0.0)
        best = k_none
        for m2 in shells:
            n_vec = shell_vectors(shape, d, m2).shape[0]
            for j, f in enumerate(fs):
                kl = resonance_kl(shape, d, m2, f)
                best = max(best, kl)
                row = {"shape": shape, "d_over_D": d, "no_line": k_none, "shell": m2, "n_vectors": n_vec, "f": f, "kl": kl,
                       "computed": True}
                prow = PAPER_RESONANCE["rows"].get(d) if shape == "cubic" else None
                if prow is not None and m2 in (9, 16) and f in (0.25, 0.5):
                    col = (1 if m2 == 9 else 3) + (0 if f == 0.25 else 1)
                    row["paper"] = {"no_line": prow[0], "kl": prow[col], "gain": prow[5]}
                rows.append(row)
        for r in rows:
            if r["d_over_D"] == d:
                r["gain"] = best / k_none
    return rows


# --- benchmark ------------------------------------------------------------------------------

def woodbury_benchmark(ranks=(6, 8, 12), n_orient: int = 20, shape: str = "cubic", d_over_D: float = 1.2,
                       seed: int = 0) -> dict:
    """Time dense Cholesky vs Woodbury for a line term rotated over random orientations.

    C0 is the lattice-frame torus covariance (factorized and inverted once, not timed); for each
    orientation R the shell is rotated, U(R) rebuilt (not timed), and ln N(x; 0, C0 + U W U^T)
    evaluated by dense Cholesky (O(n^3)) and by Woodbury (O(n^2 N)).
    Shells: cubic |m|^2 = 1 (6 vectors), 3 (8), 2 (12).
    """
    tf = load_transfer()
    rng = np.random.default_rng(seed)
    c0 = cov.torus_covariance(shape, d_over_D, None)
    x = np.linalg.cholesky(c0) @ rng.standard_normal(c0.shape[0])
    base = Factorized.of(c0)
    by_rank = {6: 1.0, 8: 3.0, 12: 2.0}
    rots = [H.random_rotation(rng) for _ in range(n_orient)]
    out = {"n": int(c0.shape[0]), "ranks": list(ranks), "n_orientations": n_orient,
           "direct_ms": [], "woodbury_ms": [], "max_rel_diff": [], "shape": shape, "d_over_D": d_over_D}
    for r in ranks:
        m2 = by_rank[r]
        kv = shell_vectors(shape, d_over_D, m2)
        assert kv.shape[0] == r
        k0 = float(np.linalg.norm(kv, axis=1).mean())
        a = line_amplitude(k0, 0.5)
        w = np.full(r, 16.0 * np.pi**2 * a / r)
        us = [cov.mode_vectors(kv @ R.T, tf, LMIN, LMAX, k_norm=k0).T for R in rots]
        t = time.perf_counter()
        direct = [gaussian_loglike(x, c0 + (u * w) @ u.T) for u in us]
        t_direct = (time.perf_counter() - t) / n_orient
        t = time.perf_counter()
        y0 = base.solve(x)
        wood = [woodbury_loglike(x, base, u, w, y0) for u in us]
        t_wood = (time.perf_counter() - t) / n_orient
        out["direct_ms"].append(1e3 * t_direct)
        out["woodbury_ms"].append(1e3 * t_wood)
        out["max_rel_diff"].append(float(np.max(np.abs(np.subtract(direct, wood)) / np.abs(direct))))
    return out


def resonance_bundle(benchmark_orientations: int = 20) -> dict:
    """Pieces of ``resonance.json`` this module owns (data contract in docs/ARCHITECTURE.md)."""
    shells = {f"{s}_{d:g}": [{"m2": sh["m2"], "k_chi": sh["k_chi"], "n_vectors": sh["n_vectors"], "index": sh["index"],
                              "grid_index": sh["grid_index"]} for sh in list_shells(s, d)]
              for s in PAPER_SEARCH["grid"]["shapes"] for d in PAPER_SEARCH["grid"]["d_over_D"]}
    return {
        "kl_table": resonance_table(),
        "shells": shells,
        "n_shells": sum(len(v) for v in shells.values()),
        "woodbury_benchmark": woodbury_benchmark(n_orient=benchmark_orientations),
        "single_field": PAPER_SINGLE_FIELD,
        "exits_3d": PAPER_EXITS_3D["rows"],
        "ceiling": PAPER_CEILING,
        "search_paper": PAPER_SEARCH,
        "source": "[TOP] Secs. 11.6-11.8",
    }
