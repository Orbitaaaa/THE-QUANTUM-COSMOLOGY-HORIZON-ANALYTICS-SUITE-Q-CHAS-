"""The quantum atmosphere: where entanglement across a horizon-like sphere lives, [BH] Secs. 10, 15.

Free massless scalar on Srednicki's radial lattice (spacing a, N sites), decomposed in
partial waves. For each l the ground state is Gaussian with X = K^{-1/2}/2 and
P = K^{1/2}/2; entropies follow from the symplectic eigenvalues of X_A P_A.
We measure the mutual information between the ball (sites 1..n) and the exterior beyond
a gap of g cutoff lengths, as a share of the total entanglement 2S.
"""

from __future__ import annotations

import numpy as np


def coupling_matrix(n_sites: int, ell: int) -> np.ndarray:
    """Srednicki's K for partial wave ell (units a = 1)."""
    j = np.arange(1, n_sites + 1, dtype=float)
    diag = ((j + 0.5) ** 2 + (j - 0.5) ** 2 + ell * (ell + 1)) / j**2
    diag[0] = (1.5**2 + ell * (ell + 1)) / 1.0
    off = -((j[:-1] + 0.5) ** 2) / (j[:-1] * j[1:])
    return np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)


def _entropy_from_xp(x: np.ndarray, p: np.ndarray) -> float:
    nu = np.sqrt(np.clip(np.linalg.eigvals(x @ p).real, 0.25, None))
    a, b = nu + 0.5, nu - 0.5
    with np.errstate(divide="ignore", invalid="ignore"):
        s = a * np.log(a) - np.where(b > 1e-300, b * np.log(np.where(b > 1e-300, b, 1.0)), 0.0)
    return float(np.sum(s))


def partial_wave_terms(n_sites: int, n_ball: int, gaps, ell: int) -> tuple[float, np.ndarray]:
    """S(ball) and I(ball : exterior beyond gap g) for one partial wave (without 2l+1)."""
    k = coupling_matrix(n_sites, ell)
    w, v = np.linalg.eigh(k)
    x = (v / np.sqrt(w)) @ v.T * 0.5
    p = (v * np.sqrt(w)) @ v.T * 0.5
    ball = np.arange(n_ball)
    s_ball = _entropy_from_xp(x[np.ix_(ball, ball)], p[np.ix_(ball, ball)])
    mi = []
    for g in gaps:
        ext = np.arange(n_ball + g, n_sites)
        both = np.concatenate([ball, ext])
        s_ext = _entropy_from_xp(x[np.ix_(ext, ext)], p[np.ix_(ext, ext)])
        s_both = _entropy_from_xp(x[np.ix_(both, both)], p[np.ix_(both, both)])
        mi.append(s_ball + s_ext - s_both)
    return s_ball, np.array(mi)


def atmosphere_profile(n_sites: int = 60, n_ball: int = 20, gaps=(0, 1, 2, 3, 5, 8, 12),
                       ell_max: int = 3000) -> dict:
    """Share of the entanglement 2S reaching beyond each gap, summed over partial waves.

    Every partial wave up to ``ell_max`` is summed with weight 2l+1; the remaining tail is
    added by fitting a power law to the last 20 % of terms, as in [BH] Sec. 15. Also
    returns the area-law coefficient S/(n + 1/2)^2 (Srednicki: 0.30).
    """
    gaps = list(gaps)
    ells = np.arange(ell_max + 1)
    s_terms = np.empty(len(ells))
    mi_terms = np.empty((len(ells), len(gaps)))
    for i, ell in enumerate(ells):
        s_terms[i], mi_terms[i] = partial_wave_terms(n_sites, n_ball, gaps, int(ell))
        s_terms[i] *= 2 * ell + 1
        mi_terms[i] *= 2 * ell + 1

    def tail(y: np.ndarray) -> float:
        sel = slice(int(0.8 * len(ells)), None)
        ok = y[sel] > 0
        if ok.sum() < 5:
            return 0.0
        slope, icpt = np.polyfit(np.log(ells[sel][ok]), np.log(y[sel][ok]), 1)
        if slope >= -1.05:
            return 0.0
        # integral of exp(icpt) l^slope from ell_max + 1/2 to infinity
        return float(np.exp(icpt) * (ell_max + 0.5) ** (slope + 1) / -(slope + 1))

    s_tot = s_terms.sum() + tail(s_terms)
    mi_tot = mi_terms.sum(axis=0) + np.array([tail(mi_terms[:, k]) for k in range(len(gaps))])
    share = mi_tot / (2.0 * s_tot)
    return {
        "n_sites": n_sites, "n_ball": n_ball, "ell_max": ell_max,
        "gaps": gaps, "share": share.tolist(), "entropy": float(s_tot),
        "area_coeff": float(s_tot / (n_ball + 0.5) ** 2),
    }


#: [BH] Sec. 10 table (60 sites, ball of 20, l <= 3000)
PAPER_SHARE = {0: 1.0, 1: 0.1551, 2: 0.0330, 3: 0.0149, 5: 0.0055, 8: 0.0022, 12: 0.0010}
