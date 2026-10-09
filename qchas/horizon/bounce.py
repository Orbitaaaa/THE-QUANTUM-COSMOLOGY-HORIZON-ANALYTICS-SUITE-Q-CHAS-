"""The quantum Big Bang of a closed radiation universe, [QG] Theorems Q, Q', Lambda (Secs. 4-6, 21).

With the radiation fluid as clock (Schutz) and units ell^2 = 8G/pi (ell = 1.596 l_P), the
Wheeler-DeWitt equation is the Schrodinger equation

    i dPsi/dT = (-(1/24) d^2/da^2 + 6 a^2) Psi  on L^2(0, inf),  Psi(0) = 0.

This is a harmonic oscillator of mass 12 and frequency 1 on the half line. With the
Dirichlet condition its eigenfunctions are the odd full-line levels, E_j = 2j + 3/2, so
|Psi|^2 is periodic with period pi. A packet released at rest at a_max with the ground-state
width (|Psi|^2 of width 1/sqrt(24)) is a coherent state g whose odd part solves the problem
exactly:

    Psi(a, T) = [g(a; x, p) - g(-a; x, p)] / sqrt(1 - exp(-12 a_max^2)),   a > 0,
    x = a_max cos T,  p = -12 a_max sin T,  |g|^2 = sqrt(12/pi) exp(-12 (a - x)^2),

which gives the closed form

    |Psi|^2 = N^2 sqrt(12/pi) [exp(-12 (a - |x|)^2) (1 - exp(-24 a |x|))^2
                              + 4 exp(-12 (a^2 + x^2)) sin^2(p a)].

Expectation values are integrals of this density (composite Gauss-Legendre, resolving the
interference fringes of period pi/|p|). The finite-difference spectral method of [QG]
Sec. 21 is implemented as a cross-check and to reproduce the printed tables.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh_tridiagonal, eigvalsh_tridiagonal
from scipy.optimize import brentq, minimize_scalar
from scipy.special import dawsn

from ..constants import BOUNCE_AMIN_ISOTROPIC, BOUNCE_AMIN_MIXMASTER, ELL_OVER_LP

MASS = 12.0  # p_a^2 / 24 = p^2 / (2 m)
OMEGA = 1.0  # 6 a^2 = m omega^2 a^2 / 2
_C2 = math.sqrt(12.0 / math.pi)  # peak of the coherent-state density |g|^2

#: [QG] Sec. 21 Methods, "Quantum universe": finite differences on (0, 40), 8000 points
PAPER_SPECTRUM_GRID = (40.0, 8000)
#: [QG] Sec. 21 Methods, "Gravitons": exact spectral evolution on (0, 24), 6000 points
PAPER_EVOLUTION_GRID = (24.0, 6000)

#: [QG] Theorem Lambda table (p. 8): a_max -> (Delta T, Delta T * a_max, peak <a^-2>, <a>_min)
PAPER_LAMBDA_TABLE = {
    2: (0.25200, 0.5040, 262.0, 0.1644),
    4: (0.12267, 0.4907, 533.0, 0.1625),
    8: (0.06083, 0.4867, 1035.0, 0.1600),
    16: (0.03048, 0.4877, 1917.0, 0.1518),
}
#: [QG] Theorem N discussion (Sec. 6): isotropic <a>_min for the mixmaster comparison runs
PAPER_ISOTROPIC_AMIN_THEOREM_N = {1.5: 0.166, 2: 0.165}
#: [QG] Theorem Q'(ii): peak <a^-2> in units of ell^-2
PAPER_A_INV2_PEAK_Q_PRIME = {2: 262.0, 4: 533.0, 8: 1035.0}
#: [QG] Theorem Lambda(i): Delta T * a_max -> 0.49 for large universes
PAPER_FWHM_PRODUCT_LIMIT = 0.49

#: [QG] Theorem N(iii) and the tables on pp. 10-11 (quantum mixmaster universe, a_max = 1.5)
MIXMASTER = {
    "a_min": BOUNCE_AMIN_MIXMASTER,
    "a_min_range": [0.256, 0.270],
    "a_mean_at_pi": 1.468,
    "a_max": 1.5,
    "norm_drift": 2.7e-13,
    "p_below_1e-3_max": 1e-4,
    "a_inv2_bound": 4798.0,
    "runs": [
        # run, <a>_min, <a> at T = pi, max <a^-2>, max P(a < 1e-3), <beta^2> at T = pi, norm drift
        {"run": "base", "a_min": 0.2573, "a_at_pi": 1.4681, "max_a_inv2": 476.5, "max_p_small": 2.1e-5, "beta2_at_pi": 0.0378, "norm_drift": 1.4e-13},
        {"run": "dt/4", "a_min": 0.2572, "a_at_pi": 1.4680, "max_a_inv2": 476.6, "max_p_small": 2.1e-5, "beta2_at_pi": 0.0378, "norm_drift": 2.9e-13},
        {"run": "36 channels", "a_min": 0.2583, "a_at_pi": 1.4805, "max_a_inv2": 516.4, "max_p_small": 2.8e-5, "beta2_at_pi": 0.0445, "norm_drift": 1.7e-13},
        {"run": "48 channels", "a_min": 0.2701, "a_at_pi": 1.5394, "max_a_inv2": 587.3, "max_p_small": 3.8e-5, "beta2_at_pi": 0.0556, "norm_drift": 2.0e-13},
        {"run": "64 channels", "a_min": 0.2628, "a_at_pi": 1.5160, "max_a_inv2": 650.2, "max_p_small": 4.5e-5, "beta2_at_pi": 0.0492, "norm_drift": 1.5e-14},
        {"run": "80 channels", "a_min": 0.2652, "a_at_pi": 1.5635, "max_a_inv2": 691.0, "max_p_small": 5.3e-5, "beta2_at_pi": 0.0646, "norm_drift": 2.8e-13},
        {"run": "120 channels", "a_min": 0.2578, "a_at_pi": 1.4971, "max_a_inv2": 751.0, "max_p_small": 6.0e-5, "beta2_at_pi": 0.0444, "norm_drift": 4.7e-14},
        {"run": "160 channels", "a_min": 0.2555, "a_at_pi": 1.4676, "max_a_inv2": 776.3, "max_p_small": 6.4e-5, "beta2_at_pi": 0.0868, "norm_drift": 2.7e-13},
        {"run": "a_min = 1e-6", "a_min": 0.2581, "a_at_pi": 1.4598, "max_a_inv2": 491.0, "max_p_small": 2.3e-5, "beta2_at_pi": 0.0287, "norm_drift": 1.3e-13},
        {"run": "1.5x grid", "a_min": 0.2601, "a_at_pi": 1.4821, "max_a_inv2": 413.6, "max_p_small": 1.5e-5, "beta2_at_pi": 0.0300, "norm_drift": 8.2e-14},
        {"run": "a_max = 2", "a_min": 0.2497, "a_at_pi": 1.9674, "max_a_inv2": 553.0, "max_p_small": 2.4e-5, "beta2_at_pi": 0.0233, "norm_drift": 6.9e-14},
    ],
    "saturation": [
        # channels, peak <H_beta>, bound 96 (E + <H_beta>), peak <a^-2>, <a>_min, max P(a < 1e-3), E
        {"channels": 24, "peak_H_beta": 46.3, "bound": 5389.0, "peak_a_inv2": 476.0, "a_min": 0.2573, "max_p_small": 2.1e-5, "E": 9.8751},
        {"channels": 48, "peak_H_beta": 50.7, "bound": 5817.0, "peak_a_inv2": 587.0, "a_min": 0.2701, "max_p_small": 3.8e-5, "E": 9.8554},
        {"channels": 80, "peak_H_beta": 49.2, "bound": 5664.0, "peak_a_inv2": 691.0, "a_min": 0.2652, "max_p_small": 5.3e-5, "E": 9.8487},
        {"channels": 120, "peak_H_beta": 43.2, "bound": 5091.0, "peak_a_inv2": 751.0, "a_min": 0.2578, "max_p_small": 6.0e-5, "E": 9.8462},
        {"channels": 160, "peak_H_beta": 40.1, "bound": 4798.0, "peak_a_inv2": 776.0, "a_min": 0.2555, "max_p_small": 6.4e-5, "E": 9.8452},
    ],
    "hardy_threshold": 1.0 / 96.0,
    "source": "[QG] Theorem N, Theorem H and tables pp. 10-11",
}

#: [QG] Sec. 19.3: Planck 2018 CMB-only spectra against three bounce-imprint templates
PLANCK_BOUNCE_TEST = {
    "data": "Planck 2018 CMB-only TT 2<=l<=2000, TE/EE 30<=l<=2000 (5989 points)",
    "chi2_no_feature": 6166.7,
    "ln_1e10_As": 3.044,
    "ns": 0.9671,
    "templates": [
        {"name": "cutoff", "form": "1 - exp(-(k/k_c)^3.35)", "extra_params": 1, "dchi2_max": 1.2, "p_chance": 0.21},
        {"name": "step", "form": "1 + A exp(-(k/k_B)^2)", "extra_params": 2, "dchi2_max": 1.6, "p_chance": 0.42},
        {"name": "bogoliubov", "form": "|alpha_k + beta_k|^2, |beta_k| = b exp(-(k/k_B)^2)", "extra_params": 3, "dchi2_max": 2.6, "p_chance": 0.34},
    ],
    "dchi2_cutoff_max": 1.2,
    "kc_best_Mpc": 2.8e-4,
    "kc_excluded_above": 4.9e-4,
    "kc_exclusion_dchi2": 3.84,
    "n_null_skies": 300,
    "verdict": "no bounce feature; the framework predicts beta_k = 0 at every scale",
    "source": "[QG] Sec. 19.3",
}


# --- exact solution --------------------------------------------------------------------

def exact_levels(n: int) -> np.ndarray:
    """Dirichlet spectrum E_j = 2j + 3/2 (odd full-line levels), [QG] Theorem Q(ii)."""
    return 2.0 * np.arange(n) + 1.5


def classical_a(T, a_max: float):
    """Classical solution a_max |cos T|, which reaches the singularity a = 0 at T = pi/2."""
    return a_max * np.abs(np.cos(np.asarray(T, dtype=float)))


def density(a, T: float, a_max: float) -> np.ndarray:
    """Exact |Psi(a, T)|^2 on a > 0 for the packet released at rest at a_max (see module doc)."""
    a = np.asarray(a, dtype=float)
    x = a_max * math.cos(T)
    p = -MASS * a_max * math.sin(T)
    ax = abs(x)
    norm2 = 1.0 / -math.expm1(-MASS * a_max**2)
    return norm2 * _C2 * (np.exp(-12.0 * (a - ax) ** 2) * np.expm1(-24.0 * a * ax) ** 2
                          + 4.0 * np.exp(-12.0 * (a * a + x * x)) * np.sin(p * a) ** 2)


def _quadrature(a_max: float, order: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """Composite Gauss-Legendre nodes on [0, a_max + 3] resolving fringes of period pi/(12 a_max)."""
    length = a_max + 3.0
    width = min(0.05, 0.5 * math.pi / (MASS * max(a_max, 1.0)))
    panels = int(math.ceil(length / width))
    xg, wg = leggauss(order)
    edges = np.linspace(0.0, length, panels + 1)
    h = np.diff(edges)
    a = (edges[:-1, None] + 0.5 * (xg[None, :] + 1.0) * h[:, None]).ravel()
    w = (0.5 * wg[None, :] * h[:, None]).ravel()
    return a, w


def moments(T, a_max: float) -> dict:
    """Exact norm, <a> and <a^-2> at the radiation times T (scalar or array)."""
    t = np.atleast_1d(np.asarray(T, dtype=float))
    a, w = _quadrature(a_max)
    out = {"norm": np.empty(t.size), "a_mean": np.empty(t.size), "a_inv2": np.empty(t.size)}
    for i, ti in enumerate(t):
        d = w * density(a, float(ti), a_max)
        out["norm"][i] = d.sum()
        out["a_mean"][i] = np.dot(d, a)
        out["a_inv2"][i] = np.dot(d, 1.0 / (a * a))
    return out


def peak_a_inv2_closed_form(a_max: float) -> float:
    """<a^-2> at the bounce T = pi/2 in closed form.

    4 N^2 C^2 int_0^inf e^{-12 a^2} sin^2(p a) / a^2 da = 2 pi N^2 C^2 [p erf(p/sqrt12)
    + sqrt(12/pi)(e^{-p^2/12} - 1)], with p = 12 a_max and C^2 = sqrt(12/pi).
    """
    p = MASS * a_max
    norm2 = 1.0 / -math.expm1(-MASS * a_max**2)
    return 2.0 * math.pi * norm2 * _C2 * (p * math.erf(p / math.sqrt(12.0)) + _C2 * math.expm1(-p * p / 12.0))


def a_min_closed_form(a_max: float) -> float:
    """<a> at the bounce in closed form: 2 N^2 C^2 (p / 12^{3/2}) F(p / sqrt 12), F = Dawson.

    For large universes it tends to the universal value sqrt(12/pi)/12 = 0.16287 ell.
    """
    p = MASS * a_max
    norm2 = 1.0 / -math.expm1(-MASS * a_max**2)
    return 2.0 * norm2 * _C2 * p / 12.0**1.5 * float(dawsn(p / math.sqrt(12.0)))


#: universal large-universe limits of the exact solution
A_MIN_LIMIT = _C2 / 12.0  # 0.16287 ell
FWHM_PRODUCT_LIMIT = 2.0 * math.sqrt(math.log(2.0) / 12.0)  # 0.4807: e^{-12 x^2} = 1/2 at x = a_max dT/2


def _fwhm(f, t_peak: float, peak: float, span: float) -> tuple[float, float]:
    lo = brentq(lambda t: f(t) - 0.5 * peak, t_peak - span, t_peak, xtol=1e-13)
    hi = brentq(lambda t: f(t) - 0.5 * peak, t_peak, t_peak + span, xtol=1e-13)
    return lo, hi


def time_grid(a_max: float, n_uniform: int = 241, n_bounce: int = 161) -> np.ndarray:
    """One period [0, pi], refined around the bounce where <a^-2> is narrow."""
    half = 2.0 / max(a_max, 1.0)
    t = np.concatenate([np.linspace(0.0, math.pi, n_uniform),
                        np.linspace(math.pi / 2 - half, math.pi / 2 + half, n_bounce)])
    return np.unique(np.round(t, 12))


def bounce_run(a_max: float, T=None) -> dict:
    """Exact <a>(T), <a^-2>(T) over one period and the bounce summary for one universe."""
    t = time_grid(a_max) if T is None else np.asarray(T, dtype=float)
    m = moments(t, a_max)
    f = lambda tt: float(moments(tt, a_max)["a_inv2"][0])  # noqa: E731
    t_peak = math.pi / 2  # |Psi|^2 is symmetric about the bounce
    peak = f(t_peak)
    lo, hi = _fwhm(f, t_peak, peak, min(1.2, 2.0 / a_max))
    a_min = float(moments(t_peak, a_max)["a_mean"][0])
    return {
        "a_max": a_max, "T": t.tolist(), "a_mean": m["a_mean"].tolist(), "a_inv2": m["a_inv2"].tolist(),
        "norm_max_dev": float(np.max(np.abs(m["norm"] - 1.0))),
        "t_bounce": t_peak, "a_min": a_min, "peak_a_inv2": peak, "fwhm_T": hi - lo,
        "fwhm_times_amax": (hi - lo) * a_max, "period": math.pi, "computed": True, "method": "exact",
    }


# --- finite-difference spectral method of [QG] Sec. 21 ----------------------------------

def fd_grid(length: float, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interior points a_i = i h, h = L/(n+1), and the tridiagonal Hamiltonian (Dirichlet at 0, L)."""
    h = length / (n + 1)
    a = h * np.arange(1, n + 1)
    diag = 2.0 / (24.0 * h * h) + 6.0 * a * a
    off = np.full(n - 1, -1.0 / (24.0 * h * h))
    return a, diag, off


def fd_levels(length: float = PAPER_SPECTRUM_GRID[0], n: int = PAPER_SPECTRUM_GRID[1],
              n_eig: int = 2000) -> np.ndarray:
    """Lowest ``n_eig`` finite-difference eigenvalues (bisection)."""
    _, d, e = fd_grid(length, n)
    return eigvalsh_tridiagonal(d, e, select="i", select_range=(0, n_eig - 1))


def spectrum_check(length: float = PAPER_SPECTRUM_GRID[0], n: int = PAPER_SPECTRUM_GRID[1],
                   n_eig: int = 2000, j_show=(0, 1, 2, 3, 4, 5, 10, 20, 50, 100, 200, 500, 1000, 1999)) -> dict:
    """Finite-difference levels against 2j + 3/2 ([QG] Sec. 21: "to 5e-5" holds for the ground level)."""
    w = fd_levels(length, n, n_eig)
    j = [k for k in j_show if k < n_eig]
    exact = exact_levels(n_eig)
    return {"grid": {"L": length, "n": n, "h": length / (n + 1)}, "j": j,
            "E_numeric": [float(w[k]) for k in j], "E_exact": [float(exact[k]) for k in j],
            "abs_error": [float(w[k] - exact[k]) for k in j], "source": "[QG] Sec. 21 Methods"}


class FDEvolution:
    """Spectral evolution of the packet on a finite-difference grid, as in [QG] Sec. 21.

    Only the eigenpairs carrying the packet are computed (MRRR on an index window), and
    the window is checked to hold the packet's norm to 1e-12.
    """

    def __init__(self, a_max: float, length: float = PAPER_EVOLUTION_GRID[0],
                 n: int = PAPER_EVOLUTION_GRID[1], width_sigmas: float = 12.0):
        self.a_max, self.length, self.n = a_max, length, n
        a, d, e = fd_grid(length, n)
        j0 = int(round(0.5 * (MASS * a_max**2 / 2.0)))  # mean odd level: n_full ~ 6 a_max^2
        jw = int(width_sigmas * math.sqrt(6.0) * a_max / 2.0) + 20
        lo, hi = max(0, j0 - jw), min(n - 1, j0 + jw)
        w, v = eigh_tridiagonal(d, e, select="i", select_range=(lo, hi), lapack_driver="stemr")
        psi0 = np.exp(-6.0 * (a - a_max) ** 2)
        psi0 /= np.linalg.norm(psi0)
        c = v.T @ psi0
        self.captured = float(np.dot(c, c))
        if self.captured < 1.0 - 1e-12:
            raise RuntimeError(f"level window misses the packet: captured {self.captured}")
        keep = np.abs(c) > 1e-10
        self.energies, self.coeffs = w[keep], c[keep]
        self.a, self.h = a, length / (n + 1)
        self._vecs = v[:, keep]
        self._a_inv2 = self._vecs.T @ (self._vecs / (a * a)[:, None])
        self._a = self._vecs.T @ (self._vecs * a[:, None])

    def _expect(self, T, op: np.ndarray) -> np.ndarray:
        t = np.atleast_1d(np.asarray(T, dtype=float))
        out = np.empty(t.size)
        for s in range(0, t.size, 256):
            amp = self.coeffs[:, None] * np.exp(-1j * np.outer(self.energies, t[s:s + 256]))
            out[s:s + 256] = np.real(np.einsum("it,ij,jt->t", amp.conj(), op, amp))
        return out

    def density(self, T: float) -> np.ndarray:
        """|Psi(a_i, T)|^2 on the grid points ``self.a`` (normalised as a density in a)."""
        psi = self._vecs @ (self.coeffs * np.exp(-1j * self.energies * T))
        return np.abs(psi) ** 2 / self.h

    def a_mean(self, T) -> np.ndarray:
        return self._expect(T, self._a)

    def a_inv2(self, T) -> np.ndarray:
        return self._expect(T, self._a_inv2)

    def summary(self) -> dict:
        """Peak <a^-2>, <a>_min and the FWHM, located on a fine grid and refined."""
        span = 0.3
        t = math.pi / 2 + np.linspace(-0.05, span, 3501)
        f = lambda tt: float(self.a_inv2(tt)[0])  # noqa: E731
        g = lambda tt: float(self.a_mean(tt)[0])  # noqa: E731
        k = int(np.argmax(self.a_inv2(t)))
        r = minimize_scalar(lambda tt: -f(tt), bounds=(t[max(k - 2, 0)], t[min(k + 2, t.size - 1)]),
                            method="bounded", options={"xatol": 1e-12})
        k2 = int(np.argmin(self.a_mean(t)))
        r2 = minimize_scalar(g, bounds=(t[max(k2 - 2, 0)], t[min(k2 + 2, t.size - 1)]),
                             method="bounded", options={"xatol": 1e-12})
        peak = -float(r.fun)
        lo, hi = _fwhm(f, float(r.x), peak, min(1.2, 1.0 / self.a_max))
        return {"a_max": self.a_max, "grid": {"L": self.length, "n": self.n}, "t_bounce": float(r.x),
                "bounce_delay": float(r.x) - math.pi / 2, "peak_a_inv2": peak, "a_min": float(r2.fun),
                "fwhm_T": hi - lo, "fwhm_times_amax": (hi - lo) * self.a_max,
                "levels_used": int(self.energies.size), "computed": True, "method": "finite differences"}


# --- wave functions, gravitons, export -------------------------------------------------

def snapshots(a_max: float = 6.0, times=(math.pi / 2 - 0.3, math.pi / 2, math.pi / 2 + 0.3),
              a_hi: float = 2.5, n_a: int = 1201) -> dict:
    """Exact |Psi(a, T)|^2 just before, at and after the bounce ([QG] Figure 1, right).

    The paper's panel shows fringes of spacing ~0.044 = pi/(12 a_max) and packets at
    ~a_max sin 0.3 = 1.77, i.e. a_max = 6, the larger universe of its left panel.
    """
    a = np.linspace(0.0, a_hi, n_a)
    return {"a_max": a_max, "a": a.tolist(), "T": [float(t) for t in times],
            "psi2": [density(a, float(t), a_max).tolist() for t in times],
            "computed": True, "source": "[QG] Figure 1 (right)"}


def snapshots_fd(a_max: float = 6.0, grid: tuple[float, int] = PAPER_SPECTRUM_GRID,
                 times=(math.pi / 2 - 0.3, math.pi / 2, math.pi / 2 + 0.3), a_hi: float = 2.5,
                 n_a: int = 1201) -> dict:
    """The same snapshots from the finite-difference evolution on the paper's grid.

    On (0, 40) with 8000 points the level spacing at a_max = 6 is ~1 % low, which delays the
    bounce and separates the T = pi/2 -+ 0.3 packets (1.82 vs 1.68) as in [QG] Figure 1;
    the exact packets coincide at a_max sin 0.3 = 1.77.
    """
    fd = FDEvolution(a_max, *grid)
    a = np.linspace(0.0, a_hi, n_a)
    psi2 = [np.interp(a, np.concatenate([[0.0], fd.a]), np.concatenate([[0.0], fd.density(float(t))]))
            for t in times]
    return {"a_max": a_max, "grid": {"L": grid[0], "n": grid[1]}, "a": a.tolist(), "T": [float(t) for t in times],
            "psi2": [p.tolist() for p in psi2], "computed": True, "method": "finite differences",
            "source": "[QG] Figure 1 (right), Sec. 21"}


def _sig(x, digits: int = 7):
    """Round floats (recursively through lists and dicts) to ``digits`` significant figures."""
    if isinstance(x, float):
        return float(f"{x:.{digits}g}")
    if isinstance(x, list):
        return [_sig(v, digits) for v in x]
    if isinstance(x, dict):
        return {k: _sig(v, digits) for k, v in x.items()}
    return x


def graviton_h2(a_inv2_peak: float, modes=(3, 4, 5, 10), occupation: int = 0) -> list[float]:
    """[QG] Theorem Q'(ii): <h_n^2> = <v_n^2> <a^-2> with <v_n^2> = (N + 1/2)/n."""
    return [(occupation + 0.5) / n * a_inv2_peak for n in modes]


def to_json(a_values=(1.5, 2, 3, 4, 6, 8, 16), fd_values=(2, 4, 8, 16), fast: bool = False) -> dict:
    """bounce.json of the data contract (docs/ARCHITECTURE.md)."""
    if fast:
        fd_values = tuple(v for v in fd_values if v <= 8)
    runs = []
    fd = {float(v): FDEvolution(float(v)).summary() for v in fd_values}
    for a_max in a_values:
        run = bounce_run(float(a_max))
        run["closed_form"] = {"a_min": a_min_closed_form(float(a_max)),
                              "peak_a_inv2": peak_a_inv2_closed_form(float(a_max))}
        key = int(a_max) if float(a_max).is_integer() else a_max
        if key in PAPER_LAMBDA_TABLE:
            dt, prod, peak, amin = PAPER_LAMBDA_TABLE[key]
            run["paper"] = {"fwhm_T": dt, "fwhm_times_amax": prod, "peak_a_inv2": peak, "a_min": amin,
                            "source": "[QG] Theorem Lambda table"}
        if key in BOUNCE_AMIN_ISOTROPIC:
            run.setdefault("paper", {})["a_min_theorem_Q"] = BOUNCE_AMIN_ISOTROPIC[key]
        if a_max in PAPER_ISOTROPIC_AMIN_THEOREM_N:
            run.setdefault("paper", {})["a_min_theorem_N"] = PAPER_ISOTROPIC_AMIN_THEOREM_N[a_max]
        if float(a_max) in fd:
            run["paper_grid"] = fd[float(a_max)]
        runs.append(run)
    spec = spectrum_check()
    spec_lo = spectrum_check(*PAPER_EVOLUTION_GRID, n_eig=1500)
    peaks = {f"{r['a_max']:g}": r["peak_a_inv2"] for r in runs}
    modes = [3, 4, 5, 10]
    ref = next(r for r in runs if r["a_max"] == 2.0)
    return {
        "runs": _sig(runs),
        "spectrum_check": spec,
        "spectrum_check_evolution_grid": spec_lo,
        "snapshots": _sig(snapshots(6.0), 6),
        "snapshots_paper_grid": _sig(snapshots_fd(6.0), 6),
        "snapshots_a_max_2": _sig(snapshots(2.0), 6),
        "graviton": {
            "n": modes, "occupation": 0, "a_max": 2.0, "h2_peak": graviton_h2(ref["peak_a_inv2"], modes),
            "h2_peak_by_a_max": {k: graviton_h2(v, modes) for k, v in peaks.items()},
            "paper_a_inv2_peak": {str(k): v for k, v in PAPER_A_INV2_PEAK_Q_PRIME.items()},
            "formula": "<h_n^2> = (N + 1/2)/n <a^-2>", "source": "[QG] Theorem Q'",
        },
        "limits": {"a_min": A_MIN_LIMIT, "fwhm_times_amax": FWHM_PRODUCT_LIMIT,
                   "paper_fwhm_times_amax": PAPER_FWHM_PRODUCT_LIMIT, "computed": True},
        "mixmaster": MIXMASTER,
        "planck_bounce_test": PLANCK_BOUNCE_TEST,
        "units": {"ell_over_lP": ELL_OVER_LP, "a_min_lP": A_MIN_LIMIT * ELL_OVER_LP,
                  "definition": "ell^2 = 8 G / pi", "source": "[QG] Sec. 4"},
        "claims": CLAIMS,
        "notes": NOTES,
    }


NOTES = [
    "Runs are the exact solution (odd coherent-state combination). 'paper_grid' repeats the "
    "computation with the finite differences of [QG] Sec. 21 on (0, 24) with 6000 points, which "
    "reproduces every printed peak <a^-2> and <a>_min of the Theorem Lambda table.",
    "The printed peaks (262, 533, 1035, 1917) are below the exact values (270.7, 565.4, 1154.9, 2333.7): "
    "the grid sum misses an O(h) boundary term at a = 0, where |Psi|^2/a^2 is largest, and second-order "
    "dispersion slows the packet at large a_max (bounce delayed by 0.063 at a_max = 16).",
    "Exact <a>_min tends to the universal sqrt(12/pi)/12 = 0.1629 ell (0.260 l_P); the printed decrease "
    "to 0.1518 at a_max = 16 is the same discretisation effect. Delta T * a_max tends to 0.481.",
    "The paper's Figure 1 (right) is a_max = 6 on the (0, 40), 8000-point grid: see snapshots_paper_grid.",
]


CLAIMS = [
    {"statement": "The Wheeler-DeWitt operator is self-adjoint with Psi(0) = 0; every self-adjoint extension gives unitary evolution through a = 0",
     "status": "Theorem", "source": "[QG] Theorem Q(i)"},
    {"statement": "With Psi(0) = 0 the spectrum is 2j + 3/2 and the universe is exactly periodic with period pi",
     "status": "Theorem", "source": "[QG] Theorem Q(ii)"},
    {"statement": "A packet released at rest bounces at a universal <a>_min ~ 0.16 ell ~ 0.26 l_P",
     "status": "Computation", "source": "[QG] Theorem Q(iii)"},
    {"statement": "Every graviton mode is an exact oscillator of integer frequency n through the bounce; h_n stays finite",
     "status": "Theorem; <a^-2> computed", "source": "[QG] Theorem Q'"},
    {"statement": "The bounce lasts Delta T ~ 0.49 ell / a_max of conformal time; gradients are suppressed as (0.49 ell / lambda_max)^2",
     "status": "Theorem; Delta T computed", "source": "[QG] Theorem Lambda"},
    {"statement": "The quantum mixmaster universe bounces at <a>_min = 0.256 (+-3 %) and re-expands",
     "status": "Theorem (i-ii); computation (iii)", "source": "[QG] Theorem N"},
    {"statement": "Planck 2018 shows no bounce feature: cutoff, step and Bogoliubov templates improve chi^2 by at most 1.2, 1.6 and 2.6",
     "status": "Analysis of public data", "source": "[QG] Sec. 19.3"},
]
