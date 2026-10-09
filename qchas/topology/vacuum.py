"""Theorem B / B-loc made visible: the local vacuum variance on E1 vs E2, [TOP] Sec. 5.3-5.4, Fig. 1.

For a free scalar of mass m in its ground state on R x R^3/Gamma the equal-time two-point
function is the image sum W(x, y) = sum_gamma W0(|x - gamma y|), with the Minkowski function
W0(r) = m K1(m r) / (4 pi^2 r). The renormalized local variance is

    <phi^2>(x) - <phi^2>_Minkowski = sum_{gamma != 1} W0(|x - gamma x|),

finite because the action is free. Group elements are enumerated as (A, b_A + t) with t in the
translation subgroup T, |t| <= cutoff (a fixed set, independent of x). On E1 every term is
x-independent, so the excess is constant (Theorem B-loc); on E2 the half-turn screws make it
largest on the rotation axes.

Smoothing, as in [TOP] Fig. 1, convolves the field with a Gaussian of width sigma; the
smoothed excess is the same image sum with W0 replaced by its 3D convolution with a Gaussian
of standard deviation sqrt(2) sigma (u - A v is Gaussian with covariance 2 sigma^2 I).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.special import k1

from .manifolds import DeckGroup, deck_group

#: default mass in units of 1/L ([TOP] Fig. 1 does not state m; we take m L = 1)
M_L_DEFAULT = 1.0
#: [TOP] Fig. 1 smoothing kernel width in lattice units
SMOOTHING_SIGMA_PAPER = 0.35


def w0(r, m: float) -> np.ndarray:
    """Minkowski equal-time two-point function of a free scalar of mass m."""
    r = np.asarray(r, dtype=float)
    return m * k1(m * r) / (4.0 * math.pi**2 * r)


def smoothed_w0_table(m: float, sigma: float, r_max: float, n: int = 600) -> CubicSpline:
    """Spline (in log) of W0 convolved with an isotropic 3D Gaussian of std sqrt(2) sigma.

    Radial convolution: (G_s * f)(r) = 1/(sqrt(2 pi) s r) int_0^inf rho f(rho)
    [e^{-(r-rho)^2/2s^2} - e^{-(r+rho)^2/2s^2}] d rho, s = sqrt(2) sigma.
    """
    s = math.sqrt(2.0) * sigma
    r_grid = np.linspace(max(1e-3, 0.05 * sigma), r_max, n)

    def integrand(rho, r):
        if rho == 0.0:
            return 2.0 * r / s**2 * math.exp(-r * r / (2 * s * s)) / (4.0 * math.pi**2)
        g = math.exp(-((r - rho) ** 2) / (2 * s * s)) - math.exp(-((r + rho) ** 2) / (2 * s * s))
        return rho * float(w0(rho, m)) * g

    vals = []
    for r in r_grid:
        lo, hi = max(0.0, r - 12 * s), r + 12 * s
        pts = [p for p in (r,) if lo < p < hi]
        v, _ = quad(integrand, lo, hi, args=(r,), points=pts or None, limit=200, epsabs=0, epsrel=1e-11)
        vals.append(v / (math.sqrt(2 * math.pi) * s * r))
    return CubicSpline(r_grid, np.log(np.asarray(vals)))


def variance_excess(group: DeckGroup, points: np.ndarray, m: float, cutoff: float,
                    kernel=None, chunk: int = 256) -> np.ndarray:
    """sum_{gamma != 1} K(|gamma x - x|) at each point (rows of ``points``), |t| <= cutoff.

    ``kernel`` defaults to W0 of mass m; pass a callable r -> K(r) for a smoothed kernel.
    """
    K = (lambda r: w0(r, m)) if kernel is None else kernel
    P = np.atleast_2d(np.asarray(points, dtype=float))
    ts = group.translation_vectors(cutoff, include_zero=True)
    Tt = group.translations.T
    out = np.zeros(len(P))
    for A, bA in group.reps:
        if np.allclose(A, np.eye(3)):
            t = ts[np.linalg.norm(ts, axis=1) > 0]
            out += float(np.sum(K(np.linalg.norm(t + bA, axis=1))))
            continue
        # the t-sum is T-periodic in D = (A - I) x + b_A: reduce D to the cell around 0
        D = P @ (A - np.eye(3)).T + bA
        c = np.linalg.solve(Tt, D.T).T
        D = (c - np.round(c)) @ Tt.T
        U, inv = np.unique(np.round(D, 12), axis=0, return_inverse=True)
        vals = np.empty(len(U))
        for i in range(0, len(U), chunk):
            d = np.linalg.norm(U[i:i + chunk, None, :] + ts[None, :, :], axis=-1)
            vals[i:i + chunk] = K(d).sum(axis=1)
        out += vals[inv.ravel()]
    return out


def slice_grid(label: str = "E2", m_L: float = M_L_DEFAULT, L: float = 1.0,
               extent: tuple[float, float] = (-0.5, 1.5), n: int = 41, z: float = 0.0,
               cutoff_mL: float = 16.0, smooth_sigma: float | None = None) -> dict:
    """Excess variance on the plane z = const (perpendicular to the screw axis) for E1-E6.

    Coordinates in units of L. ``cutoff_mL`` bounds |t| m (tail ~ e^{-cutoff_mL}).
    """
    g = deck_group(label, L)
    m = m_L / L
    xs = np.linspace(extent[0], extent[1], n) * L
    X, Y = np.meshgrid(xs, xs, indexing="xy")
    P = np.stack([X.ravel(), Y.ravel(), np.full(X.size, z * L)], axis=1)
    cutoff = cutoff_mL / m
    kernel = None
    if smooth_sigma is not None:
        spl = smoothed_w0_table(m, smooth_sigma * L, cutoff + 4.0 * L)
        kernel = lambda r: np.exp(spl(np.asarray(r)))
    ex = variance_excess(g, P, m, cutoff, kernel).reshape(X.shape)
    return {
        "manifold": label, "m_L": m_L, "L": L, "x": (xs / L).tolist(), "y": (xs / L).tolist(), "z": z,
        "excess": np.round(ex, 9).tolist(), "excess_minus_min": np.round(ex - ex.min(), 9).tolist(),
        "min": float(ex.min()), "max": float(ex.max()), "smoothing_sigma": smooth_sigma,
        "cutoff_mL": cutoff_mL, "computed": True, "source": "[TOP] Theorems B, B-loc, Fig. 1",
    }


def rotation_axes_xy(label: str, extent: tuple[float, float] = (-0.5, 1.5), L: float = 1.0) -> list[list[float]]:
    """(x, y) of the screw axes crossing the slice (E2: x, y in (L/2) Z)."""
    if label != "E2":
        raise ValueError("axes listed for E2 only")
    v = np.arange(math.ceil(extent[0] * 2 - 1e-9), math.floor(extent[1] * 2 + 1e-9) + 1) * 0.5 * L
    return [[float(a), float(b)] for a in v for b in v]


def to_json(n: int = 41) -> dict:
    """Heat-map grids for the dashboard: E1 and E2, raw and smoothed as in [TOP] Fig. 1."""
    return {
        "E1": slice_grid("E1", n=n),
        "E2": slice_grid("E2", n=n),
        "E1_smoothed": slice_grid("E1", n=n, smooth_sigma=SMOOTHING_SIGMA_PAPER),
        "E2_smoothed": slice_grid("E2", n=n, smooth_sigma=SMOOTHING_SIGMA_PAPER),
        "E2_axes": rotation_axes_xy("E2"),
        "note": "m L = 1 (not stated in [TOP]); slice z = 0 perpendicular to the screw axis",
    }
