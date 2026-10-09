"""Lattice geometry and the Casimir functional Z(2), [TOP] Sec. 6, Theorem E, Prop. E.2, Fig. 2.

Lattice bases are 3x3 arrays whose **rows** are basis vectors; a lattice vector is
``g = n @ B`` for an integer row vector ``n``. The dual lattice (with the 2 pi of physics) has
basis ``2 pi inv(B).T``; the Ewald sum below uses ``inv(B).T`` (no 2 pi).

The Epstein zeta function Z(2) = sum_{g != 0} |g|^-4 is evaluated at unit covolume with the
Ewald representation of [TOP] Prop. E.2 (split at t = 1)::

    Z = pi^2 [ 3/2 + sum_{g != 0} phi1(|g|^2) + sum_{k != 0} phi2(|k|^2) ],
    phi1(x) = int_1^inf t e^{-pi x t} dt = e^{-pi x} (1/(pi x) + 1/(pi x)^2),
    phi2(y) = int_1^inf u^{-3/2} e^{-pi y u} du = 2 e^{-pi y} - 2 pi sqrt(y) erfc(sqrt(pi y)).

For covolume V the bracket becomes ``sum phi1 + (2 + sum phi2)/V - 1/2``; this general form
is what we differentiate with respect to a Gram matrix for the shape flow.

Shape derivatives use the volume-preserving deformation ``Lambda -> e^{S/2} Lambda`` with S
symmetric traceless acting on lattice vectors as column vectors (equivalently ``B -> B e^{S/2}``
on the row basis, Gram ``h -> B e^S B^T``), in coordinates with respect to a
Frobenius-orthonormal basis of the 5-dimensional space of traceless symmetric matrices. With
this convention the analytic Ewald Hessian reproduces the certified eigenvalues of Prop. E.2
and the linearized rate ``kappa c_E * 0.566`` of Prop. H.3, with one exception: the BCC
triple is 3.39504, not the printed 3.3953. Every term obeys sum_a d^2|g|^-4/dS_a^2 =
(2/3)|g|^-4, so tr Hess = 2 Z / 3 exactly; the printed BCC spectrum violates this by 8e-4
while ours (and the printed SC and FCC spectra) satisfy it.
"""

from __future__ import annotations

import itertools
import math
from functools import lru_cache

import numpy as np
from scipy.special import erfc

from ..constants import HESSIAN_CERTIFIED, Z2_CERTIFIED

PI = math.pi
#: Ewald cutoff on |g| and |k| at unit covolume ([TOP] Prop. E.2 sums |g|, |k| <= 4.5 exactly)
EWALD_RADIUS = 5.0

#: certified gaps, [TOP] Prop. E.2(c)
GAPS_CERTIFIED = {"BCC-FCC": 0.009748764250334, "SC-BCC": 0.56043571417364}

SHAPES = ("SC", "BCC", "FCC")
_ALIASES = {"SC": "SC", "CUBIC": "SC", "SIMPLE CUBIC": "SC", "BCC": "BCC", "FCC": "FCC"}


# --- bases ------------------------------------------------------------------------------

def canonical_shape(shape: str) -> str:
    """Map 'cubic'/'SC'/'bcc'/... to one of 'SC', 'BCC', 'FCC'."""
    try:
        return _ALIASES[shape.strip().upper()]
    except KeyError as exc:
        raise ValueError(f"unknown lattice shape {shape!r}") from exc


def cubic_basis(shape: str) -> np.ndarray:
    """Primitive basis (rows) of SC, BCC or FCC in the conventional cube frame (cube side 1)."""
    s = canonical_shape(shape)
    if s == "SC":
        return np.eye(3)
    if s == "BCC":
        return 0.5 * np.array([[-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]])
    return 0.5 * np.array([[0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]])


def rhombohedral_basis(alpha_deg: float) -> np.ndarray:
    """Three unit vectors with equal pairwise angle alpha, symmetric about (1,1,1).

    alpha = 60 deg is FCC, 90 deg SC, arccos(-1/3) = 109.47 deg BCC; valid for 0 < alpha < 120.
    """
    c = math.cos(math.radians(alpha_deg))
    if not -0.5 < c < 1.0:
        raise ValueError("rhombohedral angle must lie in (0, 120) degrees")
    x = math.sqrt(1.0 - c)
    y = (math.sqrt(1.0 + 2.0 * c) - x) / 3.0
    return x * np.eye(3) + y * np.ones((3, 3))


def basis_from_gram(h: np.ndarray) -> np.ndarray:
    """A basis (rows) with Gram matrix ``h = B B^T`` (Cholesky, lower-triangular rows)."""
    h = np.asarray(h, dtype=float)
    return np.linalg.cholesky(0.5 * (h + h.T))


def gram(B: np.ndarray) -> np.ndarray:
    B = np.asarray(B, dtype=float)
    return B @ B.T


def covolume(B: np.ndarray) -> float:
    return abs(float(np.linalg.det(np.asarray(B, dtype=float))))


def dual_basis(B: np.ndarray, two_pi: bool = True) -> np.ndarray:
    """Dual-lattice basis (rows): ``2 pi inv(B).T`` (or ``inv(B).T`` with ``two_pi=False``)."""
    D = np.linalg.inv(np.asarray(B, dtype=float)).T
    return 2.0 * PI * D if two_pi else D


def normalize(B: np.ndarray, to: str = "covolume") -> np.ndarray:
    """Rescale to unit covolume (``to='covolume'``) or unit shortest vector (``to='shortest'``)."""
    B = np.asarray(B, dtype=float)
    if to == "covolume":
        return B / covolume(B) ** (1.0 / 3.0)
    if to == "shortest":
        return B / shortest_vector(B)[1]
    raise ValueError("to must be 'covolume' or 'shortest'")


# --- reduction and enumeration ----------------------------------------------------------

def _gso(B: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = B.shape[0]
    Bs = np.zeros_like(B)
    mu = np.zeros((n, n))
    for i in range(n):
        Bs[i] = B[i]
        for j in range(i):
            mu[i, j] = B[i] @ Bs[j] / (Bs[j] @ Bs[j])
            Bs[i] = Bs[i] - mu[i, j] * Bs[j]
    return Bs, mu


def lll_reduce(B: np.ndarray, delta: float = 0.99) -> tuple[np.ndarray, np.ndarray]:
    """LLL reduction of a row basis. Returns ``(B_red, U)`` with ``B_red = U @ B``, U unimodular."""
    B = np.array(B, dtype=float)
    n = B.shape[0]
    U = np.eye(n, dtype=np.int64)
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            _, mu = _gso(B)
            q = int(round(mu[k, j]))
            if q:
                B[k] -= q * B[j]
                U[k] -= q * U[j]
        Bs, mu = _gso(B)
        if Bs[k] @ Bs[k] >= (delta - mu[k, k - 1] ** 2) * (Bs[k - 1] @ Bs[k - 1]):
            k += 1
        else:
            B[[k - 1, k]] = B[[k, k - 1]]
            U[[k - 1, k]] = U[[k, k - 1]]
            k = max(k - 1, 1)
    return B, U


def minkowski_reduce(B: np.ndarray, rtol: float = 1e-12) -> tuple[np.ndarray, np.ndarray]:
    """Greedy (Minkowski, dimension 3) reduction: LLL, then b_k -> b_k + sum c_j b_j, c in {0, +-1}.

    In dimension <= 4 the greedy fixed point is Minkowski-reduced. Returns ``(B_red, U)``
    with ``B_red = U @ B``, rows sorted by length.
    """
    B, U = lll_reduce(B)
    coeffs = [c for c in itertools.product((-1, 0, 1), repeat=2) if c != (0, 0)]
    while True:
        order = np.argsort(np.einsum("ij,ij->i", B, B), kind="stable")
        B, U = B[order], U[order]
        changed = False
        for k in range(3):
            j1, j2 = [j for j in range(3) if j != k]
            best = B[k] @ B[k]
            pick = None
            for c1, c2 in coeffs:
                v = B[k] + c1 * B[j1] + c2 * B[j2]
                if v @ v < best * (1.0 - rtol):
                    best, pick = v @ v, (c1, c2)
            if pick is not None:
                B[k] = B[k] + pick[0] * B[j1] + pick[1] * B[j2]
                U[k] = U[k] + pick[0] * U[j1] + pick[1] * U[j2]
                changed = True
                break
        if not changed:
            return B, U


@lru_cache(maxsize=256)
def _int_box(n0: int, n1: int, n2: int) -> np.ndarray:
    """Nonzero integer vectors of the box |n_i| <= N_i, one of each +-pair (lexicographic half)."""
    r = [np.arange(-n, n + 1) for n in (n0, n1, n2)]
    g = np.stack(np.meshgrid(*r, indexing="ij"), axis=-1).reshape(-1, 3)
    first = np.where(g[:, 0] != 0, g[:, 0], np.where(g[:, 1] != 0, g[:, 1], g[:, 2]))
    return np.ascontiguousarray(g[first > 0]).astype(float)


def _box_for(radius: float, inv_diag: np.ndarray) -> tuple[int, int, int]:
    return tuple(int(math.floor(radius * math.sqrt(max(d, 0.0)) + 1e-9)) for d in inv_diag)


def lattice_vectors(B: np.ndarray, radius: float, include_zero: bool = False,
                    return_coeffs: bool = False):
    """All lattice vectors with ``0 < |g| <= radius`` (both signs), sorted by length.

    With ``return_coeffs`` also returns the integer coefficient rows ``n`` (``g = n @ B``).
    """
    B = np.asarray(B, dtype=float)
    hinv = np.linalg.inv(B @ B.T)
    half = _int_box(*_box_for(radius, np.diag(hinv)))
    n = np.concatenate([half, -half]) if half.size else np.zeros((0, 3))
    g = n @ B
    r2 = np.einsum("ij,ij->i", g, g)
    keep = r2 <= radius**2 * (1.0 + 1e-12)
    n, g, r2 = n[keep], g[keep], r2[keep]
    order = np.argsort(r2, kind="stable")
    n, g = n[order], g[order]
    if include_zero:
        n = np.concatenate([np.zeros((1, 3)), n])
        g = np.concatenate([np.zeros((1, 3)), g])
    return (g, n.astype(np.int64)) if return_coeffs else g


def shortest_vector(B: np.ndarray) -> tuple[np.ndarray, float]:
    """A shortest nonzero lattice vector and its length."""
    Bm, _ = minkowski_reduce(B)
    return Bm[0], float(np.linalg.norm(Bm[0]))


def shortest_vectors(B: np.ndarray, rtol: float = 1e-9) -> np.ndarray:
    """All shortest nonzero lattice vectors (kissing configuration), both signs."""
    _, lam = shortest_vector(B)
    return lattice_vectors(B, lam * (1.0 + rtol))


def shells(B: np.ndarray, radius: float, rtol: float = 1e-9) -> list[tuple[float, np.ndarray]]:
    """Lattice vectors with ``0 < |g| <= radius`` grouped into shells of equal length."""
    g = lattice_vectors(B, radius)
    r = np.linalg.norm(g, axis=1)
    out: list[tuple[float, np.ndarray]] = []
    start = 0
    for i in range(1, len(r) + 1):
        if i == len(r) or r[i] > r[start] * (1.0 + rtol):
            out.append((float(r[start]), g[start:i]))
            start = i
    return out


# --- Ewald kernels ----------------------------------------------------------------------

def phi1(x):
    """int_1^inf t e^{-pi x t} dt."""
    a = PI * np.asarray(x, dtype=float)
    return np.exp(-a) * (1.0 / a + 1.0 / a**2)


def dphi1(x):
    a = PI * np.asarray(x, dtype=float)
    return -PI * np.exp(-a) * (1.0 / a + 2.0 / a**2 + 2.0 / a**3)


def d2phi1(x):
    a = PI * np.asarray(x, dtype=float)
    return PI**2 * np.exp(-a) * (1.0 / a + 3.0 / a**2 + 6.0 / a**3 + 6.0 / a**4)


def phi2(y):
    """int_1^inf u^{-3/2} e^{-pi y u} du = 2 e^{-pi y} - 2 pi sqrt(y) erfc(sqrt(pi y))."""
    y = np.asarray(y, dtype=float)
    return 2.0 * np.exp(-PI * y) - 2.0 * PI * np.sqrt(y) * erfc(np.sqrt(PI * y))


def dphi2(y):
    y = np.asarray(y, dtype=float)
    return -PI * erfc(np.sqrt(PI * y)) / np.sqrt(y)


def d2phi2(y):
    y = np.asarray(y, dtype=float)
    a = PI * y
    return PI**2 * (np.exp(-a) / a + 0.5 * math.sqrt(PI) * a**-1.5 * erfc(np.sqrt(a)))


# --- Z(2) -------------------------------------------------------------------------------

def epstein_z2_gram(h: np.ndarray, radius: float = EWALD_RADIUS,
                    grad: bool = False):
    """Z(h) = sum_{n != 0} (n^T h n)^-2 for a Gram matrix h of any covolume (Ewald, t = 1).

    The split is tuned to covolume ~ 1 (``radius`` refers to unit covolume). With ``grad``
    returns ``(Z, M)`` where ``M = dZ/dh`` is the symmetric gradient (``dZ = tr(M dh)``).
    """
    h = np.asarray(h, dtype=float)
    h = 0.5 * (h + h.T)
    hinv = np.linalg.inv(h)
    hinv = 0.5 * (hinv + hinv.T)
    V = math.sqrt(np.linalg.det(h))
    N = _int_box(*_box_for(radius, np.diag(hinv)))
    K = _int_box(*_box_for(radius, np.diag(h)))
    x = np.einsum("ni,ij,nj->n", N, h, N)
    y = np.einsum("ni,ij,nj->n", K, hinv, K)
    mx = x <= radius**2
    my = y <= radius**2
    N, x, K, y = N[mx], x[mx], K[my], y[my]
    s1 = 2.0 * phi1(x).sum()
    s2 = 2.0 * phi2(y).sum()
    Z = PI**2 * (s1 + (2.0 + s2) / V - 0.5)
    if not grad:
        return Z
    A = 2.0 * np.einsum("n,ni,nj->ij", dphi1(x), N, N)
    C = 2.0 * np.einsum("n,ni,nj->ij", dphi2(y), K, K)
    M = PI**2 * (A - hinv @ C @ hinv / V - (2.0 + s2) / (2.0 * V) * hinv)
    return Z, 0.5 * (M + M.T)


def epstein_z2(B: np.ndarray, unit_covolume: bool = True, radius: float = EWALD_RADIUS) -> float:
    """Z(2) = sum_{g != 0} |g|^-4 of the lattice with row basis B.

    ``unit_covolume=True`` (default) rescales to covolume 1 first, as in [TOP] Prop. E.2.
    """
    B = np.asarray(B, dtype=float)
    V = covolume(B)
    Bu, _ = lll_reduce(B / V ** (1.0 / 3.0))
    z = epstein_z2_gram(Bu @ Bu.T, radius)
    return z if unit_covolume else z * V ** (-4.0 / 3.0)


def epstein_z2_direct(B: np.ndarray, radius: float) -> float:
    """Brute-force partial sum of |g|^-4 over 0 < |g| <= radius (for tests; converges slowly)."""
    g = lattice_vectors(np.asarray(B, dtype=float), radius)
    return float(np.sum(np.einsum("ij,ij->i", g, g) ** -2.0))


# --- shape derivatives ------------------------------------------------------------------

def traceless_symmetric_basis() -> np.ndarray:
    """Frobenius-orthonormal basis (5, 3, 3) of symmetric traceless 3x3 matrices.

    Order: two diagonal (E_g) directions, then three off-diagonal (T_2g) shears.
    """
    E = np.zeros((5, 3, 3))
    E[0] = np.diag([1.0, -1.0, 0.0]) / math.sqrt(2.0)
    E[1] = np.diag([1.0, 1.0, -2.0]) / math.sqrt(6.0)
    for a, (i, j) in enumerate([(0, 1), (0, 2), (1, 2)], start=2):
        E[a, i, j] = E[a, j, i] = 1.0 / math.sqrt(2.0)
    return E


def deform(B: np.ndarray, S: np.ndarray) -> np.ndarray:
    """Apply ``Lambda -> e^{S/2} Lambda`` (vectors as columns) to a row basis: ``B e^{S/2}``."""
    S = np.asarray(S, dtype=float)
    w, Q = np.linalg.eigh(0.5 * (S + S.T))
    return np.asarray(B, dtype=float) @ (Q * np.exp(0.5 * w)) @ Q.T


def _shape_terms(B: np.ndarray, radius: float):
    Bu = normalize(np.asarray(B, dtype=float), "covolume")
    Bu = lll_reduce(Bu)[0]
    g = lattice_vectors(Bu, radius)
    k = lattice_vectors(dual_basis(Bu, two_pi=False), radius)
    return Bu, g, k


def shape_gradient(B: np.ndarray, radius: float = EWALD_RADIUS) -> np.ndarray:
    """dZ/dS_a at S = 0 for ``Lambda -> e^{S/2} Lambda`` at unit covolume (5-vector)."""
    _, g, k = _shape_terms(B, radius)
    E = traceless_symmetric_basis()
    gEg = np.einsum("ni,aij,nj->na", g, E, g)
    kEk = np.einsum("ni,aij,nj->na", k, E, k)
    x = np.einsum("ij,ij->i", g, g)
    y = np.einsum("ij,ij->i", k, k)
    return PI**2 * (dphi1(x) @ gEg - dphi2(y) @ kEk)


def shape_hessian(B: np.ndarray, radius: float = EWALD_RADIUS) -> np.ndarray:
    """Analytic Ewald Hessian d^2Z/dS_a dS_b at S = 0 (5x5), unit covolume, [TOP] Prop. E.2."""
    _, g, k = _shape_terms(B, radius)
    E = traceless_symmetric_basis()
    EE = 0.5 * (np.einsum("aij,bjk->abik", E, E) + np.einsum("bij,ajk->abik", E, E))
    out = np.zeros((5, 5))
    for v, f1, f2, sgn in ((g, dphi1, d2phi1, 1.0), (k, dphi2, d2phi2, -1.0)):
        r2 = np.einsum("ij,ij->i", v, v)
        vEv = sgn * np.einsum("ni,aij,nj->na", v, E, v)
        vEEv = np.einsum("ni,abij,nj->nab", v, EE, v)
        out += np.einsum("n,na,nb->ab", f2(r2), vEv, vEv) + np.einsum("n,nab->ab", f1(r2), vEEv)
    return PI**2 * 0.5 * (out + out.T)


def shape_hessian_fd(B: np.ndarray, step: float = 1e-3) -> np.ndarray:
    """Central-difference Hessian of S -> Z(e^{S/2} Lambda) (check of the analytic one)."""
    E = traceless_symmetric_basis()
    Bu = normalize(np.asarray(B, dtype=float), "covolume")
    f = lambda s: epstein_z2(deform(Bu, np.tensordot(s, E, axes=1)))
    H = np.zeros((5, 5))
    z0 = f(np.zeros(5))
    for a in range(5):
        ea = np.eye(5)[a] * step
        H[a, a] = (f(ea) - 2.0 * z0 + f(-ea)) / step**2
        for b in range(a + 1, 5):
            eb = np.eye(5)[b] * step
            H[a, b] = H[b, a] = (f(ea + eb) - f(ea - eb) - f(-ea + eb) + f(-ea - eb)) / (4 * step**2)
    return H


def hessian_eigenvalues(B: np.ndarray) -> np.ndarray:
    return np.linalg.eigvalsh(shape_hessian(B))


# --- results for the dashboard ----------------------------------------------------------

def rhombohedral_curve(alpha_deg=None) -> dict:
    """Z(2) at unit covolume along the rhombohedral family ([TOP] Fig. 2).

    Default grid 40..119.5 deg in 0.5-deg steps; alpha = 120 deg is the degenerate (planar)
    limit where Z diverges.
    """
    a = np.arange(40.0, 119.75, 0.5) if alpha_deg is None else np.asarray(alpha_deg, dtype=float)
    z = np.array([epstein_z2(rhombohedral_basis(x)) for x in a])
    return {"alpha_deg": a.tolist(), "Z": z.tolist()}


def critical_points() -> dict:
    """Z, gradient norm and Hessian spectrum at SC, BCC, FCC with the certified values."""
    out = {}
    for s in SHAPES:
        B = cubic_basis(s)
        ev = hessian_eigenvalues(B)
        out[s] = {
            "Z": float(epstein_z2(B)),
            "Z_paper": Z2_CERTIFIED[s],
            "gradient_norm": float(np.linalg.norm(shape_gradient(B))),
            "hessian": [float(v) for v in ev],
            "hessian_paper": list(HESSIAN_CERTIFIED[s]),
            "computed": True,
            "source": "[TOP] Prop. E.2",
        }
    return out


def to_json(n_alpha: int | None = None) -> dict:
    """The ``bases``, ``rhombohedral`` and ``critical`` pieces of ``lattice.json``.

    Bases are in the conventional cube frame, scaled so the shortest vector has length 1.
    """
    alpha = None if n_alpha is None else np.linspace(40.0, 119.5, n_alpha)
    crit = critical_points()
    gaps = {
        "BCC-FCC": crit["BCC"]["Z"] - crit["FCC"]["Z"],
        "SC-BCC": crit["SC"]["Z"] - crit["BCC"]["Z"],
    }
    return {
        "bases": {s: normalize(cubic_basis(s), "shortest").tolist() for s in SHAPES},
        "rhombohedral": {**rhombohedral_curve(alpha), "angles_deg": {
            "FCC": 60.0, "SC": 90.0, "BCC": math.degrees(math.acos(-1.0 / 3.0))}},
        "critical": crit,
        "gaps": {"computed": gaps, "paper": dict(GAPS_CERTIFIED), "source": "[TOP] Prop. E.2(c)"},
    }
