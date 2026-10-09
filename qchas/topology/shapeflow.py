"""Casimir shape flow of a 3-torus, [TOP] Sec. 6.3, Theorems H.1-H.2, Prop. H.3, Sec. 14.

On unimodular Gram matrices h (det h = 1) the radiation-era shape dynamics of Theorem H.2 is

    du/dtau = -u - kappa [M h]_0,    dh/dtau = h u,    M = dE/dh,  E(h) = Z(h),

with [X]_0 the traceless part (c_E absorbed into kappa). We integrate the equivalent symmetric
form X = h u = dh/dtau:  dX/dtau = X h^-1 X - X - kappa (h M h - tr(Mh)/3 h).
Z and M come from the analytic Ewald sum of :mod:`qchas.topology.lattice`. The basis is
Minkowski-reduced whenever it becomes skewed (h -> U h U^T, X -> U X U^T, U in GL(3, Z)).
``mode='gradient'`` integrates the first-order gradient flow dh/dtau = -h [M h]_0 for the
invariant metric tr u^2 (kappa only rescales its time). The paper's basin count, 190 FCC and
50 BCC of 240 Haar-random starts, is quoted for "gradient flow" ([TOP] Prop. E.2 remark,
Prop. H.3, Sec. 14), so :func:`basin_statistics` defaults to that mode; the damped Theorem H.2
flow (``mode='inertial'``) depends on kappa, which the paper does not state for this run.

Initial lattices are Haar-random ([TOP] Sec. 14 'Shape flow'): h = k diag(e^{2t}) k^T with k
Haar on SO(3) and t in the Cartan plane (sum t = 0) with density prod_{i<j} |sinh(t_i - t_j)|,
truncated to the Riemannian ball |t| <= T_MAX (the paper states no truncation). The push-forward
of this K-bi-invariant ball to SL(3,R)/SL(3,Z) tends to the Haar measure as T_MAX grows
(Eskin-McMullen, [TOP] ref. 34); :func:`siegel_check` tests it against Siegel's mean value
theorem E #{g != 0 : |g| <= r} = 4 pi r^3 / 3. With 8000 samples the counts for r = 0.6-1.2 are
1.5-4.5 % low at T_MAX = 5 (the missing cusp of very skewed lattices) and agree within 1 %
(1 sigma) at T_MAX = 7 and 10. Basin fractions move with it: the gradient-flow FCC share is
84 +- 2 % at T_MAX = 5 (360 flows) and 79.8 +- 1.8 % at T_MAX = 10 (seeds 2026 and 7, 194 and
189 of 240), against the paper's 190 of 240. We use T_MAX = 10; the start bases are reduced as
bases, not as Gram matrices, which would square their condition number (~e^{2 sqrt2 T_MAX}).
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation

from ..constants import GPC, N_BOSONIC_DOF, N_FERMIONIC_DOF
from .lattice import (basis_from_gram, cubic_basis, epstein_z2, epstein_z2_gram, hessian_eigenvalues,
                      lattice_vectors, minkowski_reduce)

#: Ewald cutoff used inside the flow (unit covolume); terms beyond are < 1e-17 relative
FLOW_RADIUS = 3.8
#: restart from a reduced basis when prod|b_i|/covolume exceeds this (reduced bases have <= sqrt 2)
SKEW_MAX = 2.5
#: default truncation of the Haar sampler (Riemannian radius |t| of the KAK ball)
T_MAX = 10.0
#: [TOP] Prop. E.2 remark and Prop. H.3: 240 Haar-random flows end at FCC 190 times, BCC 50 times
PAPER_COUNTS = {"FCC": 190, "BCC": 50, "other": 0}
#: quoted from [TOP] Prop. H.3
RELAXATION_PAPER = {"kappa_cE_rate_fcc": 0.566, "TL_today": 1e30, "kappa_cE_coeff": 1.04,
                    "TL4_per_dtau": 0.59, "source": "[TOP] Prop. H.3"}
G_STAR_SM = 106.75
T0_EV = 2.35e-4  # [TOP] Prop. H.3
L_TODAY_GPC = 27.3  # [TOP] Prop. H.3 (L >= 0.985 D)


def traceless(X: np.ndarray) -> np.ndarray:
    return X - np.trace(X) / 3.0 * np.eye(3)


def _sym(X: np.ndarray) -> np.ndarray:
    return 0.5 * (X + X.T)


# --- Haar-random initial lattices -------------------------------------------------------

_EA = np.array([1.0, -1.0, 0.0]) / math.sqrt(2.0)
_EB = np.array([1.0, 1.0, -2.0]) / math.sqrt(6.0)


def kak_weight(t: np.ndarray) -> np.ndarray:
    """prod_{i<j} |sinh(t_i - t_j)| for rows t of the Cartan plane."""
    t = np.atleast_2d(t)
    return np.abs(np.sinh(t[:, 0] - t[:, 1]) * np.sinh(t[:, 0] - t[:, 2]) * np.sinh(t[:, 1] - t[:, 2]))


def sample_kak_t(rng: np.random.Generator, n: int, t_max: float = T_MAX) -> np.ndarray:
    """n samples of t (sum t = 0, |t| <= t_max) with density prod |sinh(t_i - t_j)|."""
    ang = np.linspace(0.0, 2.0 * math.pi, 7201)
    wmax = kak_weight(t_max * (np.cos(ang)[:, None] * _EA + np.sin(ang)[:, None] * _EB)).max() * 1.0001
    out = np.zeros((0, 3))
    while len(out) < n:
        m = 4 * (n - len(out)) + 64
        ab = rng.uniform(-t_max, t_max, size=(m, 2))
        ab = ab[np.einsum("ij,ij->i", ab, ab) <= t_max**2]
        t = ab[:, :1] * _EA + ab[:, 1:] * _EB
        acc = rng.uniform(0.0, wmax, size=len(t)) < kak_weight(t)
        out = np.concatenate([out, t[acc]])
    return out[:n]


def orthogonality_defect(h: np.ndarray) -> float:
    """prod_i |b_i| / covolume (1 for an orthogonal basis); large when the basis is skewed."""
    return float(np.sqrt(np.prod(np.diag(h)) / np.linalg.det(h)))


def reduce_gram(h: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Minkowski-reduce a Gram matrix: returns (U h U^T, U) with U in GL(3, Z)."""
    _, U = minkowski_reduce(basis_from_gram(h))
    Uf = U.astype(float)
    return _sym(Uf @ h @ Uf.T), U


def haar_random_grams(rng: np.random.Generator, n: int, t_max: float = T_MAX) -> list[np.ndarray]:
    """n reduced unimodular Gram matrices k diag(e^{2t}) k^T, KAK-Haar with |t| <= t_max.

    The basis B = k diag(e^t) (rows) is Minkowski-reduced as a basis and squared afterwards:
    reducing the Gram matrix itself would square its condition number (~e^{2 sqrt2 t_max}).
    """
    ts = sample_kak_t(rng, n, t_max)
    ks = Rotation.random(n, random_state=rng).as_matrix()
    out = []
    for t, k in zip(ts, ks):
        B = k * np.exp(t)[None, :]
        _, U = minkowski_reduce(B)
        Br = U.astype(float) @ B
        h = _sym(Br @ Br.T)
        out.append(h / np.linalg.det(h) ** (1.0 / 3.0))
    return out


def siegel_check(rng: np.random.Generator, n: int = 2000, t_max: float = T_MAX,
                 radii=(0.8, 1.0, 1.2)) -> dict:
    """Mean number of nonzero lattice points in balls vs Siegel's (4 pi/3) r^3 (Haar test)."""
    hs = haar_random_grams(rng, n, t_max)
    counts = {r: [] for r in radii}
    for h in hs:
        B = basis_from_gram(h)
        g = lattice_vectors(B, max(radii))
        rr = np.linalg.norm(g, axis=1)
        for r in radii:
            counts[r].append(int(np.sum(rr <= r)))
    return {
        "t_max": t_max, "n": n,
        "radius": list(radii),
        "mean_count": [float(np.mean(counts[r])) for r in radii],
        "stderr": [float(np.std(counts[r]) / math.sqrt(n)) for r in radii],
        "siegel": [4.0 * math.pi / 3.0 * r**3 for r in radii],
    }


# --- the flow ---------------------------------------------------------------------------

_BOX_LIMIT = 15**3


def energy(h: np.ndarray, grad: bool = False):
    """Z(h) (and M = dZ/dh) by Ewald, reducing a skewed basis first: M = U^T M_red U."""
    d = np.diag(h)
    di = np.diag(np.linalg.inv(h))
    box = np.prod(2 * np.floor(FLOW_RADIUS * np.sqrt(di)) + 1) + np.prod(2 * np.floor(FLOW_RADIUS * np.sqrt(d)) + 1)
    if box <= _BOX_LIMIT:
        return epstein_z2_gram(h, FLOW_RADIUS, grad=grad)
    hr, U = reduce_gram(h)
    if not grad:
        return epstein_z2_gram(hr, FLOW_RADIUS)
    z, Mr = epstein_z2_gram(hr, FLOW_RADIUS, grad=True)
    Uf = U.astype(float)
    return z, _sym(Uf.T @ Mr @ Uf)


def _rhs_factory(kappa: float, mode: str):
    def rhs(_tau, y):
        h = _sym(y[:9].reshape(3, 3))
        _, M = energy(h, grad=True)
        hMh = _sym(h @ M @ h)
        force = hMh - np.trace(M @ h) / 3.0 * h
        if mode == "gradient":
            return (-force).ravel()
        X = _sym(y[9:].reshape(3, 3))
        dX = X @ np.linalg.solve(h, X) - X - kappa * force
        return np.concatenate([X.ravel(), _sym(dX).ravel()])
    return rhs


@dataclass
class FlowResult:
    tau: np.ndarray
    Z: np.ndarray
    h_end: np.ndarray
    Z_end: float
    end: str
    converged: bool
    grad_norm: float
    u_norm: float
    hessian_end: list[float]
    n_reductions: int
    n_rhs: int
    seconds: float
    meta: dict = field(default_factory=dict)


@lru_cache(maxsize=1)
def _reference_points() -> dict[str, tuple[float, np.ndarray]]:
    return {s: (epstein_z2(cubic_basis(s)), np.sort(hessian_eigenvalues(cubic_basis(s))))
            for s in ("FCC", "BCC", "SC")}


def classify(h: np.ndarray, z_tol: float = 1e-5, hess_tol: float = 2e-3) -> tuple[str, list[float]]:
    """Name the critical point reached: 'FCC', 'BCC', 'SC' (by Z and Hessian) or 'other'."""
    z = epstein_z2_gram(h)
    ev = np.sort(hessian_eigenvalues(basis_from_gram(h)))
    for name, (z_ref, ev_ref) in _reference_points().items():
        if abs(z - z_ref) < z_tol and np.max(np.abs(ev - ev_ref)) < hess_tol:
            return name, ev.tolist()
    return "other", ev.tolist()


def _skew_event(_tau, y) -> float:
    return orthogonality_defect(_sym(y[:9].reshape(3, 3))) - SKEW_MAX


_skew_event.terminal = True
_skew_event.direction = 1


def run_flow(h0: np.ndarray, kappa: float = 1.0, mode: str = "inertial", tau_max: float = 400.0,
             chunk: float = 4.0, tol: float = 1e-7, rtol: float = 1e-7, atol: float = 1e-10,
             record_every: float = 0.25) -> FlowResult:
    """Integrate the Theorem H.2 flow from a unimodular Gram matrix h0 (u = 0 at tau = 0).

    Integration stops whenever the basis becomes skewed (orthogonality defect > SKEW_MAX) and
    restarts from the Minkowski-reduced basis; after every segment h is rescaled to det h = 1
    and u = h^-1 X is projected to its traceless part (both exact invariants of the flow).
    """
    t_start = time.perf_counter()
    h = _sym(np.asarray(h0, dtype=float))
    h = h / np.linalg.det(h) ** (1.0 / 3.0)
    h, _ = reduce_gram(h)
    X = np.zeros((3, 3))
    rhs = _rhs_factory(kappa, mode)
    taus, zs = [0.0], [energy(h)]
    tau, n_red, n_rhs, converged = 0.0, 0, 0, False
    gnorm = unorm = math.inf
    while tau < tau_max:
        t1 = min(tau + chunk, tau_max)
        k0 = math.floor(tau / record_every + 1e-9) + 1
        t_eval = np.arange(k0, math.floor(t1 / record_every + 1e-9) + 1) * record_every
        y0 = h.ravel() if mode == "gradient" else np.concatenate([h.ravel(), X.ravel()])
        sol = solve_ivp(rhs, (tau, t1), y0, method="DOP853", rtol=rtol, atol=atol,
                        t_eval=t_eval, events=_skew_event)
        n_rhs += sol.nfev
        if sol.status == -1:
            break
        for ti, yi in zip(sol.t, np.reshape(sol.y, (len(y0), -1)).T):
            taus.append(float(ti))
            zs.append(float(energy(_sym(yi[:9].reshape(3, 3)))))
        if sol.status == 1:
            tau, y = float(sol.t_events[0][0]), sol.y_events[0][0]
        else:
            tau, y = t1, sol.y[:, -1]
        h = _sym(y[:9].reshape(3, 3))
        X = np.zeros((3, 3)) if mode == "gradient" else _sym(y[9:].reshape(3, 3))
        s = np.linalg.det(h) ** (1.0 / 3.0)
        h, X = h / s, X / s
        X = _sym(X - np.trace(np.linalg.solve(h, X)) / 3.0 * h)
        hr, U = reduce_gram(h)
        if orthogonality_defect(h) > orthogonality_defect(hr) * (1.0 + 1e-6):
            Uf = U.astype(float)
            h, X = hr, _sym(Uf @ X @ Uf.T)
            n_red += 1
        if sol.status == 1:
            continue
        _, M = energy(h, grad=True)
        gnorm = float(np.linalg.norm(traceless(M @ h)))
        unorm = float(np.linalg.norm(np.linalg.solve(h, X)))
        if gnorm < tol and unorm < tol:
            converged = True
            break
    end, ev = classify(h)
    return FlowResult(np.array(taus), np.array(zs), h, float(epstein_z2_gram(h)), end, converged,
                      gnorm, unorm, ev, n_red, n_rhs, time.perf_counter() - t_start,
                      {"kappa": kappa, "mode": mode})


def basin_statistics(n: int = 240, kappa: float = 1.0, mode: str = "gradient", seed: int = 2026,
                     t_max: float = T_MAX, n_traj: int = 24, traj_points: int = 200, **kw) -> dict:
    """Run n Haar-random flows; endpoint counts, fractions and down-sampled (tau, Z) trajectories.

    Defaults follow [TOP]: 240 starts, first-order gradient flow (``kappa`` is then irrelevant
    and reported as None).
    """
    rng = np.random.default_rng(seed)
    hs = haar_random_grams(rng, n, t_max)
    counts = {"FCC": 0, "BCC": 0, "SC": 0, "other": 0}
    trajectories, unconverged = [], 0
    t0 = time.perf_counter()
    for i, h0 in enumerate(hs):
        r = run_flow(h0, kappa=kappa, mode=mode, **kw)
        counts[r.end] += 1
        unconverged += int(not r.converged)
        if i < n_traj:
            idx = np.unique(np.linspace(0, len(r.tau) - 1, min(traj_points, len(r.tau))).astype(int))
            trajectories.append({"tau": r.tau[idx].round(4).tolist(), "Z": r.Z[idx].round(8).tolist(),
                                 "end": r.end, "converged": r.converged})
    if counts["SC"] == 0:
        counts.pop("SC")
    f = counts["FCC"] / n
    return {
        "n": n, "kappa": None if mode == "gradient" else kappa, "mode": mode, "seed": seed, "t_max": t_max,
        "counts": counts, "fcc_fraction": f, "fcc_fraction_err": math.sqrt(f * (1 - f) / n),
        "paper_counts": {"FCC": PAPER_COUNTS["FCC"], "BCC": PAPER_COUNTS["BCC"]},
        "paper_fcc_fraction": 190 / 240,
        "unconverged": unconverged, "seconds": time.perf_counter() - t0,
        "trajectories": trajectories, "computed": True,
    }


# --- Prop. H.3 numbers ------------------------------------------------------------------

def relaxation_numbers() -> dict:
    """Quoted Prop. H.3 values plus our recomputation of the closed-form pieces.

    kappa c_E = 180 (n_F - n_B) / (pi^4 g_* (TL)^4); slowest rate near FCC kappa c_E mu_min with
    mu_min the smallest FCC Hessian eigenvalue; TL today from T0 = 2.35e-4 eV and L = 27.3 Gpc.
    """
    hbar_c_eV_m = 1.973_269_804e-7
    coeff = 180.0 * (N_FERMIONIC_DOF - N_BOSONIC_DOF) / (math.pi**4 * G_STAR_SM)
    mu_min = float(hessian_eigenvalues(cubic_basis("FCC"))[0])
    L_eV = L_TODAY_GPC * GPC / hbar_c_eV_m
    return {
        **RELAXATION_PAPER,
        "ours": {"kappa_cE_coeff": coeff, "fcc_mu_min": mu_min,
                 "TL4_per_dtau": coeff * mu_min, "TL_today": T0_EV * L_eV, "L_today_inv_eV": L_eV},
        "computed": True,
    }


def to_json(stats: dict | None = None, **kw) -> dict:
    """The ``flows`` and ``relaxation`` pieces of ``lattice.json`` (runs the basin flows if needed)."""
    s = basin_statistics(**kw) if stats is None else stats
    flows = {k: s[k] for k in ("n", "counts", "paper_counts", "trajectories")}
    flows.update({k: s[k] for k in ("kappa", "mode", "t_max", "seed", "fcc_fraction",
                                    "fcc_fraction_err", "paper_fcc_fraction", "unconverged")})
    flows.update({"computed": True, "source": "[TOP] Theorem H.2, Prop. H.3, Sec. 14"})
    return {"flows": flows, "relaxation": relaxation_numbers()}
