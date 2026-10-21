"""Phase V: the reduced quantum system, its domain of validity, sectors and information.

The reduced theory (QUANTIZATION.md, PHASE_III) uses the wall radius R, its momentum P and the
asymptotic Killing time t_inf as clock, with H = M.  Its Hamiltonian is a real function only on the
static-static domain f_in(R) > 0, f_out(R) > 0 (Phase II).  This module quantifies, branch by branch:

* where each classical branch leaves the static-static domain;
* how much asymptotic time it takes to get there, finite (an edge of the reduced phase space
  reached in finite time, where a quantum boundary condition is needed) or infinite (a genuine
  asymptotic channel);
* how the canonical momentum behaves there;
* the WKB validity parameter of the area-shift constraint, eps = hbar G |d theta/du|, in Planck units;
* information bookkeeping for the transition (channel entropy vs Bekenstein-Hawking entropy) and the
  competition with Hawking evaporation.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from .model import ThinWallModel

M_PL_KG = 2.176e-8          # Planck mass
T_PL_S = 5.391e-44          # Planck time
L_PL_M = 1.616e-35          # Planck length


def static_domain(model: ThinWallModel) -> dict:
    """Radii bounding the static-static domain along the parent (bounded) and daughter (bounce)
    branches of a regime-A model: parent R in (r_h, R1], daughter R in [R2, R_exit) with
    R_exit the first radius beyond R2 where f_in or f_out vanishes."""
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]
    r_h = model.horizons_out()[0]
    cands = []
    if model.chi > 0:
        cands.append(1 / model.chi)
    for rh in model.horizons_out():
        if rh > R2:
            cands.append(rh)
    R_exit = min(cands) if cands else np.inf
    return dict(R1=float(R1), R2=float(R2), parent_edge=float(r_h), daughter_edge=float(R_exit),
                parent_edge_kind="exterior black-hole horizon f_out = 0",
                daughter_edge_kind="interior de Sitter horizon f_in = 0" if model.chi > 0 and R_exit == 1 / model.chi
                else "exterior horizon")


def _dt_dR(model: ThinWallModel, R):
    """|d t_inf / dR| along a branch: (|A_out| / |f_out|) / sqrt(-V)."""
    return np.abs(model.A_out(R)) / (np.abs(model.f_out(R)) * np.sqrt(np.maximum(-model.V(R), 1e-300)))


def asymptotic_time(model: ThinWallModel, branch: str, cutoffs=(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)) -> dict:
    """Asymptotic time from the turning point to (edge - delta * scale) for decreasing delta.
    Divergence (log growth) signals an asymptotic channel; convergence an edge reached in finite time."""
    d = static_domain(model)
    out = []
    if branch == "parent":
        a, edge = d["R1"], d["parent_edge"]
        scale = a - edge
        for c in cutoffs:
            b = edge + c * scale
            # substitution R = a - (a - b) sin^2(s/2) removes the turning-point square root
            val = _quad_turning(model, a, b)
            out.append((c, val))
    else:
        a, edge = d["R2"], d["daughter_edge"]
        scale = edge - a
        for c in cutoffs:
            b = edge - c * scale
            out.append((c, _quad_turning(model, a, b)))
    vals = np.array([v for _, v in out])
    increments = np.diff(vals)
    return dict(branch=branch, cutoffs=list(cutoffs), t=vals.tolist(),
                divergent=bool(np.all(increments > 0.5 * increments[0]) if len(increments) else False),
                increments=increments.tolist())


def _quad_turning(model, a, b):
    def g(s):
        R = a + (b - a) * np.sin(s / 2) ** 2
        dR = (b - a) * np.sin(s / 2) * np.cos(s / 2)
        return _dt_dR(model, R) * abs(dR)
    return quad(g, 0, np.pi, limit=400, epsabs=1e-12, epsrel=1e-10)[0]


def momentum_near_edge(model: ThinWallModel, branch: str, cutoffs=(1e-2, 1e-4, 1e-6)) -> list:
    """Lorentzian canonical momentum P = (R/G)[artanh(Rdot/A_out) - artanh(Rdot/A_in)] approaching
    the edge of the static domain (real-angle form; it diverges logarithmically at a horizon)."""
    d = static_domain(model)
    res = []
    for c in cutoffs:
        if branch == "parent":
            R = d["parent_edge"] + c * (d["R1"] - d["parent_edge"])
        else:
            R = d["daughter_edge"] - c * (d["daughter_edge"] - d["R2"])
        rd = np.sqrt(max(-model.V(R), 0.0))
        P = R / model.G * (np.arctanh(rd / model.A_out(R)) - np.arctanh(rd / model.A_in(R)))
        res.append(dict(cutoff=c, R=float(R), P=float(P)))
    return res


def wkb_parameter(model: ThinWallModel, chi_planck: float, n: int = 4001) -> dict:
    """eps(u) = hbar G |d theta/du| on the barrier, theta = arccos C (the shift of u by i hbar G
    must change theta little).  The model is dimensionless (chi = G = 1); in physical units
    u_phys = u / chi^2 and hbar G = l_P^2 = 1, so eps_phys = chi_planck^2 * eps_dimless.
    The turning-point neighbourhoods (where WKB is replaced by Airy-type connection) are
    excluded: 2% of the barrier at each end."""
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]
    R = np.linspace(R1, R2, n)
    fi, fo = model.f_in(R), model.f_out(R)
    C = (fi + fo - model.kappa**2 * R**2) / (2 * np.sqrt(fi * fo))
    th = np.arccos(np.clip(C, -1, 1))
    u = R**2 / 2
    dth = np.gradient(th, u)
    k = int(0.02 * n)
    eps_dimless = float(np.max(np.abs(dth[k:-k])))
    return dict(eps_dimensionless=eps_dimless, eps_physical=eps_dimless * chi_planck**2,
                barrier_area_planck=float(4 * np.pi * R2**2 / chi_planck**2))


def binary_entropy(p: float) -> float:
    """h(p) in nats, stable for tiny p: h ~ p (1 - ln p)."""
    if p <= 0:
        return 0.0
    if p < 1e-8:
        return p * (1 - np.log(p))
    return float(-p * np.log(p) - (1 - p) * np.log1p(-p))


def information_budget(B: float, M_kg: float) -> dict:
    """Channel entropy of the transition (parent vs daughter branch, probability e^{-B}) and the
    Bekenstein-Hawking entropy of the exterior black hole, in nats; ln of the daughter probability."""
    m = M_kg / M_PL_KG
    S_BH = 4 * np.pi * m**2
    lnp = -B
    S_ch = binary_entropy(np.exp(lnp)) if lnp > -700 else float(np.exp(lnp + np.log(1 + B)) if lnp > -745 else 0.0)
    return dict(ln_p_daughter=lnp, log10_p_daughter=lnp / np.log(10), channel_entropy_nats=S_ch,
                log10_channel_entropy=(np.log10(1 + B) + lnp / np.log(10)), S_BH_nats=S_BH)


def evaporation_competition(B: float, M_kg: float) -> dict:
    """Hawking evaporation time t_ev = 5120 pi G^2 M^3 / hbar c^4 against the dynamical time 2GM/c.
    If repeated attempts (a black-hole -> white-hole bounce at r = 0) were allowed, at most
    N = t_ev / (2GM) attempts fit before the black hole is gone; log10(N e^{-B}) bounds the total
    daughter probability."""
    m = M_kg / M_PL_KG
    t_ev = 5120 * np.pi * m**3 * T_PL_S
    t_dyn = 2 * m * T_PL_S
    N = t_ev / t_dyn
    return dict(t_evaporation_s=float(t_ev), t_dynamical_s=float(t_dyn), max_attempts=float(N),
                semiclassical_black_hole=bool(m > 10.0),
                log10_total_probability_bound=float(np.log10(max(N, 1.0)) - B / np.log(10)))
