"""Causal structure, trapped surfaces and the Farhi-Guth obstruction (NO_GO_THEOREMS.md).

Round spheres of area radius r in a spherically symmetric spacetime have
null expansions theta_pm proportional to k_pm^a grad_a r.  With
g^{ab} grad_a r grad_b r = f(r):

  f > 0          normal (untrapped) sphere,
  f < 0, grad r past-directed   (r decreases to the future)  -> trapped,
  f < 0, grad r future-directed (r increases to the future)  -> anti-trapped,
  f = 0          marginal (horizon / extremal sphere of a time-symmetric slice).

Energy conditions for the thin-wall model (distributional):
  vacuum energy  T_ab = -rho g_ab            : T(k,k) = 0      (NEC saturated)
  domain wall    T_ab = -sigma h_ab delta(l) : T(k,k) = sigma (n.k)^2 delta >= 0
so NEC holds for sigma >= 0 while SEC is violated by both (the source of
inflation and of wall repulsion).  Penrose's theorem therefore applies.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np

from .model import ThinWallModel
from .branches import classify

TRAPPED = {"II", "BH", "P", "dS-past", "cosmo-past"}
ANTITRAPPED = {"IV", "WH", "F", "dS-future", "cosmo-future"}


def sphere_type(region_name: str) -> str:
    if region_name in TRAPPED:
        return "trapped"
    if region_name in ANTITRAPPED:
        return "anti-trapped"
    return "normal"


def misner_sharp_in(model: ThinWallModel, r):
    """m = (r/2G)(1 - f) = (4 pi/3) rho_in r^3."""
    return np.asarray(r) / (2 * model.G) * (1 - model.f_in(r))


def misner_sharp_out(model: ThinWallModel, r):
    return np.asarray(r) / (2 * model.G) * (1 - model.f_out(r))


def nec_wall(model: ThinWallModel, n_samples: int = 721) -> float:
    """min over null directions of T(k,k)/delta for S_ab = -sigma h_ab.

    In an orthonormal frame (u, n, e) at the wall a null vector is
    k = (1, cos b, sin b); h(k,k) = -1 + sin^2 b = -cos^2 b, hence
    T(k,k) = sigma cos^2 b >= 0 with equality for k tangent to the wall."""
    b = np.linspace(0, 2 * np.pi, n_samples)
    h_kk = -1.0 + np.sin(b) ** 2
    return float(np.min(-model.sigma * h_kk))


@dataclass
class ExtremalSphere:
    kind: str     # "max" or "min"
    r: float
    area: float
    where: str


def time_symmetric_slice(model: ThinWallModel, R_t: float) -> dict:
    """Area profile of the K = 0 slice through a turning point R_t.

    Moving outward from the interior pole: interior r grows from 0; if
    eps_in = -1 it passes the de Sitter equator (maximal sphere, r = 1/chi)
    before reaching the wall.  Beyond the wall r grows if eps_out = +1, or
    decreases to a minimal sphere (Einstein-Rosen bridge, r = 2GM) or to the
    antipodal pole (M = 0, closed slice) if eps_out = -1.
    """
    ein = int(np.sign(model.A_in(R_t)))
    eout = int(np.sign(model.A_out(R_t)))
    ext: List[ExtremalSphere] = []
    if ein < 0 and model.chi > 0:
        rc = 1 / model.chi
        ext.append(ExtremalSphere("max", rc, 4 * np.pi * rc**2, "interior de Sitter equator"))
    if ein > 0 and eout < 0:
        ext.append(ExtremalSphere("max", R_t, 4 * np.pi * R_t**2, "wall"))
    if ein < 0 and eout > 0:
        ext.append(ExtremalSphere("min", R_t, 4 * np.pi * R_t**2, "wall"))
    hz = model.horizons_out()
    closed = False
    if eout < 0:
        if model.M > 0 and model.eta == 0:
            rb = hz[0]
            ext.append(ExtremalSphere("min", rb, 4 * np.pi * rb**2, "Einstein-Rosen bridge"))
        elif model.M == 0:
            closed = True
    else:
        if model.M == 0 and model.eta > 0:
            rc = 1 / model.eta
            ext.append(ExtremalSphere("max", rc, 4 * np.pi * rc**2, "exterior de Sitter equator"))
            closed = True
    has_min = any(e.kind == "min" and e.where == "Einstein-Rosen bridge" for e in ext)
    return dict(R_t=R_t, eps_in=ein, eps_out=eout, extremal=ext, closed=closed,
                regular_laboratory_data=(not has_min and not closed and eout > 0 and ein > 0),
                contains_marginally_trapped_bridge=has_min)


def embedding_profile(model: ThinWallModel, R_t: float, n: int = 400, r_far: float = None):
    """Embedding diagram (r, z) of the time-symmetric slice through R_t
    (Schwarzschild or Minkowski exterior): dz/dl = sqrt(1 - f) > 0 along the
    slice, dr/dl = eps sqrt(f).  z is continuous across the wall; dz/dr kinks.
    Interior: a cap of the round 3-sphere of radius 1/chi; exterior: the
    Flamm paraboloid z_F(r) = 2 sqrt(r_b (r - r_b)) (r_b = 2GM)."""
    ein = int(np.sign(model.A_in(R_t)))
    eout = int(np.sign(model.A_out(R_t)))
    if model.chi > 0:
        rc = 1 / model.chi
        aw = np.arcsin(min(1.0, model.chi * R_t))
        if ein < 0:
            aw = np.pi - aw
        a = np.linspace(0, aw, n)
        rin, zin = rc * np.sin(a), rc * (1 - np.cos(a))
    else:
        rin, zin = np.linspace(0, R_t, n), np.zeros(n)
    zw = float(zin[-1])
    rb = 2 * model.G * model.M
    r_far = r_far or 4 * max(R_t, rb)

    def zF(r):
        return 2 * np.sqrt(rb * np.maximum(r - rb, 0.0))

    if eout > 0:
        er = np.linspace(R_t, r_far, n)
        ez = zw + zF(er) - zF(R_t)
    else:
        r1 = np.linspace(R_t, rb, n)
        z1 = zw + zF(R_t) - zF(r1)
        r2 = np.linspace(rb, r_far, n)
        z2 = zw + zF(R_t) + zF(r2)
        er, ez = np.concatenate([r1, r2]), np.concatenate([z1, z2])
    return dict(interior=(rin, zin), exterior=(er, ez), eps_in=ein, eps_out=eout)


def farhi_guth_audit(model: ThinWallModel) -> list:
    """For every branch: does the interior reach anti-trapped spheres, is the
    spacetime past-singular, and what does the time-symmetric slice contain."""
    rows = []
    for b in classify(model):
        names_in = b.names_in()
        names_out = b.names_out()
        anti_in = any(sphere_type(n) == "anti-trapped" for n in names_in)
        anti_out = any(sphere_type(n) == "anti-trapped" for n in names_out)
        slices = [time_symmetric_slice(model, R) for R in b.turning_points]
        rows.append(dict(kind=b.kind,
                         interior_reaches_antitrapped=anti_in,
                         exterior_wall_in_antitrapped_region=anti_out,
                         future_inflating=bool(b.flags["future_inflating"]),
                         emerges_from_singularity=bool(b.flags["emerges_from_singularity"]),
                         exterior_contains_past_singularity=b.flags["exterior_contains_past_singularity"],
                         slices=[dict(R_t=s["R_t"], eps_in=s["eps_in"], eps_out=s["eps_out"],
                                      extremal=[(e.kind, e.where, e.area) for e in s["extremal"]],
                                      regular_laboratory_data=s["regular_laboratory_data"])
                                 for s in slices]))
    return rows
