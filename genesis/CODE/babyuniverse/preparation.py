"""Round 6d (P-R6j): the full route from ordinary initial data, in the characteristic solver of ``throat_evolution``.

Data at a stage before any black hole exists:
* the corner sits inside a false-vacuum ball (pure de Sitter, no mass inside: normal, r < l_H);
* along the outgoing ray: de Sitter, then an **imploding shell of ordinary radiation** (classical null fluid with T_vv > 0, the stand-in for a
  collapsing star; NEC-satisfying) of local energy E, then the region outside the shell, then a **light** wall (exterior normal: no
  Einstein-Rosen bridge on the data), then the true vacuum extending outward (asymptotically flat side);
* along the ingoing ray: de Sitter towards the centre (the run ends before the ray reaches r = 0), with an optional outgoing flux N_u (the
  ghost-condensate NEC phase, negative) emitted inside the ball.
E is tuned so that the region outside the shell has the Schwarzschild-de Sitter mass parameter M = M_N (1 - delta): when the shell implodes it
leaves a black hole inside the false vacuum with a near-Nariai static band around it.  Every sphere on the data surface is normal and r increases
outward along the outgoing ray, so no white hole or bridge is required; the NEC holds everywhere except in the optional N_u.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from . import throat_evolution as TE
from .throat_mol import _line


def _ray_mass(a, b, info, h, G, Q2, k, pot):
    """Misner-Sharp mass along the outgoing data ray u = 0: r_u from the transport equations dZ/dv = F (exact on the data ray)."""
    ra, sa, pa = a
    RU, SU, PU = _line(np.asarray(ra, float), np.asarray(sa, float), np.asarray(pa, float), len(ra), float(info["ru_b"][0]), 0.0, 0.0, h, G, Q2, k,
                       pot.lam, pot.v, pot.rho_f)
    RV = np.gradient(ra, h)
    return 0.5 * ra * (1 + 4 * RU * RV * np.exp(-2 * np.asarray(sa))), RU, RV


def build(E, l_H=30.0, r0=15.0, eps=4.0, vphi=0.0015, h=0.04, Lu=30.0, Lv=1000.0, v_shell=3.0, w_shell=0.5, wall_after=None,
          wall_frac=0.6, delta_target=4e-3, tu_func=None, c=0.0):
    """Characteristic data for the route; returns (a, b, info, pot, meta).  wall_after: v-distance from the shell to the wall (if None, placed
    at wall_frac of the way from r_N to r_c of the target mass)."""
    Lam = 3 / l_H**2
    rho_f = Lam / (8 * np.pi)
    lam = 4 * rho_f / (eps * vphi**4)
    w = np.sqrt(2) / (np.sqrt(lam) * vphi)
    pot = TE.Potential(lam, vphi, rho_f)
    rN = 1 / np.sqrt(Lam)
    MN = 1 / (3 * np.sqrt(Lam))
    f0 = 1 - Lam * r0**2 / 3
    rv0 = np.sqrt(f0) / 2
    M = MN * (1 - delta_target)
    fM = lambda r: 1 - 2 * M / r - Lam * r**2 / 3
    rc = brentq(fM, rN, 3 * rN)
    Nv = lambda v: E / (np.sqrt(np.pi) * w_shell) * np.exp(-((np.asarray(v, float) - v_shell) / w_shell) ** 2)
    # first pass without the wall to locate r(v) beyond the shell, then put the wall where r reaches the target radius
    rw_target = rN + wall_frac * (rc - rN)
    a0, b0, _ = TE.characteristic_data(0.0, c, r0, -rv0, rv0, h, int(Lv / h), int(Lu / h), phi_v=lambda v: -vphi + 0 * np.asarray(v, float),
                                       phi_u=lambda u: -vphi + 0 * np.asarray(u, float), t_u=tu_func, N_v=Nv)
    ra = a0[0]
    j = np.where(ra >= rw_target)[0]
    v_wall = (j[0] * h) if len(j) else Lv * 0.9
    phi_v = lambda v: vphi * np.tanh((np.asarray(v, float) - v_wall) / w)
    a, b, info = TE.characteristic_data(0.0, c, r0, -rv0, rv0, h, int(Lv / h), int(Lu / h), phi_v=phi_v,
                                        phi_u=lambda u: -vphi + 0 * np.asarray(u, float), t_u=tu_func, N_v=Nv)
    m, RU, RV = _ray_mass(a, b, info, h, 1.0, 0.0, c / (12 * np.pi), pot)
    vgrid = np.arange(len(a[0])) * h
    sel = (vgrid > v_shell + 4 * w_shell) & (vgrid < v_wall - 4 * w)
    Cpar = 2 * m - Lam * a[0] ** 3 / 3
    M_out = float(np.median(Cpar[sel]) / 2) if sel.any() else np.nan
    meta = dict(M_N=MN, r_N=rN, v_wall=float(v_wall), M_beyond_shell=M_out, delta=float(1 - M_out / MN), rv_min=float(np.min(RV)),
                normal_everywhere=bool(np.all(RU < 0) and np.all(RV > 0)), r_c=float(rc), wall_found=bool(len(j)))
    return a, b, info, pot, meta


def tune(delta_target=4e-3, lo=6.0, hi=8.0, n=30, **kw):
    """Largest shell energy E (closest to Nariai) for which the data stay normal everywhere with the wall inside the band and
    M <= M_N (1 - delta_target).  With the shell at v_shell = 4.5 the normality constraint binds first, at delta ~ 0.009."""
    for _ in range(n):
        E = 0.5 * (lo + hi)
        meta = build(E, delta_target=delta_target, **kw)[4]
        ok = meta["normal_everywhere"] and meta["wall_found"] and np.isfinite(meta["delta"]) and meta["delta"] >= delta_target
        if ok:
            lo = E
        else:
            hi = E
    return lo, build(lo, delta_target=delta_target, **kw)[4]


def run(E, N_amp=0.0, u_p=1.0, w_p=0.1, c=0.0, **kw):
    tu = (lambda u: -N_amp / (np.sqrt(np.pi) * w_p) * np.exp(-((np.asarray(u, float) - u_p) / w_p) ** 2)) if N_amp else None
    a, b, info, pot, meta = build(E, tu_func=tu, c=c, **kw)
    res = TE.evolve(a, b, 0.0, c, pot, kw.get("h", 0.04), r_inf=1e4, phi_max=10 * pot.v)
    cl = TE.classify(res)
    g = TE.inflated(res, cl, meta["r_N"])
    return res, cl, g, meta
