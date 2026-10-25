"""Round 6d: generalised entropy of the simulated Nariai-throat daughter (classical sources: S_out = 0, so S_gen = A/4G = pi r^2 in Planck units).

* ``parent_event_horizon``  The parent's event horizon is the last outgoing ray u = u_EH that still reaches the end of the domain (v = L_v) in the
                            untrapped true vacuum with r_v > 0; S_gen along it is pi r(u_EH, v)^2, which must not decrease in v.
* ``daughter_horizons``     For points p = (u*, v*) deep in the inflating daughter (anti-trapped false vacuum, r > factor r_N), the boundary of the
                            past I^-(p) is {v = v*, u <= u*} (generators d/du) and {u = u*, v <= v*} (generators d/dv).  The GSL deficit of each
                            piece is the largest decrease pi (r(x1)^2 - r(x2)^2) with x1 < x2 along the generators.  An observer reaching the
                            daughter's future infinity has a horizon that is the limit of these.
"""
from __future__ import annotations

import numpy as np


def _deficit(r):
    r = np.asarray(r, float)
    r = r[np.isfinite(r)]
    if len(r) < 2:
        return 0.0
    runmax = np.maximum.accumulate(r)
    return float(np.max(np.pi * (runmax**2 - r**2)))


def parent_event_horizon(o: dict) -> dict:
    R, C, P = o["res"]["r"], o["cls"], o["res"]["phi"]
    h = o["res"]["h"]
    Nv = R.shape[1] - 1
    esc = []
    for i in range(R.shape[0]):
        row = R[i]
        if np.isfinite(row[Nv]) and np.isfinite(row[Nv - 1]) and row[Nv] > row[Nv - 1] and P[i, Nv] > 0:
            esc.append(i)
        else:
            break
    if not esc:
        return dict(u_EH=None)
    iEH = esc[-1]
    r = R[iEH]
    return dict(u_EH=float(iEH * h), deficit=_deficit(r), S_start=float(np.pi * r[0] ** 2), S_end=float(np.pi * np.nanmax(r) ** 2))


def daughter_horizons(o: dict, factor: float = 5.0, n_points: int = 12) -> dict:
    R, C, P = o["res"]["r"], o["cls"], o["res"]["phi"]
    h = o["res"]["h"]
    rN = o["r_N"]
    dau = (C == 1) & (P < 0) & (R > factor * rN) & np.isfinite(R)
    idx = np.argwhere(dau)
    if len(idx) == 0:
        return dict(points=0)
    pick = idx[np.linspace(0, len(idx) - 1, n_points).astype(int)]
    out = []
    for i, j in pick:
        d_in = _deficit(R[: i + 1, j])          # v = v*, u increasing
        d_out = _deficit(R[i, : j + 1])         # u = u*, v increasing
        out.append(dict(u=float(i * h), v=float(j * h), deficit_ingoing_piece=d_in, deficit_outgoing_piece=d_out,
                        S_at_corner_of_piece=float(np.pi * R[0, j] ** 2)))
    return dict(points=len(out), max_deficit=float(max(max(x["deficit_ingoing_piece"], x["deficit_outgoing_piece"]) for x in out)),
                S_N=float(np.pi * rN**2), rows=out)


def ending_observer_deficits(o: dict, region: np.ndarray, n_points: int = 12) -> dict:
    """Generalised-entropy decrease (S_out = 0) on the past boundaries of observers whose histories end at the last computed point of a
    row (next to a singularity or crunch): for endpoints p = (u*, v*) in `region` at the end of their rows, the pieces {v = v*, u <= u*}
    and {u = u*, v <= v*}.  Returns the largest decrease and the area at the start of each piece."""
    R = o["res"]["r"]
    h = o["res"]["h"]
    ends = []
    for i in range(R.shape[0]):
        js = np.where(np.isfinite(R[i]))[0]
        if len(js) and region[i, js.max()]:
            ends.append((i, js.max()))
    if not ends:
        return dict(points=0)
    pick = [ends[k] for k in np.linspace(0, len(ends) - 1, min(n_points, len(ends))).astype(int)]
    rows = []
    for i, j in pick:
        rows.append(dict(u=float(i * h), v=float(j * h), r_end=float(R[i, j]), deficit_ingoing=_deficit(R[: i + 1, j]),
                         deficit_outgoing=_deficit(R[i, : j + 1]), S_max_on_boundary=float(np.pi * np.nanmax(np.concatenate([R[: i + 1, j], R[i, : j + 1]])) ** 2)))
    return dict(points=len(rows), max_deficit=float(max(max(x["deficit_ingoing"], x["deficit_outgoing"]) for x in rows)), rows=rows)
