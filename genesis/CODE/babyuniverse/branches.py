"""Classification of classical thin-wall solution branches.

Two ingredients (GLOBAL_GEOMETRY.md):

1. *Branches* of Rdot^2 + V(R) = 0: since V -> -infinity at R -> 0 (M > 0) and
   at R -> infinity, V has 0 or 2 positive roots (or a double root).  Each
   connected component of {V <= 0} traversed in proper time is a branch.

2. *Chart-free itinerary*: the causal region occupied by the wall on each side
   is determined by three locally computable data,
       sign f(R)   (static or dynamical region),
       eps(R)      (which static region, where f > 0),
       sign Rdot   (future or past dynamical region, where f < 0),
   because eps can only change where A = 0, which forces f = -Rdot^2 <= 0,
   i.e. inside a dynamical region.  The itinerary therefore needs no
   coordinates; it is cross-checked against explicit Kruskal / global de
   Sitter integrations in ``geometry.py`` and the tests.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np

from .model import ThinWallModel

INF = float("inf")


@dataclass
class Segment:
    side: str          # "out" or "in"
    label: str         # chart-free label, e.g. "static[1]+", "dyn[0]future"
    name: str          # conventional name: I, II, III, IV (Schwarzschild); N, S, F, P (dS) ...
    R_start: float
    R_end: float
    direction: int     # sign of Rdot on the segment
    eps: int           # orientation sign on the segment


@dataclass
class Branch:
    kind: str
    legs: List[Tuple[float, float, int]]
    turning_points: List[float]
    itinerary_out: List[Segment] = field(default_factory=list)
    itinerary_in: List[Segment] = field(default_factory=list)
    flags: dict = field(default_factory=dict)

    def names_out(self) -> List[str]:
        return [s.name for s in self.itinerary_out]

    def names_in(self) -> List[str]:
        return [s.name for s in self.itinerary_in]

    def to_dict(self) -> dict:
        d = asdict(self)
        return d


# ----------------------------------------------------------------- region names
def _exterior_namer(model: ThinWallModel) -> Callable[[int, bool, int, int], str]:
    hz = model.horizons_out()
    GM = model.G * model.M

    def name(interval: int, static: bool, eps: int, direction: int) -> str:
        if model.eta == 0 and GM > 0:            # Schwarzschild
            if static:
                return "I" if eps > 0 else "III"
            return "II" if direction < 0 else "IV"
        if model.eta == 0 and GM == 0:           # Minkowski
            return "Mink(i0)" if eps > 0 else "Mink(ball)"
        if GM == 0:                              # pure de Sitter exterior
            if static:
                return "dS-static(+)" if eps > 0 else "dS-static(-)"
            return "dS-future" if direction > 0 else "dS-past"
        # Schwarzschild-de Sitter, horizons r_b < r_c (if two)
        if len(hz) == 2:
            if interval == 0:
                return "BH" if direction < 0 else "WH"
            if interval == 1:
                return "SdS-static(+)" if eps > 0 else "SdS-static(-)"
            return "cosmo-future" if direction > 0 else "cosmo-past"
        # over-extremal SdS: f < 0 everywhere
        return "Kantowski-Sachs(" + ("expanding" if direction > 0 else "contracting") + ")"

    return name


def _interior_namer(model: ThinWallModel) -> Callable[[int, bool, int, int], str]:
    def name(interval: int, static: bool, eps: int, direction: int) -> str:
        if model.chi == 0:
            return "Mink-ball" if eps > 0 else "Mink-inverted"
        if static:
            return "N" if eps > 0 else "S"
        return "F" if direction > 0 else "P"

    return name


def _itinerary(side: str, f: Callable, A: Callable, horizons: Sequence[float],
               flip: Optional[float], legs, namer) -> List[Segment]:
    segs: List[Segment] = []
    hz = sorted(h for h in horizons if h is not None and h > 0)
    for (Ra, Rb, d) in legs:
        lo, hi = min(Ra, Rb), max(Ra, Rb)
        cuts = [c for c in hz + ([flip] if flip else []) if lo < c < hi]
        pts = sorted(set([lo, hi] + cuts))
        if d < 0:
            pts = pts[::-1]
        for a, b in zip(pts[:-1], pts[1:]):
            if np.isinf(a) or np.isinf(b):
                fin = b if np.isinf(a) else a
                Rm = max(2 * fin, (max(hz) * 4 if hz else 1.0), 10.0 * fin if fin > 0 else 10.0)
            else:
                Rm = 0.5 * (a + b)
            fv = float(f(Rm))
            static = fv > 0
            eps = int(np.sign(A(Rm)))
            interval = int(sum(1 for h in hz if h < Rm))
            label = (f"static[{interval}]{'+' if eps > 0 else '-'}" if static
                     else f"dyn[{interval}]{'future' if (d < 0) == (interval == 0) else 'past'}")
            nm = namer(interval, static, eps, d)
            if segs and segs[-1].name == nm and segs[-1].label == label and segs[-1].direction == d:
                segs[-1].R_end = b
            else:
                segs.append(Segment(side, label, nm, a, b, d, eps))
    return segs


def classify(model: ThinWallModel) -> List[Branch]:
    """All classical branches for the given parameters, with itineraries."""
    tps = model.turning_points()
    branches: List[Branch] = []
    if model.M == 0:
        R0 = tps[0]
        branches.append(Branch("bounce-M0", [(INF, R0, -1), (R0, INF, +1)], [R0]))
    else:
        Rc = model.M_critical
        if len(tps) >= 2 and abs(model.M - Rc) > 1e-12 * Rc:
            R1, R2 = tps[0], tps[-1]
            branches.append(Branch("bounded", [(0.0, R1, +1), (R1, 0.0, -1)], [R1]))
            branches.append(Branch("bounce", [(INF, R2, -1), (R2, INF, +1)], [R2]))
        elif len(tps) == 0 or (len(tps) == 1 and model.M > Rc):
            branches.append(Branch("monotonic-expanding", [(0.0, INF, +1)], []))
            branches.append(Branch("monotonic-collapsing", [(INF, 0.0, -1)], []))
        else:
            Rs = tps[0]
            branches.append(Branch("static-unstable", [(Rs, Rs, 0)], [Rs]))
    nout = _exterior_namer(model)
    nin = _interior_namer(model)
    for b in branches:
        if b.kind == "static-unstable":
            Rs = b.turning_points[0]
            eo, ei = int(np.sign(model.A_out(Rs))), int(np.sign(model.A_in(Rs)))
            b.itinerary_out = [Segment("out", "static", nout(1, True, eo, 0), Rs, Rs, 0, eo)]
            b.itinerary_in = [Segment("in", "static", nin(0, True, ei, 0), Rs, Rs, 0, ei)]
        else:
            b.itinerary_out = _itinerary("out", model.f_out, model.A_out, model.horizons_out(),
                                         model.flip_radius_out(), b.legs, nout)
            b.itinerary_in = _itinerary("in", model.f_in, model.A_in, [model.horizon_in()],
                                        model.flip_radius_in(), b.legs, nin)
        b.flags = _flags(model, b)
    return branches


def _flags(model: ThinWallModel, b: Branch) -> dict:
    legs = b.legs
    fut_inf = bool(legs) and np.isinf(legs[-1][1])
    past_inf = bool(legs) and np.isinf(legs[0][0])
    from_zero = bool(legs) and legs[0][0] == 0.0 and model.M > 0
    to_zero = bool(legs) and legs[-1][1] == 0.0 and model.M > 0
    names = b.names_out()
    visible_names = {"I", "Mink(i0)", "SdS-static(+)", "dS-static(+)"}
    visible = any(n in visible_names for n in names)
    final_eps_out = b.itinerary_out[-1].eps if b.itinerary_out else None
    eps_at_tp = {}
    for R in b.turning_points:
        eps_at_tp[float(R)] = dict(eps_in=int(np.sign(model.A_in(R))), eps_out=int(np.sign(model.A_out(R))))
    # eternal exterior white hole: present whenever the exterior region (the side
    # into which n points) contains the past singularity of a Schwarzschild(-dS)
    # white hole.  If the wall never enters the white-hole region, the whole WH
    # region lies to the exterior side iff eps_out = -1 throughout (wall in III).
    exterior_contains_WH = None
    if model.M > 0:
        if from_zero:
            exterior_contains_WH = True
        else:
            exterior_contains_WH = all(s.eps < 0 for s in b.itinerary_out
                                       if s.name in ("III", "SdS-static(-)"))
    return dict(future_inflating=fut_inf, past_from_infinity=past_inf,
                emerges_from_singularity=from_zero, ends_in_singularity=to_zero,
                visible_from_parent_infinity=visible, final_eps_out=final_eps_out,
                eps_at_turning_points=eps_at_tp,
                exterior_contains_past_singularity=exterior_contains_WH,
                daughter_causally_disconnected=bool(fut_inf and final_eps_out is not None
                                                    and final_eps_out < 0 and model.eta == 0
                                                    and model.M > 0))


def bgg_regime(model: ThinWallModel) -> str:
    """Regime label for the Schwarzschild-exterior model (BGG 1987)."""
    if model.M == 0:
        return "M=0"
    E, gam = model.E, model.gamma
    gmin, g1 = model.g_min(gam), 2 - 2 * gam
    if abs(E - gmin) < 1e-12 * gmin:
        return "critical"
    if E < gmin:
        return "C: monotonic (no barrier)"
    if E > g1:
        return "A: bounded branch visible (eps_out(R1)=+1)"
    return "B: bounded branch hidden (eps_out(R1)=-1)"


def summarize(model: ThinWallModel) -> dict:
    br = classify(model)
    return dict(params=model.describe(), regime=bgg_regime(model) if model.eta == 0 else "SdS",
                branches=[dict(kind=b.kind, turning_points=b.turning_points,
                               itinerary_out=[(s.name, s.R_start, s.R_end, s.direction, s.eps)
                                              for s in b.itinerary_out],
                               itinerary_in=[(s.name, s.R_start, s.R_end, s.direction, s.eps)
                                             for s in b.itinerary_in],
                               flags=b.flags) for b in br])
