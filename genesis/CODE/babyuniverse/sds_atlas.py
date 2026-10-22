"""Stage 1 task 10 (P-G1): a conformal atlas for the Schwarzschild-de Sitter exterior and thin-wall worldlines in it.

f(r) = 1 - 2GM/r - eta^2 r^2 with horizons r_b < r_c (and a negative root r_n), surface gravities
kappa_b = f'(r_b)/2, kappa_c = -f'(r_c)/2, tortoise coordinate r* = sum_i A_i ln|r - r_i| with A_i = 1/f'(r_i).
Null coordinates u = t - r*, v = t + r* in every region (same formula on both sides of each horizon).

Compactified null coordinates: sigma = T + x, delta = T - x.  Psi_v(v) is logistic with rate kappa_b for v < 0 and kappa_c
for v > 0, Psi_u(u) with kappa_c for u < 0 and kappa_b for u > 0 (Kruskal-regular at each horizon, C^0 at v = 0, u = 0).
Blocks (chosen so that every horizon is crossed continuously):

    I    (Psi_v,          Psi_u)          static, eps = +1 side           III  (-2 - Psi_v,  2 - Psi_u)   static, beyond the bridge
    BH   (Psi_v,          2 - Psi_u)      r < r_b, future                 WH   (-2 - Psi_v,  Psi_u)       r < r_b, past
    CF   (2 - Psi_v,      Psi_u)          r > r_c, future (beside I)      CP   (Psi_v,      -2 - Psi_u)   r > r_c, past (beside I)
    CF3  (-2 - Psi_v,     4 + Psi_u)      r > r_c, future (beside III)    CP3  (-4 + Psi_v,  2 - Psi_u)   r > r_c, past (beside III)

III is the point reflection of I through the black-hole bifurcation sphere, (x, T) -> (-2 - x, -T), which is the Kruskal map
(U, V) -> (-U, -V); its Killing time is past-directed.  A static wall segment with eps_out = +1 is drawn in I, with
eps_out = -1 in III.  Integration constants of t are fixed at each horizon crossing by continuity of whichever null
coordinate stays finite there (Kruskal continuity).  Correctness check: along every wall worldline, in proper-time
order, sigma and delta must both increase (timelike, future-directed) -- ``causal_check``.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import quad

from .model import ThinWallModel
from . import branches as Br


class SdSExterior:
    def __init__(self, M: float, eta: float, G: float = 1.0):
        roots = np.roots([-eta**2, 0.0, 1.0, -2 * G * M])
        real = sorted(r.real for r in roots if abs(r.imag) < 1e-12)
        if len(real) != 3 or real[1] <= 0:
            raise ValueError("need two horizons (sub-Nariai SdS)")
        self.r_n, self.r_b, self.r_c = real
        self.M, self.eta, self.G = M, eta, G
        fp = lambda r: 2 * G * M / r**2 - 2 * eta**2 * r
        self.kb, self.kc = fp(self.r_b) / 2, -fp(self.r_c) / 2
        self.A = {k: 1 / fp(r) for k, r in (("n", self.r_n), ("b", self.r_b), ("c", self.r_c))}

    def f(self, r):
        return 1 - 2 * self.G * self.M / r - self.eta**2 * r**2

    def rstar(self, r):
        return (self.A["n"] * np.log(abs(r - self.r_n)) + self.A["b"] * np.log(abs(r - self.r_b))
                + self.A["c"] * np.log(abs(r - self.r_c)))

    # compactification
    def Psi_v(self, v):
        v = np.asarray(v, dtype=float)
        k = np.where(v < 0, self.kb, self.kc)
        return np.tanh(0.5 * k * v)

    def Psi_u(self, u):
        u = np.asarray(u, dtype=float)
        k = np.where(u < 0, self.kc, self.kb)
        return np.tanh(0.5 * k * u)

    BLOCKS = {
        "I": lambda pv, pu: (pv, pu), "III": lambda pv, pu: (-2 - pv, 2 - pu),
        "BH": lambda pv, pu: (pv, 2 - pu), "WH": lambda pv, pu: (-2 - pv, pu),
        "CF": lambda pv, pu: (2 - pv, pu), "CP": lambda pv, pu: (pv, -2 - pu),
        "CF3": lambda pv, pu: (-2 - pv, 4 + pu), "CP3": lambda pv, pu: (-4 + pv, 2 - pu),
    }

    def to_diagram(self, block: str, t, r):
        t, r = np.asarray(t, float), np.asarray(r, float)
        rs = np.vectorize(self.rstar)(r)
        u, v = t - rs, t + rs
        s, d = self.BLOCKS[block](self.Psi_v(v), self.Psi_u(u))
        return (s - d) / 2, (s + d) / 2, s, d


def _block_for(seg_name: str, static_side: str) -> str:
    """Map an itinerary region name to a block, given the side (I or III) of the adjacent static segment."""
    if seg_name == "SdS-static(+)":
        return "I"
    if seg_name == "SdS-static(-)":
        return "III"
    if seg_name in ("BH", "WH"):
        return seg_name
    if seg_name == "cosmo-future":
        return "CF" if static_side == "I" else "CF3"
    if seg_name == "cosmo-past":
        return "CP" if static_side == "I" else "CP3"
    raise ValueError(seg_name)


def wall_worldline(model: ThinWallModel, branch, n: int = 400, R_far: float = 6.0, delta: float = 1e-7) -> list:
    """Wall worldline of one branch on the exterior SdS atlas: list of pieces dict(block, R, t, x, T, sigma, delta)
    in proper-time order."""
    ext = SdSExterior(model.M, model.eta, model.G)
    segs = branch.itinerary_out
    statics = [s for s in segs if s.name.startswith("SdS-static")]
    side_default = "I" if (statics and statics[0].eps > 0) else "III"

    def dtdR(R, d):
        rd = np.sqrt(max(-float(model.V(R)), 1e-300))
        return float(model.A_out(R) / model.f_out(R)) / (d * rd)

    horizons = [ext.r_b, ext.r_c]
    pieces = []
    # proper-time ordered list of (Ra, Rb, d, name, eps)
    for s in segs:
        Ra, Rb = s.R_start, s.R_end
        Ra = R_far * ext.r_c if np.isinf(Ra) else (1e-3 * ext.r_b if Ra == 0 else Ra)
        Rb = R_far * ext.r_c if np.isinf(Rb) else (1e-3 * ext.r_b if Rb == 0 else Rb)
        pieces.append(dict(Ra=Ra, Rb=Rb, d=s.direction, name=s.name, eps=s.eps))
    # anchor: the static piece that ends or starts at a turning point gets t = 0 there
    tp = branch.turning_points[0] if branch.turning_points else None
    out = [None] * len(pieces)

    def sample(p, t_start, R_start_eff, R_end_eff):
        Rs = R_start_eff + (R_end_eff - R_start_eff) * (1 - np.cos(np.linspace(0, np.pi, n))) / 2
        ts = [t_start]
        for a, b in zip(Rs[:-1], Rs[1:]):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ts.append(ts[-1] + quad(lambda R: dtdR(R, p["d"]), a, b, limit=200)[0])
        return Rs, np.array(ts)

    def shrink(Ra, Rb):
        """move endpoints off horizons by delta (relative)"""
        def off(R, toward):
            for h in horizons:
                if abs(R - h) < 1e-9 * h:
                    return h + np.sign(toward - h) * delta * h
            return R
        return off(Ra, Rb), off(Rb, Ra)

    # find anchor piece
    k0 = next((i for i, p in enumerate(pieces) if tp is not None and (abs(p["Ra"] - tp) < 1e-9 * tp or abs(p["Rb"] - tp) < 1e-9 * tp)), 0)
    side = {}
    for i, p in enumerate(pieces):
        if p["name"].startswith("SdS-static"):
            side[i] = "I" if p["eps"] > 0 else "III"

    def static_side_near(i):
        for j in sorted(side, key=lambda j: abs(j - i)):
            return side[j]
        return side_default

    # anchor piece: start at its turning-point end with t = 0
    p = pieces[k0]
    Ra, Rb = shrink(p["Ra"], p["Rb"])
    if tp is not None and abs(p["Rb"] - tp) < 1e-9 * tp:          # piece ends at the turning point: integrate backwards
        Rs, ts = sample(dict(p, d=p["d"]), 0.0, Rb, Ra)
        Rs, ts = Rs[::-1], ts[::-1]
    else:
        Rs, ts = sample(p, 0.0, Ra, Rb)
    out[k0] = dict(block=_block_for(p["name"], static_side_near(k0)), R=Rs, t=ts)

    def continue_to(i, j, forward):
        prev = out[i]
        Rh_prev = prev["R"][-1] if forward else prev["R"][0]
        th_prev = prev["t"][-1] if forward else prev["t"][0]
        h = min(horizons, key=lambda hh: abs(hh - Rh_prev))
        if abs(Rh_prev - h) > 10 * delta * h:              # join at a turning point or a flip radius: t is continuous
            q = pieces[j]
            Ra2, Rb2 = shrink(q["Ra"], q["Rb"])
            R_entry, R_exit = (Ra2, Rb2) if forward else (Rb2, Ra2)
            Rs2, ts2 = sample(q, th_prev, R_entry, R_exit)
            if not forward:
                Rs2, ts2 = Rs2[::-1], ts2[::-1]
            out[j] = dict(block=_block_for(q["name"], static_side_near(j)), R=Rs2, t=ts2, matched="t")
            return
        # which null coordinate is finite at the crossing: compare at delta and delta/100 on the known side
        rs1 = ext.rstar(Rh_prev)
        Rq = h + (Rh_prev - h) / 100
        # t at Rq by extrapolating the log divergence: integrate from Rh_prev to Rq
        pp = pieces[i]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tq = th_prev + quad(lambda R: dtdR(R, pp["d"]), Rh_prev, Rq, limit=200)[0]
        rs2 = ext.rstar(Rq)
        du, dv = abs((tq - rs2) - (th_prev - rs1)), abs((tq + rs2) - (th_prev + rs1))
        use_v = dv < du
        c_val = (tq + rs2) if use_v else (tq - rs2)
        q = pieces[j]
        Ra2, Rb2 = shrink(q["Ra"], q["Rb"])
        R_entry = Ra2 if forward else Rb2
        R_exit = Rb2 if forward else Ra2
        R_entry_q = h + (R_entry - h) / 100
        rs_e = ext.rstar(R_entry_q)
        t_entry_q = (c_val - rs_e) if use_v else (c_val + rs_e)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            t_entry = t_entry_q + quad(lambda R: dtdR(R, q["d"]), R_entry_q, R_entry, limit=200)[0]
        Rs2, ts2 = sample(q, t_entry, R_entry, R_exit)
        if not forward:
            Rs2, ts2 = Rs2[::-1], ts2[::-1]
        out[j] = dict(block=_block_for(q["name"], static_side_near(j)), R=Rs2, t=ts2, matched="v" if use_v else "u")

    for j in range(k0 + 1, len(pieces)):
        continue_to(j - 1, j, True)
    for j in range(k0 - 1, -1, -1):
        continue_to(j + 1, j, False)
    for o in out:
        x, T, s, d = ext.to_diagram(o["block"], o["t"], o["R"])
        o.update(x=x, T=T, sigma=s, delta=d)
    return out


def causal_check(pieces: list) -> dict:
    """In proper-time order the worldline must have d sigma > 0 and d delta > 0 everywhere (timelike, future-directed),
    and consecutive pieces must join continuously."""
    worst = np.inf
    jumps = []
    for k, p in enumerate(pieces):
        ds, dd = np.diff(p["sigma"]), np.diff(p["delta"])
        worst = min(worst, float(np.min(np.minimum(ds, dd) / np.maximum(1e-300, np.hypot(ds, dd) + 1e-300))) if len(ds) else np.inf)
        if k > 0:
            q = pieces[k - 1]
            jumps.append(float(np.hypot(p["sigma"][0] - q["sigma"][-1], p["delta"][0] - q["delta"][-1])))
    most_negative = min(float(np.min(np.minimum(np.diff(p["sigma"]), np.diff(p["delta"])))) for p in pieces)
    # increments below ~1e-12 in absolute value are rounding where tanh saturates (r -> 0, r -> infinity)
    return dict(min_normalized_increment=worst, most_negative_increment=most_negative,
                timelike_future=bool(most_negative > -1e-12), max_join_gap=max(jumps) if jumps else 0.0)


def branch_worldlines(model: ThinWallModel) -> dict:
    out = {}
    for b in Br.classify(model):
        if b.kind in ("static-unstable",):
            continue
        out[b.kind] = wall_worldline(model, b)
    return out
