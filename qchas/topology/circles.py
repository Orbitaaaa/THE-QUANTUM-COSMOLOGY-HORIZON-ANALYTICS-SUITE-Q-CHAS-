"""Matched circles and ghost images on a 3-torus, [TOP] Sec. 7, Theorem F, Corollary F.1.

Lengths are in units of the last-scattering diameter D = 2 chi_*, so chi_* = 1/2. A torus is
given by a row basis (any scale; it is rescaled so its shortest vector has length d = d/D) and
an orientation R (rotation, lattice frame -> Galactic frame). For every lattice vector g with
0 < |g| < 2 chi_* there is a circle pair centred on -/+ R g_hat with angular radius
arccos(|g| / 2 chi_*), matched with zero twist by n' = n + g / chi_* (Corollary F.1).
Sky directions use Galactic (l, b) with unit vector (cos b cos l, cos b sin l, sin b).
"""

from __future__ import annotations

import itertools
import math

import numpy as np

from ..constants import D_LSS_GPC, FCC_PREREGISTERED_AXES_GAL
from .lattice import canonical_shape, cubic_basis, lattice_vectors, shortest_vector

CHI_STAR = 0.5  # in units of D


def lb_to_vec(l_deg, b_deg) -> np.ndarray:
    l, b = np.radians(l_deg), np.radians(b_deg)
    return np.stack([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], axis=-1)


def vec_to_lb(v) -> tuple[np.ndarray, np.ndarray]:
    v = np.asarray(v, dtype=float)
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    l = np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360.0
    b = np.degrees(np.arcsin(np.clip(v[..., 2], -1.0, 1.0)))
    return l, b


def _as_direction(d) -> np.ndarray:
    d = np.asarray(d, dtype=float)
    if d.shape == (2,):
        return lb_to_vec(d[0], d[1])
    if d.shape == (3,):
        return d / np.linalg.norm(d)
    raise ValueError("direction must be (l, b) in degrees or a 3-vector")


def scaled_basis(basis: np.ndarray, d_over_D: float) -> np.ndarray:
    """Rescale a row basis so its shortest lattice vector has length d/D."""
    if not (math.isfinite(d_over_D) and d_over_D > 0.0):
        raise ValueError("d_over_D must be a positive number")
    B = np.asarray(basis, dtype=float)
    return B * (d_over_D / shortest_vector(B)[1])


def _rotation(R) -> np.ndarray:
    return np.eye(3) if R is None else np.asarray(R, dtype=float)


def circle_pairs(basis: np.ndarray, d_over_D: float, R=None, chi_star: float = CHI_STAR) -> list[dict]:
    """One entry per +-g with 0 < |g| < 2 chi_*: centres -/+R g_hat, radius arccos(|g|/2chi_*)."""
    Rm = _rotation(R)
    B = scaled_basis(basis, d_over_D)
    g, n = lattice_vectors(B, 2.0 * chi_star * (1.0 - 1e-12), return_coeffs=True)
    out = []
    seen = set()
    for gi, ni in zip(g, n):
        key = tuple(ni)
        if tuple(-ni) in seen:
            continue
        seen.add(key)
        length = float(np.linalg.norm(gi))
        if length >= 2.0 * chi_star:
            continue
        gs = Rm @ gi
        u = gs / length
        lp, bp = vec_to_lb(u)
        lm, bm = vec_to_lb(-u)
        out.append({
            "n": [int(x) for x in ni], "g_lattice": gi.tolist(), "g_sky": gs.tolist(),
            "length_over_D": length,
            "radius_deg": math.degrees(math.acos(length / (2.0 * chi_star))),
            "centre_minus": (-u).tolist(), "centre_plus": u.tolist(),
            "centre_minus_lb": [float(lm), float(bm)], "centre_plus_lb": [float(lp), float(bp)],
        })
    out.sort(key=lambda c: c["length_over_D"])
    return out


def matched_point(n, g, chi_star: float = CHI_STAR, atol: float = 1e-9) -> np.ndarray:
    """Partner direction n' = n + g/chi_* of a direction n on the circle centred at -g_hat."""
    n = _as_direction(n)
    g = np.asarray(g, dtype=float)
    m = n + g / chi_star
    if abs(np.linalg.norm(m) - 1.0) > atol:
        raise ValueError("n does not lie on the circle C_g (|n + g/chi_*| != 1)")
    return m


def circle_points(g, n_points: int = 120, chi_star: float = CHI_STAR) -> tuple[np.ndarray, np.ndarray]:
    """Points on the circle centred at -g_hat and their matched partners (same azimuth phi)."""
    g = np.asarray(g, dtype=float)
    length = np.linalg.norm(g)
    u = g / length
    e1 = np.cross(u, [0.0, 0.0, 1.0] if abs(u[2]) < 0.9 else [1.0, 0.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(u, e1)
    ca = length / (2.0 * chi_star)
    sa = math.sqrt(max(0.0, 1.0 - ca * ca))
    phi = np.linspace(0.0, 2.0 * math.pi, n_points, endpoint=False)
    ring = sa * (np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2)
    return -ca * u + ring, ca * u + ring


def ghost_images(source_direction, source_distance_over_D: float, basis: np.ndarray,
                 d_over_D: float, R=None, horizon_over_D: float = CHI_STAR,
                 max_distance_over_D: float | None = None) -> list[dict]:
    """Every image x + g of a source at x (direction, distance in D) out to ``max_distance``.

    Default ``max_distance_over_D`` is twice the horizon. Each entry gives the observed
    Galactic (l, b), distance (D and Gpc) and whether the image lies inside the horizon; an
    image at the observer (source distance 0) has ``l = b = None``.
    """
    if not source_distance_over_D >= 0.0:
        raise ValueError("source_distance_over_D must be >= 0")
    Rm = _rotation(R)
    x = _as_direction(source_direction) * float(source_distance_over_D)
    rmax = 2.0 * horizon_over_D if max_distance_over_D is None else float(max_distance_over_D)
    B = scaled_basis(basis, d_over_D)
    g, n = lattice_vectors(B, rmax + np.linalg.norm(x), include_zero=True, return_coeffs=True)
    y = x + g @ Rm.T
    r = np.linalg.norm(y, axis=1)
    keep = r <= rmax
    out = []
    for yi, ri, ni in zip(y[keep], r[keep], n[keep]):
        l, b = (float(a) for a in vec_to_lb(yi)) if ri > 1e-12 else (None, None)
        out.append({
            "n": [int(v) for v in ni], "is_source": bool(not np.any(ni)),
            "l": l, "b": b, "distance_over_D": float(ri),
            "distance_Gpc": float(ri * D_LSS_GPC), "inside_horizon": bool(ri <= horizon_over_D),
        })
    out.sort(key=lambda e: e["distance_over_D"])
    return out


def _acute_angles(axes: np.ndarray) -> list[float]:
    angs = set()
    for i, j in itertools.combinations(range(len(axes)), 2):
        c = abs(float(axes[i] @ axes[j]))
        angs.add(round(math.degrees(math.acos(min(1.0, c))), 1))
    return sorted(angs)


def circle_table(shape: str, d_over_D: float) -> dict:
    """Circle pattern of a cubic-family torus in the lattice frame ([TOP] Sec. 7 table, Fig. 4)."""
    s = canonical_shape(shape)
    pairs = circle_pairs(cubic_basis(s), d_over_D)
    shells: list[dict] = []
    for p in pairs:
        if shells and abs(shells[-1]["length_over_D"] - p["length_over_D"]) < 1e-9:
            shells[-1]["n_pairs"] += 1
            shells[-1]["centres"].append(p["centre_plus"])
        else:
            shells.append({"length_over_D": p["length_over_D"], "radius_deg": p["radius_deg"],
                           "n_pairs": 1, "centres": [p["centre_plus"]]})
    first = np.array(shells[0]["centres"]) if shells else np.zeros((0, 3))
    return {
        "shape": s, "d_over_D": d_over_D, "n_pairs": len(pairs),
        "radius_deg": shells[0]["radius_deg"] if shells else None,
        "angles_between_centres_deg": _acute_angles(first),
        "shells": [{**sh, "angles_between_centres_deg": _acute_angles(np.array(sh["centres"]))}
                   for sh in shells],
        "computed": True,
    }


#: [TOP] Sec. 7 table (circle pairs for 0.985 D <= d < D)
CIRCLE_TABLE_PAPER = [
    {"shape": "SC", "pairs": 3, "centres": "±e1, ±e2, ±e3", "radius": "<= 9.9°", "angles_deg": [90.0]},
    {"shape": "BCC", "pairs": 4, "centres": "(±1, ±1, ±1)/√3", "radius": "<= 9.9°", "angles_deg": [70.5]},
    {"shape": "FCC", "pairs": 6, "centres": "(±1, ±1, 0)/√2 and permutations", "radius": "<= 9.9°",
     "angles_deg": [60.0, 90.0]},
    {"shape": "any, d >= D", "pairs": 0, "centres": "no circles; topology visible only statistically (Sec. 8)",
     "radius": None, "angles_deg": []},
]
#: [TOP] Fig. 4 captions at d/D = 0.8
FIG4_PAPER = {"d_over_D": 0.8, "radius_deg": 36.9, "pairs": {"SC": 3, "BCC": 7, "FCC": 6}}


def fcc_axes_lattice_frame() -> np.ndarray:
    """The six <110> directions of FCC (one of each antipodal pair), unit vectors."""
    ax = []
    for i, j in itertools.combinations(range(3), 2):
        for s in (1.0, -1.0):
            v = np.zeros(3)
            v[i], v[j] = 1.0, s
            ax.append(v / math.sqrt(2.0))
    return np.array(ax)


def fcc_preregistered_rotation(axes_gal=FCC_PREREGISTERED_AXES_GAL) -> dict:
    """Best-fit R in SO(3) mapping the FCC <110> axes onto the [TOP] Sec. 9 Galactic axes.

    Axes are unsigned: we minimize sum |R a_i - s_i t_pi(i)|^2 over assignments pi, signs s
    and R (Kabsch). R is unique up to the 24 rotations of the cubic group.
    """
    A = fcc_axes_lattice_frame()
    T = lb_to_vec(*np.array(axes_gal, dtype=float).T)
    perms = np.array(list(itertools.permutations(range(6))))
    signs = np.array(list(itertools.product((1.0, -1.0), repeat=6)))
    Tp = T[perms]  # (720, 6, 3)
    Ts = Tp[:, None, :, :] * signs[None, :, :, None]  # (720, 64, 6, 3)
    H = np.einsum("ik,psil->pskl", A, Ts).reshape(-1, 3, 3)  # sum_i a_i t_i^T
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.transpose(0, 2, 1) @ U.transpose(0, 2, 1)))
    score = S[:, 0] + S[:, 1] + d * S[:, 2]
    best = int(np.argmax(score))
    D = np.diag([1.0, 1.0, d[best]])
    Rm = Vt[best].T @ D @ U[best].T
    p, s = divmod(best, len(signs))
    target = Ts[p, s]
    res = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", A @ Rm.T, target), -1.0, 1.0)))
    return {
        "R": Rm.tolist(), "assignment": perms[p].tolist(), "signs": signs[s].tolist(),
        "residual_rms_deg": float(np.sqrt(np.mean(res**2))), "residual_max_deg": float(res.max()),
        "axes_gal_paper": [list(a) for a in axes_gal], "source": "[TOP] Sec. 9 (best fit, not significant)",
    }


def to_json() -> dict:
    """Circle tables (d/D = 0.985 and 0.8) and the pre-registered FCC orientation."""
    tables = [circle_table(s, d) for d in (0.985, 0.8) for s in ("SC", "BCC", "FCC")]
    for t in tables:
        for sh in t["shells"]:
            sh.pop("centres")
    return {
        "chi_star_over_D": CHI_STAR, "D_Gpc": D_LSS_GPC,
        "tables": tables, "paper_table": CIRCLE_TABLE_PAPER, "paper_fig4": FIG4_PAPER,
        "fcc_preregistered": fcc_preregistered_rotation(), "source": "[TOP] Sec. 7",
    }
