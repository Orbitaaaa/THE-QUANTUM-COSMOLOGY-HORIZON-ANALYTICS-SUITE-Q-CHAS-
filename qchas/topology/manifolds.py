"""The 18 flat 3-manifolds, principles P1-P4 and the Main Theorem, [TOP] Secs. 1-5, 9, 13, 14.

Each flat space is Sigma = R^3 / Gamma with gamma(x) = A x + b, A in O(3) ([TOP] Sec. 1).
The principle statuses follow from three flags by the proved theorems:

* P1 (uniform chiral physics)  <=> orientable                       (Theorem A)
* P2 (finite thermodynamics)   <=> closed                           (Theorem C')
* P3 (local Copernican)        <=> Gamma consists of translations   (Theorem B-loc)
* P4 (primordial SUSY)         <=> translations + periodic spin     (Theorem D);
  non-orientable spaces carry no spin structure at all.

The derived table is checked against the quoted [TOP] Sec. 4 table. The module also ships
explicit deck groups for E1-E6 ([TOP] Sec. 14 'Group checks') and the quoted evidence
tables of [TOP] Secs. 9, 9.2-9.4 and the status table of Sec. 13.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

import numpy as np

ALLOWED, EXCLUDED, NO_SPIN = "allowed", "excluded", "no spin structure"
PRINCIPLE_KEYS = ("P1", "P2", "P3", "P4")


@dataclass(frozen=True)
class Manifold:
    label: str
    name: str
    closed: bool
    orientable: bool
    holonomy: str
    holonomy_order: int | None
    translations_only: bool
    note: str = ""

    def principles(self) -> dict[str, str]:
        """Status of P1-P4 derived from Theorems A, C', B-loc and D."""
        return {
            "P1": ALLOWED if self.orientable else EXCLUDED,
            "P2": ALLOWED if self.closed else EXCLUDED,
            "P3": ALLOWED if self.translations_only else EXCLUDED,
            "P4": NO_SPIN if not self.orientable else (ALLOWED if self.translations_only else EXCLUDED),
        }

    @property
    def periodic_spin_possible(self) -> bool:
        return self.principles()["P4"] == ALLOWED

    def to_json(self) -> dict:
        return {
            "label": self.label, "name": self.name, "closed": self.closed,
            "orientable": self.orientable, "holonomy": self.holonomy,
            "holonomy_order": self.holonomy_order,
            "translations_only": self.translations_only,
            "periodic_spin_possible": self.periodic_spin_possible,
            "principles": self.principles(), "note": self.note,
        }


_R = "contains a reflection"
#: [TOP] Sec. 1 table, grouped rows expanded; individual holonomy groups of E7-E10, E13-E15,
#: E17 follow the Conway-Rossetti / COMPACT classification ([TOP] refs. 1-3, 21).
MANIFOLDS: tuple[Manifold, ...] = (
    Manifold("E1", "3-torus", True, True, "trivial", 1, True),
    Manifold("E2", "half-turn space", True, True, "Z2", 2, False),
    Manifold("E3", "quarter-turn space", True, True, "Z4", 4, False),
    Manifold("E4", "third-turn space", True, True, "Z3", 3, False),
    Manifold("E5", "sixth-turn space", True, True, "Z6", 6, False),
    Manifold("E6", "Hantzsche-Wendt space", True, True, "Z2 x Z2", 4, False),
    Manifold("E7", "Klein space", True, False, f"Z2, {_R}", 2, False),
    Manifold("E8", "Klein space with horizontal flip", True, False, f"Z2, {_R}", 2, False),
    Manifold("E9", "Klein space with vertical flip", True, False, f"Z2 x Z2, {_R}", 4, False),
    Manifold("E10", "Klein space with half-turn", True, False, f"Z2 x Z2, {_R}", 4, False),
    Manifold("E11", "chimney space T^2 x R", False, True, "trivial", 1, True),
    Manifold("E12", "chimney space with half-turn", False, True, "Z2", 2, False),
    Manifold("E13", "chimney space with vertical flip", False, False, f"Z2, {_R}", 2, False),
    Manifold("E14", "chimney space with horizontal flip", False, False, f"Z2, {_R}", 2, False),
    Manifold("E15", "chimney space with half-turn and flip", False, False, f"Z2 x Z2, {_R}", 4, False),
    Manifold("E16", "slab space S^1 x R^2", False, True, "trivial", 1, True,
             "an optional screw angle adds a rotation to the metric, not to the topology; "
             "P3 and P4 hold only without screw"),
    Manifold("E17", "slab space with flip", False, False, f"Z2, {_R}", 2, False),
    Manifold("E18", "covering space R^3", False, True, "trivial", 1, True),
)
BY_LABEL = {m.label: m for m in MANIFOLDS}

#: [TOP] Sec. 4 status table exactly as printed (grouped rows).
PRINCIPLE_TABLE_PAPER: dict[str, tuple[str, str, str, str]] = {
    "E1": (ALLOWED, ALLOWED, ALLOWED, ALLOWED),
    "E2-E5": (ALLOWED, ALLOWED, EXCLUDED, EXCLUDED),
    "E6": (ALLOWED, ALLOWED, EXCLUDED, EXCLUDED),
    "E7-E10": (EXCLUDED, ALLOWED, EXCLUDED, NO_SPIN),
    "E11": (ALLOWED, EXCLUDED, ALLOWED, ALLOWED),
    "E12": (ALLOWED, EXCLUDED, EXCLUDED, EXCLUDED),
    "E13-E15": (EXCLUDED, EXCLUDED, EXCLUDED, NO_SPIN),
    "E16": (ALLOWED, EXCLUDED, ALLOWED, ALLOWED),
    "E17": (EXCLUDED, EXCLUDED, EXCLUDED, NO_SPIN),
    "E18": (ALLOWED, EXCLUDED, ALLOWED, ALLOWED),
}


def expand_label_range(key: str) -> list[str]:
    """'E2-E5' -> ['E2', 'E3', 'E4', 'E5']."""
    if "-" not in key:
        return [key]
    a, b = key.split("-")
    return [f"E{i}" for i in range(int(a[1:]), int(b[1:]) + 1)]


@dataclass(frozen=True)
class Principle:
    key: str
    name: str
    status: str
    statement: str
    theorem: str

    def to_json(self) -> dict:
        return {"key": self.key, "name": self.name, "status": self.status,
                "statement": self.statement, "theorem": self.theorem, "source": "[TOP] Sec. 3"}


PRINCIPLES: tuple[Principle, ...] = (
    Principle("P1", "Uniform chiral physics", "experimentally grounded",
              "The Standard Model, with its measured CP violation, holds throughout N. The CKM "
              "Jarlskog invariant is J = (3.12 +0.13/-0.12) x 10^-5 (PDG 2024), nonzero at about 25 sigma.",
              "Theorem A"),
    Principle("P2", "Finite primordial thermodynamics", "hypothesis",
              "A free quantum field in thermal equilibrium at the Big Bang had a finite partition "
              "function at every temperature. Theorem C' proves this is equivalent to each of: finite "
              "spatial volume; closed space; and Hartle-Hawking no-boundary creation (space is the "
              "boundary of a compact 4-geometry).",
              "Theorem C'"),
    Principle("P3", "Local Copernican principle", "hypothesis, minimal form",
              "Every local observable of the primordial quantum field has the same expectation value "
              "at every point of space. This constrains only what an observer can measure in an "
              "arbitrarily small neighbourhood of themselves. It is strictly weaker than asking the "
              "whole two-point function to be statistically homogeneous.",
              "Theorem B-loc"),
    Principle("P4", "Primordial supersymmetry", "speculative, optional",
              "At creation the state preserves an untwisted global supersymmetry: some supercharge "
              "has a covariantly constant spinor parameter on Sigma.",
              "Theorem D"),
)


def satisfies(label: str, principles) -> bool:
    st = BY_LABEL[label].principles()
    return all(st[p] == ALLOWED for p in principles)


def survivors(principles) -> list[str]:
    """Labels of the manifolds compatible with every principle in ``principles`` (Main Theorem)."""
    ps = list(principles)
    for p in ps:
        if p not in PRINCIPLE_KEYS:
            raise ValueError(f"unknown principle {p!r}")
    return [m.label for m in MANIFOLDS if satisfies(m.label, ps)]


def failed_principles(label: str) -> list[str]:
    st = BY_LABEL[label].principles()
    return [p for p in PRINCIPLE_KEYS if st[p] != ALLOWED]


def minimal_selecting_sets() -> list[list[str]]:
    """Inclusion-minimal sets of principles whose survivors are exactly {E1}."""
    found: list[set[str]] = []
    for r in range(1, 5):
        for c in itertools.combinations(PRINCIPLE_KEYS, r):
            if survivors(c) == ["E1"] and not any(f <= set(c) for f in found):
                found.append(set(c))
    return [sorted(f) for f in found]


# --- deck groups of E1-E6 ---------------------------------------------------------------

def rot_z(theta: float) -> np.ndarray:
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def compose(g1, g2):
    """(A1, b1) o (A2, b2): x -> A1 (A2 x + b2) + b1."""
    A1, b1 = g1
    A2, b2 = g2
    return A1 @ A2, A1 @ b2 + b1


def inverse(g):
    A, b = g
    return A.T, -A.T @ b


@dataclass
class DeckGroup:
    """A flat-space group Gamma: generators (A, b), translation lattice T (rows) and holonomy reps."""

    label: str
    L: float
    generators: list[tuple[np.ndarray, np.ndarray]]
    translations: np.ndarray
    reps: list[tuple[np.ndarray, np.ndarray]] = field(default_factory=list)

    def __post_init__(self):
        if not self.reps:
            self.reps = _holonomy_reps(self.generators, self.translations)

    @property
    def holonomy(self) -> list[np.ndarray]:
        return [A for A, _ in self.reps]

    @property
    def order(self) -> int:
        return len(self.reps)

    def reduce(self, b: np.ndarray) -> np.ndarray:
        """Translation part reduced into the fundamental cell of T (coordinates in [0, 1))."""
        c = np.linalg.solve(self.translations.T, b)
        c = c - np.floor(c + 1e-9)
        return self.translations.T @ c

    def decompose(self, A: np.ndarray, b: np.ndarray) -> tuple[int, np.ndarray]:
        """Index of the holonomy rep with linear part A and the coefficients of b - b_A in T."""
        for i, (Ai, bi) in enumerate(self.reps):
            if np.allclose(Ai, A, atol=1e-9):
                return i, np.linalg.solve(self.translations.T, b - bi)
        raise ValueError("linear part not in the holonomy group")

    def contains(self, A: np.ndarray, b: np.ndarray, atol: float = 1e-8) -> bool:
        try:
            _, c = self.decompose(A, b)
        except ValueError:
            return False
        return bool(np.allclose(c, np.round(c), atol=atol))

    def translation_vectors(self, radius: float, include_zero: bool = True) -> np.ndarray:
        """Vectors t of the translation subgroup with |t| <= radius."""
        from .lattice import lattice_vectors
        return lattice_vectors(self.translations, radius, include_zero=include_zero)

    def elements(self, radius: float) -> list[tuple[np.ndarray, np.ndarray]]:
        """Elements (A, b_A + t) with |t| <= radius (identity included)."""
        ts = self.translation_vectors(radius)
        return [(A, bA + t) for A, bA in self.reps for t in ts]

    def words(self, length: int) -> list[tuple[np.ndarray, np.ndarray]]:
        """All products of at most ``length`` generators and inverses."""
        letters = list(self.generators) + [inverse(g) for g in self.generators]
        out = [(np.eye(3), np.zeros(3))]
        frontier = list(out)
        for _ in range(length):
            frontier = [compose(w, s) for w in frontier for s in letters]
            out.extend(frontier)
        return out


def _holonomy_reps(generators, T) -> list[tuple[np.ndarray, np.ndarray]]:
    Tt = np.asarray(T, dtype=float).T

    def red(b):
        c = np.linalg.solve(Tt, b)
        return Tt @ (c - np.floor(c + 1e-9))

    reps = [(np.eye(3), np.zeros(3))]
    frontier = list(reps)
    while frontier:
        new = []
        for g in frontier:
            for s in generators:
                A, b = compose(g, s)
                if not any(np.allclose(A, Ar, atol=1e-9) for Ar, _ in reps):
                    reps.append((A, red(b)))
                    new.append((A, red(b)))
        frontier = new
    return reps


def deck_group(label: str, L: float = 1.0) -> DeckGroup:
    """Standard generators of E1-E6 of side L ([TOP] Sec. 14 'Group checks').

    E1 translations of a cube; E2/E3 screws about z over a square xy lattice; E4/E5 screws over
    a hexagonal xy lattice; E6 Hantzsche-Wendt with three half-turn screws about x, y, z.
    """
    e = np.eye(3) * L
    hexa = np.array([[L, 0.0, 0.0], [-0.5 * L, 0.5 * math.sqrt(3.0) * L, 0.0], [0.0, 0.0, L]])
    z = lambda frac: np.array([0.0, 0.0, frac * L])
    if label == "E1":
        gens = [(np.eye(3), e[i]) for i in range(3)]
        return DeckGroup(label, L, gens, e.copy())
    if label in ("E2", "E3"):
        n = 2 if label == "E2" else 4
        gens = [(rot_z(2 * math.pi / n), z(1.0 / n)), (np.eye(3), e[0]), (np.eye(3), e[1])]
        return DeckGroup(label, L, gens, e.copy())
    if label in ("E4", "E5"):
        n = 3 if label == "E4" else 6
        gens = [(rot_z(2 * math.pi / n), z(1.0 / n)), (np.eye(3), hexa[0]), (np.eye(3), hexa[1])]
        return DeckGroup(label, L, gens, hexa)
    if label == "E6":
        h = 0.5 * L
        gens = [
            (np.diag([1.0, -1.0, -1.0]), np.array([h, h, 0.0])),
            (np.diag([-1.0, 1.0, -1.0]), np.array([0.0, h, h])),
            (np.diag([-1.0, -1.0, 1.0]), np.array([h, 0.0, h])),
        ]
        return DeckGroup(label, L, gens, e.copy())
    raise ValueError("explicit deck groups are provided for E1-E6")


def fixed_point_free(group: DeckGroup, radius: float = 3.0) -> bool:
    """True iff no element (A, b_A + t) != 1 with |t| <= radius has a fixed point.

    (A - I) x = -b is solvable iff b has no component in Fix(A) = ker(A - I).
    """
    for A, b in group.elements(radius * group.L):
        if np.allclose(A, np.eye(3)) and np.allclose(b, 0.0):
            continue
        _, s, vt = np.linalg.svd(A - np.eye(3))
        fix = vt[s < 1e-9]
        if fix.shape[0] == 0 or np.linalg.norm(fix @ b) < 1e-9:
            return False
    return True


def min_displacement(group: DeckGroup, points: np.ndarray, radius: float = 3.0) -> float:
    """min over sample points x and elements gamma != 1 of |gamma x - x|."""
    pts = np.atleast_2d(points)
    best = math.inf
    for A, b in group.elements(radius * group.L):
        if np.allclose(A, np.eye(3)) and np.allclose(b, 0.0):
            continue
        d = np.linalg.norm(pts @ (A - np.eye(3)).T + b, axis=1)
        best = min(best, float(d.min()))
    return best


def displacement_is_constant(group: DeckGroup, points: np.ndarray, radius: float = 2.0) -> bool:
    """Whether every |gamma x - x| (truncated set of gamma) is independent of x."""
    pts = np.atleast_2d(points)
    for A, b in group.elements(radius * group.L):
        d = np.linalg.norm(pts @ (A - np.eye(3)).T + b, axis=1)
        if np.ptp(d) > 1e-9:
            return False
    return True


def _su2_lift(A: np.ndarray) -> np.ndarray:
    """One of the two SU(2) lifts of a rotation A in SO(3)."""
    w, V = np.linalg.eig(A)
    axis = np.real(V[:, np.argmin(np.abs(w - 1.0))])
    axis /= np.linalg.norm(axis)
    c = np.clip((np.trace(A) - 1.0) / 2.0, -1.0, 1.0)
    th = math.acos(c)
    if np.linalg.norm(A - A.T) > 1e-12:
        sgn = np.sign(np.dot([A[2, 1] - A[1, 2], A[0, 2] - A[2, 0], A[1, 0] - A[0, 1]], axis))
        axis = axis * (sgn if sgn != 0 else 1.0)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    n_sigma = axis[0] * sx + axis[1] * sy + axis[2] * sz
    return math.cos(th / 2) * np.eye(2) - 1j * math.sin(th / 2) * n_sigma


def parallel_spinor_dimension(group: DeckGroup) -> int:
    """max over sign choices of the generator lifts of dim{psi0 : eps(gamma) psi0 = psi0} ([TOP] Thm D)."""
    if any(np.linalg.det(A) < 0 for A, _ in group.generators):
        return 0
    lifts = [_su2_lift(A) for A, _ in group.generators]
    best = 0
    for signs in itertools.product((1.0, -1.0), repeat=len(lifts)):
        M = np.concatenate([s * U - np.eye(2) for s, U in zip(signs, lifts)])
        sv = np.linalg.svd(M, compute_uv=False)
        best = max(best, int(np.sum(sv < 1e-10)))
    return best


# --- quoted tables ----------------------------------------------------------------------

CURVED = [
    {"name": "Poincaré dodecahedral S^3/I*", "order": 120, "shortest_loop": "π/5",
     "omega_k_range": [-0.0104, -0.0010], "lnB_vs_S3": -0.05,
     "note": "ln B (N <= 100) = -0.04; near the circle edge (left/mirror) ln B = -4.77/-3.74 at Omega_k = -0.0104",
     "source": "[TOP] Sec. 9.3"},
    {"name": "Binary octahedral S^3/O*", "order": 48, "shortest_loop": "π/4",
     "omega_k_range": [-0.0183, -0.0010], "lnB_vs_S3": -0.05,
     "note": "ln B (N <= 100) = -0.04; near the circle edge (left/mirror) ln B = -1.99/-3.94 at Omega_k = -0.0183",
     "source": "[TOP] Sec. 9.3"},
    {"name": "Binary tetrahedral S^3/T*", "order": 24, "shortest_loop": "π/3",
     "omega_k_range": [-0.0200, -0.0010], "lnB_vs_S3": 0.01,
     "note": "ln B (N <= 100) = +0.01; all KL <= 0.2 (undetectable)",
     "source": "[TOP] Sec. 9.3"},
]

#: Weeks manifold systole in curvature units, [TOP] Sec. 9.4 (SnapPy m003(-3,1))
WEEKS_SYSTOLE = 0.58460
WEEKS_VOLUME = 0.9427073628
WEEKS_OMEGA_K_THRESHOLD = 0.0080

_W = "no circles for any observer"
_C = "circles possible near a systolic geodesic"
WEEKS = [
    {"omega_k": ok, "Rc_Gpc": rc, "chi_over_Rc": x, "shortest_loop_over_D": s, "volume_ratio": v,
     "circles": c, "source": "[TOP] Prop. 9.4"}
    for ok, rc, x, s, v, c in [
        (0.0005, 198.3, 0.070, 4.18, 659.0, _W), (0.0010, 139.7, 0.099, 2.95, 230.8, _W),
        (0.0020, 98.1, 0.141, 2.07, 80.1, _W), (0.0045, 64.2, 0.215, 1.36, 22.6, _W),
        (0.0070, 50.4, 0.272, 1.08, 11.1, _W), (0.0100, 41.1, 0.331, 0.88, 6.1, _C),
        (0.0150, 32.0, 0.420, 0.70, 2.9, _C), (0.0200, 26.3, 0.506, 0.58, 1.7, _C),
    ]
]

EVIDENCE_TORUS = [
    {"shape": s, "d_over_D": d, "max_dlnL": m, "lnB": b, "null95": n, "p": p, "source": "[TOP] Table 9a"}
    for s, d, m, b, n, p in [
        ("FCC", 1.05, 7.20, -1.20, 1.12, 0.44), ("FCC", 1.1, 5.83, 0.07, 0.61, 0.31),
        ("FCC", 1.2, 2.82, 0.15, 0.27, 0.19), ("FCC", 1.3, 1.62, 0.13, 0.15, 0.11),
        ("cubic", 1.05, 7.84, 0.25, 0.68, 0.15), ("cubic", 1.1, 3.93, 0.14, 0.41, 0.25),
        ("cubic", 1.2, 2.03, 0.15, 0.21, 0.12), ("cubic", 1.3, 1.18, 0.17, 0.13, 0.03),
        ("BCC", 1.05, 9.70, 1.12, 0.85, 0.04), ("BCC", 1.1, 6.33, 0.19, 0.52, 0.18),
        ("BCC", 1.2, 2.56, 0.15, 0.24, 0.15), ("BCC", 1.3, 1.32, 0.12, 0.13, 0.11),
    ]
]

_X, _A = "excluded (95%)", "allowed"
EXCLUSION_TORUS = [
    {"shape": s, "d_over_D": d, "lnB_planck": b, "median_if_true": m, "q05_if_true": q,
     "p_le_planck": p, "sims": n, "verdict": v, "source": "[TOP] Table 9b"}
    for s, d, b, m, q, p, n, v in [
        ("FCC", 0.985, -22.14, 17.99, 3.07, "<=0.016", 60, _X),
        ("FCC", 1.0, -21.03, 15.46, 0.06, "<=0.016", 60, _X),
        ("FCC", 1.02, -12.56, 14.40, -2.30, "<=0.016", 60, _X),
        ("FCC", 1.03, -5.99, 4.98, -2.88, "<=0.016", 60, _X),
        ("FCC", 1.04, -2.60, 2.84, -2.07, "0.034", 260, _X),
        ("FCC", 1.05, -1.20, 1.35, -1.71, "0.092", 260, _A),
        ("FCC", 1.1, 0.07, 0.17, -0.72, "0.410", 60, _A),
        ("FCC", 1.15, 0.16, 0.11, -0.29, "0.557", 60, _A),
        ("FCC", 1.2, 0.15, 0.04, -0.29, "0.770", 60, _A),
        ("cubic", 0.985, -5.81, 6.23, -1.08, "0.033", 60, _X),
        ("cubic", 1.0, -2.39, 3.45, -3.78, "0.115", 60, _A),
        ("cubic", 1.05, 0.25, -0.03, -0.85, "0.689", 60, _A),
        ("cubic", 1.1, 0.14, 0.02, -0.26, "0.656", 60, _A),
        ("BCC", 0.985, -8.72, 10.55, -0.98, "<=0.016", 60, _X),
        ("BCC", 1.0, -4.22, 10.63, 1.05, "0.033", 60, _X),
        ("BCC", 1.05, 1.12, 0.33, -1.25, "0.689", 60, _A),
        ("BCC", 1.1, 0.19, 0.08, -0.35, "0.672", 60, _A),
    ]
]

OTHER_FLAT = [
    {"space": s, "orientable": o, "samples": 48, "lnB": b, "largest": g, "source": "[TOP] Sec. 9.2"}
    for s, o, b, g in [
        ("E2 half-turn", True, "+0.20 ± 0.03", 0.89), ("E3 quarter-turn", True, "-0.05 ± 0.05", 0.79),
        ("E4 third-turn", True, "-0.00 ± 0.03", 0.22), ("E5 sixth-turn", True, "-0.73 ± 0.20", 1.32),
        ("E6 Hantzsche-Wendt", True, "+0.23 ± 0.06", 1.07), ("E7 Klein", False, "+0.12 ± 0.04", 0.81),
        ("E8 Klein, horizontal flip", False, "+0.21 ± 0.05", 0.64),
        ("E9 Klein, vertical flip", False, "+0.04 ± 0.03", 0.39),
        ("E10 Klein, half-turn", False, "+0.17 ± 0.04", 0.75),
    ]
]

_CLAIM_ROWS = [
    ("18 flat 3-manifolds, 10 closed; local physics cannot distinguish them (Lemma 1)",
     "Theorem, unconditional"),
    ("Chirality with J ≠ 0 forces orientability (Theorem A)",
     "Theorem, from P1 (experimentally grounded)"),
    ("Homogeneous vacuum on a closed flat space ⇔ torus (Theorem B)", "Theorem, unconditional"),
    ("Locally homogeneous vacuum on any flat space ⇔ translations only (Theorem B-loc)",
     "Theorem, unconditional"),
    ("Finite partition function ⇔ finite volume ⇔ closed ⇔ no-boundary (Theorem C′)",
     "Theorem, unconditional"),
    ("Parallel spinor on a flat space ⇔ translations + periodic spin structure (Theorem D)",
     "Theorem, unconditional"),
    ("Selection of T³ (Main Theorem)",
     "Theorem, from P2 with P3 or P4 (hypotheses); P2 and P3 are logically independent and both needed"),
    ("No test on the ℓ ≤ 30 sky can separate a torus with d ≥ 1.5D from infinite space by more than "
     "0.15 in power minus size (Theorem P)", "Theorem, with exactly computed divergences"),
    ("No evidence for a torus of allowed size; exclusion d_FCC > 1.04D, d_BCC > 1.00D, "
     "d_cubic > 0.985D (Section 9)", "Data analysis, calibrated on simulations, pre-registered where stated"),
    ("Cubic lattices are critical for every rotation-invariant functional (Theorem E)",
     "Theorem, unconditional"),
    ("FCC and BCC are strict local minima, SC a saddle, Z_FCC < Z_BCC (Proposition E.2)",
     "Theorem, computer-assisted with certified error bounds"),
    ("Shape equations and convergence to critical points (Theorems H.1, H.2)",
     "Theorem (H.2 for a radiation background)"),
    ("FCC is the actual shape of our torus",
     "Conditional: needs a pre-inflationary Casimir era and a start in the FCC basin "
     "(Proposition H.3); not assumed in the data analysis"),
    ("FCC is the global minimum of Z(2)",
     "Not used; announced in an unrefereed preprint, consistent with our numerics"),
    ("Orientation cannot be predicted by rotation-invariant laws (Proposition G)", "Theorem, unconditional"),
    ("Circle positions, radii and matching (Theorem F, Corollary F.1)", "Theorem, exact given (R, d, χ*)"),
    ("Low-ℓ power, S_1/2, KL values (Section 8)",
     "Numerical, linear theory, Planck 2018 parameters, converged to 0.1%"),
    ("No evidence for any of E2–E10 in the Planck sky, ln B ∈ [−0.73, +0.23] (Section 9.2)",
     "Pre-registered data analysis"),
    ("No evidence for the dodecahedral, octahedral or tetrahedral spherical spaces; near-edge "
     "dodecahedral and octahedral configurations disfavoured (Section 9.3)", "Pre-registered data analysis"),
    ("Weeks manifold: no circles unless Ω_k ≥ 0.0080 (Proposition 9.4); statistical signal below that "
     "estimated to be undetectable", "Theorem (circles, given measured curvature); estimate (statistics)"),
    ("Everything any observer can learn is bounded by the relative entropy on the causal ball (Theorem U)",
     "Theorem, from Einstein causality, the time-slice axiom and Uhlmann monotonicity"),
    ("Perfect maps identify the topology with certainty iff the causal ball overlaps its images; "
     "otherwise the information is finite, unless the spectrum has lines (Theorem V)",
     "Theorem (Feldman–Hájek)"),
    ("Ceiling values for a perfect 3D map (Section 11.2)",
     "Numerical; every entry is a rigorous lower bound, the true ceiling is finite by Theorem V"),
    ("Laboratory bound 0.28(ℓ/L)^8 (Proposition E)",
     "Exact lattice computation, continuum-extrapolated; exponent explained by the first law of entanglement"),
    ("Defects cannot beat their seeds; T-duality adds nothing (Propositions D, T)",
     "Theorem (D(i), D(ii), T); estimate (D(iii))"),
    ("No topological resonance in the Planck sky (p = 0.78 and 0.53, Section 11.6)",
     "Pre-registered data analysis, 160 models, calibrated on simulations"),
    ("Single-field inflation limits resonant lines to f < 0.05 on most shells; multi-field lines mapped "
     "in 3D decide tori up to 3D (Section 11.8)",
     "Theorem (growth bounds) with exact computation; 3D values are rigorous lower bounds"),
    ("Beyond causality, only a messenger with v > d/2η₀ helps; instantaneous response does not "
     "(Proposition C)", "Theorem"),
    ("Size d and orientation R of our torus", "Must be measured; currently d ≥ 0.985D"),
]
CLAIMS = [{"statement": s, "status": t, "source": "[TOP] Sec. 13"} for s, t in _CLAIM_ROWS]


def to_json() -> dict:
    """The full ``topology.json`` bundle of docs/ARCHITECTURE.md."""
    return {
        "manifolds": [m.to_json() for m in MANIFOLDS],
        "principles": [p.to_json() for p in PRINCIPLES],
        "curved": [dict(r) for r in CURVED],
        "weeks": [dict(r) for r in WEEKS],
        "evidence_torus": [dict(r) for r in EVIDENCE_TORUS],
        "exclusion_torus": [dict(r) for r in EXCLUSION_TORUS],
        "other_flat": [dict(r) for r in OTHER_FLAT],
        "claims": [dict(r) for r in CLAIMS],
    }
