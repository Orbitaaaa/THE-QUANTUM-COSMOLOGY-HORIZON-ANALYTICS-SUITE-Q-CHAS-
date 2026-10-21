"""Project Baby Universe -- Genesis Branch: computational research package.

The package implements, and only implements, equations that are derived in
the accompanying documents (``genesis/*.md``):

* ``model``       thin-wall de Sitter / Schwarzschild-de Sitter junction model,
                  effective potential, orientation signs, turning points.
* ``symbolic``    sympy derivations: extrinsic curvature, Israel conditions,
                  squaring, mass function, reduced Hamiltonian, Euclidean
                  continuation, closed forms at M = 0.
* ``branches``    classification of classical solution branches and of the
                  causal regions visited by the wall (chart-free itinerary).
* ``geometry``    explicit global trajectories in Kruskal and global de Sitter
                  null coordinates; Penrose compactification.
* ``causal``      trapped / anti-trapped spheres, Misner-Sharp mass, minimal
                  surfaces of time-symmetric data, singularity-theorem checks.
* ``tunnelling``  canonical WKB exponent with explicit branch prescriptions,
                  Euclidean wall trajectories, multiplicity (pseudomanifold)
                  analysis, closed-form compact-instanton actions.
* ``thickwall``   scalar-field walls: time-symmetric constraint solutions and
                  O(4) Coleman-De Luccia bounces with gravity.
* ``quantum``     reduced-phase-space Hamiltonian, Hamilton-equation checks,
                  comparison of inequivalent quantizations.

Units: G = c = hbar = 1 unless a function takes G explicitly.  Signature
(-,+,+,+).  The orientation signs eps_in, eps_out follow the convention of
JUNCTION_CONDITIONS.md: eps = sign of n^r, n pointing from interior to exterior.
"""

from .model import ThinWallModel  # noqa: F401

__all__ = ["ThinWallModel"]
__version__ = "0.1.0"
