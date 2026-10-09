"""The papers' own status tables, quoted verbatim in substance: [BH] Sec. 16 and [QG] Sec. 20.

The topology paper's table ([TOP] Sec. 13) lives in :mod:`qchas.topology.manifolds`.
The dashboard's Claims Ledger renders all three, so a reader always sees whether a
statement is a theorem, a computation, a data analysis or a hypothesis.
"""

from __future__ import annotations

BH_CLAIMS = [
    ("Interior and radiation are commuting, entangled subalgebras of one irreducible algebra (Theorem A)",
     "Theorem (Einstein causality; Verch-Werner; Reeh-Schlieder)"),
    ("Field theory on a fixed background: no information out, monotone entanglement, mixed final state (Theorem B)",
     "Theorem, for every Einstein-causal QFT and every Reeh-Schlieder state"),
    ("4D Hawking radiation is a density matrix with S = beta S_BH, beta = 1.855; 2D radiation is disjoint from Fock space; soft sectors are central (Theorem C)",
     "Theorem (Bures-Kakutani) with computed greybody factors"),
    ("Gravitational dressing turns the type III_1 exterior algebra into a type II_1 factor (energy window) with maximal entropy A/4G; its entropy is the generalized entropy (Theorem G)",
     "Theorem (Takesaki; crossed-product trace; first law)"),
    ("Exterior capacity is compatible with Hawking's curve (Theorem I)", "Theorem"),
    ("Page curve <=> H*, no independent interior (Theorem J)", "Theorem, given unitarity"),
    ("H* for the eternal black hole: the dressed algebras of the two sides are each other's commutants (Theorem K)",
     "Theorem (Takesaki commutation)"),
    ("H* for an evaporating black hole",
     "Derived by the gravitational path integral (replica wormholes, island rule of Sec. 9)"),
    ("The island is forced once S_H > A/4G, with minimal size S_H - A/4G (Theorem D)",
     "Theorem, given capacity (Theorem G) and H*"),
    ("S(R) <= min{S_H, A/4G}; equality within e^{-|S_H - A/4G|} for scrambling; dS/du = dA/(4G du) after the Page time (Theorem E)",
     "Theorem, given capacity (Theorem G) and H* (2-design scrambling for the lower bound)"),
    ("Interior volume grows to ~M0^5 while the information capacity vanishes (Theorem F)",
     "Theorem (Christodoulou-Rovelli volume)"),
    ("Island location, scrambling-time lag, universal atmosphere profile in the continuum limit, exact unitary simulation",
     "Computation"),
]

QG_CLAIMS = [
    ("Wheeler-DeWitt constraint with a clock: Schrodinger dynamics in clock time, positive inner product, no ghosts (Theorem T)", "Theorem"),
    ("Universe with Lambda > 0 and an observer: type II_1 algebra; unique trace = no-boundary state; maximal entropy 3 pi / G Lambda (Theorem C)",
     "Theorem (crossed product; CLPW)"),
    ("Closed radiation universe: unitary evolution through the classical singularity; cyclic; universal bounce at ~0.26 l_P (Theorem Q)",
     "Theorem (i-ii); computation (iii)"),
    ("All linear gravitational waves through the Big Bang: exact factorization, integer frequencies, unitary, finite metric perturbation (Theorem Q')",
     "Theorem; <a^-2> computed"),
    ("Scalar sound waves of the radiation decouple exactly through the Big Bang (Theorem B)", "Theorem (verified symbolically)"),
    ("Nonlinear (mixmaster) Big Bang: unitary; every anisotropy channel falls below Hardy's threshold; bounce computed (Theorem N)",
     "Theorem (i-ii); exact computation (iii): bounce scale converged to +-3 % and P(a < 1e-3) < 1e-4, while the peak of <a^-2> grows slowly with the channel number"),
    ("Localization at the bounce: exact causal decoupling beyond 2 DeltaT, gradient corrections (l/lambda_max)^2, DeltaT a_max -> 0.49 (Theorem Lambda)",
     "Theorem; DeltaT computed"),
    ("Curvature bound <a^-2> <= 96 (E + <H_beta>); saturation of metric fluctuations (Theorem H)",
     "Theorem; anisotropy energy converged in channel number"),
    ("Local constraints of 4D gravity: no anomaly; strong Lie algebra (G, H_a) without structure functions; dust clock; unitary dust-time evolution (Theorem L)",
     "Theorem; (ii) verified exactly (spherical gravity + scalar) and in full 3+1 to 1e-18"),
    ("Lie algebra with Standard-Model matter; the clock equation; every scalar clock gives commuting local Hamiltonians (Theorem W)",
     "Theorem; symbolic proof of (ii); verified in full 3+1"),
    ("Gravity + non-abelian gauge fields: one Lie algebra off the Gauss surface; anomaly-free quantum Diff x G (Theorem Y)",
     "Theorem; verified in full 3+1 (SU(2)) to 1e-17"),
    ("Type III_1 to type II for every region of every globally hyperbolic spacetime (Theorem G_inf)",
     "Theorem (Verch; Takesaki; Jensen-Sorce-Speranza)"),
    ("Discrete geometry with minimal area 5.97 l_P^2; finite horizon Hilbert spaces; finite Standard-Model Hamiltonians (Theorem U)",
     "Theorem (uniqueness of the representation; LQG spectra; Thiemann), with gamma fixed by A/4G"),
    ("The whole Standard Model with gravity and the clocks: one Lie algebra off every constraint surface; all quantum anomalies vanish (Theorem Y_SM)",
     "Theorem; verified in full 3+1 to 4.8e-16; anomalies exact and in Lean 4"),
    ("Smooth metric from the discrete representation; Einstein dynamics at low energy (Theorem C)",
     "Theorem; metric computed exactly to j = 48; long-time peakedness proved in the cosmological sectors"),
    ("Einstein's equations from entanglement equilibrium, G = 1/(4 eta); maximal vacuum entanglement and the local first law from the crossed product (Theorems E, M, Corollary E')",
     "Theorem, for conformal matter at leading order in the ball size; eta computed; lattice check of the modular Hamiltonian"),
    ("H* in non-perturbative JT gravity; replica wormholes = statistics of the factorized Hilbert space (Theorem J)",
     "Theorem; gravity-QM agreement computed to +-0.003"),
    ("Complete evaporation with quantum geometry: S(R) -> 0 (Theorem J(iii))", "Theorem; computed through the Schwarzian regime"),
    ("Near-extremal 4D Reissner-Nordstrom black holes (Corollary 4D)", "Theorem, within the Schwarzian (JT) description of the throat"),
    ("H* for every 2D dilaton gravity, including the s-wave sector of 4D Schwarzschild black holes (Theorem P)",
     "Theorem (planar replica sum; all orders in the matrix completion); gravity-QM agreement computed to +-0.00016"),
    ("Complete evaporation with Schwarzschild thermodynamics (negative heat capacity): S(R) -> 0 (Theorem S)", "Theorem; computed to the endpoint"),
    ("Recovery of the interior from the radiation after the Page time, infidelity <= 2^{1+k-s} (Theorem R)",
     "Theorem, for 2-design scrambling; tested with the ETH dynamics"),
    ("Kerr-Newman exterior: type II algebra, anomaly-free constraint group R x SU(2) x U(1) (Theorem K)",
     "Theorem (crossed product; Kudler-Flam-Leutheusser-Satishchandran; Whitehead, Bargmann)"),
    ("Replica wormholes of Kerr-Newman black holes = rotation-covariant matrix completion; H* (Theorem P_K)",
     "Theorem (planar replica sum); gravity-QM agreement computed to +-0.0023"),
    ("A rotating black hole spins down and evaporates completely, S(R) -> 0 (Theorem S_K)", "Theorem; computed to the endpoint"),
    ("Exact Kerr greybody factors in the completion; Page spin-down; scalar fixed point a* = 0.555; generalized second law along the trajectory (Theorem K')",
     "Theorem; Teukolsky computation validated to 5e-6"),
    ("Predictions: gravity-mediated entanglement, no echoes, a transparent bounce, no remnants",
     "Computed; testable with tabletop interferometry, LIGO-Virgo-KAGRA, ET, CE, LISA, CMB and gamma-ray data"),
    ("Existing data: no echoes in nine LIGO-Virgo mergers (R < 0.28, CLs 90 %); Kerr ringdowns in GW150914 and GW250114; no bounce feature in Planck 2018; branch-phase law on superposed matter to 4e-4",
     "Our analyses of public data; consistent with every prediction. Two-body entanglement and primordial-black-hole bursts await experiments"),
]


def horizon_claims() -> list[dict]:
    rows = [{"statement": s, "status": st, "source": "[BH] Sec. 16"} for s, st in BH_CLAIMS]
    rows += [{"statement": s, "status": st, "source": "[QG] Sec. 20"} for s, st in QG_CLAIMS]
    return rows
