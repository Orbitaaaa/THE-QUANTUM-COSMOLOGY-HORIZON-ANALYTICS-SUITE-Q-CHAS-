"""Roadmap Stage 4 (tasks 31-38): alternative mechanisms, each with its quantitative cost.

* ``SeededModel``                 a false-vacuum bubble around a seed black hole (interior Schwarzschild-de Sitter).
                                  The exponent drops by at most ~30%, and the interior time-symmetric slice beyond the seed's
                                  bridge already contains a de Sitter equator: the initial data already contain the daughter
                                  (FAILED_APPROACHES F15).
* ``qei_max_antitrapped_radius``  quantum-energy-inequality bound on the size of a sphere that NEC-violating quantum fields
  ``species_required``            can anti-trap (Prop. IV-1 + a QEI of Ford-Roman type); the number of species needed for an
                                  inflating daughter is of order its de Sitter entropy, which pushes the species cutoff
                                  m_P / sqrt(N) below the daughter's Hubble rate.
* ``ds_parent_upward_exponent``   the de Sitter (closed-parent) loophole for our own universe: upward tunnelling costs about
                                  the parent's de Sitter entropy, ~3e122.
* ``fmm_junction``                limiting curvature (Frolov-Markov-Mukhanov): radius where the Schwarzschild interior reaches
                                  the curvature cap; every black hole above ~m_P would contain a Planck-curvature de Sitter phase.
* ``mechanism_ledger``            the summary table of Part E.

Planck units unless stated (G = hbar = c = 1).
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

M_PL_GEV = 1.22e19
FORD_ROMAN_C = 3 / (32 * np.pi**2)      # <rho> >= -3/(32 pi^2 tau^4) for a Lorentzian sampling of width tau


# ------------------------------------------------------------------ seeded bubble (pre-existing throat)
class SeededModel:
    """Interior f_in = 1 - 2 G m / r - chi^2 r^2 (seed mass m inside the bubble), exterior f_out = 1 - 2 G M / r.
    The junction formulas A_in/out = (f_in - f_out +- kappa^2 R^2)/(2 kappa R) and V = f_out - A_out^2 hold for any static f."""

    def __init__(self, chi: float, kappa: float, M: float, m: float, G: float = 1.0):
        self.chi, self.kappa, self.M, self.m, self.G = chi, kappa, M, m, G

    def f_in(self, R):
        return 1 - 2 * self.G * self.m / R - self.chi**2 * R**2

    def f_out(self, R):
        return 1 - 2 * self.G * self.M / R

    def A_out(self, R):
        return (self.f_in(R) - self.f_out(R) - self.kappa**2 * R**2) / (2 * self.kappa * R)

    def A_in(self, R):
        return (self.f_in(R) - self.f_out(R) + self.kappa**2 * R**2) / (2 * self.kappa * R)

    def V(self, R):
        return self.f_out(R) - self.A_out(R) ** 2

    def turning_points(self):
        R = np.geomspace(1e-5 / self.chi, 10 / self.chi, 40001)
        v = self.V(R)
        idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
        return [brentq(self.V, R[i], R[i + 1], xtol=1e-15) for i in idx]

    def B(self) -> float:
        tp = self.turning_points()
        if len(tp) < 2:
            return np.nan
        R1, R2 = tp[0], tp[1]

        def th(R):
            fi, fo = self.f_in(R), self.f_out(R)
            V = max(self.V(R), 0.0)
            return np.arctan2(self.kappa * R * np.sqrt(V / (fi * fo)), (fi + fo - self.kappa**2 * R**2) / (2 * np.sqrt(fi * fo)))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            v = quad(lambda R: R * th(R), R1, R2, limit=400, epsabs=0.0, epsrel=1e-10)[0]
        return 2 * v / self.G

    def interior_horizons(self):
        """Positive roots of f_in: seed horizon r_b and the interior cosmological horizon r_c (if m < 1/(3 sqrt3 G chi))."""
        roots = np.roots([-self.chi**2, 0.0, 1.0, -2 * self.G * self.m])
        return sorted(float(r.real) for r in roots if abs(r.imag) < 1e-12 and r.real > 0)

    def initial_data_contains_equator(self) -> bool:
        """With m > 0 the interior side of the wall continues through the seed bridge (minimal sphere r_b) into a static
        SdS patch whose time-symmetric slice reaches the cosmological bifurcation sphere r_c, a maximal sphere: a de Sitter
        equator, i.e. the seed of an inflating region, is already present in the initial data."""
        return self.m > 0 and len(self.interior_horizons()) == 2


# ------------------------------------------------------------------ quantum energy inequalities
def qei_max_antitrapped_radius(N: float = 1.0, C: float = FORD_ROMAN_C) -> float:
    """Prop. IV-1 along an ingoing null ray from flat data (r_u / alpha^2 = -1/2) gives
    int r T_uu du <= -1/(8 pi G) for an anti-trapped sphere.  With |T_uu| <= C N / (4 tau^4) sampled over tau ~ Delta u ~ 2 r:
    C N / (2 r^2) >= 1/(8 pi)  =>  r <= sqrt(4 pi C N) (Planck lengths).  Order-of-magnitude (constants O(1))."""
    return float(np.sqrt(4 * np.pi * C * N))


def species_required(E_vac_GeV: float, C: float = FORD_ROMAN_C) -> dict:
    """Species needed for QEI-bounded fields to anti-trap a de Sitter sphere of radius 1/chi (the minimum for inflation)."""
    rho = (E_vac_GeV / M_PL_GEV) ** 4
    chi = np.sqrt(8 * np.pi * rho / 3)
    N = 1 / (4 * np.pi * C * chi**2)
    S_dS = np.pi / chi**2
    cutoff = 1 / np.sqrt(N)                       # species scale m_P / sqrt(N), Planck units
    return dict(chi=float(chi), N_required=float(N), S_dS=float(S_dS), N_over_SdS=float(N / S_dS),
                species_cutoff_over_H=float(cutoff / chi))


# ------------------------------------------------------------------ de Sitter / closed parent
def ds_parent_upward_exponent(H0_km_s_Mpc: float = 67.4, Omega_L: float = 0.685) -> dict:
    """Upward (Lee-Weinberg) tunnelling from our Lambda vacuum costs ~ S_dS(parent) - S_dS(false vacuum) ~ S_dS(parent)."""
    H0 = H0_km_s_Mpc * 1e3 / 3.0857e22 * 5.391e-44          # Planck units
    eta = H0 * np.sqrt(Omega_L)
    return dict(eta=float(eta), S_dS_parent=float(np.pi / eta**2), horizon_radius_m=float(1.616e-35 / eta))


# ------------------------------------------------------------------ limiting curvature
def fmm_junction(M_planck: float, k: float = 1.0) -> dict:
    """Schwarzschild interior reaches Kretschmann 48 M^2 / r^6 = k at r_* = (48 M^2 / k)^(1/6); inside the horizon iff
    M > (3/(4k))^(1/4).  The de Sitter phase at the same curvature has 24 chi^4 = k."""
    r_star = (48 * M_planck**2 / k) ** (1 / 6)
    return dict(r_star=float(r_star), inside_horizon=bool(r_star < 2 * M_planck), M_min=float((3 / (4 * k)) ** 0.25),
                chi_daughter=float((k / 24) ** 0.25))


# ------------------------------------------------------------------ ledger
def mechanism_ledger() -> list:
    """Part E summary: which assumption each route relaxes, what it costs, and the verdict."""
    sp = species_required(1e17)
    dsp = ds_parent_upward_exponent()
    return [
        dict(mechanism="FGG tunnelling, laboratory data (thin barrier)", relaxes="classical GR (quantum tunnelling)",
             cost="B = B_CdL Phi((R_CdL/R_M)^2) >= B_CdL; B_lab = 4 pi^2 sigma^2 R_M^2 / rho",
             verdict="allowed but suppressed; P ~ exp(-1.32 sqrt(B_CdL) (M/m_P)^(2/3) (m_P^4/rho)^(1/6))"),
        dict(mechanism="FGG tunnelling, heavy walls (kappa > chi)", relaxes="classical GR",
             cost="B >= 0.564 B(M=0) (numerical inf, A/B boundary, gamma -> -1)", verdict="suppressed"),
        dict(mechanism="seed black hole inside the bubble", relaxes="regular initial data",
             cost="B reduced by <= ~30%; initial data already contain a de Sitter equator behind the seed throat",
             verdict="not a creation mechanism (F15)"),
        dict(mechanism="pre-existing throat / white hole (data with r_u > 0)", relaxes="past completeness",
             cost="requires a past singularity (Cor. IV-2, NG2)", verdict="not producible from laboratory data"),
        dict(mechanism="QEI-bounded NEC violation by quantum fields", relaxes="NEC",
             cost=f"N >~ {sp['N_over_SdS']:.2f} S_dS species; species cutoff = {sp['species_cutoff_over_H']:.2f} H at 1e17 GeV",
             verdict="self-defeating: effective theory breaks down below the daughter's Hubble scale"),
        dict(mechanism="closed / de Sitter parent (compact Cauchy surfaces)", relaxes="non-compact Cauchy surface",
             cost=f"bubble beyond the parent horizon; upward tunnelling ~ exp(-{dsp['S_dS_parent']:.1e})",
             verdict="cosmological only (inflationary era: Garriga-Vilenkin-Zhang)"),
        dict(mechanism="bubbles nucleated during inflation (GVZ 2016)", relaxes="laboratory initial data",
             cost="M > M_cr ~ 1/(2 G L); no tunnelling suppression beyond the nucleation rate",
             verdict="viable in the early universe; predicts a PBH population; not a laboratory route"),
        dict(mechanism="monopoles / charged exteriors (Sakai et al. 2006)", relaxes="global hyperbolicity",
             cost="tunnelling, wall collision, or continuation beyond a Cauchy horizon (mass inflation)",
             verdict="no predictable classical route; core inflation needs eta_sym ~ 0.3 m_P"),
        dict(mechanism="aspherical collapse / thick walls", relaxes="symmetry, thin wall",
             cost="Farhi-Guth/Penrose are symmetry-free; thick-wall CdL converge to thin-wall values",
             verdict="no escape classically; quantum: spherical saddle expected dominant (open)"),
        dict(mechanism="limiting curvature (FMM 1990)", relaxes="GR at Planck curvature",
             cost=f"any black hole with M > {fmm_junction(1.0)['M_min']:.2f} m_P; daughter at Planck curvature",
             verdict="hypothesis-dependent; daughter invisible from outside"),
    ]


# ------------------------------------------------------------------ de Sitter parent: the non-flip window
def sds_nonflip_masses(kappa: float, eta: float, chi: float = 1.0, n: int = 300) -> list:
    """Masses for which a laboratory-type transition (wall at R1 in the static region, eps_out(R1) = +1) in a
    Schwarzschild-de Sitter parent lands on a daughter with eps_out(R2) = +1: no exterior flip, so the embedding
    obstruction NG3 does not apply.  Empty for eta = 0 (Prop. G1)."""
    from .model import ThinWallModel
    out = []
    for M in np.geomspace(1e-5, 1.0, n) / chi:
        m = ThinWallModel(chi, kappa, M, eta)
        tps, hz = m.turning_points(), m.horizons_out()
        if len(tps) < 2 or len(hz) < 2:
            continue
        R1, R2 = tps[0], tps[-1]
        if not (m.A_out(R1) > 0 and R1 > hz[0] and m.f_out(R1) > 0):
            continue
        if m.A_out(R2) > 0:
            out.append(float(M))
    return out


def sds_nonflip_threshold(kappa: float, chi: float = 1.0, iters: int = 25) -> float:
    """Smallest eta/chi at which the non-flip window is non-empty (bisection; Proposition-N)."""
    lo, hi = 0.0, 0.999 * np.sqrt(chi**2 + kappa**2)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if sds_nonflip_masses(kappa, mid, chi, 250):
            hi = mid
        else:
            lo = mid
    return float(hi / chi)
