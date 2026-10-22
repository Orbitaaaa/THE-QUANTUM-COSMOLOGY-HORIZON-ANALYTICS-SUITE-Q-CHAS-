"""Roadmap Stage 4 (task 30) and Stage 1 (P-T5): analytic suppression bounds.

Everything here is for laboratory transitions: Schwarzschild exterior (eta = 0), regime A.  Notation as in
model.py; gamma, E are the dimensionless parameters, E_AB = 2 - 2 gamma is the A/B boundary.

* ``sine_bound``      Lemma S1:  B_an >= B_sq >= B_mem = (2 kappa/G) int R^2 sqrt(V/(f_in f_out)) dR.
* ``elementary_bound`` Theorem S2: B >= 2 pi kappa^2 R1^2 f_out(R1) / (G (L^2 + 4GM/R1^3) sqrt(f_in(R1) f_out(R2))).
* ``Phi``             the universal thin-barrier function of Theorem S3, B = B(0) Phi(E) (1 + O(kappa^2/chi^2)).
* ``M_sd``            Theorem S4: the transition is double-flip (FGG-unambiguous) iff kappa < chi and M < M_sd.
* ``omega_condition`` Theorem S5: a pointwise sufficient condition for Delta t_E^W < 0 (hence dB/dM > 0 in the
                      double-flip class, B >= B_CdL).
* ``laboratory_law``  Corollary S6: B_lab = 4 pi^2 sigma^2 R_M^2 / rho  (R_M >> R_CdL, kappa << chi).

See ROADMAP_STAGE1_2_4.md, Part A.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from .model import ThinWallModel
from . import tunnelling as T


# ------------------------------------------------------------------ Lemma S1 / Theorem S2
def sine_bound(model: ThinWallModel) -> float:
    """B_mem = (2/G) int R sin(theta) dR with sin(theta) = kappa R sqrt(V/(f_in f_out)).
    theta >= sin(theta) and arccos|C| >= sqrt(1 - C^2) = sin(theta) give B_an >= B_sq >= B_mem."""
    Ra, Rb = T.barrier(model)

    def f(R):
        fi, fo = float(model.f_in(R)), float(model.f_out(R))
        V = max(float(model.V(R)), 0.0)
        return R * R * np.sqrt(V / (fi * fo))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        v = quad(f, Ra, Rb, limit=600, epsabs=0.0, epsrel=1e-11)[0]
    return 2 * model.kappa * v / model.G


def elementary_bound(model: ThinWallModel) -> float:
    """Theorem S2 (eta = 0, regime A).  On [R1, R'] with A_out(R') = -sqrt(f_out(R1)):
    R^2 >= R1^2, V >= f_out(R1) - A_out^2, f_in f_out <= f_in(R1) f_out(R2), |A_out'| <= (L^2 + 4GM/R1^3)/(2 kappa);
    the w = A_out integral of sqrt(f1 - w^2) is pi f1 / 2."""
    if model.eta != 0:
        raise ValueError("elementary_bound needs a Schwarzschild exterior")
    R1, R2 = T.barrier(model)
    GM = model.G * model.M
    return float(2 * np.pi * model.kappa**2 * R1**2 * model.f_out(R1)
                 / (model.G * (model.L2 + 4 * GM / R1**3) * np.sqrt(model.f_in(R1) * model.f_out(R2))))


# ------------------------------------------------------------------ Theorem S3: thin barrier
def _g1(z):
    return (z**3 - 1) ** 2 / z**4


def Phi(E: float) -> float:
    """Phi(E) = (16 / (pi E^2)) int z^2 sqrt(E - (z^3-1)^2/z^4) dz over {g_1 < E}.
    Computed in the variable s = (z^3-1)/z^2 (monotone, g_1 = s^2), z^2 dz = z^5/(z^3+2) ds,
    with s = sqrt(E) sin(phi): no endpoint singularities."""
    rE = np.sqrt(E)

    def z_of_s(s):
        # z^3 - s z^2 - 1 = 0 has exactly one positive root
        lo, hi = 1e-12, max(2.0, abs(s) + 2.0)
        return brentq(lambda z: z**3 - s * z**2 - 1, lo, hi, xtol=1e-15, rtol=1e-15, maxiter=200)

    def integrand(phi):
        s = rE * np.sin(phi)
        z = z_of_s(s)
        return z**5 / (z**3 + 2) * E * np.cos(phi) ** 2

    v = quad(integrand, -np.pi / 2, np.pi / 2, limit=400, epsabs=0.0, epsrel=1e-12)[0]
    return 16 / (np.pi * E**2) * v


def Phi_small_E(E):
    """Phi(E) = 8/(3E) + 2/3 + O(E)   (h(s) = z^5/(z^3+2) = (1 + 4s/3 + s^2)/3 + O(s^3))."""
    return 8 / (3 * np.asarray(E)) + 2 / 3


def Phi_large_E(E):
    """Phi(E) = 1 + (16/pi) E^(-3/2) + o(E^(-3/2))   (first-order effect of the -2/z term of g)."""
    return 1 + 16 / np.pi * np.asarray(E) ** -1.5


def thin_barrier_B(model: ThinWallModel) -> float:
    """Leading thin-barrier prediction B(0) Phi(E), with B(0) the exact closed form T1."""
    B0 = T.B_canonical(ThinWallModel(model.chi, model.kappa, 0.0, model.eta, model.G))
    return B0 * Phi(model.E)


def R_flip(model: ThinWallModel) -> float:
    return (2 * model.G * model.M / model.L2) ** (1 / 3)


def thin_barrier_asymptote(model: ThinWallModel) -> float:
    """E -> 0 end of Theorem S3: B -> (2 pi / 3) (kappa/chi)^2 R_f^2 / G  = (2/3)(kappa/chi)^2 A_f / 4G."""
    return 2 * np.pi / 3 * (model.kappa / model.chi) ** 2 * R_flip(model) ** 2 / model.G


# ------------------------------------------------------------------ Theorem S4: flip class
def M_sd(model: ThinWallModel) -> float:
    """Single/double-flip boundary: R2 = 1/chi, i.e. 2 G M chi = 1 - kappa^2/chi^2 (eta = 0).
    Double flip iff kappa < chi and M < M_sd; returns 0 if kappa >= chi."""
    if model.chi <= model.kappa:
        return 0.0
    return (model.chi**2 - model.kappa**2) / (2 * model.G * model.chi**3)


def E_sd(gamma: float) -> float:
    """Dimensionless form of M_sd: E_sd = (1 - gamma^2) / gamma^(2/3) (gamma > 0)."""
    return (1 - gamma**2) / gamma ** (2 / 3) if gamma > 0 else np.inf


def flip_class_predicted(model: ThinWallModel) -> str:
    return "double" if model.M < M_sd(model) else "single"


# ------------------------------------------------------------------ Theorem S5: sign of Delta t_E
def _dimless_barrier(gamma: float, E: float):
    g = lambda z: z**2 - 2 * gamma / z + 1 / z**4
    zs = ThinWallModel.z_star(gamma)
    z1 = brentq(lambda z: g(z) - E, 1e-10, zs, xtol=1e-15)
    z2 = brentq(lambda z: g(z) - E, zs, 10 * np.sqrt(E) + 10, xtol=1e-15)
    return z1, z2


def omega(z):
    """omega(z) = z^2/12 + 7/(6 z): Delta t_E^W < 0 whenever E / E_AB >= max(omega(z1), omega(z2))."""
    return z**2 / 12 + 7 / (6 * z)


def omega_condition(gamma: float, E: float) -> bool:
    z1, z2 = _dimless_barrier(gamma, E)
    return E / (2 - 2 * gamma) >= max(omega(z1), omega(z2))


def pointwise_D_min(gamma: float, E: float, n: int = 4001) -> float:
    """min over the barrier of D = 6z(W - s^2)W - b(z^3+2)(W/2 - s^2) (b = E_AB, W = E - b/z, s = (z^3-1)/z^2).
    D >= 0 on the barrier implies Delta t_E^W < 0 (proof of Theorem S5)."""
    z1, z2 = _dimless_barrier(gamma, E)
    b = 2 - 2 * gamma
    z = np.linspace(z1, z2, n)[1:-1]
    W = E - b / z
    s = (z**3 - 1) / z**2
    return float((6 * z * (W - s**2) * W - b * (z**3 + 2) * (W / 2 - s**2)).min())


def delta_tE_dimless(gamma: float, E: float) -> float:
    """J(gamma, E) = int (1 - z^3) / (z^2 f_out sqrt(E - g)) dz, proportional (positive factor) to Delta t_E^W.
    Computed in s = sqrt(E) sin-type variables to remove the turning-point singularities."""
    z1, z2 = _dimless_barrier(gamma, E)
    b = 2 - 2 * gamma
    g = lambda z: z**2 - 2 * gamma / z + 1 / z**4

    def f(t):
        z = z1 + (z2 - z1) * (1 - np.cos(t)) / 2
        dz = (z2 - z1) * np.sin(t) / 2
        r = E - g(z)
        return (1 - z**3) / (z**2 * (1 - b / (E * z))) * dz / np.sqrt(max(r, 1e-300))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return quad(f, 0, np.pi, limit=800)[0]


# ------------------------------------------------------------------ Corollary S6: laboratory law
def laboratory_law(E_vac_GeV: float, E_wall_GeV: float, M_planck: float) -> dict:
    """B_lab = 4 pi^2 sigma^2 R_M^2 / rho with R_M = (3M / 4 pi rho)^(1/3)  (Planck units),
    and the equivalent form 4 pi sqrt(2/27) sqrt(B_flat) (3M/4pi)^(2/3) rho^(-1/6) with
    B_flat = 27 pi^2 sigma^4 / (2 rho^3)."""
    from .observability import planck_params
    p = planck_params(E_vac_GeV, E_wall_GeV)
    rho, sigma = p["rho"], p["sigma"]
    RM = (3 * M_planck / (4 * np.pi * rho)) ** (1 / 3)
    R0 = 3 * sigma / rho
    B_flat = 27 * np.pi**2 * sigma**4 / (2 * rho**3)
    B_law = 4 * np.pi**2 * sigma**2 * RM**2 / rho
    B_phi = B_flat * Phi((R0 / RM) ** 2)
    return dict(R_M_lP=float(RM), R_CdL_lP=float(R0), x=float(RM / R0), B_flat=float(B_flat),
                B_law=float(B_law), B_Phi=float(B_phi), kappa_over_chi=float(p["k"]),
                B_law_identity=float(4 * np.pi * np.sqrt(2 / 27) * np.sqrt(B_flat)
                                     * (3 * M_planck / (4 * np.pi)) ** (2 / 3) * rho ** (-1 / 6)))


def minimum_exponent(B_survive: float, M_planck: float, rho_planck: float) -> float:
    """Smallest thin-barrier laboratory exponent compatible with metastability B_CdL >= B_survive,
    a black hole of at least M_planck Planck masses and vacuum energy density at most rho_planck:
    B_min = 4 pi sqrt(2 B_survive / 27) (3 M / 4 pi)^(2/3) rho^(-1/6)."""
    return float(4 * np.pi * np.sqrt(2 * B_survive / 27) * (3 * M_planck / (4 * np.pi)) ** (2 / 3)
                 * rho_planck ** (-1 / 6))
