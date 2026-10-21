"""Thin-wall junction model: de Sitter interior, Schwarzschild-de Sitter exterior.

Geometry (JUNCTION_CONDITIONS.md, CLASSICAL_MODEL.md)
-----------------------------------------------------
Interior  ds^2 = -f_in dt^2 + dr^2/f_in + r^2 dOmega^2,   f_in  = 1 - chi^2 r^2
Exterior  ds^2 = -f_out dt^2 + dr^2/f_out + r^2 dOmega^2, f_out = 1 - 2GM/r - eta^2 r^2
Wall      S_ab = -sigma h_ab,  kappa = 4 pi G sigma > 0.

Israel junction condition (theta-theta component, derived in symbolic.py):

    eps_in sqrt(f_in + Rdot^2) - eps_out sqrt(f_out + Rdot^2) = kappa R .

On every solution the signed square roots are rational functions of R:

    A_in(R)  := eps_in  sqrt(f_in  + Rdot^2) = (f_in - f_out + kappa^2 R^2) / (2 kappa R)
    A_out(R) := eps_out sqrt(f_out + Rdot^2) = (f_in - f_out - kappa^2 R^2) / (2 kappa R)

and Rdot^2 + V(R) = 0 with V = f_out - A_out^2 = f_in - A_in^2.  Because kappa > 0
the signs eps = sign(A) are *determined* by R: the squared equation does not
admit spurious sign assignments once sigma > 0 is imposed (the same V arises
for sigma < 0 with both signs reversed; see JUNCTION_CONDITIONS.md Sec. 5).

Dimensionless form (valid for all chi, eta, kappa, M > 0):

    R = (2GM/L^2)^(1/3) z,  L^4 = (chi^2 - eta^2 + kappa^2)^2 + 4 kappa^2 eta^2,
    V = 1 - g(z)/E,  g(z) = z^2 - 2 gamma / z + 1/z^4,
    gamma = (chi^2 - eta^2 - kappa^2)/L^2,  E = 4 kappa^2 / (L^(8/3) (2GM)^(2/3)),
    eps_out = sign(1 - beta_out z^3),  beta_out = (chi^2 + kappa^2 - eta^2)/L^2,
    eps_in  = sign(1 - gamma z^3).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class ThinWallModel:
    """Spherical thin domain wall separating dS (inside) from SdS (outside).

    Parameters are the interior Hubble rate ``chi``, the wall parameter
    ``kappa = 4 pi G sigma``, the exterior mass ``M`` and the exterior Hubble
    rate ``eta``.  ``G`` is Newton's constant (default 1).
    """

    chi: float
    kappa: float
    M: float
    eta: float = 0.0
    G: float = 1.0

    def __post_init__(self):
        if self.kappa <= 0:
            raise ValueError("wall tension must be positive (kappa > 0)")
        if self.chi < 0 or self.eta < 0 or self.M < 0:
            raise ValueError("chi, eta, M must be non-negative")

    # ------------------------------------------------------------------ physics
    @property
    def sigma(self) -> float:
        return self.kappa / (4 * np.pi * self.G)

    @property
    def rho_in(self) -> float:
        return 3 * self.chi**2 / (8 * np.pi * self.G)

    @property
    def rho_out(self) -> float:
        return 3 * self.eta**2 / (8 * np.pi * self.G)

    @property
    def Lambda_in(self) -> float:
        return 3 * self.chi**2

    @property
    def Lambda_out(self) -> float:
        return 3 * self.eta**2

    def f_in(self, R):
        R = np.asarray(R, dtype=float)
        return 1.0 - self.chi**2 * R**2

    def f_out(self, R):
        R = np.asarray(R, dtype=float)
        return 1.0 - 2 * self.G * self.M / R - self.eta**2 * R**2

    def df_in(self, R):
        return -2 * self.chi**2 * np.asarray(R, dtype=float)

    def df_out(self, R):
        R = np.asarray(R, dtype=float)
        return 2 * self.G * self.M / R**2 - 2 * self.eta**2 * R

    # ----------------------------------------------------- on-shell square roots
    def A_in(self, R):
        """eps_in sqrt(f_in + Rdot^2) on shell (rational in R)."""
        R = np.asarray(R, dtype=float)
        k = self.kappa
        return (self.f_in(R) - self.f_out(R) + k**2 * R**2) / (2 * k * R)

    def A_out(self, R):
        """eps_out sqrt(f_out + Rdot^2) on shell (rational in R)."""
        R = np.asarray(R, dtype=float)
        k = self.kappa
        return (self.f_in(R) - self.f_out(R) - k**2 * R**2) / (2 * k * R)

    def eps_in(self, R):
        return np.sign(self.A_in(R))

    def eps_out(self, R):
        return np.sign(self.A_out(R))

    # -------------------------------------------------------- effective potential
    @property
    def L4(self) -> float:
        c2, e2, k2 = self.chi**2, self.eta**2, self.kappa**2
        return (c2 - e2 + k2) ** 2 + 4 * k2 * e2

    @property
    def L2(self) -> float:
        return float(np.sqrt(self.L4))

    def V(self, R):
        """Effective potential, Rdot^2 + V(R) = 0, explicit polynomial form."""
        R = np.asarray(R, dtype=float)
        c2, e2, k2 = self.chi**2, self.eta**2, self.kappa**2
        GM = self.G * self.M
        return (1.0 - self.L4 / (4 * k2) * R**2
                + GM * (c2 - e2 - k2) / (k2 * R) - GM**2 / (k2 * R**4))

    def V_squared_form(self, R):
        """Effective potential from its definition V = f_out - A_out^2."""
        return self.f_out(R) - self.A_out(R) ** 2

    def dV(self, R):
        R = np.asarray(R, dtype=float)
        c2, e2, k2 = self.chi**2, self.eta**2, self.kappa**2
        GM = self.G * self.M
        return (-self.L4 / (2 * k2) * R - GM * (c2 - e2 - k2) / (k2 * R**2)
                + 4 * GM**2 / (k2 * R**5))

    def junction_residual(self, R, Rdot, eps_in, eps_out):
        """Residual of the *unsquared* Israel equation."""
        a = np.sqrt(np.maximum(self.f_in(R) + Rdot**2, 0.0))
        b = np.sqrt(np.maximum(self.f_out(R) + Rdot**2, 0.0))
        return eps_in * a - eps_out * b - self.kappa * R

    def mass_function(self, R, Rdot, eps_in):
        """Exterior mass implied by (R, Rdot) and eps_in (BGG mass formula).

        M = (4 pi/3)(rho_in - rho_out) R^3 + eps_in 4 pi sigma R^2 sqrt(f_in+Rdot^2)
            - 8 pi^2 G sigma^2 R^3        (derived in symbolic.mass_formula)
        """
        G = self.G
        return ((self.chi**2 - self.eta**2) * R**3 / (2 * G)
                + eps_in * self.kappa * R**2 / G * np.sqrt(self.f_in(R) + Rdot**2)
                - self.kappa**2 * R**3 / (2 * G))

    # ------------------------------------------------------- dimensionless form
    @property
    def gamma(self) -> float:
        return (self.chi**2 - self.eta**2 - self.kappa**2) / self.L2

    @property
    def beta_out(self) -> float:
        return (self.chi**2 + self.kappa**2 - self.eta**2) / self.L2

    @property
    def E(self) -> float:
        if self.M == 0:
            return np.inf
        return 4 * self.kappa**2 / (self.L2 ** (4.0 / 3.0) * (2 * self.G * self.M) ** (2.0 / 3.0))

    @property
    def R_unit(self) -> float:
        """Length unit (2GM/L^2)^(1/3) of the dimensionless variable z."""
        return (2 * self.G * self.M / self.L2) ** (1.0 / 3.0)

    @staticmethod
    def g(z, gamma):
        z = np.asarray(z, dtype=float)
        return z**2 - 2 * gamma / z + 1.0 / z**4

    @staticmethod
    def z_star(gamma) -> float:
        """Minimum of g: z^6 + gamma z^3 - 2 = 0."""
        return ((-gamma + np.sqrt(gamma**2 + 8)) / 2) ** (1.0 / 3.0)

    @classmethod
    def g_min(cls, gamma) -> float:
        return float(cls.g(cls.z_star(gamma), gamma))

    @property
    def M_critical(self) -> float:
        """Mass at which the barrier disappears (E = g_min): above it the
        branch is monotonic, below it bounded + bounce branches exist."""
        gm = self.g_min(self.gamma)
        return 4 * self.kappa**3 / (self.G * gm**1.5 * self.L4)

    # ----------------------------------------------------------- special radii
    def turning_points(self) -> List[float]:
        """Positive roots of V(R) = 0, refined by bracketing, ascending."""
        if self.M == 0:
            return [2 * self.kappa / self.L2]
        c2, e2, k2 = self.chi**2, self.eta**2, self.kappa**2
        GM = self.G * self.M
        coeffs = [-self.L4 / 4, 0.0, k2, GM * (c2 - e2 - k2), 0.0, 0.0, -GM**2]
        roots = np.roots(coeffs)
        cand = sorted(r.real for r in roots if abs(r.imag) < 1e-7 * max(1, abs(r)) and r.real > 0)
        out: List[float] = []
        for r in cand:
            lo, hi = r * (1 - 1e-6), r * (1 + 1e-6)
            try:
                if np.sign(self.V(lo)) != np.sign(self.V(hi)):
                    r = brentq(self.V, lo, hi, xtol=1e-15 * r, rtol=1e-14, maxiter=200)
            except ValueError:
                pass
            if not out or abs(r - out[-1]) > 1e-9 * r:
                out.append(float(r))
        return out

    def horizons_out(self) -> List[float]:
        """Positive roots of f_out: black-hole and cosmological horizons."""
        GM = self.G * self.M
        if self.eta == 0:
            return [2 * GM] if GM > 0 else []
        if GM == 0:
            return [1.0 / self.eta]
        roots = np.roots([-self.eta**2, 0.0, 1.0, -2 * GM])
        return sorted(float(r.real) for r in roots if abs(r.imag) < 1e-10 and r.real > 0)

    def horizon_in(self) -> Optional[float]:
        return 1.0 / self.chi if self.chi > 0 else None

    def flip_radius_in(self) -> Optional[float]:
        """Radius where A_in changes sign (eps_in flips), if any."""
        d = self.chi**2 - self.eta**2 - self.kappa**2
        if self.M == 0 or d <= 0:
            return None
        return (2 * self.G * self.M / d) ** (1.0 / 3.0)

    def flip_radius_out(self) -> Optional[float]:
        """Radius where A_out changes sign (eps_out flips), if any."""
        d = self.chi**2 + self.kappa**2 - self.eta**2
        if self.M == 0 or d <= 0:
            return None
        return (2 * self.G * self.M / d) ** (1.0 / 3.0)

    # ------------------------------------------------------------ conveniences
    def with_mass(self, M: float) -> "ThinWallModel":
        return ThinWallModel(self.chi, self.kappa, M, self.eta, self.G)

    def describe(self) -> dict:
        d = dict(chi=self.chi, kappa=self.kappa, M=self.M, eta=self.eta, G=self.G,
                 sigma=self.sigma, rho_in=self.rho_in, rho_out=self.rho_out,
                 L4=self.L4, gamma=self.gamma, beta_out=self.beta_out,
                 turning_points=self.turning_points(), horizons_out=self.horizons_out(),
                 horizon_in=self.horizon_in(), flip_in=self.flip_radius_in(),
                 flip_out=self.flip_radius_out())
        if self.M > 0:
            d.update(E=self.E, g_min=self.g_min(self.gamma), M_critical=self.M_critical,
                     z_star=self.z_star(self.gamma), R_unit=self.R_unit)
        return d


def model_from_dimensionless(gamma: float, E: float, kappa: float = 1.0,
                             G: float = 1.0) -> ThinWallModel:
    """Schwarzschild-exterior (eta = 0) model with prescribed (gamma, E).

    For eta = 0: gamma = (chi^2-kappa^2)/(chi^2+kappa^2), so
    chi^2 = kappa^2 (1+gamma)/(1-gamma); E fixes M.
    """
    if not -1 < gamma < 1:
        raise ValueError("need -1 < gamma < 1")
    chi = kappa * np.sqrt((1 + gamma) / (1 - gamma))
    L2 = chi**2 + kappa**2
    twoGM = (4 * kappa**2 / (E * L2 ** (4.0 / 3.0))) ** 1.5
    return ThinWallModel(chi=chi, kappa=kappa, M=twoGM / (2 * G), eta=0.0, G=G)
