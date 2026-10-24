"""Round 6c (P-R7): 3+1 greybody factors for the Round 7 flux.

Regge-Wheeler equation for a massless field of spin s in Schwarzschild (G = 1),
    d^2 psi / dr*^2 + (omega^2 - V_l) psi = 0,   V_l = f [ l(l+1)/r^2 + 2M(1 - s^2)/r^3 ],   f = 1 - 2M/r,
with purely ingoing data at the horizon, psi = e^{-i omega r*}.  The Wronskian flux Im(psi* dpsi/dr*) = -omega is exact; at large r the
incident amplitude is read off adiabatically, |alpha|^2 = |sqrt(k) psi + i psi'/sqrt(k)|^2 / 4 with k = sqrt(omega^2 - V), and
Gamma_l(omega) = omega / |alpha|^2 (the adiabatic error falls like V'/k^2 and is below 1e-6 at the matching radius used).

* ``transmission``      Gamma_l(omega); checks: Gamma_0 -> 16 M^2 omega^2 (sigma_abs -> A_H, Das-Gibbons-Mathur) and
                        sum_l (2l+1) Gamma_l -> 27 M^2 omega^2 (geometric optics, sigma -> 27 pi M^2).
* ``hawking_rates``     3+1 Hawking power and number rate, sum_l (2l+1) int d omega/2pi {omega, 1} Gamma_l / (e^{omega/T} - 1), T = 1/8 pi M.
* ``branch_difference`` Proposition 6c-12.  Round 7's s-wave transient is a deficit of about two thermal units, int (F - F_H) du
                        = -D F_H/kappa with D ~ 2, i.e. a delay of the switch-on by D/kappa = 4 D M, set by the wall's geometry and
                        frequency-independent to leading order.  In 3+1 the same delay removes N_dot * 4 D M quanta: IR-finite (Gamma ~ omega^2
                        at low frequency), unlike the 1+1 count, which diverges logarithmically.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad, solve_ivp


def _rstar(r, M):
    return r + 2 * M * np.log(r / (2 * M) - 1)


def transmission(l: int, omega: float, M: float = 1.0, s: int = 0) -> float:
    """Gamma_l(omega) for a massless spin-s field (s = 0, 1, 2) on Schwarzschild of mass M."""
    V = lambda r: (1 - 2 * M / r) * (l * (l + 1) / r**2 + 2 * M * (1 - s * s) / r**3)
    r0 = 2 * M * (1 + 1e-7)
    x0 = _rstar(r0, M)
    psi0 = np.exp(-1j * omega * x0)
    k0 = np.sqrt(omega**2 - V(r0) + 0j)
    y0 = [psi0.real, psi0.imag, (-1j * k0 * psi0).real, (-1j * k0 * psi0).imag]

    def rhs(r, y):
        f = 1 - 2 * M / r
        q = (V(r) - omega**2) / f
        return [y[2] / f, y[3] / f, q * y[0], q * y[1]]

    r1 = max(400 * M, 200 / omega, 30 * (l + 1) / omega)
    sol = solve_ivp(rhs, (r0, r1), y0, method="DOP853", rtol=1e-10, atol=1e-13, max_step=min(0.1 / omega, 10 * M))
    psi = sol.y[0, -1] + 1j * sol.y[1, -1]
    dpsi = sol.y[2, -1] + 1j * sol.y[3, -1]
    k = np.sqrt(omega**2 - V(r1))
    alpha2 = abs(np.sqrt(k) * psi + 1j * dpsi / np.sqrt(k)) ** 2 / 4
    return float(omega / alpha2)


def hawking_rates(M: float = 1.0, s: int = 0, l_max: int = 4, n_omega: int = 70) -> dict:
    """3+1 Hawking power P and number rate N_dot for one massless spin-s degree of freedom (per helicity), with the l-decomposition."""
    T = 1 / (8 * np.pi * M)
    om = np.linspace(2e-3 / M, 14 * T, n_omega)
    out, P, N = {}, 0.0, 0.0
    for l in range(max(s, 0), l_max + 1):
        G = np.array([transmission(l, w, M, s) for w in om])
        n = 1 / np.expm1(om / T)
        Pl = (2 * l + 1) * np.trapezoid(om * G * n, om) / (2 * np.pi)
        Nl = (2 * l + 1) * np.trapezoid(G * n, om) / (2 * np.pi)
        out[l] = dict(P=float(Pl), N=float(Nl))
        P += Pl
        N += Nl
    return dict(P=float(P), N_dot=float(N), P_M2=float(P * M**2), per_l=out, T=float(T),
                P_swave_1p1=float(np.pi * T**2 / 12))


def branch_difference(M: float = 1.0, D: float = 2.0, rates: dict = None) -> dict:
    """Proposition 6c-12: quanta and energy removed in 3+1 by a switch-on delay of D/kappa = 4 D M (Round 7: D ~ 1.9-2.2)."""
    rates = hawking_rates(M) if rates is None else rates
    delay = 4 * D * M
    return dict(delay=float(delay), quanta=float(rates["N_dot"] * delay), energy=float(rates["P"] * delay),
                quanta_swave=float(rates["per_l"][0]["N"] * delay))
