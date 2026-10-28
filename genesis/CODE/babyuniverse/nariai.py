"""Round 6c: the Nariai-throat daughter, quantified (the answers table of ROADMAP_ROUND6C).

Configuration (``throat_evolution.nariai_false_vacuum``): a Schwarzschild-de Sitter hole of mass M = M_N (1 - delta), M_N = r_N / 3G,
r_N = 1/sqrt(Lambda) = l_H / sqrt 3, inside a false-vacuum layer whose outer wall lies in the static band (r_b, r_c); true vacuum outside.
The data are normal.  An outgoing negative-energy flux N(u) on the ingoing ray v = 0 tips it.

* ``weighted_deficit``      Proposition 6c-12 (universal tipping law) [Numerical].  All 25 tipping runs (Gaussian pulses of 8 amplitudes and
                            4 emission times, negative/positive pairs with 5 separations, 6 extended top-hat sources) are separated by one
                            criterion: the daughter grows iff
                                W = -int N(u) e^{-u/tau} du > W_th,   tau = 3.2 r_N (2.8-3.9),   W_th = 1.85 (1.73-1.98) x 10^-3 r_N / G,
                            in the static (geodesic, f'(r_N) = 0) frame of the corner, for delta <= 1e-8.  The law is linear response
                            of an instability: what counts is the exponentially weighted NEC deficit, not its total.  Classical scale
                            invariance makes W_th and tau proportional to r_N.
* ``invariant_cost``        The gauge-invariant size of the pulse is the Misner-Sharp jump: m = (r/2)(1 + 4 r_u r_v e^{-2s}) and
                            Delta r_u = (G/r) int|N| du give dm = 2 r_v e^{-2s} int|N| du = sqrt(f_0) W (r_v = sqrt(f_0)/2 at the
                            corner, f_0 = 2 delta / 3): (2.5-5)e-7 M_N at delta = 1e-8, and ~(0.3-2)e-2 sqrt(delta) M_N for
                            delta = 1e-10 to 1e-3.  The local-frame number W = 0.55% M_N is the figure quoted in the Round 6c
                            checkpoint ("0.6-0.9% of M_N", which was this local energy, not a mass).  Flat space needs local energy >= r/2G
                            just to flip a sphere (Theorem G1), i.e. 270 times more than W at the same radius.
* ``ghost_condensate_window``  Proposition 6c-14 [Proposition, conditional on the ghost-condensate EFT].  A ghost condensate at rest in the
                            corner frame with rho + p < 0 for a time T gives N = pi r^2 |rho + p| on the ray, so tipping needs
                            |rho + p| >= W_th / (pi r_N^2 tau (1 - e^{-T/tau})) = 4.6e-3 rho_f (T >> tau).  Creminelli-Luty-Nicolis-
                            Senatore (hep-th/0606090) show the softly broken ghost condensate violates the NEC stably up to
                            |dH/dt| ~ H^2, i.e. |rho + p| <~ H^2/(4 pi G): with H = 1/l_H (1/r_N) the ceiling is 145 (430) times the need.
* ``partner_vacuum``        The Nariai condition 9 G^2 M^2 Lambda = 1 ties the hole to the false vacuum: rho_f = Lambda/8 pi G = 1/(72 pi G^3 M^2)
                            (Planck units 1/(72 pi M^2)), i.e. rho_f^{1/4} = 0.329 GeV (M_sun / M)^{1/2}; a GUT-scale vacuum (10^16 GeV)^4
                            pairs with a 2.2 g hole, an electroweak one (100 GeV)^4 with 3.6 Earth masses.  (Corrected in Round 6j: the
                            formula had 24 pi, a factor 3 in rho_f; accessibility.py always used the correct 1/(8 pi r_N^2), M_N = r_N/3.)
* ``preparation``           Proposition 6c-15 (what the initial geometry needs) [Numerical + Theorem pieces]: the outer wall inside the
                            static band, |1 - M/M_N| = delta <~ 1e-4 (W_th grows only like delta^0.15, the band width like sqrt delta), and
                            false vacuum continuing inside the black-hole horizon to u >~ 3 u_hole along the ingoing ray (inner wall at
                            u = 100, 150, 195: grows and is hidden; at u = 60, 40, 25, 12, 6: no growth even at 4x the pulse).  The hole must
                            have formed inside the false vacuum or the vacuum must be lifted in situ across the horizon.
"""

from __future__ import annotations

import numpy as np

TAU_OVER_RN = 3.2           # fitted window 48-68 in u-units at r_N = 17.32 (2.8-3.9 r_N)
WTH_OVER_RN = 1.85e-3       # fitted window 0.0300-0.0342 at r_N = 17.32

M_SUN_KG = 1.98892e30
M_PLANCK_KG = 2.176434e-8
M_PLANCK_GEV = 1.220890e19


def weighted_deficit(N, tau: float, u: np.ndarray = None) -> float:
    """W = -int N(u) e^{-u/tau} du for a flux profile N(u) (callable) on the ingoing data ray."""
    u = np.linspace(-10.0, 40 * tau, 200001) if u is None else u
    return float(-np.trapezoid(N(u) * np.exp(-u / tau), u))


def tipping_law_fit(profiles, grew, taus=None, u=None) -> dict:
    """Range of (tau, W_th) for which W(tau) > W_th reproduces every outcome.  profiles: callables N(u); grew: booleans."""
    taus = np.arange(10.0, 150.0, 1.0) if taus is None else taus
    u = np.linspace(-10.0, 600.0, 122001) if u is None else u
    g = np.asarray(grew, bool)
    ok = []
    for tau in taus:
        W = np.array([weighted_deficit(f, tau, u) for f in profiles])
        lo, hi = W[~g].max() if (~g).any() else -np.inf, W[g].min() if g.any() else np.inf
        if hi > lo:
            ok.append((float(tau), float(lo), float(hi)))
    return dict(windows=ok, tau_range=(ok[0][0], ok[-1][0]) if ok else None)


def invariant_cost(W_local: float, delta: float) -> float:
    """Misner-Sharp mass jump sqrt(f_0) W of a pulse with local-frame energy W (f_0 = 2 delta/3 + O(delta^2))."""
    return float(np.sqrt(2 * delta / 3) * W_local)


def flat_space_flip(r: float) -> float:
    """Local energy r/2G that a pulse must carry to anti-trap a sphere of radius r with f = 1 (Theorem G1)."""
    return float(r / 2)


def ghost_condensate_window(r_N: float = 1.0, T_over_tau: float = np.inf) -> dict:
    """Proposition 6c-14: required |rho + p| / rho_f against the CLNS stability ceiling H^2/(4 pi G) (G = 1)."""
    tau = TAU_OVER_RN * r_N
    W = WTH_OVER_RN * r_N
    rho_f = 1 / (8 * np.pi * r_N**2)
    need = W / (np.pi * r_N**2 * tau * (1 - np.exp(-T_over_tau)))
    ceil_lH = (1 / (3 * r_N**2)) / (4 * np.pi)
    ceil_rN = (1 / r_N**2) / (4 * np.pi)
    return dict(need_over_rho_f=float(need / rho_f), ceiling_over_need_lH=float(ceil_lH / need),
                ceiling_over_need_rN=float(ceil_rN / need))


def partner_vacuum(M_kg: float) -> float:
    """rho_f^{1/4} in GeV for which a hole of mass M is at the Nariai point: rho_f = 1/(72 pi M^2) (Planck units)."""
    m = M_kg / M_PLANCK_KG
    return float((1 / (72 * np.pi * m * m)) ** 0.25 * M_PLANCK_GEV)


def partner_mass(rho_quarter_GeV: float) -> float:
    """Inverse of partner_vacuum: the Nariai hole mass in kg for a false vacuum rho_f = (rho_quarter_GeV)^4."""
    rho = (rho_quarter_GeV / M_PLANCK_GEV) ** 4
    return float(M_PLANCK_KG / np.sqrt(72 * np.pi * rho))
