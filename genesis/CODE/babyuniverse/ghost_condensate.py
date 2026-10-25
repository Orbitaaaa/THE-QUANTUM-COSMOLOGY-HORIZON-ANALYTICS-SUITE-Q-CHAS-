"""Round 6d: the ghost condensate in the near-Nariai throat (background, NEC phase, stability of perturbations).

EFT (Arkani-Hamed, Cheng, Luty, Mukohyama 2004; Creminelli, Luty, Nicolis, Senatore 2006):  L = M^4 P(X) - (alpha / 2M^2)(box phi)^2,
X = -(d phi)^2 / M^4, condensate at P'(1) = 0.  Units G = c = hbar = 1 (m_P = 1), M_Pl^2 = 1/(8 pi).

Background in the Nariai region.  Near the Nariai radius the geometry is dS_2(r_N) x S^2(r_N); in flat dS_2 slicing
ds^2 = -dt^2 + e^{2t/r_N} dx^2 + r_N^2 dOmega^2, and dt is a geodesic congruence through the corner (f'(r_N) = 0 there).  phi = M^2 t has X = 1.
* ``nec_phase``          Proposition 6d-1 [Theorem in the EFT].  Shift symmetry gives d/dt (a P'(X) phi_dot) = 0 with a = e^{t/r_N} (only the
                         throat direction expands), so a displaced condensate relaxes as P' = eps(t) = eps_0 e^{-t/r_N}, and
                         rho + p = 2 X P' M^4 = 2 eps_0 M^4 e^{-t/r_N}: a NEC-violating phase (eps_0 < 0) switches itself off on the dS_2 Hubble time.
                         In the corner frame it gives the null fluxes N_u = N_v = pi r^2 (rho + p), so the tipping law's weighted deficit is
                         W = 2 pi r_N^2 |eps_0| M^4 / (1/r_N + 1/tau).
* ``mode_growth``        Perturbations pi with L_2 = A pi_dot^2 - eps (grad pi)^2 - (alpha/2M^2)(lap pi)^2, A = 2P'' + eps, plus the Jeans term of
                         the condensate (ACLM: omega^2 contains -(alpha M^2 / 2 M_Pl^2) k^2).  On dS_2 x S^2 a mode (k, l) has
                         k_ph^2 = k^2 e^{-2t/r_N} + l(l+1)/r_N^2 and obeys
                             pi_tt + H pi_t + [eps(t) k_ph^2 + (alpha/2M^2) k_ph^4 - (alpha M^2/2M_Pl^2) A k_ph^2] / A pi = 0,   H = 1/r_N.
                         Integrated numerically; returns the largest amplification over a (k, l) grid.
* ``window``             Proposition 6d-2 [Proposition, conditional on the EFT].  Gradient modes grow at most at |eps| M / sqrt(2 alpha A); with the
                         natural decay their total growth is |eps_0| M r_N / sqrt(2 alpha A).  The Jeans rate is 2 pi sqrt(alpha) M^3 (Planck units).
                         Tipping (W >= W_th) with at most one e-fold of any instability over 3 tau requires
                             1.9e-4 / sqrt(alpha P'') <= M^3 r_N <= 1 / (6 pi sqrt(alpha) tau/r_N),
                         a window of width ~ 80 sqrt(P'') in M^3, inside the EFT (M r_N >> 1, |eps_0| << 1) once r_N >~ 10^3 l_P.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

M_PL2 = 1 / (8 * np.pi)          # reduced Planck mass squared in G = 1 units
TAU_OVER_RN = 3.2
WTH_OVER_RN = 1.85e-3


def nec_phase(eps0: float, M: float, r_N: float, t):
    """rho + p (t) = 2 eps0 M^4 exp(-t/r_N) of the relaxing condensate (Proposition 6d-1)."""
    return 2 * eps0 * M**4 * np.exp(-np.asarray(t, float) / r_N)


def weighted_deficit(eps0: float, M: float, r_N: float, tau_over_rN: float = TAU_OVER_RN) -> float:
    """W of the natural NEC phase (one null direction): int pi r^2 |rho + p| e^{-t/tau} dt = 2 pi r_N^2 |eps0| M^4 / (1/r_N + 1/tau)."""
    tau = tau_over_rN * r_N
    return float(2 * np.pi * r_N**2 * abs(eps0) * M**4 / (1 / r_N + 1 / tau))


def eps0_for_tipping(M: float, r_N: float, margin: float = 1.0) -> float:
    """|eps0| for which the natural phase reaches W = margin * W_th (W_th = 1.85e-3 r_N)."""
    tau = TAU_OVER_RN * r_N
    return float(margin * WTH_OVER_RN * r_N * (1 / r_N + 1 / tau) / (2 * np.pi * r_N**2 * M**4))


def mode_growth(M: float, eps0: float, r_N: float, alpha: float = 1.0, Ppp: float = 1.0, t_end_over_rN: float = 10.0,
                ks=None, ls=None, jeans: bool = True, eps_t=None, r_t=None) -> dict:
    """Largest amplification max_t |pi|/|pi_0| (energy-normalised) over modes (k, l) for t in [0, t_end]; pi_0 = 1, pi_t0 = 0.
    eps_t(t): optional NEC-phase history (default eps0 e^{-t/r_N}); r_t(t): optional sphere radius for the l modes (default r_N)."""
    H = 1 / r_N
    A = 2 * Ppp
    C = alpha / (2 * M * M)
    J = alpha * M * M / (2 * M_PL2) * A if jeans else 0.0
    kmax = np.sqrt(max(abs(eps0), 1e-300) / C) * 3 + np.sqrt(abs(J) / C) * 3 + 10 / r_N
    ks = np.geomspace(1e-2 / r_N, kmax, 60) if ks is None else ks
    ls = np.unique(np.concatenate([np.arange(0, 6), np.geomspace(6, max(7, kmax * r_N), 25).astype(int)])) if ls is None else ls
    best = (1.0, None)
    for l in ls:
        for k in ks:
            def rhs(t, y):
                rr = r_N if r_t is None else r_t(t)
                kp2 = k * k * np.exp(-2 * t * H) + l * (l + 1) / (rr * rr)
                eps = eps0 * np.exp(-t * H) if eps_t is None else eps_t(t)
                w2 = (eps * kp2 + C * kp2 * kp2 - J * kp2) / A
                return [y[1], -H * y[1] - w2 * y[0]]
            sol = solve_ivp(rhs, (0, t_end_over_rN * r_N), [1.0, 0.0], rtol=1e-6, atol=1e-9, max_step=r_N / 20)
            g = float(np.max(np.abs(sol.y[0])))
            if g > best[0]:
                best = (g, (float(k * r_N), int(l)))
    return dict(max_amplification=best[0], at_k_rN_l=best[1])


def window(r_N: float, alpha: float = 1.0, Ppp: float = 1.0, efolds: float = 1.0, t_over_tau: float = 3.0) -> dict:
    """Proposition 6d-2: the range of M^3 r_N (Planck units) where the NEC phase tips the throat and every instability stays below
    `efolds` e-folds: gradient (needs |eps0| M r_N / sqrt(2 alpha A) <= efolds with |eps0| from W_th) and Jeans (2 pi sqrt(alpha) M^3 T)."""
    A = 2 * Ppp
    tau = TAU_OVER_RN * r_N
    # |eps0| M^4 = K / r_N^2 with K = W_th (1 + r_N/tau) / (2 pi)
    K = WTH_OVER_RN * (1 + 1 / TAU_OVER_RN) / (2 * np.pi)
    lo = K * np.sqrt(2 * alpha * A) ** -1 / efolds          # M^3 r_N >= K / (efolds sqrt(2 alpha A))  ... from |eps0| M r_N <= efolds sqrt(2 alpha A)
    lo = K / (efolds * np.sqrt(2 * alpha * A))
    hi = efolds / (2 * np.pi * np.sqrt(alpha) * t_over_tau * TAU_OVER_RN)
    return dict(M3rN_min=float(lo), M3rN_max=float(hi), width=float(hi / lo), exists=bool(hi > lo))
