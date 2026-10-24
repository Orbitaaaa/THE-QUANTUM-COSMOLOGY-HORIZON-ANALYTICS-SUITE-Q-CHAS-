"""Round 6c, Section II: the quantum physics of the magnetic throat (P-R6g, P-R6h).

* ``lll_spectrum``          Proposition 6c-6 (exact reduction).  The Dirac operator on a sphere of radius a threading q flux quanta has
                            exactly |q| zero modes (index theorem) and nonzero eigenvalues +-sqrt(n (n + |q|))/a with degeneracy |q| + 2n.
                            A massless 4D Dirac fermion on M_2 x S^2 is therefore |q| massless 2D Dirac fermions (c = |q| per flavour) plus a
                            tower with 2D masses >= sqrt(1 + |q|)/a.  A 4D mass m_f passes unchanged to the 2D modes.
* ``strong_qnec_curved``    Proposition 6c-7 (strong QNEC on a curved 2D background).  For the conformal vacuum of a coordinate X(u) on
                            ds^2 = -e^{2 sigma} du dv, with S_out = (c/12)(2 sigma - ln X') + (endpoint-independent terms) and the Polyakov
                            stress tensor T_uu = (c/12 pi)(sigma_uu - sigma_u^2) + (c/24 pi)(y'^2/2 - y''), y = ln X', the strong QNEC
                            T_ll = (1/2 pi)(S_ll + (6/c) S_l^2) holds with EQUALITY in the affine parameter l (dl = e^{2 sigma} du): the
                            Weyl anomaly is part of the saturation.  Verified symbolically.
* ``schwarzian_fluctuation``  Proposition 6c-8 (where semiclassics ends).  Schwarzian Z ~ beta^{-3/2} exp(2 pi^2 C / beta), C = Q^3: exactly
                            <E> = 2 pi^2 C T^2 + 3T/2 and Var E = 4 pi^2 C T^3 + 3 T^2/2, so dE/E = sqrt(2 T_brk/T + 1.5 T_brk^2/T^2).  The
                            anti-trapping budget at depth rho is E sinh^2(rho/Q); quantum energy fluctuations exceed it for
                            sinh^2(rho/Q) < dE/E.  Below T_brk the classical trapped/normal status of the deep throat is undefined.
* ``weighted_anec_bound``   Proposition 6c-10 (exponentially weighted ANEC, P-A1 with exact constants).  For a unitary 2D CFT in a state
                            that is the vacuum before u_0, the Fewster-Hollands null QEI  int T f^2 >= -(c/12 pi) int f'^2  with
                            f^2 = e^{-(u-u_0)/tau} (switched on adiabatically inside the vacuum region) gives
                                int T_uu e^{-(u-u_0)/tau} du >= -c/(48 pi tau),
                            and the conformal-map states of ``throat_evolution.conformal_pulse`` with A = 1/tau, L -> infinity saturate it.
* ``quantum_tipping_c``     Corollary 6c-11 (a QFT cannot tip a Nariai throat inside the effective theory).  The throat responds to a
                            source through W = int (-T) e^{-u/tau} (Proposition 6c-12: tau = 3.2 r_N, W_th = 1.85e-3 r_N/G); a CFT supplies
                            at most c/(48 pi tau).  Tipping needs c >= 48 pi tau W_th = 0.89 (r_N/l_P)^2: the species length sqrt(c) l_P
                            reaches 0.95 r_N, where the semiclassical description of the throat ends.  The simulations agree: no
                            conformal-map pulse with c <= 500 = 1.7 r_N^2 tips the l_H = 30 throat, and the Bunch-Davies back-reaction of
                            the same c fields raises the threshold (classical pulse W = 0.035 grows at c <= 150, not at c = 300) and never
                            tips the throat by itself (c = 30-1000).
* ``nariai_fluctuation``    Assessment.  Quantum fluctuations of the Nariai radion are delta r / r ~ l_P / r_N; the classical tipping
                            pulse dominates them when f_th >> l_P/r_N, i.e. r_N >> 170 l_P, and the one-loop (trace-anomaly) correction
                            k / r_N^2 = G c / (12 pi r_N^2) is below 1% for c < 0.38 (r_N/l_P)^2.
* ``monopole_charge``       Magnetic charge in flux quanta: Q = q l_P / (2 sqrt(alpha)), so q = 2 sqrt(alpha) Q / l_P (Dirac quantisation,
                            Gaussian units).  A throat of radius Q ~ r carries q ~ 0.17 r/l_P quanta and c_LLL = q N_f lowest-Landau-level
                            2D fermions (Maldacena 2020): linear in r, far below the (r/l_P)^2 a quantum tip would need.

Units G = c = hbar = k_B = 1.
"""

from __future__ import annotations

import numpy as np
import sympy as sp


def lll_spectrum(q: int, a: float = 1.0, n_max: int = 5) -> dict:
    """Dirac spectrum on S^2 (radius a) with q flux quanta: zero modes and the first levels (eigenvalue, degeneracy per sign)."""
    q = abs(int(q))
    levels = [(float(np.sqrt(n * (n + q)) / a), q + 2 * n) for n in range(1, n_max + 1)]
    return dict(zero_modes=q, levels=levels, gap=float(np.sqrt(1 + q) / a))


def lll_central_charge(q: int, N_f: int = 1) -> int:
    """c of the massless 2D sector: |q| 2D Dirac fermions (c = 1 each) per 4D Dirac flavour."""
    return abs(int(q)) * int(N_f)


def pulse_inside_lll(q: int, eta: float, a: float = 1.0) -> bool:
    """A pulse of local scale a/eta excites only the lowest Landau level if eta/a << gap = sqrt(1 + q)/a."""
    return bool(eta < 0.1 * np.sqrt(1 + abs(q)))


def strong_qnec_curved() -> dict:
    """Symbolic proof that conformal vacua saturate the strong 2D QNEC on ds^2 = -e^{2 sigma(u,v)} du dv, Weyl anomaly included."""
    u, c = sp.symbols("u c", positive=True)
    sig = sp.Function("sigma")(u)
    y = sp.Function("y")(u)                      # y = ln X'(u)
    T_uu = (c / (12 * sp.pi)) * (sp.diff(sig, u, 2) - sp.diff(sig, u) ** 2) \
        + (c / (24 * sp.pi)) * (sp.diff(y, u) ** 2 / 2 - sp.diff(y, u, 2))
    S = (c / 12) * (2 * sig - y)                 # u-dependent part of S_out for the half-line beyond the point
    e2s = sp.exp(2 * sig)
    # affine parameter l with dl = e^{2 sigma} du along the null line v = const
    S_l = S.diff(u) / e2s
    S_ll = S_l.diff(u) / e2s
    T_ll = T_uu / e2s**2
    residual = sp.simplify(T_ll - (S_ll + (6 / c) * S_l**2) / (2 * sp.pi))
    return dict(residual=residual, saturated=bool(residual == 0))


def schwarzian_fluctuation(T: float, Q: float) -> dict:
    """Exact Schwarzian canonical moments (C = Q^3) and the depth below which energy fluctuations exceed the anti-trapping budget."""
    C = Q**3
    T_brk = 1 / (2 * np.pi**2 * C)
    E_cl = 2 * np.pi**2 * C * T**2
    E_mean = E_cl + 1.5 * T
    var = 4 * np.pi**2 * C * T**3 + 1.5 * T**2
    rel = np.sqrt(var) / E_cl
    rho_q = Q * np.arcsinh(np.sqrt(rel))
    return dict(T_brk=float(T_brk), E_cl=float(E_cl), E_mean=float(E_mean), dE_over_E=float(rel),
                dE_over_E_formula=float(np.sqrt(2 * T_brk / T + 1.5 * (T_brk / T) ** 2)), depth_quantum=float(rho_q))


def schwarzian_numeric(T: float, Q: float, n: int = 400000) -> dict:
    """Check of the closed-form moments against the density of states rho(E) ~ sinh(2 pi sqrt(2 C E)) (numerical quadrature)."""
    C = Q**3
    beta = 1 / T
    E_cl = 2 * np.pi**2 * C * T**2
    Emax = 50 * (E_cl + 10 * T)
    E = np.linspace(0, Emax, n)[1:]
    logw = 2 * np.pi * np.sqrt(2 * C * E) - beta * E
    w = np.exp(logw - logw.max()) * (1 - np.exp(-4 * np.pi * np.sqrt(2 * C * E)))
    m1 = np.sum(E * w) / np.sum(w)
    m2 = np.sum(E * E * w) / np.sum(w)
    return dict(E_mean=float(m1), var=float(m2 - m1 * m1))


def weighted_anec_bound(c: float, tau: float) -> float:
    """Proposition 6c-10: the sharp lower bound -c/(48 pi tau) on int T_uu e^{-(u-u0)/tau} du for a 2D CFT state equal to the vacuum before u0."""
    return float(-c / (48 * np.pi * tau))


def weighted_anec_of_pulse(c: float, A: float, tau: float, w: float = 0.05, L: float = None, W2: float = None, n: int = 400001) -> float:
    """int t_u e^{-(u - u_p)/tau} du for the conformal-map pulse (rise width w at u_p = 0, fall at L of width W2)."""
    from .throat_evolution import conformal_pulse
    L = 40 * tau if L is None else L
    W2 = 5 * tau if W2 is None else W2
    tu, _ = conformal_pulse(c, A, 0.0, w, L, W2)
    u = np.linspace(-20 * w, L + 10 * W2, n)
    return float(np.trapezoid(tu(u) * np.exp(-u / tau), u))


def quantum_tipping_c(r_N: float, W_over_rN: float = 1.85e-3, tau_over_rN: float = 3.2) -> dict:
    """Corollary 6c-11: smallest central charge whose best weighted pulse c/(48 pi tau) reaches W_th (G = l_P = 1)."""
    tau = tau_over_rN * r_N
    c = 48 * np.pi * tau * W_over_rN * r_N
    return dict(c_min=float(c), c_over_rN2=float(c / r_N**2), species_length_over_rN=float(np.sqrt(c) / r_N))


def nariai_fluctuation(r_N: float, f_th: float = 0.006, c: float = 0.0) -> dict:
    """Assessment: radion fluctuation l_P/r_N against the tipping fraction, and the one-loop parameter G c / (12 pi r_N^2)."""
    fluct = 1 / r_N
    return dict(fluctuation=float(fluct), classical_dominates=bool(f_th > 10 * fluct), one_loop=float(c / (12 * np.pi * r_N**2)))


def monopole_charge(Q_over_lP: float, alpha: float = 1 / 137.036) -> float:
    """Number of Dirac flux quanta q = 2 sqrt(alpha) Q / l_P carried by a geometric magnetic charge Q."""
    return float(2 * np.sqrt(alpha) * Q_over_lP)
