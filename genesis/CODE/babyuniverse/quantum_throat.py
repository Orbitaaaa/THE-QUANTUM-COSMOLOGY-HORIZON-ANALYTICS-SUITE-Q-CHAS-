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
