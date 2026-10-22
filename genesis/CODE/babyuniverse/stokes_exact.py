"""Roadmap Stage 1, P-Q2: exact (numerical) Stokes multiplier of the squared constraint at M = 0, kappa > chi.

Momentum (Laplace) representation psi(u) = int exp(-u theta / h) phi(theta) d theta, h = hbar G, P_u = i theta / G,
u -> h d/dtheta.  The squared constraint g(u) cosh^2(G P_u) - N(u)^2 = 0 with g = 4(1 - 2 chi^2 u), N = 2(1 - a u),
a = chi^2 + kappa^2 (M = 0, eta = 0) becomes, for the ordering family "g to the left with weight lam" (lam = 1: g left,
lam = 0: g right, lam = 1/2: Weyl-symmetric), the exact second-order ODE

    a^2 h^2 phi'' - 2 h (a - chi^2 cos^2 theta) phi' + (sin^2 theta - 2 lam chi^2 h sin 2theta + beta h) phi = 0 .

beta = 0 is the squared constraint; beta != 0 adds a constant subprincipal term and is used only as a control
(Weber index nu = -1/2 + beta / (2 r_q), r_q = sqrt(kappa^4 - chi^4)).

The two characteristic roots u_{1,2}(theta) = [(a - chi^2 c^2) +- c r(theta)] / a^2 cross at theta_q = pi/2 (C = 0: the
coalescence point of PHASE_III Sec. 7).  WKB solutions phi_i = A_i exp(S_i/h)(1 + h B_i), S_i(theta_q) = 0,
A_i = |theta - theta_q|^{alpha_i} * Atilde_i with Atilde_i(theta_q) = 1 (alpha_1 = nu, alpha_2 = -1 - nu).  The exact
solution that is phi_1 (recessive) on the left contains s * phi_2 on the right; Weber theory predicts

    s = sqrt(2 pi) / Gamma(-nu) * (2 r_q / (a^2 h))^(-nu - 1/2)   ->   sqrt(2)  for nu = -1/2 .

The Laplace transform maps the Weber index nu to -1 - nu (u-space recessive <-> theta-space dominant), so nu = -1/2 is
self-dual and the theta-space multiplier tests the u-space one (PHASE_III Prop. III-8).
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import gamma as _Gamma


class SquaredConstraintODE:
    def __init__(self, chi: float, kappa: float, lam: float = 0.5, beta: float = 0.0):
        if kappa <= chi:
            raise ValueError("need kappa > chi (single-flip class at M = 0)")
        self.chi, self.kappa, self.lam, self.beta = chi, kappa, lam, beta
        self.a = chi**2 + kappa**2
        self.rq = np.sqrt(kappa**4 - chi**4)
        self.thq = np.pi / 2
        self.nu = -0.5 + beta / (2 * self.rq)

    # characteristic data
    def r(self, t):
        return np.sqrt(self.kappa**4 - self.chi**4 * np.sin(t) ** 2)

    def u(self, i, t):
        c = np.cos(t)
        return ((self.a - self.chi**2 * c**2) + (1 if i == 1 else -1) * c * self.r(t)) / self.a**2

    def du(self, i, t):
        c, s = np.cos(t), np.sin(t)
        dr = -self.chi**4 * s * c / self.r(t)
        return (2 * self.chi**2 * c * s + (1 if i == 1 else -1) * (-s * self.r(t) + c * dr)) / self.a**2

    def delta(self, i, t):
        return (1 if i == 1 else -1) * np.cos(t) * self.r(t)

    def tau(self, i, t):
        return self.a**2 * self.du(i, t) - 2 * self.lam * self.chi**2 * np.sin(2 * t) + self.beta

    def ell(self, i, t):                      # A_i'/A_i from the transport equation 2 delta A' + tau A = 0
        return -self.tau(i, t) / (2 * self.delta(i, t))

    def alpha(self, i):
        return self.nu if i == 1 else -1 - self.nu

    def alpha_from_transport(self, i, eps: float = 1e-6) -> float:
        """Local exponent of A_i at theta_q read off the transport equation: lim x * ell_i."""
        return 0.5 * (eps * self.ell(i, self.thq + eps) - eps * self.ell(i, self.thq - eps))

    # WKB pieces
    def S(self, i, t):
        return quad(lambda x: self.u(i, x), self.thq, t, epsabs=0.0, epsrel=1e-13, limit=400)[0]

    def ln_Atilde(self, i, t):
        al, q = self.alpha(i), self.thq

        def f(x):
            y = x - q
            if abs(y) < 1e-4:
                return 0.5 * ((self.ell(i, q + 1e-4) - al / 1e-4) + (self.ell(i, q - 1e-4) + al / 1e-4))
            return self.ell(i, x) - al / y
        return quad(f, q, t, epsabs=1e-15, epsrel=1e-13, limit=400)[0]

    def _dell(self, i, t, e=1e-5):
        return (self.ell(i, t + e) - self.ell(i, t - e)) / (2 * e)

    def _Bprime(self, i, t):
        return -self.a**2 * (self._dell(i, t) + self.ell(i, t) ** 2) / (2 * self.delta(i, t))

    def B(self, i, t, tref):
        return quad(lambda x: self._Bprime(i, x), tref, t, epsabs=1e-14, epsrel=1e-12, limit=400)[0]

    def wkb(self, i, t, h, order, tref):
        """(w, w') with w = phi_i exp(-S_1/h)."""
        A = np.exp(self.ln_Atilde(i, t)) * abs(t - self.thq) ** self.alpha(i)
        Bv = self.B(i, t, tref) if order > 0 else 0.0
        dB = self._Bprime(i, t) if order > 0 else 0.0
        E = np.exp((self.S(i, t) - self.S(1, t)) / h)
        w = A * (1 + h * Bv) * E
        dlog = self.u(i, t) / h + self.ell(i, t) + h * dB / (1 + h * Bv)
        return w, w * dlog - self.u(1, t) / h * w

    def multiplier(self, h: float, half_window: float = 1.0, order: int = 1, tref: float = 1.0) -> float:
        q = self.thq
        tL, tR = q - half_window, q + half_window
        a2 = self.a**2

        def rhs(t, y):
            return [y[1], -(2 * self.delta(1, t) * y[1] + self.tau(1, t) * y[0]) / (a2 * h)]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            w0, dw0 = self.wkb(1, tL, h, order, q - tref)
            sol = solve_ivp(rhs, (tL, tR), [w0, dw0], method="DOP853", rtol=1e-13, atol=1e-300)
            w, dw = sol.y[0, -1], sol.y[1, -1]
            w1, dw1 = self.wkb(1, tR, h, order, q + tref)
            w2, dw2 = self.wkb(2, tR, h, order, q + tref)
        return float((w1 * dw - dw1 * w) / (w1 * dw2 - dw1 * w2))

    def predicted(self, h: float) -> float:
        nu = self.nu
        if abs(nu - round(nu)) < 1e-14 and round(nu) >= 0:
            return 0.0
        return float(np.sqrt(2 * np.pi) / _Gamma(-nu) * (2 * self.rq / (self.a**2 * h)) ** (-nu - 0.5))


def convergence_table(chi=1.0, kappa=1.5, hs=(0.01, 0.005, 0.0025, 0.00125, 0.000625), lams=(1.0, 0.5, 0.0)) -> dict:
    out = {}
    for lam in lams:
        ode = SquaredConstraintODE(chi, kappa, lam)
        errs = [ode.multiplier(h) / np.sqrt(2) - 1 for h in hs]
        # fit err = h (A ln h + B) on the two smallest h
        y = np.array(errs[-2:]) / np.array(hs[-2:])
        A = (y[1] - y[0]) / (np.log(hs[-1]) - np.log(hs[-2]))
        out[str(lam)] = dict(h=list(hs), rel_err=errs, log_coeff=float(A), const=float(y[1] - A * np.log(hs[-1])))
    return out
