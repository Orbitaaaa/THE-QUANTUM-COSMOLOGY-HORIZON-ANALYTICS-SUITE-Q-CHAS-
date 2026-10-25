"""Round 7 (problem 3): creating the false-vacuum ball by supercooling, starting from the true vacuum.

The route of Round 6d starts from a false-vacuum ball that is sub-horizon (wall at r ~ 0.6 l_H), momentarily static at its corner, and empty
of mass inside (pure de Sitter).  Here that state is produced from ordinary matter: a hot, expanding ball of radiation in which the scalar
starts in the TRUE vacuum phi = +v.  The particles that make up the radiation get their mass from phi, m_i^2 = g_i^2 (phi + v)^2, so they are
massless at phi = -v; the leading thermal correction to the free energy (high-T expansion, (T^2/24) sum m_i^2) then favours phi = -v.

Homogeneous interior (closed FRW, k = +1; the ball is the comoving cap chi < chi_b, matched later to the exterior by the full route code):
    F(phi, T) = -(pi^2/30) g T^4 + (c_T/2) T^2 (phi + v)^2 + U(phi),         U = Round 5 tilted quartic (false vacuum -v, true +v),
    rho = phi_dot^2/2 + F - T F_T,   p = phi_dot^2/2 - F,   s = -F_T,
    phi_dd + (3H + Gamma) phi_dot + F_phi = 0,    d(s a^3)/dt = a^3 Gamma phi_dot^2 / T   (dissipation of phi oscillations into the bath),
    a_dd / a = -(4 pi G/3)(rho + 3p),            constraint H^2 = (8 pi G/3) rho - 1/a^2 (monitored).
Energy conservation rho_dot = -3H (rho + p) follows identically.  A sphere at chi is anti-trapped iff a H tan(chi) > 1 (trapped iff < -1),
so the ball is normal throughout iff |a_dot| tan(chi_b) < 1.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp


class Fireball:
    def __init__(self, l_H=30.0, eps=2.0, vphi=0.0015, g_star=10.0, alpha_pt=1.0, Gamma_over_m=0.3, G=1.0):
        """alpha_pt = rho_f / rho_rad(T_c) fixes c_T (the strength of the phase transition); Gamma in units of the false-vacuum mass."""
        self.G, self.v, self.g = G, vphi, g_star
        self.rho_f = 3 / (8 * np.pi * G * l_H**2)
        self.lam = 4 * self.rho_f / (eps * vphi**4)
        self.l_H = l_H
        # T_c: F(+v) = F(-v)  ->  2 c_T T_c^2 v^2 = rho_f ;  rho_rad(T_c) = (pi^2/30) g T_c^4 = rho_f / alpha_pt
        self.Tc = (self.rho_f / alpha_pt / (np.pi**2 / 30 * g_star)) ** 0.25
        self.cT = self.rho_f / (2 * self.Tc**2 * vphi**2)
        self.m_f = np.sqrt(max(2 * self.lam * vphi**2 - 1.5 * self.rho_f / vphi**2, 0.0))
        self.Gamma = Gamma_over_m * self.m_f

    def U(self, p):
        x = p / self.v
        return self.lam / 4 * (p * p - self.v**2) ** 2 + self.rho_f * (0.5 - 0.75 * x + 0.25 * x**3)

    def dU(self, p):
        x = p / self.v
        return self.lam * p * (p * p - self.v**2) - self.rho_f * (0.75 - 0.75 * x * x) / self.v

    def entropy(self, T, p):
        return 2 * np.pi**2 / 15 * self.g * T**3 - self.cT * T * (p + self.v) ** 2

    def temperature(self, s, p):
        """Invert s(T, phi) for the branch continuously connected to pure radiation (largest root)."""
        b = 2 * np.pi**2 / 15 * self.g
        c = self.cT * (p + self.v) ** 2
        T = max((max(s, 0.0) / b) ** (1 / 3), np.sqrt(c / b) * 1.0001)
        for _ in range(60):
            f = b * T**3 - c * T - s
            T_new = T - f / (3 * b * T * T - c)
            if abs(T_new - T) < 1e-14 * max(T, 1e-30):
                T = T_new
                break
            T = max(T_new, np.sqrt(c / (3 * b)) * 1.0001)
        return T

    def thermo(self, p, pd, T):
        F = -np.pi**2 / 30 * self.g * T**4 + 0.5 * self.cT * T * T * (p + self.v) ** 2 + self.U(p)
        FT = -4 * np.pi**2 / 30 * self.g * T**3 + self.cT * T * (p + self.v) ** 2
        rho = 0.5 * pd * pd + F - T * FT
        pr = 0.5 * pd * pd - F
        Fp = self.cT * T * T * (p + self.v) + self.dU(p)
        return rho, pr, Fp

    def rhs(self, t, y):
        a, ad, p, pd, S = y
        T = self.temperature(S / a**3, p)
        rho, pr, Fp = self.thermo(p, pd, T)
        H = ad / a
        add = -4 * np.pi * self.G / 3 * (rho + 3 * pr) * a
        pdd = -(3 * H + self.Gamma) * pd - Fp
        Sd = a**3 * self.Gamma * pd * pd / T
        return [ad, add, pd, pdd, Sd]

    def run(self, T_i, a_i, t_max=400.0, phi_i=None, expanding=True, max_step=None):
        """Start hot (T_i) with phi in the TRUE vacuum (phi_i = +v by default) at scale factor a_i; a_dot from the Friedmann constraint."""
        p = self.v if phi_i is None else phi_i
        S = self.entropy(T_i, p) * a_i**3
        rho = self.thermo(p, 0.0, T_i)[0]
        H2 = 8 * np.pi * self.G / 3 * rho - 1 / a_i**2
        if H2 <= 0:
            raise ValueError("no expanding solution at these initial conditions")
        ad = a_i * np.sqrt(H2) * (1 if expanding else -1)
        turn = lambda t, y: y[1]
        turn.terminal, turn.direction = True, -1
        ms = max_step or 0.05 / max(self.m_f, 1e-9)
        sol = solve_ivp(self.rhs, (0, t_max), [a_i, ad, p, 0.0, S], method="DOP853", rtol=1e-9, atol=1e-14, events=turn, max_step=ms,
                        dense_output=False)
        return sol

    def summary(self, sol, sin_chi_b=None):
        a, ad, p, pd, S = sol.y
        T = np.array([self.temperature(S[k] / a[k] ** 3, p[k]) for k in range(len(a))])
        rho_rad = np.pi**2 / 10 * self.g * T**4 - 0.5 * self.cT * T * T * (p + self.v) ** 2
        rho = np.array([self.thermo(p[k], pd[k], T[k])[0] for k in range(len(a))])
        cons = np.abs(ad**2 / a**2 - (8 * np.pi * self.G / 3 * rho - 1 / a**2)) / (8 * np.pi * self.G / 3 * rho)
        k = -1
        out = dict(turnaround=bool(sol.status == 1), t_turn=float(sol.t[k]), a_max=float(a[k]), a_max_over_lH=float(a[k] / self.l_H),
                   x_end=float(p[k] / self.v), kinetic_over_rho_f=float(0.5 * pd[k] ** 2 / self.rho_f), rad_over_rho_f=float(rho_rad[k] / self.rho_f),
                   T_end_over_Tc=float(T[k] / self.Tc), max_adot=float(np.max(np.abs(ad))), constraint_max=float(np.max(cons)),
                   x_min_during=float(p.min() / self.v), frac_time_true_basin=float(np.mean(p / self.v > -0.375)))
        if sin_chi_b is not None:
            tanb = sin_chi_b / np.sqrt(1 - sin_chi_b**2)
            out.update(sin_chi_b=sin_chi_b, R_ball_turn_over_lH=float(a[k] * sin_chi_b / self.l_H), normal_throughout=bool(np.max(np.abs(ad)) * tanb < 1),
                       max_adot_tanchib=float(np.max(np.abs(ad)) * tanb))
        return out


def static_hot_ball(fb: Fireball, R_over_lH=0.6, T_over_Tc=1.5, t_restore=6.0, t_cool_over_R=1.0, n_cool=8.0):
    """Stages A and B of the creation (problem 3), in the ball's rest frame.
    A (restoration): a momentarily static hot ball of radius R, phi = +v (true vacuum), T = T_over_Tc * T_c.  Normal iff 2m/r < 1, i.e.
      (8 pi G/3) rho R^2 < 1.  The field obeys phi_dd + Gamma phi_dot + F_phi(phi, T) = 0 for t_restore (<< R, so the geometry is frozen).
    B (adiabatic emptying): the radiation free-streams out on the light-crossing time, T(t) = T_A exp(-t / (t_cool_over_R R)), for n_cool
      e-folds; the field equation is integrated through it.  Returns where the field ends, its kinetic energy, and the normality margin."""
    R = R_over_lH * fb.l_H
    T_A = T_over_Tc * fb.Tc
    rho = fb.thermo(fb.v, 0.0, T_A)[0]
    compact = 8 * np.pi * fb.G / 3 * rho * R**2
    f = lambda t, y, T: [y[1], -fb.Gamma * y[1] - (fb.cT * T * T * (y[0] + fb.v) + fb.dU(y[0]))]
    ms = 0.05 / max(fb.m_f, 1e-9)
    a = solve_ivp(lambda t, y: f(t, y, T_A), (0, t_restore), [fb.v, 0.0], method="DOP853", rtol=1e-10, atol=1e-14, max_step=ms)
    tc = t_cool_over_R * R
    b = solve_ivp(lambda t, y: f(t, y, T_A * np.exp(-t / tc)), (0, n_cool * tc), a.y[:, -1], method="DOP853", rtol=1e-10, atol=1e-14,
                  max_step=ms)
    x = b.y[0] / fb.v
    return dict(R_over_lH=R_over_lH, T_over_Tc=T_over_Tc, compactness_2m_over_r=float(compact), normal=bool(compact < 1),
                x_after_restoration=float(a.y[0, -1] / fb.v), x_end=float(x[-1]), max_x_during_cooling=float(x.max()),
                kinetic_end_over_rho_f=float(0.5 * b.y[1, -1] ** 2 / fb.rho_f), T_end_over_Tc=float(np.exp(-n_cool) * T_over_Tc),
                adiabaticity=float(1 / (tc * np.sqrt(fb.m_f**2 + fb.cT * T_A**2))), rho_rad_over_rho_f=float((rho - fb.rho_f) / fb.rho_f),
                E_rad=float(4 * np.pi / 3 * R**3 * (rho - fb.U(-fb.v))))
