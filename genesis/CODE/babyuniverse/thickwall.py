"""Scalar-field (thick) walls coupled to gravity (CLASSICAL_MODEL.md Sec. 7).

Matter: S_m = int sqrt(-g) [ -(1/2)(grad phi)^2 - U(phi) ], with the quartic
potential

    U(phi) = (lam/4)(phi^2 - v^2)^2 + rho_f * (1 - s(phi/v)),
    s(x)   = 1/2 + (3/4) x - (1/4) x^3,

which keeps the extrema at phi = -v (false vacuum, U = rho_f) and phi = +v
(true vacuum, U = 0) for every rho_f (s'(+-1) = 0).  The flat-space kink of
the symmetric part has width w = sqrt(2)/(sqrt(lam) v) and tension
sigma_0 = (2 sqrt 2/3) sqrt(lam) v^3 = (4/3) v^2 / w.

(a) Time-symmetric initial data (K_ij = 0, pi_phi = 0): the only constraint
is the Hamiltonian constraint 3R = 16 pi G rho, which for
dl^2 + r(l)^2 dOmega^2 reads

    2 r r'' = 1 - r'^2 - 8 pi G r^2 rho(l),   rho = phi'^2/2 + U(phi).

Any profile phi(l) gives admissible data.  Outside the wall r'^2 = 1 - 2GM/r,
so r' -> eps_out sqrt(f_out) and the thin-wall junction at Rdot = 0,
r'_+ - r'_- = -4 pi G sigma R, is recovered as w -> 0 (tested).

(b) O(4) Coleman-De Luccia bounce with gravity, ds^2 = dxi^2 + a(xi)^2 dOmega_3^2:
    phi'' + 3 (a'/a) phi' = U'(phi),   a'' = -(8 pi G/3) a (phi'^2 + U),
    a'^2 = 1 + (8 pi G/3) a^2 (phi'^2/2 - U)   (constraint, monitored),
on-shell action S_E = 4 pi^2 int dxi [a^3 U - 3 a/(8 pi G)].
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad, solve_ivp


@dataclass(frozen=True)
class QuarticPotential:
    lam: float
    v: float
    rho_f: float

    def s(self, x):
        return 0.5 + 0.75 * x - 0.25 * x**3

    def ds(self, x):
        return 0.75 - 0.75 * x**2

    def U(self, phi):
        x = phi / self.v
        return self.lam / 4 * (phi**2 - self.v**2) ** 2 + self.rho_f * (1 - self.s(x))

    def dU(self, phi):
        x = phi / self.v
        return self.lam * phi * (phi**2 - self.v**2) - self.rho_f * self.ds(x) / self.v

    @property
    def width(self):
        return np.sqrt(2) / (np.sqrt(self.lam) * self.v)

    @property
    def sigma0(self):
        return 2 * np.sqrt(2) / 3 * np.sqrt(self.lam) * self.v**3

    def sigma_thin(self):
        """Thin-wall tension int_{-v}^{v} sqrt(2 (U - U_true)) dphi (true side)."""
        return quad(lambda p: np.sqrt(max(2 * (self.U(p) - 0.0), 0.0)), -self.v, self.v, limit=200)[0]

    @classmethod
    def from_width_tension(cls, w, sigma, rho_f):
        v2 = 3 * sigma * w / 4
        lam = 2 / (w**2 * v2)
        return cls(lam=lam, v=np.sqrt(v2), rho_f=rho_f)


# ------------------------------------------------------- (a) initial data
def kink_profile(pot: QuarticPotential, l_w: float):
    w = pot.width

    def phi(l):
        return pot.v * np.tanh((l - l_w) / w)

    def dphi(l):
        x = np.abs((l - l_w) / w)
        e = np.exp(-2 * x)
        return pot.v / w * 4 * e / (1 + e) ** 2      # = v/w sech^2, overflow-safe

    return phi, dphi


def time_symmetric_data(pot: QuarticPotential, l_w: float, G: float = 1.0,
                        l_max: float = None, n: int = 4000):
    """Integrate the Hamiltonian constraint outward from a regular centre."""
    phi, dphi = kink_profile(pot, l_w)
    w = pot.width
    l_max = l_max or (l_w + 60 * w + 6 * l_w)

    def rho(l):
        return 0.5 * dphi(l) ** 2 + pot.U(phi(l))

    def rhs(l, y):
        r, rp = y
        return [rp, (1 - rp**2 - 8 * np.pi * G * r**2 * rho(l)) / (2 * r)]

    l0 = 1e-6 * max(l_w, 1.0)
    rho_c = rho(0.0)
    r0 = l0 - (4 * np.pi * G / 9) * rho_c * l0**3
    rp0 = 1 - (4 * np.pi * G / 3) * rho_c * l0**2
    ls = np.unique(np.concatenate([np.linspace(l0, l_max, n),
                                   np.linspace(max(l0, l_w - 16 * w), min(l_max, l_w + 16 * w), 4001)]))
    sol = solve_ivp(rhs, (l0, l_max), [r0, rp0], t_eval=ls, method="DOP853",
                    rtol=1e-11, atol=1e-13, max_step=w / 8)
    l, (r, rp) = sol.t, sol.y
    rh = np.array([rho(x) for x in l])
    mass = r / (2 * G) * (1 - rp**2)          # Misner-Sharp (time-symmetric slice)
    return dict(l=l, r=r, rp=rp, rho=rh, mass=mass, phi=phi(l))


def thin_wall_comparison(pot: QuarticPotential, l_w: float, G: float = 1.0, k_out: float = 8.0):
    """Compare thick time-symmetric data with the thin-wall junction.

    The one-sided slopes are extrapolated to the wall centre with the exact
    vacuum solutions, A_- = eps_in sqrt(1 - chi^2 R^2) and
    A_+ = eps_out sqrt(1 - 2 G M_far / R), where M_far is the asymptotic
    Misner-Sharp mass and eps = sign r' on either side of the energy lump;
    the thin-wall prediction is A_+ - A_- = -4 pi G sigma_eff R with
    sigma_eff the background-subtracted proper surface energy density.
    All discrepancies are O(w/R)."""
    d = time_symmetric_data(pot, l_w, G)
    w = pot.width
    l, r, rp, rho = d["l"], d["r"], d["rp"], d["rho"]
    i_w = np.searchsorted(l, l_w)
    i_m = np.searchsorted(l, l_w - k_out * w)
    i_p = np.searchsorted(l, l_w + k_out * w)
    R = r[i_w]
    bg = np.where(l < l_w, pot.rho_f, 0.0)
    sel = slice(i_m, i_p + 1)
    sigma_eff = np.trapezoid((rho - bg)[sel], l[sel])
    chi = np.sqrt(8 * np.pi * G * pot.rho_f / 3)
    M_far = d["mass"][-1]
    eps_in = int(np.sign(rp[i_m]))
    eps_out = int(np.sign(rp[i_p]))
    A_minus = eps_in * np.sqrt(max(1 - chi**2 * R**2, 0.0))
    A_plus = eps_out * np.sqrt(max(1 - 2 * G * M_far / R, 0.0))
    kappa = 4 * np.pi * G * sigma_eff
    M_thin = chi**2 * R**3 / (2 * G) + kappa * R**2 / G * A_minus - kappa**2 * R**3 / (2 * G)
    i_min = i_p + int(np.argmin(r[i_p:]))
    return dict(w=w, R=R, A_minus=A_minus, A_plus=A_plus, jump=A_plus - A_minus,
                jump_thin=-kappa * R, sigma_eff=sigma_eff, sigma0=pot.sigma0,
                M_far=M_far, M_thin=M_thin, eps_in=eps_in, eps_out=eps_out,
                r_min_outside=float(r[i_min]), bridge_thin=2 * G * M_far,
                data=d)


# ----------------------------------------------------------- (b) CdL bounce
def cdl_bounce(pot: QuarticPotential, G: float = 1.0, tol: float = 1e-10, max_iter: int = 200,
               eps0: float = 1e-6):
    """O(4) bounce from the true vacuum (centre, phi ~ +v) to the false vacuum
    (far pole), found by bisection on t = ln(1/|phi(0) - v|).

    Thin walls require phi(0) exponentially close to +v (|phi(0)-v| ~ e^{-m R0}
    with m = sqrt(U''(v))), far below double precision.  We therefore start
    the integration at xi_s where the regular linearised solution
        phi - v = -|delta| 2 I_1(m xi)/(m xi),   a = xi   (true vacuum is flat)
    has grown to |phi - v| = eps0 * v, so delta is never formed explicitly:
        ln|delta| + ln(2 I_1(m xi_s)/(m xi_s)) = ln(eps0 v).
    Overshoot  = field passes the false vacuum before a returns to zero;
    undershoot = field turns back (phi' changes sign).
    Returns the profile and B = S_E[bounce] - S_E[false vacuum]."""
    from scipy.special import ive
    from scipy.optimize import brentq as _brentq
    v = pot.v
    m = np.sqrt(abs(2 * pot.lam * v**2 - pot.rho_f * (-1.5 * 1.0) / v**2 * 0 + 0.0))
    # exact curvature at +v: U'' = 2 lam v^2 + lam (3v^2 - v^2) ... use numerical second derivative
    h = 1e-6 * v
    m = np.sqrt((pot.U(v + h) - 2 * pot.U(v) + pot.U(v - h)) / h**2)
    Hf = np.sqrt(8 * np.pi * G * pot.rho_f / 3)

    def log_growth(x):  # ln(2 I_1(x)/x), stable for large x (ive = I e^{-x})
        return np.log(2 * ive(1, x) / x) + x

    def start(t):
        target = np.log(eps0 * v) + t          # ln(2I1/x) = ln(eps0 v) - ln|delta|
        if log_growth(1e-8) >= target:
            xs = 1e-8
        else:
            xs = _brentq(lambda x: log_growth(x) - target, 1e-8, 1e6)
        xi = xs / m
        dev = eps0 * v
        # derivative of 2 I_1(x)/x is 2 I_2(x)/x  ->  dphi/dxi = -|delta| m 2 I_2/x
        dphi = -dev * m * ive(2, xs) / ive(1, xs)
        return xi, [v - dev, dphi, xi, 1.0]

    def rhs(x, y):
        p, dp, a, da = y
        if a <= 1e-300:
            return [dp, 0.0, da, 0.0]
        return [dp, pot.dU(p) - 3 * da / a * dp, da, -(8 * np.pi * G / 3) * a * (dp**2 + pot.U(p))]

    def ev_turn(x, y):
        return y[1]
    ev_turn.terminal = True
    ev_turn.direction = 1

    def ev_a0(x, y):
        return y[2]
    ev_a0.terminal = True
    ev_a0.direction = -1

    def ev_over(x, y):
        return y[0] + 1.5 * v
    ev_over.terminal = True

    def shoot(t, dense=False):
        x0, y0 = start(t)
        sol = solve_ivp(rhs, (x0, x0 + 1.5 * np.pi / Hf), y0, method="DOP853", rtol=1e-11,
                        atol=1e-14 * max(1.0, v), events=[ev_turn, ev_a0, ev_over],
                        dense_output=dense, max_step=pot.width / 6)
        if len(sol.t_events[0]):
            return -1, sol
        if len(sol.t_events[2]):
            return +1, sol
        if len(sol.t_events[1]):
            return (+1 if sol.y[1, -1] < 0 else -1), sol
        return -1, sol

    t_lo, t_hi = 0.0, 5.0 * m / Hf + 50.0
    s_lo, _ = shoot(t_lo)
    s_hi, _ = shoot(t_hi)
    if s_lo == s_hi:
        raise RuntimeError(f"could not bracket the bounce (signs {s_lo}, {s_hi})")
    for _ in range(max_iter):
        t_mid = 0.5 * (t_lo + t_hi)
        s_mid, _ = shoot(t_mid)
        if s_mid == s_lo:
            t_lo = t_mid
        else:
            t_hi = t_mid
        if t_hi - t_lo < tol * max(1.0, t_hi):
            break
    s_f, sol = shoot(0.5 * (t_lo + t_hi))
    x0, _ = start(0.5 * (t_lo + t_hi))
    x = sol.t
    p, dp, a, da = sol.y
    cons = da**2 - 1 - (8 * np.pi * G / 3) * a**2 * (0.5 * dp**2 - pot.U(p))
    # The far side sits in the false vacuum; near the far pole the Hubble term
    # turns into anti-friction and any residual error grows, so we cut the
    # numerical solution where it is closest to (phi, phi') = (-v, 0) after the
    # wall and complete it with the exact false-vacuum cap
    #   a = sin(H s)/H,  s = distance to the pole,  s_c = atan2(H a, -a')/H.
    after = np.where(p < 0)[0]
    i0w = after[0] if len(after) else len(x) - 1
    dist = np.abs(p + v) / v + np.abs(dp) * pot.width / v
    i_c = i0w + int(np.argmin(dist[i0w:]))
    xs_, ps_, as_ = x[: i_c + 1], p[: i_c + 1], a[: i_c + 1]
    integrand = as_**3 * pot.U(ps_) - 3 * as_ / (8 * np.pi * G)
    S = 4 * np.pi**2 * np.trapezoid(integrand, xs_)
    S += 4 * np.pi**2 * (-3 / (8 * np.pi * G)) * x0**2 / 2          # inner flat ball
    Hc = np.sqrt(8 * np.pi * G * pot.U(-v) / 3)
    s_c = np.arctan2(Hc * a[i_c], -da[i_c]) / Hc
    c = np.cos(Hc * s_c)
    S += 4 * np.pi**2 * (pot.U(-v) / Hc**4 * (2 / 3 - c + c**3 / 3)
                         - 3 / (8 * np.pi * G * Hc**2) * (1 - c))
    a_end = a[i_c]
    S_false = -3.0 / (8 * G**2 * pot.rho_f)
    # wall position: where phi crosses the barrier top region (phi = 0)
    i0 = int(np.argmin(np.abs(p)))
    return dict(x=x, phi=p, a=a, i_cut=int(i_c), x_start=float(x0), xi_pole=float(x[i_c] + s_c),
                constraint_max=float(np.max(np.abs(cons[: i_c + 1]))), S_bounce=S,
                S_false=S_false, B=S - S_false, a_wall=float(a[i0]), a_end=float(a_end),
                bisection_width=t_hi - t_lo, outcome=s_f)
