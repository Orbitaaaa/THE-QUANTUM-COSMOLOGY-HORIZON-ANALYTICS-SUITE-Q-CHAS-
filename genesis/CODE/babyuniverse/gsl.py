"""Round 6b: genesis and the generalized second law (GENESIS -- a serious attempt to get around the GSL).

Correction F21 (accepted): violating the NEC does not violate the GSL.  The GSL constrains causal horizons, and quantum matter may violate
the NEC while respecting the QNEC and the GSL (Gao-Jafferis-Wall).  This module computes the quantities behind the Round 6b results
(ROADMAP_ROUND6B_GSL.md); the theorems themselves are proved there.

* ``wall_dichotomy``        Theorem GSL-3 (genesis dichotomy).  With the GSL on causal horizons of complete observers, global
                            hyperbolicity and normal, weakly gravitating laboratory data, a complete daughter observer must keep a
                            non-compact past (the neck never closes) and cannot end in a de Sitter phase.  An eternal daughter and a
                            daughter determined by the laboratory data exclude each other.
* ``qnec_cost``             Theorem GSL-4 (QNEC cost of genesis).  Semiclassical Einstein equations, spherical symmetry and the QNEC
                            give r'' <= -(G hbar / 2 pi r) S_out'' on an ingoing null ray, so flipping theta_- from theta_0 < 0 needs
                                -int (r_0/r) S_out'' dlambda  >=  (A_0 / 4 G hbar) |theta_0|.
* ``cft_flip``              Proposition GSL-5 (two dimensions, exact).  For 2D CFT matter (strong QNEC
                            T >= (hbar/2pi)(S'' + 6 S'^2/c)) the quantum expansion obeys Theta' <= -(6/c) S'^2, and an anti-trapped
                            episode that starts quantum-normal lasts an affine length < c / (6 |Theta_0|).
* ``NearExtremalThroat``    Proposition GSL-6 (throat discount).  At proper depth rho in a near-extremal throat of radius Q and
                            temperature T, sqrt f = 2 pi Q T sinh(rho/Q); a QNEC-saturating pulse of local scale eta/Q anti-traps the
                            spheres iff T sinh(rho/Q) < (c eta / 12) T_brk, T_brk = 1/(2 pi^2 Q^3).
* ``vaidya_flip``           Proposition GSL-7 (exact metric).  Outgoing Bonnor-Vaidya-Reissner-Nordstrom with the CFT flux: a sphere
                            is anti-trapped once the mass has risen by delta m > r f(r)/2 -- the budget of Theorem G1, which vanishes
                            at an extremal horizon.
* ``nariai_parameter``      Proposition GSL-8 (growth).  A T-region of Einstein-Maxwell-Lambda (charge Q, mass parameter C = 2M)
                            expands without bound iff C > max_{r >= r_0} (r + Q^2/r - Lambda r^3/3); an extremal-throat bud
                            (C = 2Q, r_0 = Q) grows iff Lambda Q^2 > 3/16 (0.2212 if the vacuum energy is converted locally,
                            ``switch_on_threshold``); a static throat in a fixed vacuum never grows; none exists for Lambda Q^2 > 1/4.
* ``cauchy_route_window``   Proposition GSL-9 (Cauchy-horizon route near extremality).  The barrier-free mass window of the charged
                            bubble (Proposition G4) survives as Q/M -> 1 but lies inside [1/(sqrt L (1+s)), 1/(sqrt L (1-s))],
                            L = chi^2 + kappa^2, s = sqrt(1 - Q^2/M^2) (the outer flip must lie between the horizons), with the upper
                            edge attained: its width is <= r_+ - r_- = 2 kappa_- r_-^2, so weakening the inner-horizon blueshift
                            closes the window.
* ``hiding_depth``          Proposition GSL-10 (deep throats hide): energy dM kept inside an extremal throat puts every region deeper
                            than rho ~ Q (1 - ln(2 dM/Q)/2) behind the event horizon (exactly 1.98 Q for the bud energy Q/32).
* ``electroweak_corona_bound``  Remark GSL-11.  Self-consistent bound on Lambda Q^2 in the electroweak corona of a magnetic black hole.

Units: G = c = hbar = k_B = 1 unless a function says SI or GeV.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from .mechanisms import C_SI, G_SI, GEV4_J_M3, HBAR_SI, L_P_SI, M_PL_GEV, ChargedBubble

K_B_SI = 1.380649e-23
T_P_K = 1.416784e32                          # Planck temperature
M_P_KG = 2.176434e-8
ALPHA_EM = 1 / 137.035999


# ------------------------------------------------------------------ Theorem GSL-3: the genesis dichotomy
def wall_dichotomy(observer_complete: bool, globally_hyperbolic: bool = True, final_phase: str = "de_sitter",
                   Lambda_final: float = 0.0, Lambda_parent: float = 0.0, neck_persists: bool = True,
                   past_planckian: bool = False) -> dict:
    """Verdict of Theorem GSL-3 for a daughter formed from normal, weakly gravitating laboratory data Sigma_0.

    final_phase: 'de_sitter' (Lambda_final > 0 at late times for the observer), 'flat', 'anti_de_sitter' or 'crunch'.
    Returns the verdict and the premise that decides it."""
    if past_planckian:
        return dict(allowed=None, verdict="outside semiclassical control: the observer's past on Sigma_0 is sub-ell_* (escape E3)",
                    premise="semiclassical GSL")
    if not globally_hyperbolic:
        return dict(allowed=True, verdict="not excluded by the GSL; the daughter lies beyond the Cauchy horizon of Sigma_0 and is not "
                    "determined by the laboratory data (escape E2)", premise="global hyperbolicity")
    if not observer_complete or final_phase == "crunch":
        return dict(allowed=True, verdict="the GSL (causal horizons of complete observers) does not constrain this daughter; its formation "
                    "still needs NEC violation along ingoing rays (Farhi-Guth), priced by Theorem GSL-4", premise="observer completeness")
    if not neck_persists:
        return dict(allowed=False, verdict="excluded: a complete observer with a compact past on Sigma_0 violates the GSL there "
                    "(Proposition GSL-2; Wall's disconnected baby)", premise="GSL on Sigma_0")
    if final_phase == "de_sitter" and Lambda_final > Lambda_parent:
        why = ("excluded: dependent (non-compact past) observers cannot end in a compact-horizon phase (Lemma GSL-1)"
               if Lambda_parent == 0 else "excluded: the GSL on the observer's horizon needs A_daughter >= A_parent (Corollary GSL-3b)")
        return dict(allowed=False, verdict=why, premise="global hyperbolicity + GSL")
    return dict(allowed=True, verdict="allowed by the GSL: transient, dependent daughter (neck never closes, Lambda_final <= Lambda_parent)",
                premise="none violated")


# ------------------------------------------------------------------ Theorem GSL-4: QNEC cost
def qnec_cost(r0: float, theta0: float, G: float = 1.0, hbar: float = 1.0) -> float:
    """(A_0 / 4 G hbar) |theta_0| = pi r_0^2 |theta_0| / (G hbar): the minimum of -int (r_0/r) S_out'' dlambda for a flip."""
    return float(np.pi * r0**2 * abs(theta0) / (G * hbar))


def static_frame_theta(f: float, r: float) -> float:
    """|theta_-| = sqrt(2 f)/r for the ingoing null normal k = (u - n)/sqrt2 of a static observer u (k.u = -1/sqrt2)."""
    return float(np.sqrt(2 * f) / r)


def qnec_raychaudhuri(lam, S, T_extra, r0: float, rp0: float, G: float = 1.0, hbar: float = 1.0) -> dict:
    """Integrate r'' = -4 pi G r T_kk along an affinely parametrised ingoing ray, with
        T_kk = hbar S''(lam) / (2 pi A) + T_extra(lam),   A = 4 pi r^2,
    i.e. a state that satisfies the (integrated) QNEC with slack T_extra >= 0.  Returns r'(end) - r'(0) and the bound
    -(G hbar / 2 pi) int S''/r dlambda of Theorem GSL-4 (equality for T_extra = 0)."""
    lam = np.asarray(lam, float)
    Sd2 = np.gradient(np.gradient(np.asarray(S, float), lam), lam)
    Tx = np.asarray(T_extra(lam) if callable(T_extra) else T_extra, float) * np.ones_like(lam)

    def rhs(x, y):
        s2 = np.interp(x, lam, Sd2)
        tx = np.interp(x, lam, Tx)
        T = hbar * s2 / (2 * np.pi * 4 * np.pi * y[0] ** 2) + tx
        return [y[1], -4 * np.pi * G * y[0] * T, -G * hbar * s2 / (2 * np.pi * y[0])]

    sol = solve_ivp(rhs, (lam[0], lam[-1]), [r0, rp0, 0.0], rtol=1e-10, atol=1e-13, max_step=(lam[-1] - lam[0]) / 400)
    return dict(delta_rp=float(sol.y[1, -1] - rp0), bound=float(sol.y[2, -1]), r_end=float(sol.y[0, -1]))


# ------------------------------------------------------------------ Proposition GSL-8: product solutions and growth
def product_solutions(Q: float, Lam: float) -> dict:
    """Einstein-Maxwell-Lambda product solutions M_2 x S^2 (sphere radius a, magnetic charge Q, G = 1):
        1/a^2 = Lambda + Q^2/a^4,   R(M_2) = 2 (Lambda - Q^2/a^4)   (dS_2 if > 0, AdS_2 if < 0, flat if = 0).
    Real iff Lambda Q^2 <= 1/4; at 1/4 the two branches meet at a^2 = 2 Q^2 (flat x S^2, Plebanski-Hacyan)."""
    if Lam == 0:
        return dict(exists=True, branches=[dict(a=float(Q), kind="AdS2", ell=float(Q))], x=0.0)
    x = Lam * Q**2
    disc = 1 - 4 * x
    if disc < 0:
        return dict(exists=False, branches=[], x=float(x))
    out = []
    for s in (-1, 1):
        a2 = (1 + s * np.sqrt(disc)) / (2 * Lam)
        curv = Lam - Q**2 / a2**2
        kind = "flat" if abs(curv) < 1e-14 * Lam else ("dS2" if curv > 0 else "AdS2")
        out.append(dict(a=float(np.sqrt(a2)), kind=kind, ell=float(1 / np.sqrt(abs(curv))) if kind != "flat" else np.inf))
    return dict(exists=True, branches=out, x=float(x))


def _g(r, Q, Lam):
    return r + Q**2 / r - Lam * r**3 / 3


def nariai_parameter(Q: float, Lam: float, r0: float) -> float:
    """C_N = max_{r >= r0} g(r), g = r + Q^2/r - Lambda r^3/3.  An expanding T-region (r_dot^2 = -f_C, f_C = 1 - C/r + Q^2/r^2 -
    Lambda r^2/3) starting at r0 grows without bound iff C > C_N; for r0 below the AdS_2 radius C_N is the charged-Nariai value."""
    cands = [r0]
    if Lam > 0:
        disc = 1 - 4 * Lam * Q**2
        if disc >= 0:
            for s in (-1, 1):
                rr = np.sqrt((1 + s * np.sqrt(disc)) / (2 * Lam))
                if rr >= r0:
                    cands.append(rr)
    else:
        return np.inf
    return float(max(_g(c, Q, Lam) for c in cands))


def bud_threshold(C_over_Q: float = 2.0, r0_over_Q: float = 1.0) -> float:
    """x = Lambda Q^2 above which a T-region with C = C_over_Q * Q starting at r0 = r0_over_Q * Q grows without bound (3/16 for an
    extremal-throat bud C = 2Q, r0 = Q)."""
    return float(brentq(lambda x: C_over_Q - nariai_parameter(1.0, x, r0_over_Q), 1e-6, 0.25 + 1e-9, xtol=1e-15))


def switch_on_threshold(energy: str = "injected") -> float:
    """Lambda Q^2 above which vacuum energy switched on in an extremal Lambda = 0 throat (r = Q, mass parameter 2Q) makes the region
    grow without bound.  'injected': the vacuum energy (4 pi/3) rho Q^3 is brought in from outside, C stays 2Q, threshold 3/16.
    'conserved': local conversion at fixed Misner-Sharp mass, C = 2Q - Lambda Q^3/3 = g(Q), threshold g(Q) = g(r_N), x = 0.22121.
    A static throat in a fixed vacuum never grows (C = g(a) < C_N for Lambda Q^2 < 1/4)."""
    if energy == "injected":
        return bud_threshold(2.0, 1.0)
    rN = lambda x: np.sqrt((1 + np.sqrt(1 - 4 * x)) / (2 * x))
    return float(brentq(lambda x: _g(1.0, 1.0, x) - _g(rN(x), 1.0, x), 0.2, 0.25 - 1e-12, xtol=1e-15))


def ks_growth(C: float, Q: float, Lam: float, r0: float, r_max: float, n: int = 4000) -> dict:
    """Proper-time history of the expanding T-region r_dot = sqrt(-f_C(r)) from r0: unbounded iff f_C < 0 on [r0, r_max]; returns
    the e-folds of r and whether the late expansion accelerates (r_ddot > 0)."""
    rs = np.geomspace(r0, r_max, n)
    fC = 1 - C / rs + Q**2 / rs**2 - Lam * rs**2 / 3
    grows = bool(np.all(fC < 0))
    rdd = -0.5 * np.gradient(fC, rs)                       # r_ddot = -f_C'(r)/2
    return dict(grows=grows, efolds=float(np.log(r_max / r0)) if grows else None, accelerating_late=bool(rdd[-1] > 0),
                stall_radius=None if grows else float(rs[np.argmax(fC >= 0)]))


# ------------------------------------------------------------------ Proposition GSL-6: near-extremal throat
class NearExtremalThroat:
    """Reissner-Nordstrom with charge Q and Hawking temperature T (Planck units), f = (r - r_+)(r - r_-)/r^2, r_+ r_- = Q^2,
    T = (r_+ - r_-)/(4 pi r_+^2).  Near extremality M = Q + 2 pi^2 Q^3 T^2 and the semiclassical description of the throat breaks down
    at T_brk = 1/(2 pi^2 Q^3) (Schwarzian scale)."""

    def __init__(self, Q: float, T: float):
        self.Q, self.T = Q, T
        T_max = 1 / (6 * np.sqrt(3) * np.pi * Q)            # hottest RN of charge Q (r_+ = sqrt3 Q); near-extremal branch below it
        if not 0 <= T < T_max:
            raise ValueError("need 0 <= T < 1/(6 sqrt3 pi Q)")
        F = lambda s: (s - Q**2 / s) / (4 * np.pi * s**2) - T
        self.rp = brentq(F, Q, np.sqrt(3) * Q, xtol=1e-16 * Q) if T > 0 else Q
        self.rm = Q**2 / self.rp
        self.M = 0.5 * (self.rp + self.rm)
        self.T_brk = 1 / (2 * np.pi**2 * Q**3)

    def f(self, r):
        return (r - self.rp) * (r - self.rm) / r**2

    def depth(self, r: float) -> float:
        """Proper radial distance from the horizon r_+ to r (substitution r = r_+ + s^2 removes the endpoint singularity)."""
        rp, rm = self.rp, self.rm
        s1 = np.sqrt(r - rp)
        val, _ = quad(lambda s: 2 * (rp + s**2) / np.sqrt(rp + s**2 - rm), 0, s1, limit=200, epsabs=0, epsrel=1e-12)
        return float(val)

    def r_at_depth(self, rho: float) -> float:
        hi = self.rp + 1e-30
        while self.depth(hi) < rho:
            hi = self.rp + 2 * (hi - self.rp) + 1e-12 * self.Q
        return float(brentq(lambda r: self.depth(r) - rho, self.rp, hi, xtol=1e-15 * self.Q))

    def sqrt_f_approx(self, rho: float) -> float:
        """AdS_2 approximation sqrt f = (r_+ - r_-)/(2 Q) sinh(rho/Q) = 2 pi Q T sinh(rho/Q) (leading order in T Q)."""
        return float((self.rp - self.rm) / (2 * self.Q) * np.sinh(rho / self.Q))

    def flip_condition(self, rho: float, c: float, eta: float = 1.0) -> dict:
        """Local pulse scale Y_loc = eta/Q (static observer's proper time) of a QNEC-saturating CFT state with central charge c:
        flip iff sqrt f(r) < c eta / (12 pi Q r); the AdS_2 form is T sinh(rho/Q) < (c eta / 12) T_brk."""
        r = self.r_at_depth(rho)
        sf = np.sqrt(self.f(r))
        need = c * eta / (12 * np.pi * self.Q * r)
        return dict(r=r, sqrt_f=float(sf), threshold=float(need), flips=bool(sf < need),
                    T_max_approx=float(c * eta / 12 * self.T_brk / np.sinh(rho / self.Q)),
                    episode_over_Q_max=float(c / (12 * np.pi * r * sf) / self.Q))


def flip_temperature(c: float, Q: float, rho_over_Q: float = 1.0, eta: float = 1.0) -> float:
    """T_flip = (c eta / 12) T_brk / sinh(rho/Q): above it a throat-scale pulse cannot anti-trap spheres at depth rho."""
    return float(c * eta / 12 / (2 * np.pi**2 * Q**3) / np.sinh(rho_over_Q))


# ------------------------------------------------------------------ Proposition GSL-5: 2D CFT flips (exact s-wave model)
def _pulse(kind, Y, w, lc):
    x = lambda l: np.clip((l - lc) / w, -300, 300)
    if kind == "step":
        yp = lambda l: 0.5 * Y * (1 + np.tanh(x(l)))
        ypp = lambda l: 0.5 * Y / w / np.cosh(x(l)) ** 2
    elif kind == "bump":
        yp = lambda l: Y / np.cosh(x(l)) ** 2
        ypp = lambda l: -2 * Y / w * np.tanh(x(l)) / np.cosh(x(l)) ** 2
    else:
        raise ValueError(kind)
    return yp, ypp


def cft_flip(c: float, r0: float, rp0: float, Y: float, w: float, lam_max: float, kind: str = "step", lc: float = None,
             G: float = 1.0, hbar: float = 1.0, n: int = 20001) -> dict:
    """s-wave model on an ingoing radial null ray (affine lambda): r'' = -(G/r) T2, with the conformal-vacuum CFT state
        T2 = (c hbar / 24 pi)(y'^2/2 - y''),   S_out' = -(c/12) y',   y' = Y * pulse((lambda - lc)/w),
    which saturates the strong 2D QNEC.  Returns whether theta_- = 2 r'/r turns non-negative, the anti-trapped duration, its bound
    c/(6 |Theta_0|) (Theta_0 = 2 pi r_0 r'_0/(G hbar) + S'_0), and the largest increment of the 4D quantum expansion
    Theta = theta + 4 G hbar S'/A (quantum focusing: <= 0 up to discretisation)."""
    lc = 3 * w if lc is None else lc
    yp, ypp = _pulse(kind, Y, w, lc)
    T2 = lambda l: c * hbar / (24 * np.pi) * (0.5 * yp(l) ** 2 - ypp(l))
    lam = np.linspace(0, lam_max, n)
    sol = solve_ivp(lambda l, y: [y[1], -G / y[0] * T2(l)], (0, lam_max), [r0, rp0], t_eval=lam, rtol=1e-11, atol=1e-14,
                    max_step=min(w / 20, lam_max / 2000))
    r, rp = sol.y[0], sol.y[1]
    theta = 2 * rp / r
    Sp = -(c / 12) * yp(lam)
    Theta = theta + 4 * G * hbar * Sp / (4 * np.pi * r**2)
    anti = theta >= 0
    dur = float(np.sum(anti) * (lam[1] - lam[0]))
    Theta0 = 2 * np.pi * r0 * rp0 / (G * hbar) + Sp[0]
    return dict(flips=bool(anti.any()), theta_max=float(theta.max()), duration=dur,
                duration_bound=float(c / (6 * abs(Theta0))) if Theta0 < 0 else np.inf,
                max_dTheta=float(np.max(np.diff(Theta))), lam=lam, theta=theta, r=r, Theta=Theta,
                kick_bound=float(c * Y / (24 * np.pi * r0)))


def flip_duration_bound(c: float, Theta0: float) -> float:
    """Proposition GSL-5: an anti-trapped episode starting from Theta_0 < 0 lasts less than c/(6|Theta_0|) (affine, entropy units)."""
    return float(c / (6 * abs(Theta0)))


# ------------------------------------------------------------------ Proposition GSL-7: exact charged Vaidya flip
def vaidya_flip(rp: float, rm: float, c: float, Y_u: float, w_u: float, r_grid, kind: str = "step") -> dict:
    """Outgoing Bonnor-Vaidya-Reissner-Nordstrom, ds^2 = -(1 - 2m(u)/r + Q^2/r^2) du^2 - 2 du dr + r^2 dOmega^2, an exact solution of
    Einstein-Maxwell + null fluid with T_uu = -m'(u)/(4 pi r^2); initial horizons r_+, r_- (M = (r_+ + r_-)/2, Q^2 = r_+ r_-).
    With the CFT flux m'(u) = -(c/24 pi)(y'^2/2 - y'') (pulse in u), theta_+ = 2/r and theta_- = -f(u,r)/r, so a sphere is anti-trapped
    iff m(u) - M > r f(r)/2 with f = (r - r_+)(r - r_-)/r^2 (factored: no cancellation near extremality).  Returns the peak mass
    increase, its bound c Y_u / (24 pi), the budget r f(r)/2 (Theorem G1 in mass units) and the radii that become anti-trapped."""
    yp, ypp = _pulse(kind, Y_u, w_u, 3 * w_u)
    us = np.linspace(0, 12 * w_u + 6 / max(Y_u, 1e-300) if kind == "step" else 6 * w_u, 20001)
    mdot = -(c / (24 * np.pi)) * (0.5 * yp(us) ** 2 - ypp(us))
    dm = np.concatenate([[0.0], np.cumsum(0.5 * (mdot[1:] + mdot[:-1]) * np.diff(us))])
    r = np.asarray(r_grid, float)
    budget = r * (r - rp) * (r - rm) / r**2 / 2
    return dict(delta_m_peak=float(dm.max()), delta_m_bound=float(c * Y_u / (24 * np.pi)), budget=budget,
                antitrapped=r[dm.max() > budget], u=us, delta_m=dm)


# ------------------------------------------------------------------ Proposition GSL-9: Cauchy-horizon route near extremality
def _free(chi, kappa, M, q):
    try:
        return bool(ChargedBubble(chi, kappa, M, q * M).classical_route().get("barrier_free", False))
    except ValueError:                       # brentq bracket failure at the window edge = not barrier-free
        return False


def flip_between_horizons_bounds(chi: float, kappa: float, q: float) -> tuple:
    """Masses at which the outer flip of A_out ((chi^2 + kappa^2) R^4 = 2 M R - Q^2) sits exactly on r_- and on r_+ at Q = q M:
    (chi^2 + kappa^2) M^2 (1 +- s)^2 = 1, s = sqrt(1 - q^2).  The classical route needs r_- < R_flip < r_+, i.e. M between them."""
    L, s = chi**2 + kappa**2, np.sqrt(1 - q**2)
    return float(1 / (np.sqrt(L) * (1 + s))), float(1 / (np.sqrt(L) * (1 - s)))


def cauchy_route_window(chi: float, kappa: float, q: float, M_lo: float = 0.2, M_hi: float = 3.0, n_scan: int = 141) -> dict:
    """Barrier-free mass window [M_a, M_b] of the charged bubble at fixed Q/M = q (Proposition G4 route through the Cauchy horizon),
    with the horizon data at the window centre."""
    Ms = np.linspace(M_lo, M_hi, n_scan)
    fr = [M for M in Ms if _free(chi, kappa, M, q)]
    if not fr:
        return dict(exists=False)
    mid = fr[len(fr) // 2]

    def edge(a, b, inside_at_a):
        for _ in range(48):
            m = 0.5 * (a + b)
            if _free(chi, kappa, m, q) == inside_at_a:
                a = m
            else:
                b = m
        return 0.5 * (a + b)

    Ma = edge(mid, M_lo, True)
    Mb = edge(mid, M_hi * 1.5, True)
    rp = mid * (1 + np.sqrt(1 - q**2))
    rm = mid * (1 - np.sqrt(1 - q**2))
    Mc = 0.5 * (Ma + Mb)
    rpc, rmc = Mc * (1 + np.sqrt(1 - q**2)), Mc * (1 - np.sqrt(1 - q**2))
    return dict(exists=True, M_a=float(Ma), M_b=float(Mb), width=float(Mb - Ma), M_centre=float(Mc),
                horizon_gap_centre=float(rpc - rmc), kappa_minus_rm=float((rp - rm) / (2 * rm)),
                kappa_minus=float((rpc - rmc) / (2 * rmc**2)), M_extremal_limit=float(1 / np.sqrt(chi**2 + kappa**2)))


# ------------------------------------------------------------------ Proposition GSL-10: deep throats hide
def hiding_depth(dM_over_Q: float) -> dict:
    """Extremal Reissner-Nordstrom (Q = M = 1): energy dM that stays inside moves the event horizon to r_+ - Q = dM + sqrt(2 dM + dM^2).
    Returns the proper depth below the throat mouth r = 2Q (s(r) = (r - Q) + Q ln(r - Q), depth = s(2Q) - s(r)) beyond which a region
    ends up inside the horizon, exactly and in the approximation 1 - ln(2 dM/Q)/2."""
    d = dM_over_Q
    x = d + np.sqrt(2 * d + d**2)                                  # r_+ - Q in units of Q
    s = lambda y: y + np.log(y)                                    # s(r) with y = (r - Q)/Q
    return dict(depth_exact=float(s(1.0) - s(x)), depth_approx=float(1 - 0.5 * np.log(2 * d)), horizon_shift=float(x))


# ------------------------------------------------------------------ Remark GSL-11: electroweak corona
def electroweak_corona_bound(criterion: str = "higgs", alpha: float = ALPHA_EM, lam: float = 0.129, g: float = 0.652) -> dict:
    """Largest x = Lambda_EW Q^2 compatible with a symmetric-phase corona at the extremal horizon of a magnetic black hole.
    The corona needs e B(a) >= m_th^2 at the horizon radius a (AdS_2 branch, Q^2/a^2 = y(x) = 2x/(1 - sqrt(1-4x))) and
    B(Q)^2 = 1/(4 pi G Q^2); with rho_EW = lam v^4 / 4 this gives x <= K y(x)^2, K = 8 pi alpha rho_EW / m_th^4
    (K = pi alpha/(2 lam) for m_th = m_H, 32 pi alpha lam / g^4 for m_th = m_W, the latter with the full rho_EW, an overestimate).
    Returns K and the self-consistent x_max, to be compared with 1/4 (no static throat beyond it).  The bud threshold 3/16 refers to
    vacuum energy switched on in a Lambda = 0 throat and does not apply to a corona that is part of the static solution."""
    K = np.pi * alpha / (2 * lam) if criterion == "higgs" else 32 * np.pi * alpha * lam / g**4
    y = lambda x: 2 * x / (1 - np.sqrt(1 - 4 * x)) if x > 0 else 1.0
    F = lambda x: x - K * y(x) ** 2
    x_max = brentq(F, 1e-12, 0.25) if F(0.25) > 0 else 0.25
    return dict(K=float(K), x_max=float(x_max), decompactification=0.25, reaches_decompactification=bool(F(0.25) <= 0))


# ------------------------------------------------------------------ physical scales
def magnetic_charge_radius(q: float, alpha: float = ALPHA_EM) -> float:
    """Geometric charge radius Q = q l_P / (2 sqrt alpha) (metres) of a magnetic black hole with q Dirac quanta."""
    return float(q * L_P_SI / (2 * np.sqrt(alpha)))


def antitrapped_floor(R_m: float) -> dict:
    """An anti-trapped sphere of radius R has Misner-Sharp mass > R c^2 / 2G; a homogeneous region needs rho > 3 c^4/(8 pi G R^2)."""
    m = R_m * C_SI**2 / (2 * G_SI)
    rho = 3 * C_SI**4 / (8 * np.pi * G_SI * R_m**2)
    return dict(misner_sharp_kg=float(m), rho_J_m3=float(rho), rho_quarter_GeV=float((rho / GEV4_J_M3) ** 0.25))


def throat_window_SI(Q_m: float, c: float, eta: float = 1.0, rho_over_Q: float = 1.0) -> dict:
    """Temperatures (kelvin) bounding the semiclassical flip window of a throat of radius Q: T_brk (Schwarzian) and T_flip."""
    Qp = Q_m / L_P_SI
    T_brk = 1 / (2 * np.pi**2 * Qp**3)
    T_flip = flip_temperature(c, Qp, rho_over_Q, eta)
    return dict(Q_planck=float(Qp), T_brk_K=float(T_brk * T_P_K), T_flip_K=float(T_flip * T_P_K),
                window_ratio=float(T_flip / T_brk), mass_kg=float(Q_m * C_SI**2 / G_SI))


def bud_density(Q_m: float, x: float = 3 / 16) -> dict:
    """Vacuum energy density with Lambda Q^2 = x in a throat of radius Q (metres): rho = x c^4 / (8 pi G Q^2)."""
    rho = x * C_SI**4 / (8 * np.pi * G_SI * Q_m**2)
    return dict(rho_J_m3=float(rho), rho_quarter_GeV=float((rho / GEV4_J_M3) ** 0.25))
