"""Round 6: the mechanism atlas -- which physical routes to a daughter universe exist, what each must break, what each costs.

Four results organise the search (ROADMAP_ROUND6_MECHANISM_ATLAS.md):

* ``genesis_budget``       Theorem G1 (classical budget).  From the E_uu constraint of Prop. IV-1,
                           alpha^2 d_u(r_u/alpha^2) = -4 pi G r T_uu, an anti-trapped sphere can form from laboratory data only if,
                           along some ingoing null ray (v normalised so that r_v = 1/2 where the ray leaves the initial cone),
                               2 G m / r |_0  +  8 pi G int r (-T_uu) du / alpha^2  >=  1 :
                           the data must already be near their own Schwarzschild radius, or a NEC-violating "null force"
                           int r (-T_uu) du of order c^4 / (8 pi G) = 4.8e42 N must act.
* ``GeneralShell``         Theorem G2 (no bridge crossing for thin walls, any wall matter).  At the flip A_out = 0 the potential
                           is V = f_out, so a classical flip of the exterior orientation needs f_out <= 0: it happens inside a
                           horizon.  With growing R that is a white-hole region, whose past is a singularity (Schwarzschild), a
                           Cauchy horizon (Reissner-Nordstrom) or a cosmological region (SdS).  NEC violation confined to the wall
                           cannot produce a hidden daughter; it can only produce a *visible* inflating bubble, and that needs
                           negative asymptotic tension (Proposition G3).
* ``ChargedBubble``        Proposition G4 (charge).  A Reissner-Nordstrom exterior changes the laboratory exponent by <~2% but
                           supplies a bounce inside r_-; above the barrier-free mass the wall passes the Cauchy horizon classically and
                           expands into the next universe (Sakai et al. 2006).  Threshold M_cr ~ 0.16 m_P^3 / sqrt(rho).
* ``seed_lifetime``        Proposition G5 (stabilised seeds).  A quasi-stable seed attempts the transition at frequency omega, so
                           its lifetime is e^B / omega; at the semiclassical frontier B >= 5.28 sqrt(R/l_P) this is ~1e-18 s at
                           R = 100 l_P and ~0.1 s at R = 300 l_P.
* ``genesis_interface``    Proposition G7 (the edge of a genesis seed).  A seed bounded by a thin interface in vacuum keeps its exterior
                           orientation (Theorem G2); the interface energy must turn negative at H^2 R^2 = 2M/R and the seed becomes a
                           *visible* inflating bubble that grows into the parent.  A hidden daughter needs the throat to open inside a
                           NEC-violating boundary layer.
* ``genesis_timeline``     Proposition G6 (self-amplifying NEC violation).  In a Galilean-genesis phase H = h/(-t)^3 and
                           |p| = 3hc^4/(4 pi G (ct)^4) grow without external input; a seed of size L produces anti-trapped spheres when
                           H = 1/L, and the budget of G1 is then met automatically.

Units: G = c = hbar = 1 unless a function says SI or GeV.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

C_SI, G_SI, HBAR_SI = 2.99792458e8, 6.67430e-11, 1.054571817e-34
L_P_SI, T_P_SI = 1.616255e-35, 5.391247e-44
M_PL_GEV = 1.22e19
GEV4_J_M3 = 2.085e37                     # 1 GeV^4 / (hbar c)^3 in J/m^3
HBARC_GEV_M = 1.97327e-16


# ------------------------------------------------------------------ Theorem G1: the classical budget
def nec_force_SI() -> float:
    """c^4 / (8 pi G) in newtons (= 0.159 x the Gibbons maximum force c^4/4G)."""
    return C_SI**4 / (8 * np.pi * G_SI)


def genesis_budget(compactness0: float = 0.0) -> dict:
    """Theorem G1 in Planck units: required 8 pi G int r (-T_uu) du / alpha^2 >= 1 - 2Gm/r|_0, and the same in newtons."""
    need = max(1.0 - compactness0, 0.0)
    return dict(dimensionless=need, null_force_N=need * nec_force_SI(),
                fraction_of_max_force=need * (C_SI**4 / (8 * np.pi * G_SI)) / (C_SI**4 / (4 * G_SI)))


def stress_scale(L_m: float, compactness0: float = 0.0) -> dict:
    """NEC-violating stress |T_uu| ~ (1 - C0) c^4 / (8 pi G L^2) acting over a region of size L (r ~ du ~ L), in Pa and as an
    energy-density scale |T|^(1/4) in GeV."""
    p = max(1.0 - compactness0, 0.0) * C_SI**4 / (8 * np.pi * G_SI * L_m**2)
    return dict(stress_Pa=float(p), scale_GeV=float((p / GEV4_J_M3) ** 0.25))


def casimir_shortfall(a_m: float = 1e-8) -> dict:
    """Casimir energy density between plates, |rho| = pi^2 hbar c / (720 a^4), sampled over a region of size a: the null force
    |rho| a^2 compared with the budget c^4 / (8 pi G).  Orders of magnitude only."""
    rho = np.pi**2 * HBAR_SI * C_SI / (720 * a_m**4)
    F = rho * a_m**2
    return dict(energy_density_J_m3=float(rho), null_force_N=float(F), shortfall_orders=float(np.log10(nec_force_SI() / F)))


# ------------------------------------------------------------------ Theorem G2 / Proposition G3: thin walls of any matter
class GeneralShell:
    """Thin wall with interior f_in, exterior f_out and arbitrary surface energy kappa(R) = 4 pi G sigma(R) (surface stress fixed
    by conservation, p = -sigma - R sigma'/2).  The first junction condition gives, at each R,
        A_in = (f_in - f_out + kappa^2 R^2)/(2 kappa R),  A_out = (f_in - f_out - kappa^2 R^2)/(2 kappa R),  V = f_out - A_out^2.
    NEC on the wall: sigma >= 0 and sigma + p = -R sigma'/2 >= 0."""

    def __init__(self, f_in, f_out, kappa):
        self.f_in, self.f_out, self.kappa = f_in, f_out, kappa

    def A_in(self, R):
        k = self.kappa(R)
        return (self.f_in(R) - self.f_out(R) + k**2 * R**2) / (2 * k * R)

    def A_out(self, R):
        k = self.kappa(R)
        return (self.f_in(R) - self.f_out(R) - k**2 * R**2) / (2 * k * R)

    def V(self, R):
        return self.f_out(R) - self.A_out(R) ** 2

    def nec_on_wall(self, R, h=1e-6):
        k = self.kappa(R)
        dk = (self.kappa(R * (1 + h)) - self.kappa(R * (1 - h))) / (2 * R * h)
        return bool(k >= 0 and -R * dk / 2 >= 0)

    def flip_radii(self, lo, hi, n=20000):
        R = np.geomspace(lo, hi, n)
        a = np.array([self.A_out(r) for r in R])
        idx = np.where(np.sign(a[:-1]) * np.sign(a[1:]) < 0)[0]
        out = []
        for i in idx:
            try:
                out.append(brentq(self.A_out, R[i], R[i + 1], xtol=1e-14))
            except ValueError:      # sign change through a pole of A_out (kappa = 0), not a flip
                pass
        return out

    def bridge_crossing_report(self, lo, hi) -> list:
        """At every flip of the exterior orientation: V, f_out and whether a classical crossing is possible (needs V <= 0,
        hence f_out <= 0, i.e. a horizon interior)."""
        rows = []
        for Rf in self.flip_radii(lo, hi):
            rows.append(dict(R_flip=float(Rf), V=float(self.V(Rf)), f_out=float(self.f_out(Rf)),
                             classical_crossing_possible=bool(self.V(Rf) <= 0), inside_horizon=bool(self.f_out(Rf) <= 0)))
        return rows


def visible_bubble_orientation(chi: float, kappa_inf: float, R: float = 1e3) -> int:
    """Sign of A_out at large R for interior de Sitter and Schwarzschild exterior: -sign(kappa_inf).  A visible inflating
    bubble (exterior orientation +1 as R -> infinity, no bridge) needs negative asymptotic tension (Proposition G3)."""
    sh = GeneralShell(lambda r: 1 - chi**2 * r**2, lambda r: 1 - 0.0 / r, lambda r: kappa_inf)
    return int(np.sign(sh.A_out(R)))


def upward_nucleation_exponent(chi: float, kappa: float) -> dict:
    """M = 0 bubble with interior de Sitter, flat exterior.  The barrier is [0, R_*], R_* = 2|kappa|/(chi^2 + kappa^2), and the
    canonical exponent B = 2 int_0^{R_*} R theta dR depends on kappa^2 only: a negative-tension wall nucleates a *visible* expanding
    de Sitter bubble in Minkowski space with the same exponent as the positive-tension CdL bounce (Proposition G3)."""
    k = abs(kappa)
    Rs = 2 * k / (chi**2 + k**2)

    def integrand(t):
        R = Rs * (1 - np.cos(t)) / 2
        dR = Rs * np.sin(t) / 2
        fi = 1 - chi**2 * R**2
        Ao = (fi - 1 - k**2 * R**2) / (2 * k * R) if R > 0 else 0.0
        V = 1 - Ao**2
        c = (fi + 1 - k**2 * R**2) / (2 * np.sqrt(fi))
        s = k * R * np.sqrt(max(V, 0.0) / fi)
        return R * np.arctan2(s, c) * dR

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        B = 2 * quad(integrand, 0, np.pi, limit=400)[0]
    B_flat = np.pi * k * (2 * k / (chi**2 + k**2)) ** 3 / 8
    return dict(R_star=float(Rs), B=float(B), B_thin_flat=float(B_flat))


def vacuum_survival_exponent(R_star_lP: float = 1.0, H0_inv_lP: float = 8.7e60) -> float:
    """Minimum exponent for upward nucleation not to have happened in our past light cone: B > ln(H0^-4 / R_*^4)."""
    return float(4 * np.log(H0_inv_lP / R_star_lP))


# ------------------------------------------------------------------ Proposition G4: charged exteriors
class ChargedBubble:
    """Interior de Sitter (chi), wall kappa, exterior Reissner-Nordstrom f_out = 1 - 2M/R + Q^2/R^2 (G = 1).  The junction is the
    same as for Schwarzschild (the Maxwell stress is a bulk term)."""

    def __init__(self, chi: float, kappa: float, M: float, Q: float):
        self.chi, self.kappa, self.M, self.Q = chi, kappa, M, Q
        self.shell = GeneralShell(lambda R: 1 - chi**2 * R**2, lambda R: 1 - 2 * M / R + Q**2 / R**2, lambda R: kappa)

    def horizons(self):
        d = self.M**2 - self.Q**2
        return (self.M - np.sqrt(d), self.M + np.sqrt(d)) if d > 0 else None

    def zeros(self, lo=1e-5, hi=None, n=60000):
        hi = hi if hi is not None else 5.0 / self.chi
        R = np.geomspace(lo, hi, n)
        v = np.array([self.shell.V(r) for r in R])
        idx = np.where(np.sign(v[:-1]) * np.sign(v[1:]) < 0)[0]
        return [brentq(self.shell.V, R[i], R[i + 1], xtol=1e-14) for i in idx]

    def laboratory_barrier(self):
        """(R1, R2) with V > 0 between, A_out(R1) > 0 (laboratory side) and A_out(R2) < 0 (daughter side); None if absent."""
        z = self.zeros(hi=0.99999 / self.chi)
        for a, b in zip(z[:-1], z[1:]):
            if self.shell.V(0.5 * (a + b)) > 0 and self.shell.A_out(a) > 0 > self.shell.A_out(b):
                return a, b
        return None

    def B_lab(self):
        br = self.laboratory_barrier()
        if br is None:
            return None
        a, b = br
        fi, fo, k = self.shell.f_in, self.shell.f_out, self.kappa

        def integrand(t):
            R = a + (b - a) * (1 - np.cos(t)) / 2
            dR = (b - a) * np.sin(t) / 2
            c = (fi(R) + fo(R) - k**2 * R**2) / (2 * np.sqrt(fi(R) * fo(R)))
            s = k * R * np.sqrt(max(self.shell.V(R), 0.0) / (fi(R) * fo(R)))
            return R * np.arctan2(s, c) * dR

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return float(2 * quad(integrand, 0, np.pi, limit=400)[0])

    def classical_route(self) -> dict:
        """Barrier-free case: the allowed set contains the whole range from the inner bounce R_min (inside r_-) to R = 5/chi, and the
        outer flip of A_out lies between the horizons (f_out < 0), so the wall leaves the next universe's white-hole region with
        exterior orientation -1 and expands (Sakai et al. route through the Cauchy horizon)."""
        hz = self.horizons()
        if hz is None:
            return dict(barrier_free=False)
        rminus, rplus = hz
        A, V, fo = self.shell.A_out, self.shell.V, self.shell.f_out
        # inner flip: A_out ~ (2MR - Q^2)/(2 kappa R^2) near R -> 0, so it changes sign at R_fi < r_- (r_- > Q^2/2M); there V = f_out > 0
        # and V(r_-) = -A_out^2 < 0, so the bounce is the zero of V in (R_fi, r_-).  The barrier is O(Q^2/M) thin: no scan.
        if np.sign(A(1e-6 * rminus)) == np.sign(A(rminus * (1 - 1e-12))):
            # r_- > 1/L: A_out < 0 on the whole of (0, r_-), there is no inner flip and no bounce with exterior orientation +1
            return dict(barrier_free=False, r_minus=float(rminus), r_plus=float(rplus), inner_flip=None, R_bounce=None)
        R_fi = brentq(A, 1e-6 * rminus, rminus * (1 - 1e-12), xtol=1e-16 * rminus)
        R_bounce = brentq(V, R_fi, rminus * (1 - 1e-12), xtol=1e-16 * rminus) if V(R_fi) > 0 else None
        Rs = np.geomspace(rminus * (1 + 1e-9), 5.0 / self.chi, 20000)
        vmax = max(V(r) for r in Rs)
        outer_flips = [f for f in self.shell.flip_radii(rminus * (1 + 1e-9), 5.0 / self.chi)]
        return dict(barrier_free=bool(vmax < 0), r_minus=float(rminus), r_plus=float(rplus), inner_flip=float(R_fi),
                    R_bounce=None if R_bounce is None else float(R_bounce), max_V_outside_r_minus=float(vmax),
                    outer_flip=None if not outer_flips else float(outer_flips[0]),
                    flip_inside_horizons=bool(outer_flips and rminus < outer_flips[0] < rplus and fo(outer_flips[0]) < 0))


def critical_mass(chi: float, kappa: float, Q_over_M: float = 0.5, lo: float = 0.05, hi: float = 3.0) -> float:
    """Smallest M (G = 1) at which a charged bubble is barrier-free outside r_- (bisection on ChargedBubble.classical_route)."""
    def free(M):
        return ChargedBubble(chi, kappa, M, Q_over_M * M).classical_route().get("barrier_free", False)
    if free(lo) or not free(hi):
        raise ValueError("bracket does not contain the threshold")
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if free(mid) else (mid, hi)
    return float(hi)


def classical_threshold_physical(E_vac_GeV: float, coefficient: float = 0.45) -> dict:
    """M_cr = c_cr / (G chi), chi = sqrt(8 pi G rho / 3): the black-hole mass whose radius is ~ the false-vacuum Hubble radius.
    Returns M_cr in kg and GeV and the Hubble radius in m (coefficient from critical_mass, 0.39-0.5 for kappa/chi = 1-0.3)."""
    rho = (E_vac_GeV / M_PL_GEV) ** 4                       # Planck units
    chi = np.sqrt(8 * np.pi * rho / 3)
    M = coefficient / chi                                   # Planck masses
    m_P_kg = 2.176434e-8
    return dict(M_cr_kg=float(M * m_P_kg), M_cr_GeV=float(M * M_PL_GEV), hubble_radius_m=float(L_P_SI / chi),
                energy_J=float(M * m_P_kg * C_SI**2))


# ------------------------------------------------------------------ Proposition G5: stabilised seeds
def seed_lifetime(B: float, R_lP: float) -> dict:
    """Lifetime e^B / omega of a quasi-stable seed with attempt frequency omega ~ c / R (seconds)."""
    omega = 1.0 / (R_lP * T_P_SI)
    return dict(omega_Hz=float(omega), lifetime_s=float(np.exp(B) / omega))


def frontier_seed_window(B_CdL: float = 10.0, M_planck: float = 1.0, Rs=(100, 200, 300, 350, 400, 1000)) -> list:
    """Lifetimes of stabilised seeds sitting exactly at the semiclassical frontier B = 1.671 sqrt(B_CdL M R)."""
    from .bounds import FRONTIER_CONST
    rows = []
    for R in Rs:
        B = FRONTIER_CONST * np.sqrt(B_CdL * M_planck * R)
        rows.append(dict(R_lP=R, B=float(B), **seed_lifetime(B, R)))
    return rows


# ------------------------------------------------------------------ Proposition G6: self-amplifying NEC violation
def genesis_timeline(L_m: float, h_m2: float, p0_Pa: float) -> dict:
    """Galilean-genesis scaling H = h/(-t)^3 (t < 0), dH/dt = 4 pi G |p| / c^2, so |p| = 3 h c^4 /(4 pi G (ct)^4) (SI; h in m^2, times as ct in m).
    A seed of size L produces anti-trapped spheres when H = c/L, at |t_*| = (h L)^(1/3); validity of the flat-space genesis
    background there needs H|t| = (h/L^2)^(1/3) << 1.  Starting from a seed stress p0, the elapsed time is |t_0| - |t_*| with
    |t_0| = |t_*| (p_*/p0)^(1/4).  Returns the G1 budget check |p_*| L^2 / (c^4/8 pi G)."""
    t_star = (h_m2 * L_m) ** (1 / 3)                         # metres (c t)
    p_star = 3 * h_m2 * C_SI**4 / (4 * np.pi * G_SI * t_star**4)
    t0 = t_star * (p_star / p0_Pa) ** 0.25
    return dict(t_star_s=float(t_star / C_SI), p_star_Pa=float(p_star), elapsed_s=float((t0 - t_star) / C_SI),
                amplification=float(p_star / p0_Pa), validity_H_t=float((h_m2 / L_m**2) ** (1 / 3)),
                budget_ratio=float(p_star * L_m**2 / nec_force_SI()))


# ------------------------------------------------------------------ Proposition G7: the edge of a genesis seed
def genesis_interface(h: float, r_c: float, M: float, n: int = 4000) -> dict:
    """Interior: flat FRW genesis patch H = h/(-t)^3, a = exp(h/2t^2) (a -> 1 as t -> -inf); comoving thin interface at r_c,
    R = a r_c; exterior Schwarzschild of mass M (G = c = 1).  For a comoving wall in flat FRW, A_in^2 = 1 - 2 m_in/R + Rdot^2 = 1
    exactly (m_in = R^3 H^2 / 2, Rdot = H R), so the junction A_in - A_out = kappa R gives
        visible branch (exterior orientation +1):  kappa R = 1 - sqrt(1 - 2M/R + H^2 R^2),
        hidden branch  (exterior orientation -1):  kappa R = 1 + sqrt(1 - 2M/R + H^2 R^2).
    The branches meet only where 1 - 2M/R + H^2 R^2 = 0, i.e. R < 2M (inside the horizon; Theorem G2).  Starting from regular data
    (visible branch), the interface energy must turn negative at H^2 R^2 = 2M/R, before the interior anti-traps (H R = 1), and the
    seed ends as a visible inflating bubble that grows into the parent."""
    t_end = -(h ** 0.5) * 0.999                              # stay where the genesis scaling holds roughly (H|t| < 1)
    t = -np.geomspace(1e4 * h ** (1 / 3) + 1.0, -t_end, n)
    H = h / (-t) ** 3
    R = np.exp(h / (2 * t**2)) * r_c
    disc = 1 - 2 * M / R + H**2 * R**2
    kR_vis = 1 - np.sqrt(disc)
    i_neg = int(np.argmax(kR_vis < 0)) if np.any(kR_vis < 0) else None
    i_at = int(np.argmax(H * R > 1)) if np.any(H * R > 1) else None
    return dict(min_discriminant=float(disc.min()), branches_meet=bool(disc.min() <= 0),
                HR_when_tension_turns_negative=None if i_neg is None else float(H[i_neg] * R[i_neg]),
                HR_predicted=float(np.sqrt(2 * M / r_c)), antitrapped_reached=i_at is not None,
                tension_negative_before_antitrapping=bool(i_neg is not None and (i_at is None or i_neg <= i_at)),
                final_kappaR=float(kR_vis[-1]), final_HR=float(H[-1] * R[-1]))
