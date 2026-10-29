"""Round 6k: the two doors.  Y1, black-hole interiors that turn into expanding daughters, in definite theories; Y2, black holes at collider
energies; and Y1 + Y2, a conditional machine specification.

Y1 (ROADMAP_ROUND6K.md, section 1):
* ``ecsk_bounce``            Einstein-Cartan-Sciama-Kibble gravity (torsion sourced by fermion spin; identical to GR below extreme densities).
                             Spin fluid: epsilon_eff = epsilon - kappa s^2/4, s^2 = hbar^2 c^2 n^2 / 8 (unpolarised fermions), kappa = 8 pi G/c^4,
                             so a nucleon fluid bounces at n_b = 4 m c^4/(pi G hbar^2), rho_b = 4 m^2 c^4/(pi G hbar^2) = 3.9e58 kg/m^3 [Proposition,
                             from the standard spin-fluid equations; Poplawski quotes deviations from GR above ~1e45 kg/m^3].
* ``ecsk_min_hole``          A horizon forms before the bounce only if 2GM/c^2 > (3M / 4 pi rho_b)^{1/3}: M > M_min = (3/4 pi rho_b)^{1/2}
                             (c^2/2G)^{3/2} = 4.3e10 kg for rho_b = 3.9e58 kg/m^3 (2.7e17 kg for 1e45).  Every astrophysical black hole is far
                             above it: in ECSK every one of them holds a bounce.
* ``ecsk_interior``          [Proposition-N 6k-1] Oppenheimer-Snyder interior (closed FRW dust) with the spin term: (da/dt)^2 = A/a - S/a^4 - 1
                             + a^2 / l_d^2.  It bounces at a_min (= (S/A)^{1/3} to leading order) and re-expands beyond the horizon (into the
                             white-hole region of the extended exterior: a new universe).  Without vacuum energy it returns to the original
                             size and cycles; a contracting region carrying vacuum energy with l_d < (3 sqrt3/2) A bounces and then inflates
                             for ever.  For a stellar collapse with our dark energy (l_d ~ 1e26 m >> A ~ 1e4 m) it cycles: a lasting cosmos needs
                             growth from cycle to cycle (Poplawski: particle production).
* ``fmm_D``                  [Proposition 6k-2] Limiting curvature in D = 4 + n dimensions: the Schwarzschild-Tangherlini interior
                             f = 1 - (r_s/r)^{n+1} joins de Sitter f = 1 - r^2/l^2 at r_0 = (r_s^{n+1} l^2)^{1/(n+3)}, valid iff r_s > l, and
                             r_0/l = (r_s/l)^{(n+1)/(n+3)} > 1: the de Sitter region starts larger than its Hubble length and inflates.
Y2 and Y1 + Y2 (sections 2-3):
* ``add_daughter_machine``   [Proposition 6k-3, conditional] With limiting length l = c_l / M_D, a collider black hole holds an FMM daughter iff
                             R_S > l: M_BH > M_Y = threshold mass.  For c_l = 1, M_Y = 3.9 M_D (n = 4) to 4.7 M_D (n = 6): inside the reach of
                             a 100 TeV collider for M_D <~ 10-15 TeV.  The interior reaches r_0 in a proper time ~R_S/c, shorter than the
                             Hawking lifetime ~(1/M_D)(M/M_D)^{(n+3)/(n+1)} by a factor ~5-10.  Semiclassical control (R_S >= 3 l) needs
                             M_BH ~ 1e3-1e4 M_D.
* ``collider_reach``         Scaling from the LHC exclusion edge (tau = M^2/s ~ 0.5): a 100 TeV collider reaches M_BH ~ 70 TeV > M_Y = 24-57 TeV.
* ``pbh_partner``            A 5-10 Earth-mass primordial black hole (Planet Nine hypothesis, Scholtz & Unwin 2020) has Nariai partner vacuum
                             ~60-85 GeV: within a factor 1.2-1.7 of the electroweak vacuum energy (104.5 GeV)^4, whose partner is 3.3 Earth
                             masses.  Laser-launched probes at 0.001c reach 500 AU in 7.9 years of flight (Witten 2020: about a decade).
* ``regular_core_temperature`` [Proposition-N 6k-4] The class with a static de Sitter core (Hayward-type, D dimensions) has T(M) rising, peaking and
                             falling to zero at a remnant mass: a collider signature that separates it from GR and from the FMM class (whose
                             exterior evaporation is the GR one).
"""

from __future__ import annotations

import numpy as np
from math import gamma
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G_SI = 6.67430e-11
C_SI = 2.99792458e8
HBAR_SI = 1.054571817e-34
M_N_KG = 1.67492750e-27
M_SUN = 1.98847e30
HBAR_GEV_S = 6.582119569e-25
HBARC_GEV_M = 1.97327e-16


# --------------------------------------------------------------------------- Y1: Einstein-Cartan-Sciama-Kibble


def ecsk_bounce(m_kg=M_N_KG) -> dict:
    """Bounce of a spin-1/2 fluid of particles of mass m: m c^2 n = (kappa/4)(hbar^2 c^2 n^2 / 8)  ->  n_b = 4 m c^4/(pi G hbar^2)."""
    n_b = 4 * m_kg * C_SI**4 / (np.pi * G_SI * HBAR_SI**2)
    rho_b = n_b * m_kg
    return dict(m_kg=m_kg, n_b_m3=float(n_b), rho_b_kg_m3=float(rho_b), rho_over_planck=float(rho_b / 5.155e96),
                rho_over_nuclear=float(rho_b / 2.3e17))


def ecsk_min_hole(rho_b_kg_m3) -> dict:
    """Smallest collapsing mass that forms a horizon before reaching the bounce density."""
    M = np.sqrt(3 / (4 * np.pi * rho_b_kg_m3)) * (C_SI**2 / (2 * G_SI)) ** 1.5
    return dict(rho_b=rho_b_kg_m3, M_min_kg=float(M), r_s_m=float(2 * G_SI * M / C_SI**2))


def ecsk_inflation_criterion(A=1.0) -> float:
    """Closed dust FRW with vacuum energy (adot^2 = A/a - 1 + a^2/l_d^2): no turnaround (eternal expansion) iff l_d < (3 sqrt3 / 2) A."""
    return float(1.5 * np.sqrt(3) * A)


def ecsk_interior(A=1.0, S=1e-6, l_d=None, T=60.0) -> dict:
    """Closed FRW dust interior with the spin term, in units of the turnaround scale: adot^2 = A/a - S/a^4 - 1 + a^2/l_d^2, i.e.
    a'' = -A/(2a^2) + 2S/a^5 + a/l_d^2.  Without vacuum energy, start at the GR turnaround (adot = 0).  With vacuum energy, start contracting
    at a = A (adot = -sqrt(V(A))): the region collapses, bounces at a_min and then either turns around (l_d > 2.6 A) or inflates."""
    def Vfun(a):
        return A / a - S / a**4 - 1 + (0 if l_d is None else a**2 / l_d**2)
    if l_d is None:
        a0, ad0 = brentq(lambda a: Vfun(a), 0.5 * A, 1.5 * A), -1e-9
    else:
        a0 = A
        v = Vfun(a0)
        if v <= 0:
            return dict(A=A, S=S, l_d=l_d, fate="no contracting solution through a = A")
        ad0 = -np.sqrt(v)
    def rhs(t, y):
        a, ad = y
        return [ad, -A / (2 * a**2) + 2 * S / a**5 + (0 if l_d is None else a / l_d**2)]
    def big(t, y):
        return y[0] - (50 * l_d if l_d else 1e9)
    big.terminal = True
    s = solve_ivp(rhs, (0, T), [a0, ad0], rtol=1e-11, atol=1e-14, max_step=0.01, events=[big])
    a = s.y[0]
    i_min = int(np.argmin(a))
    a_min = float(a[i_min])
    after = a[i_min:]
    if l_d is not None and (len(s.t_events[0]) or after[-1] > 10 * l_d):
        fate, a_max2 = "bounces, then inflates", None
    else:
        fate, a_max2 = "bounces, recollapses (cycles)", float(after.max())
    a_min_analytic = brentq(lambda x: A * x**3 - S - x**4, 1e-12, (S / A) ** (1 / 3) * 2)
    return dict(A=A, S=S, l_d=l_d, a_start=float(a0), a_min=a_min, a_min_analytic=float(a_min_analytic), a_max_after=a_max2, fate=fate,
                bounce_time=float(s.t[i_min]))


# --------------------------------------------------------------------------- Y1 in D dimensions: limiting curvature


def fmm_D(n, r_s, l=1.0) -> dict:
    """Schwarzschild-Tangherlini interior f_S = 1 - (r_s/r)^{n+1} meets de Sitter f_dS = 1 - r^2/l^2 (both T-regions) at f_S = f_dS."""
    r0 = (r_s ** (n + 1) * l**2) ** (1 / (n + 3))
    return dict(n=n, r_s=r_s, l=l, r0=float(r0), valid=bool(l < r0 < r_s), r0_over_l=float(r0 / l))


# --------------------------------------------------------------------------- Y2 and Y1 + Y2


def _rs_MD(M_over_MD, n):
    x = 8 * gamma((n + 3) / 2) / (n + 2) * M_over_MD
    return x ** (1 / (n + 1)) / np.sqrt(np.pi)


def add_daughter_machine(n=6, M_D_TeV=8.0, c_l=1.0, safety=1.0) -> dict:
    """Threshold black-hole mass for an FMM daughter (R_S > safety * l, l = c_l / M_D), its geometric parton cross-section, Hawking
    temperature, lifetime (order of magnitude, (1/M_D)(M/M_D)^{(n+3)/(n+1)}) and the interior proper time to the junction (~R_S)."""
    target = safety * c_l
    M_Y = brentq(lambda x: _rs_MD(x, n) - target, 1e-3, 1e12)
    R = _rs_MD(M_Y, n) / (M_D_TeV * 1e3)                       # GeV^-1
    T = (n + 1) / (4 * np.pi * R)                              # GeV
    tau = (1 / (M_D_TeV * 1e3)) * M_Y ** ((n + 3) / (n + 1))     # GeV^-1
    sigma_pb = np.pi * (R * HBARC_GEV_M) ** 2 / 1e-40
    return dict(n=n, M_D_TeV=M_D_TeV, c_l=c_l, safety=safety, M_Y_over_MD=float(M_Y), M_Y_TeV=float(M_Y * M_D_TeV),
                R_S_m=float(R * HBARC_GEV_M), T_GeV=float(T), lifetime_s=float(tau * HBAR_GEV_S), interior_time_s=float(R * HBAR_GEV_S),
                lifetime_over_interior=float(tau / R), sigma_parton_pb=float(sigma_pb))


def regular_core_temperature(n=6, l=1.0, n_pts=400) -> dict:
    """D-dimensional Hayward-type black hole f = 1 - mu r^2 / (r^{n+3} + mu l^2), mu = r_s^{n+1}: horizon radius r_h(mu), temperature
    T = f'(r_h)/4 pi, and the remnant (extremal) mass where T -> 0.  Mass ~ mu (units with the D-dimensional constants set to 1)."""
    def f(r, mu):
        return 1 - mu * r**2 / (r ** (n + 3) + mu * l**2)
    def fp(r, mu):
        den = r ** (n + 3) + mu * l**2
        return -(2 * mu * r * den - mu * r**2 * (n + 3) * r ** (n + 2)) / den**2
    # extremal: f = 0 and f' = 0 simultaneously
    def extremal(mu):
        rs = np.geomspace(1e-3 * l, 1e3 * l, 4000)
        return np.min(f(rs, mu))
    mu_rem = brentq(extremal, 1e-6, 1e8)
    rows = []
    for mu in np.geomspace(mu_rem * 1.0001, mu_rem * 1e4, n_pts):
        rs = np.geomspace(1e-3 * l, 1e4 * l, 20000)
        fv = f(rs, mu)
        idx = np.where(fv[:-1] * fv[1:] < 0)[0]
        if not len(idx):
            continue
        r_h = brentq(lambda r: f(r, mu), rs[idx[-1]], rs[idx[-1] + 1])        # outer horizon
        rows.append(dict(mu=float(mu), r_h=float(r_h), T=float(fp(r_h, mu) / (4 * np.pi))))
    Ts = [r["T"] for r in rows]
    k = int(np.argmax(Ts))
    return dict(n=n, l=l, mu_remnant=float(mu_rem), mu_Tmax=rows[k]["mu"], T_max=rows[k]["T"], rows=rows,
                gr_T_at_large_mu=float((n + 1) / (4 * np.pi * rows[-1]["r_h"])), T_large_mu=rows[-1]["T"])


def collider_reach(n=6, M_D_TeV=8.0, sqrt_s_TeV=100.0, tau_reach=0.5) -> dict:
    """Scaling estimate: the LHC (13 TeV, 138 fb^-1) excludes black holes up to tau = M^2/s ~ 0.4-0.75 (8.4-11.4 TeV); a collider of
    energy sqrt(s) with comparable luminosity reaches M ~ sqrt(tau_reach) sqrt(s).  Daughter threshold M_Y (R_S = 1/M_D) against it."""
    m = add_daughter_machine(n, M_D_TeV)
    M_reach = np.sqrt(tau_reach) * sqrt_s_TeV
    return dict(n=n, M_D_TeV=M_D_TeV, sqrt_s_TeV=sqrt_s_TeV, M_reach_TeV=float(M_reach), M_Y_TeV=m["M_Y_TeV"],
                reachable=bool(m["M_Y_TeV"] <= M_reach), tau_at_threshold=float((m["M_Y_TeV"] / sqrt_s_TeV) ** 2),
                sqrt_s_needed_TeV=float(m["M_Y_TeV"] / np.sqrt(tau_reach)))


def pbh_partner(M_earth=5.0) -> dict:
    """Nariai partner vacuum of a primordial black hole of M_earth Earth masses (rho_f = 1/(72 pi M^2)), and its size; compared with the
    electroweak vacuum energy V(0) = lambda v^4 / 4 = (104.5 GeV)^4 of the Standard Model Higgs."""
    from . import nariai as NA
    M = M_earth * 5.972e24
    vac = NA.partner_vacuum(M)
    r_s = 2 * G_SI * M / C_SI**2
    return dict(M_earth=M_earth, M_kg=M, r_s_m=float(r_s), r_N_m=float(1.5 * r_s), partner_vac_GeV=float(vac),
                ew_vac_GeV=float((0.129 * 246.22**4 / 4) ** 0.25), ew_partner_M_earth=float(NA.partner_mass((0.129 * 246.22**4 / 4) ** 0.25) / 5.972e24),
                travel_years_at_0p001c=float(500 * 1.495978707e11 / (1e-3 * C_SI) / 3.15576e7))


def y1_theories() -> list:
    """What each definite framework puts inside a black hole (literature; ROADMAP_ROUND6K.md section 1)."""
    return [
        dict(theory="Einstein-Cartan-Sciama-Kibble (Poplawski 2010-2025)", interior="bounce -> new closed universe beyond the horizon; "
             "cycles grow with particle production until dark energy dominates", daughter="expanding", status="classical, passes all tests"),
        dict(theory="Limiting curvature (Frolov-Markov-Mukhanov 1989/90; Barrabes-Frolov 1996; Easson-Brandenberger 2001)",
             interior="spacelike transition to de Sitter; possibly several disconnected de Sitter worlds", daughter="inflating",
             status="effective (higher-derivative constructions exist)"),
        dict(theory="Mimetic limiting curvature (Chamseddine-Mukhanov 2017, 2019)", interior="nested Schwarzschild 'dolls' / static de Sitter "
             "patch; stable remnants after evaporation", daughter="static core (not expanding)", status="effective"),
        dict(theory="Loop quantum gravity effective (Ashtekar-Olmedo-Singh 2018)", interior="transition surface -> anti-trapped (white-hole) "
             "region -> new asymptotic region", daughter="new region, not inflating", status="effective quantum geometry"),
        dict(theory="Holography (Cooper et al. 2018; Sasieta-Swingle-Vilar Lopez 2026)", interior="certain microstates: closed FRW universe "
             "behind the horizon (big bang / big crunch)", daughter="crunching", status="controlled quantum gravity (AdS, JT/SYK)"),
    ]
