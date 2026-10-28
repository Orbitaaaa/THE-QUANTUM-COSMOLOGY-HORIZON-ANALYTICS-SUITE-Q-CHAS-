"""Round 6j, Groups B-D: let an environment prepare the state, quantum spacetime without a classical throat, and a different gravitational
theory.  G = c = hbar = 1 unless stated (8 pi G = 1 in the wormhole section).

Group B (environment prepares the state):
* ``accretion_requirement``  #6 [Proposition 6j-6] Accretion restores a vacuum whose Hubble length matches the hole (the near-Nariai partner,
                             rho_f = c^6/(72 pi G^3 M^2)) only if the infalling energy density reaches ~rho_f at the horizon: 10^18.8-10^19.7
                             times an Eddington flow; neutron-star merger debris (0.1 M_sun within 3 r_s) supplies only 0.25-1.8% of it.  Marsden et al. (2025) find
                             restored shells and outward-growing (parent-eating) decay bubbles; a hidden daughter still needs the NEC escape
                             (Theorem 6h-1 holds for collapse-formed holes).
* ``bag_star``               #7 [Proposition-N 6j-7] The most compact self-bound false-vacuum star (MIT-bag equation of state, TOV) has
                             2GM/(c^2 R) = 0.54, radius 0.28 l_f, vacuum share 14%: sub-Hubble, so its collapse ends in a black hole
                             (Theorem 6h-1).  The QCD bag constant B^{1/4} = 145 MeV is exactly the Nariai partner of a 5.2 M_sun black hole.
* ``oy_channel``             #8 [Proposition-N 6j-8] Oshita-Yokoyama, both steps: a false-vacuum bubble around a hole of mass M_in is made
                             thermally (weight e^{-Delta M / T_H}, T_H = 1/8 pi M_in) and tunnels to the hidden inflating branch (B, thin
                             shell with the hole inside).  Total exponent Sigma = Sigma_hat (l_f / l_P)^2; integrated over a hole's evaporation
                             the expected number of daughters is ~1.6e4 M_in^2 e^{-Sigma}.  Sigma_hat = 0.013 (wall k_w = 0.5)
                             to ~4e-4 (k_w = 0.1): a GUT vacuum ((l_f/l_P)^2 ~ 5e10) gives Sigma >~ 2e7; O(1) needs l_f <~ 10^2-10^3 l_P.
Group C (quantum spacetime):
* ``wineglass_family``       #11 [Proposition-N 6j-9] Euclidean axion wormholes with a cosmological constant (Giddings-Strominger + Lambda):
                             a'^2 = 1 - (a^2/3)(Lambda + q^2/a^6).  The half-wormhole action S(q) -> -12 pi^2 / Lambda (no-boundary) as the
                             throat pinches off (q -> 0), reproducing the Lavrelashvili-Lehners limit.  The excess is Delta S = S(q) - S(0)
                             = 3 pi^3 a_0^2 (8 pi G = 1), independent of Lambda to 1e-5: Delta S = (3 pi^2/8)(a_0/l_P)^2 = 1.18 x (4 pi a_0^2)/4G.
                             A channel that stays connected through a throat pays action in proportion to the throat's area; the throat-free
                             limit is the disconnected (no-boundary) universe.
Group D (different gravity):
* ``fmm_junction``           #15 [Proposition 6j-10] Frolov-Markov-Mukhanov: a Schwarzschild interior joins de Sitter of Hubble length l across
                             a spacelike surface at r_0 = (2 M l^2)^{1/3} (continuous K_theta_theta) for M > l/2; the de Sitter region starts
                             larger than its Hubble length (r_0/l = (2M/l)^{1/3}) and inflates.  If quantum gravity does this, every black hole
                             contains a daughter, and making a daughter reduces to making a black hole.
* ``add_black_holes``        #14 [Proposition 6j-11, conditional] With n large extra dimensions and fundamental scale M_D, black holes form at
                             M_BH >~ few M_D: the bridge invariant of 4D gravity is replaced by a collider-scale threshold.  Allowed today:
                             M_D >~ 5-11 TeV (n = 7 to 2, LHC), semiclassical black holes excluded below ~8-11 TeV; n >= 4 escapes the strong
                             astrophysical bounds.  A 100 TeV collider reaches M_D of order 15-20 TeV (scaling estimate).  But in D dimensions
                             R_S grows only as M^{1/(n+1)}: at M_BH ~ 5 M_D, R_S M_D ~ 1 (n = 6), so collider holes are objects at the
                             fundamental scale, with no semiclassical interior; R_S = 3/M_D needs M_BH ~ 1e3 M_D (n = 4) to 1e4 M_D (n = 6),
                             i.e. 10-100 PeV collisions (cosmic rays reach ~0.75 PeV).
"""

from __future__ import annotations

import numpy as np
from math import gamma
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from .stable_shells import Shell
from . import nariai as NA

G_SI = 6.67430e-11
C_SI = 2.99792458e8
M_SUN = 1.98847e30
YEAR = 3.15576e7
GEV_J = 1.602176634e-10
HBARC_GEV_M = 1.97327e-16


# --------------------------------------------------------------------------- Group B


def accretion_requirement(M_sun=10.0, efficiency=0.1) -> dict:
    """Near-Nariai partner density rho_f = c^6/(72 pi G^3 M^2) [kg/m^3 equivalent] versus an Eddington-limited flow at the horizon
    (rho ~ Mdot_Edd / (4 pi r_s^2 c)) and merger debris (0.1 M_sun within 3 r_s)."""
    M = M_sun * M_SUN
    rho_f = C_SI**6 / (72 * np.pi * G_SI**3 * M**2)            # Lambda c^2/(8 pi G), Nariai Lambda = c^4/(9 G^2 M^2)
    r_s = 2 * G_SI * M / C_SI**2
    sigma_T, m_p = 6.6524587e-29, 1.67262192e-27
    Mdot_edd = 4 * np.pi * G_SI * M * m_p / (sigma_T * C_SI * efficiency)
    rho_edd = Mdot_edd / (4 * np.pi * r_s**2 * C_SI)
    rho_merger = 0.1 * M_SUN / (4 / 3 * np.pi * (3 * r_s) ** 3)
    return dict(M_sun=M_sun, partner_vac_GeV=float(NA.partner_vacuum(M)), rho_f_kg_m3=float(rho_f), rho_eddington_kg_m3=float(rho_edd),
                log10_eddington_shortfall=float(np.log10(rho_f / rho_edd)), rho_merger_kg_m3=float(rho_merger),
                merger_over_required=float(rho_merger / rho_f))


def qcd_bag_partner(B_quarter_GeV=0.145) -> dict:
    """The black-hole mass whose Nariai partner vacuum equals the QCD bag constant."""
    M_kg = NA.partner_mass(B_quarter_GeV)
    return dict(B_quarter_GeV=B_quarter_GeV, M_kg=float(M_kg), M_sun=float(M_kg / M_SUN))


def bag_star(n=120) -> dict:
    """TOV for the MIT-bag equation of state p = (rho - 4B)/3 (units B = 1, G = c = 1).  Returns the maximum-mass star: compactness,
    radius over the bag vacuum's Hubble length l_f = sqrt(3/(8 pi B)), and the vacuum share of the energy."""
    def star(rho_c):
        def rhs(r, y):
            m, p = y
            rho = 3 * p + 4.0
            return [4 * np.pi * r**2 * rho, -(rho + p) * (m + 4 * np.pi * r**3 * p) / (r * (r - 2 * m))]
        ev = lambda r, y: y[1]
        ev.terminal = True
        r0 = 1e-6
        p0 = (rho_c - 4) / 3
        s = solve_ivp(rhs, (r0, 10.0), [4 / 3 * np.pi * r0**3 * rho_c, p0], events=ev, rtol=1e-10, atol=1e-14, max_step=1e-3)
        R = s.t_events[0][0]
        M = s.y_events[0][0][0]
        return M, R
    best = None
    for rho_c in np.geomspace(4.05, 400, n):
        M, R = star(rho_c)
        if best is None or M > best[0]:
            best = (M, R, rho_c)
    M, R, rho_c = best
    l_f = np.sqrt(3 / (8 * np.pi * 1.0))
    return dict(M_max=float(M), R=float(R), compactness=float(2 * M / R), R_over_l_f=float(R / l_f), rho_c_over_B=float(rho_c),
                vacuum_share=float(4 / 3 * np.pi * R**3 * 1.0 / M))


def oy_channel(M_in, chi=1.0, k_w=0.5, n=80) -> dict:
    """Bubble of false vacuum (Hubble rate chi, wall 4 pi sigma = k_w) around a hole of mass M_in, exterior mass M_out = M_in + dM.
    Thermal cost 8 pi M_in dM (Boltzmann with T_H = 1/(8 pi M_in)), tunnelling B(M_out) from the bound branch to the hidden inflating branch.
    Returns the minimum over dM of the total exponent."""
    if M_in >= 1 / (3 * np.sqrt(3) * chi):
        return dict(M_in=M_in, exists=False)
    best = None
    for dM in np.geomspace(1e-4, 2.0, n) / chi:
        sh = Shell(chi=chi, M=M_in + dM, k_w=k_w, M_in=M_in)
        try:
            regs = sh.regions(2.05 * M_in + 1e-4, 40.0 / chi)
        except Exception:
            continue
        if len(regs) < 3 or not regs[-1]["allowed"] or regs[-2]["allowed"] or regs[-1]["beta_out"] >= 0:
            continue
        B = sh.B(*regs[-2]["R"])
        tot = 8 * np.pi * M_in * dM + B
        if best is None or tot < best["Sigma"]:
            best = dict(M_in=M_in, dM=float(dM), thermal=float(8 * np.pi * M_in * dM), B=float(B), Sigma=float(tot), barrier=regs[-2]["R"])
    return dict(exists=best is not None, **(best or {}))


def oy_scan(chi=1.0, k_w=0.5, fracs=(0.01, 0.03, 0.06, 0.1, 0.15)) -> list:
    """Sigma_hat versus the hole mass (in units of the false-vacuum Hubble length 1/chi; Nariai bound 0.192/chi)."""
    return [oy_channel(f / chi, chi, k_w) for f in fracs]


def oy_expected_daughters(Sigma_hat, M_in_hat, L_over_lP) -> float:
    """Expected daughters while an evaporating hole passes M_in: attempts ~ (5120 pi M_in^3)(1/M_in) ~ 1.6e4 M_in^2, rate e^{-Sigma}."""
    M_in = M_in_hat * L_over_lP
    return float(1.6e4 * M_in**2 * np.exp(-Sigma_hat * L_over_lP**2))


# --------------------------------------------------------------------------- Group C: Euclidean wineglass wormholes (8 pi G = 1)


def _wg_roots(q, Lam):
    """Neck a_0 and maximum a_m: roots of 1 - (a^2/3)(Lam + q^2/a^6) = 0, i.e. Lam a^6 - 3 a^4 + q^2 = 0 in x = a^2: Lam x^3 - 3x^2 + q^2."""
    roots = np.roots([Lam, -3.0, 0.0, q**2])
    xs = sorted(r.real for r in roots if abs(r.imag) < 1e-12 and r.real > 0)
    if len(xs) < 2:
        return None
    return np.sqrt(xs[0]), np.sqrt(xs[1])


def wineglass_action(q, Lam=0.1) -> dict:
    """Half-wormhole action from the neck to the maximum: S = 2 pi^2 int (-6a + 2 Lam a^3 + 2 q^2/a^3) da / a'.  q = 0: -12 pi^2 / Lam."""
    if q == 0:
        return dict(q=0.0, a0=0.0, am=float(np.sqrt(3 / Lam)), S=float(-12 * np.pi**2 / Lam))
    rt = _wg_roots(q, Lam)
    if rt is None:
        return dict(q=q, exists=False)
    a0, am = rt

    def ap(a):
        return np.sqrt(max(1 - a**2 / 3 * (Lam + q**2 / a**6), 0.0))

    def f(t):                                   # substitution a = a0 + (am - a0) sin^2(t) removes both endpoint singularities
        a = a0 + (am - a0) * np.sin(t) ** 2
        da = 2 * (am - a0) * np.sin(t) * np.cos(t)
        return (-6 * a + 2 * Lam * a**3 + 2 * q**2 / a**3) / max(ap(a), 1e-300) * da
    val = quad(f, 0, np.pi / 2, limit=400, epsabs=0, epsrel=1e-10)[0]
    return dict(q=q, a0=float(a0), am=float(am), S=float(2 * np.pi**2 * val))


def wineglass_family(Lam=0.1, qs=None) -> dict:
    """Action excess Delta S = S(q) - S(0) against the throat radius a_0; fit Delta S = c a_0^2 for small throats."""
    qs = np.geomspace(1e-4, 0.5, 30) if qs is None else qs
    S0 = -12 * np.pi**2 / Lam
    rows = [wineglass_action(q, Lam) for q in qs]
    rows = [dict(r, dS=r["S"] - S0) for r in rows if r.get("S") is not None]
    small = [r for r in rows if r["a0"] < 0.2 * np.sqrt(3 / Lam)]
    c = float(np.median([r["dS"] / r["a0"] ** 2 for r in small])) if small else None
    return dict(Lam=Lam, S0=float(S0), rows=rows, c_fit=c)


# --------------------------------------------------------------------------- Group D


def fmm_junction(M, l=1.0) -> dict:
    """Spacelike junction of the Schwarzschild interior (f_S = 1 - 2M/r) with de Sitter (f_dS = 1 - r^2/l^2), both T-regions, at
    f_S = f_dS: r_0 = (2 M l^2)^{1/3}.  Needs l < r_0 < 2M (M > l/2).  K_theta_theta is continuous; the jump of K_t^t gives a layer with
    tangential stress S_theta^theta = [K_t^t]/(8 pi)."""
    r0 = (2 * M * l**2) ** (1 / 3)
    ok = l < r0 < 2 * M
    f0 = 1 - 2 * M / r0
    fpS, fpD = 2 * M / r0**2, -2 * r0 / l**2
    jump = -(fpD - fpS) / (2 * np.sqrt(abs(f0))) if ok else None
    return dict(M=M, l=l, r0=float(r0), valid=bool(ok), r0_over_l=float(r0 / l), kretschmann_S=float(48 * M**2 / r0**6),
                kretschmann_dS=float(24 / l**4), K_tt_jump=None if jump is None else float(jump),
                S_thth=None if jump is None else float(jump / (8 * np.pi)), proper_time_to_r0=float(np.pi * M))


def add_schwarzschild_radius(M_BH_TeV, M_D_TeV, n) -> float:
    """Dimopoulos-Landsberg: R_S = (1/(sqrt(pi) M_D)) [8 Gamma((n+3)/2)/(n+2) * M_BH/M_D]^{1/(n+1)} (metres)."""
    x = 8 * gamma((n + 3) / 2) / (n + 2) * M_BH_TeV / M_D_TeV
    return float(HBARC_GEV_M / (np.sqrt(np.pi) * M_D_TeV * 1e3) * x ** (1 / (n + 1)))


def add_black_holes(n=6, M_D_TeV=8.0, M_BH_TeV=30.0) -> dict:
    """Geometric cross-section pi R_S^2 (pb), and the bounds that apply (LHC monojet, CMS black-hole searches; astrophysical bounds weak for
    n >= 4)."""
    R = add_schwarzschild_radius(M_BH_TeV, M_D_TeV, n)
    sigma_pb = np.pi * R**2 / 1e-40
    bounds = {2: 10.7, 3: 8.0, 4: 7.0, 5: 6.2, 6: 5.6, 7: 5.2}           # CMS 2021 monojet (approximate interpolation between the quoted ends)
    return dict(n=n, M_D_TeV=M_D_TeV, M_BH_TeV=M_BH_TeV, R_S_m=R, sigma_parton_pb=float(sigma_pb),
                M_D_lhc_bound_TeV=bounds.get(n), allowed_by_lhc=bool(M_D_TeV >= bounds.get(n, 11.2)),
                astro_bounds_weak=bool(n >= 4), threshold_over_MD=float(M_BH_TeV / M_D_TeV))


def add_semiclassical_mass(n=6, k=3.0) -> float:
    """M_BH / M_D at which R_S = k / M_D (a black hole k fundamental lengths across):
    x = (k sqrt(pi))^{n+1} (n+2) / (8 Gamma((n+3)/2))."""
    return float((k * np.sqrt(np.pi)) ** (n + 1) * (n + 2) / (8 * gamma((n + 3) / 2)))


def add_rs_over_fundamental(M_BH_TeV, M_D_TeV, n) -> float:
    """R_S M_D: the black hole's size in fundamental lengths (Dimopoulos-Landsberg)."""
    x = 8 * gamma((n + 3) / 2) / (n + 2) * M_BH_TeV / M_D_TeV
    return float(x ** (1 / (n + 1)) / np.sqrt(np.pi))


def ledger6j() -> list:
    """The fifteen candidates: what changed, what was computed, verdict."""
    return [
        ("1 Non-uniform matter", "peak density >= areal mean at the first marginal sphere (Prop. 6j-1)", "no gain: non-uniformity raises the peak"),
        ("2 Non-spherical data", "Riemannian Penrose inequality: k_A >= 1 for time-symmetric data with rho >= 0", "same invariant; dynamical case open"),
        ("3 Field-dependent G", "strong-gravity core costs compactness (3/4) ln^2 F (Prop. 6j-4)", "horizon-class for F < 0.32; no gain"),
        ("4 Higher curvature", "f(R) = scalar-tensor (Prop. 6j-4); Gauss-Bonnet (Round 6i); Weyl^2 has ghosts", "no gain"),
        ("5 Finite-lived pocket", "sub-Hubble pocket expands by ~R/l_f before collapse (Prop. 6j-5)", "no measurable spacetime creation"),
        ("6 Accretion-driven vacuum", "partner density needs merger-class inflow; NEC escape still needed (Prop. 6j-6)", "natural state preparation only"),
        ("7 Horizonless precursor", "max-compact bag star: vacuum region 0.43 l_f, collapses to a black hole (Prop.-N 6j-7)", "no daughter classically"),
        ("8 Evaporating-hole catalysis", "full two-step exponent Sigma_hat (l_f/l_P)^2 (Prop.-N 6j-8)", "O(1) only near the Planck scale"),
        ("9 Cosmological critical state", "trigger ~10% of patch vacuum energy (R137); patches decided at t ~ their Hubble time", "early universe only"),
        ("10 Quantum transition without a throat", "throat-free saddles are the disconnected ones (Prop.-N 6j-9)", "no connected throat-free channel found"),
        ("11 Wormhole / no-boundary", "Delta S grows with the throat; S -> no-boundary as q -> 0 (Prop.-N 6j-9)", "creation not attached to us"),
        ("12 Holographic baby universe", "SSV baby universe lives behind a black-hole horizon (a bridge)", "consistent with Theorem 6i-11; toy model"),
        ("13 More quantum degrees of freedom", "minimal-action path can only drop, but the throat-area cost remains", "open; no barrier removal shown"),
        ("14 Dimensional crossover", "bridges at M_BH ~ few M_D (M_D ~ 6-20 TeV allowed, n >= 4); but R_S M_D ~ 1 there (Prop. 6j-11)", "CONDITIONAL door: 100 TeV makes Planckian holes"),
        ("15 Quantum-corrected interiors", "FMM junction r_0 = (2Ml^2)^{1/3}: every hole holds an inflating region (Prop. 6j-10)", "CONDITIONAL door: daughter = any black hole"),
    ]
