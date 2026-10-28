"""Round 6f, problems 4 and 5: does nature support the route, and could it be made technologically accessible?

Nothing here is a new physical prediction.  The route's requirements were computed in Planck units (Rounds 6d-6f) with r_N as the only
scale; this module converts them to SI units and sets them against published bounds and records (LITERATURE.md, Round 6f entries).

Requirements of the route (G = c = hbar = 1, r_N = Nariai radius of the false vacuum):
    rho_f = 1 / (8 pi r_N^2)                   Nariai: Lambda = 8 pi rho_f = 1 / r_N^2
    M_N = r_N / 3                              mass of the Nariai throat; the route's ball has ADM mass 0.516 r_N (8.93 at r_N = 17.32)
    K_lo <= M^3 r_N <= K_hi                    condensate scale M: tipping and stability window (Proposition-N 6f-6, universal in r_N)
    |n| r_N^2 ~ 1.5e-3, lam_f r_N^3 ~ 4e-4     condensate charge to tip; tilt for the self-generated NEC phase (Round 6e)

Theorem 6f-7 (the invariant).  Any route that ends in a near-Nariai throat of the false vacuum must place at least the Nariai mass
M_N = r_N / 3 at energy density rho_f = 1 / (8 pi r_N^2) inside a region of size ~ r_N.  Hence
    E_req * sqrt(rho_req) >= (1/3) (8 pi)^{-1/2} E_P sqrt(rho_P) ~ 8.8e64 J^{3/2} m^{-3/2},
the same for every vacuum scale: lowering rho_f raises the mass needed as rho_f^{-1/2}.  This is the black-hole-formation invariant
M sqrt(rho) ~ (c^2/G)^{3/2}; it is a statement about gravity, not about the condensate.
"""

from __future__ import annotations

import numpy as np

# CODATA 2018 / IAU
L_P = 1.616255e-35            # m
M_P = 2.176434e-8             # kg
T_P = 5.391247e-44            # s
C_LIGHT = 2.99792458e8        # m / s
E_P = M_P * C_LIGHT**2        # J  (1.956e9)
E_P_GEV = 1.220890e19         # GeV
RHO_P = E_P / L_P**3          # J / m^3 (4.63e113)
GEV_J = 1.602176634e-10       # J
M_SUN = 1.98847e30            # kg
MPC = 3.0857e22               # m
AU = 1.495978707e11           # m
H0_SI = 67.4e3 / MPC          # s^-1 (Planck 2018)
H0_P = H0_SI * T_P            # Planck units (1.18e-61)

# Route requirements (Planck units, per unit r_N); see the module docstring
M_BALL_OVER_RN = 8.93 / np.sqrt(300.0)        # route data: ADM mass 8.93 at r_N = 17.32 (l_H = 30)
N_TIP = 1.5e-3                                 # |n_th| r_N^2 (Round 6e, R107-R109)
LAM_TIP = 4e-4                                 # lam_th r_N^3, self-generated NEC phase (Round 6e)
WINDOW_ALPHA1 = (1.9e-4, 1.7e-2)               # analytic Proposition 6d-2 window for M^3 r_N at alpha = 1 (overridden by the map)


def physical_scales(r_N: float) -> dict:
    """The route's scales in SI units for a Nariai radius r_N (Planck lengths)."""
    rho = 1 / (8 * np.pi * r_N**2)
    return dict(r_N=r_N, r_N_m=r_N * L_P, rho_f_J_m3=rho * RHO_P, rho_f_kg_m3=rho * RHO_P / C_LIGHT**2,
                vac_scale_GeV=rho**0.25 * E_P_GEV, M_N_kg=r_N / 3 * M_P, M_N_Msun=r_N / 3 * M_P / M_SUN,
                M_ball_kg=M_BALL_OVER_RN * r_N * M_P, E_N_J=r_N / 3 * E_P, crossing_time_s=r_N * T_P)


def condensate_mass_range(r_N: float, window=WINDOW_ALPHA1) -> tuple:
    """Condensate scale M (GeV) allowed by the window K_lo <= M^3 r_N <= K_hi."""
    return tuple(float((k / r_N) ** (1 / 3) * E_P_GEV) for k in window)


def r_N_min_from_condensate_bound(M_max_GeV: float, window=WINDOW_ALPHA1) -> dict:
    """Smallest throat the route can tip when the condensate scale is bounded by M <= M_max (the window's lower edge)."""
    M = M_max_GeV / E_P_GEV
    r_N = window[0] / M**3
    return physical_scales(r_N)


def condensate_needed_for_mass(M_kg: float, window=WINDOW_ALPHA1) -> float:
    """Condensate scale (GeV) needed to tip a throat of Nariai mass M_kg."""
    r_N = 3 * M_kg / M_P
    return float((window[0] / r_N) ** (1 / 3) * E_P_GEV)


def jeans_bound_GeV(alpha: float = 1.0, t_age_over_hubble: float = 1.0) -> float:
    """Linear-theory bound on the condensate scale from its Jeans instability (rate 4 pi sqrt(alpha) M^3 in Planck units, Section 6f-6)
    being slower than the Hubble rate today: 4 pi sqrt(alpha) M^3 <= H0 / t_age_over_hubble."""
    return float((H0_P / (t_age_over_hubble * 4 * np.pi * np.sqrt(alpha))) ** (1 / 3) * E_P_GEV)


def requirement_invariant(c_E: float = 1 / 3) -> float:
    """Theorem 6f-7: E_req sqrt(rho_req) in J^{3/2} m^{-3/2}, with E_req = c_E r_N E_P and rho_req = rho_P / (8 pi r_N^2)."""
    return float(c_E * (8 * np.pi) ** -0.5 * E_P * np.sqrt(RHO_P))


def figure_of_merit(E_J: float, rho_J_m3: float) -> float:
    return float(E_J * np.sqrt(rho_J_m3))


# Achieved energies and the energy densities at which they were delivered.  Orders of magnitude; sources in LITERATURE.md (Round 6f).
CAPABILITIES = [
    dict(name="LHC Pb-Pb central collision", kind="technology", E_J=208 * 5.36e3 * GEV_J, rho_J_m3=12.3 * GEV_J / 1e-45,
         note="sqrt(s_NN) = 5.36 TeV, 208 nucleon pairs; Bjorken energy density 12.3 GeV/fm^3 (ALICE, 0-5%, tau = 1 fm/c)"),
    dict(name="Record laser focus (CoReLS 2021)", kind="technology", E_J=1e2, rho_J_m3=1.1e27 / C_LIGHT,
         note="1.1e23 W/cm^2 (Yoon et al., Optica 8, 630); pulse energy ~1e2 J (order of magnitude)"),
    dict(name="NIF ignition shot (April 2025)", kind="technology", E_J=8.6e6, rho_J_m3=5e16,
         note="8.6 MJ yield from 2.08 MJ on target; hot-spot pressure ~ few x 1e16 Pa (order of magnitude)"),
    dict(name="Largest thermonuclear explosion (50 Mt)", kind="technology", E_J=2.1e17, rho_J_m3=1e17,
         note="energy 2.1e17 J; core energy density ~ a T^4 at 1e8 K ~ 1e17 J/m^3 (order of magnitude)"),
    dict(name="Core-collapse supernova", kind="nature", E_J=3e46, rho_J_m3=2.4e34,
         note="binding energy ~3e46 J at nuclear saturation density 2.7e17 kg/m^3"),
    dict(name="Binary neutron-star merger", kind="nature", E_J=2.7 * M_SUN * C_LIGHT**2, rho_J_m3=1e18 * C_LIGHT**2,
         note="2.7 solar masses at ~1e18 kg/m^3"),
]


def accessibility_table(c_E: float = 1 / 3) -> list:
    """log10 of the shortfall of each capability against the invariant of Theorem 6f-7 (negative: the capability exceeds it)."""
    req = requirement_invariant(c_E)
    return [dict(c, fom=figure_of_merit(c["E_J"], c["rho_J_m3"]), log10_gap=float(np.log10(req / figure_of_merit(c["E_J"], c["rho_J_m3"]))))
            for c in CAPABILITIES]


def gap_report(window=WINDOW_ALPHA1, M_bounds_GeV=(100.0, 1e-2), lab_masses_kg=(1.0, 1e6), vac_scale_max_GeV=1.5e16) -> dict:
    """Problems 4 and 5 in numbers: the smallest throats the observationally allowed condensate can tip, the condensate needed for
    laboratory masses, the lightest throat at the largest vacuum scale considered, and the technology gap."""
    out = dict(invariant=requirement_invariant(), capabilities=accessibility_table())
    out["r_N_min"] = {f"M<={m:g} GeV": r_N_min_from_condensate_bound(m, window) for m in M_bounds_GeV}
    out["condensate_for_lab"] = {f"{m:g} kg": condensate_needed_for_mass(m, window) for m in lab_masses_kg}
    rho = (vac_scale_max_GeV / E_P_GEV) ** 4
    out["lightest_throat"] = physical_scales((8 * np.pi * rho) ** -0.5)
    tech = [c for c in out["capabilities"] if c["kind"] == "technology"]
    out["best_technology_gap_log10"] = float(min(c["log10_gap"] for c in tech))
    return out


# ------------------------------------------------------------------ Round 6g: can any machine do it?
G_SI = 6.67430e-11
HBAR = 1.054571817e-34
PLANCK_POWER = C_LIGHT**5 / G_SI          # 3.63e52 W
PLANCK_FORCE = C_LIGHT**4 / G_SI          # 1.21e44 N
WORLD_POWER = 1.9e13                       # W, world primary energy use (~600 EJ / yr)
WORLD_ENERGY_PER_YEAR = 6.0e20             # J


def assembly_bounds(compactness=2.0 / 3.0, alpha=1.0) -> dict:
    """Proposition 6g-2 (conditional engineering bound; corrected in Round 6h).  A region of radius R holding energy E with
    k = 2GE/(Rc^4) (the near-Nariai ball: 2/3) evolves in t_dyn = alpha R / c.  Two preparation strategies cost:
      (i)  delivering E from outside within t_dyn:  P >= (k / (2 alpha)) c^5 / G  (1.2e52 W for alpha = 1; 8.1e51 W for alpha = 1.5);
      (ii) holding it static while it is assembled: a stress across a cross-section ~ (k^2 / 4) c^4 / G (1.3e43 N).
    These are NOT the only possible histories (a system can evolve dynamically, or a small trigger can act on a pre-existing near-critical
    background), and the Planck power and force are not proven limits: the maximum-tension c^4/4G (3.0e43 N, above the holding stress)
    and maximum-luminosity c^5/4G conjectures have counterexamples in GR (Jowsey & Visser 2021; Cardoso et al. 2018).  The bound says
    what external assembly costs; it does not exclude every machine.  The classical no-go that holds for any strategy is Penrose /
    Farhi-Guth (``penrose_farhi_guth``)."""
    P = compactness / (2 * alpha) * PLANCK_POWER
    F = compactness**2 / 4 * PLANCK_FORCE
    return dict(power_W=float(P), force_N=float(F), max_luminosity_W=float(PLANCK_POWER / 4), max_tension_N=float(PLANCK_FORCE / 4),
                power_over_world=float(P / WORLD_POWER), log10_power_gap=float(np.log10(P / WORLD_POWER)),
                GW150914_peak_W=3.6e49)


def quantum_requirement() -> float:
    """Energy delivered by quanta of energy E_q can be localised no better than hbar c / E_q, so the invariant of Theorem 6f-7 becomes
    E_total * E_q >= hbar c * (E sqrt(rho))_req^(2/3)  (J^2)."""
    return float(HBAR * C_LIGHT * requirement_invariant() ** (2 / 3))


def single_quantum_energy() -> float:
    """The quantum energy at which one collision (E_total ~ 2 E_q) meets the requirement: E_q = sqrt(requirement / 2) (J)."""
    return float(np.sqrt(quantum_requirement() / 2))


def accelerator_size(E_q_J, gradient_V_per_m=None, dipole_T=None) -> float:
    """Length (m) of a linear accelerator with the given gradient, or circumference of a ring with the given dipole field (protons,
    E[GeV] = 0.3 B[T] R[m]; synchrotron losses ignored), reaching the quantum energy E_q."""
    E_eV = E_q_J / 1.602176634e-19
    if gradient_V_per_m is not None:
        return float(E_eV / gradient_V_per_m)
    return float(2 * np.pi * (E_eV / 1e9) / (0.3 * dipole_T))


MACHINES = [
    dict(name="LHC (whole stored beam, 7 TeV protons, focused ideally)", E_total_J=3.6e8, E_quantum_J=7e12 * 1.602176634e-19),
    dict(name="FCC-hh design (8 GJ beam, 50 TeV protons, focused ideally)", E_total_J=8.4e9, E_quantum_J=5e13 * 1.602176634e-19),
    dict(name="NIF laser (2.1 MJ, 351 nm photons)", E_total_J=2.1e6, E_quantum_J=3.5 * 1.602176634e-19),
    dict(name="50 Mt thermonuclear explosion (keV X-rays)", E_total_J=2.1e17, E_quantum_J=1e3 * 1.602176634e-19),
    dict(name="one year of world energy in 100 TeV quanta (invented, ideal focus)", E_total_J=WORLD_ENERGY_PER_YEAR, E_quantum_J=1e14 * 1.602176634e-19),
]


def machine_audit() -> dict:
    """Every machine against the quantum form of the invariant and against the assembly power; the ideal-focus rows assume the whole
    energy can be put into one quantum-limited spot (hbar c / E_q) within one light-crossing of it, which no device can approach."""
    need = quantum_requirement()
    rows = []
    for m in MACHINES:
        prod = m["E_total_J"] * m["E_quantum_J"]
        L = HBAR * C_LIGHT / m["E_quantum_J"]
        P = m["E_total_J"] * C_LIGHT / L
        rows.append(dict(m, product_J2=float(prod), log10_shortfall=float(np.log10(need / prod)), spot_m=float(L), power_if_delivered_W=float(P)))
    Eq = single_quantum_energy()
    return dict(requirement_J2=need, machines=rows, single_collision_quantum_J=Eq, single_collision_quantum_GeV=float(Eq / 1.602176634e-10),
                linac_100GV_per_m_m=accelerator_size(Eq, gradient_V_per_m=1e11), ring_16T_m=accelerator_size(Eq, dipole_T=16.0),
                assembly=assembly_bounds())


def tev_gravity_conditional(M_D_TeV=10.1) -> dict:
    """If gravity became strong at M_D (large extra dimensions), the Planck power and force would be replaced by M_D^2 c^4/hbar-class
    values, reachable in a single collision above M_D.  CMS excludes semiclassical black holes up to 10.1 TeV (JHEP 11 (2018) 042), so a
    collider would need E_cm well above that, and the whole route would have to be rebuilt in higher dimensions."""
    E = M_D_TeV * 1e12 * 1.602176634e-19
    return dict(M_D_TeV=M_D_TeV, power_scale_W=float(E**2 / HBAR), excluded_below_TeV=10.1, fcc_hh_E_cm_TeV=100.0)



# ------------------------------------------------------------------ Round 6h: what survives, and what a background would take
def penrose_farhi_guth() -> dict:
    """Theorem 6h-1 (Penrose 1965; Farhi-Guth 1987; non-spherical form via the time-reversed Penrose theorem).  In a globally hyperbolic
    spacetime with non-compact Cauchy surfaces obeying the null energy condition, a closed anti-trapped surface implies past null
    incompleteness.  An inflating daughter region contains closed anti-trapped surfaces.  So a daughter needs (i) NEC violation, (ii) a
    past singularity (the Big Bang: the cosmological route), or (iii) a quantum transition.  A machine of NEC-respecting matter acting on
    data with a regular past has none of these, whatever its power, duration or strategy (rapid, slow, dynamical, triggered)."""
    return dict(escapes=("NEC violation (the ghost condensate; quantum fields are limited by QEIs, Round 6c)",
                         "a past singularity (cosmology, Round 6g)", "quantum tunnelling (bounded by R42)"),
                applies_to_any_assembly_strategy=True)


def background_requirements(vac_scale_GeV: float) -> dict:
    """With NEC violation available, the route still needs a near-Nariai false-vacuum ball (normal, not anti-trapped, so Theorem 6h-1
    does not forbid it).  For a vacuum scale rho_f^{1/4}: its size r_N, mass M_N = r_N c^2 / (3G), the smallest quantum energy that can
    be localised in it (hbar c / r_N), how many such quanta carry M_N c^2, the light-crossing time, and a linac at 100 GV/m reaching that
    quantum energy."""
    rho = (vac_scale_GeV / E_P_GEV) ** 4
    s = physical_scales((8 * np.pi * rho) ** -0.5)
    Eq = HBAR * C_LIGHT / s["r_N_m"]
    return dict(vac_scale_GeV=float(vac_scale_GeV), r_N_m=s["r_N_m"], M_N_kg=s["M_N_kg"], E_N_J=s["E_N_J"], quantum_J=float(Eq),
                quantum_GeV=float(Eq / GEV_J), n_quanta=float(s["E_N_J"] / Eq), crossing_time_s=s["crossing_time_s"],
                linac_100GV_m=accelerator_size(Eq, gradient_V_per_m=1e11))


def route_ledger() -> list:
    """Routes A-D of Round 6h with what each still needs (see ROADMAP_ROUND6H.md)."""
    return [
        dict(route="A: trigger a near-critical background", status="works with a small positive-energy trigger only in an anti-trapped "
             "(cosmological) background; in a normal background the trigger must violate the NEC (Theorem 6h-1, tipping law R88)",
             needs="a pre-existing near-critical false-vacuum configuration; none is accessible today"),
        dict(route="B: black-hole-assisted (Oshita-Yokoyama)", status="nucleation only once the black hole's mass is below ~ the false "
             "vacuum's Hubble scale; quantum, suppressed by a Euclidean action", needs="a black hole of about the Nariai mass class at "
             "the end of its evaporation (primordial only; making one is itself a horizon-class assembly)"),
        dict(route="C: dynamical assembly", status="cannot evade Theorem 6h-1 with NEC matter; with the condensate it still has to "
             "produce the near-Nariai ball (mass r_N c^2 / 3G), e.g. by gravitational collapse of natural matter",
             needs=">= 2.5e9 solar masses if the condensate is the universal one allowed by observation (M <= 100 GeV)"),
        dict(route="D: different gravitational physics", status="TeV-scale gravity would put horizon-class objects in reach of colliders",
             needs="M_D above the CMS bound of 10.1 TeV and the whole route rebuilt in higher dimensions"),
    ]
