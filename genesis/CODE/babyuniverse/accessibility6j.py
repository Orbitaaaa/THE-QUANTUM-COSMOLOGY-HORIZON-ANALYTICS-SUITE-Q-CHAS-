"""Round 6j, Part I and Group A: what the bridge invariant bounds (geometry, stored energy, supplied energy), and the classical loopholes
that change the gravitational assumptions (non-uniform matter, non-spherical data, a field-dependent Planck mass, higher curvature, a
finite-lived pocket).

Part I (ROADMAP_ROUND6J.md, section 1):
* ``compactness_identity``  E sqrt(rho_mean) = sqrt(3/(32 pi)) c^6/G^{3/2} k^{3/2}, k = 2GE/(c^4 R): an identity.  Theorem 6i-11 uses it at
                            k = 1, the value it takes at a marginally trapped sphere.
* ``misner_sharp_identity`` [Theorem, symbolic] d_u m = r^2 e^{-2s} (G_uv r_u - G_uu r_v) (and u <-> v) for ds^2 = -e^{2s} du dv + r^2 dOmega^2.
                            In an orthonormal frame, along a radial slice, dm/dr = 4 pi r^2 (rho + F r_t / r_n), F the outward energy flux.
                            So the density in Theorem 6i-11 is the areal-volume mean of the local energy density in the frame with no radial
                            energy flux (comoving for fluids); on a time-symmetric slice it is the ordinary energy density.
* ``which_density``         [Proposition 6j-1] At the first marginally trapped sphere (radius R) of regular data, some region has local
                            zero-flux energy density >= rho_mean = 3 c^4 / (8 pi G R^2); non-uniformity raises the peak, and negative-energy
                            components raise the positive energy needed (m = int rho_+ - int rho_- = R c^2 / 2G).
* ``three_bounds``          [Proposition 6j-2] Three different statements: the geometric identity (exact), a bound on stored energy (local
                            matter energy in the ball >= Misner-Sharp mass >= R c^4/2G for rho >= 0, since proper volume >= areal volume),
                            and the externally supplied energy E_sup = E_stored - E_natural, which is not bounded by the invariant at all.
* ``precursor_window``      [Proposition 6j-3] Supplied energy E_sup << E_stored needs a natural precursor within delta = E_sup / (M c^2) of the
                            threshold that is (a) dynamically stable (near-Nariai states are not: they decide within tau = 3.2 r_N / c),
                            (b) old enough to be found (Hawking lifetime > 13.8 Gyr needs M > 1.7e11 kg) and (c) not pushed over by natural
                            perturbations.  For E_sup = 1e20 J and M >= 1.7e11 kg, delta <= 7e-9.
* ``qei_negative_energy``   [Assessment, Ford-Roman scaling] A machine's negative energy in a region of size R over a time ~R/c is |E_-| <~
                            hbar c / R, against the tipping requirement W_th c^4/G ~ 1.85e-3 R c^4/G: short by (R / l_P)^2 / 540.
* ``penrose_form``          [Theorem, Huisken-Ilmanen 2001; Bray 2001] Non-spherical, time-symmetric data with rho >= 0: M_ADM >= sqrt(A/16 pi)
                            c^2/G for the outermost minimal surface, i.e. k_A >= 1 and the same invariant with the area radius.  With rho < 0
                            the inequality fails (Ellis wormhole: a throat with M_ADM = 0), but formation from regular data still needs
                            positive energy >= R c^4/2G inside the first marginal sphere (Proposition 6j-1).
Group A:
* ``strong_gravity_core``   [Proposition 6j-4] F(phi) R / 16 pi G: a core where F = M_P,eff^2 / M_P^2 << 1 lowers the requirement measured
                            outside by F^{3/4}, but needs the canonical excursion Delta chi = sqrt(3/(16 pi)) |ln F| M_P; by the capacitance
                            (Dirichlet) bound its gradient energy is E >= 2 pi L Delta chi^2, so the core's own compactness is >= (3/4) ln^2 F:
                            horizon-class for F < 0.32.  f(R) gravity is the same theory (F = f'(R)).
* ``pocket_expansion``      [Proposition 6j-5] A classical pocket must stay sub-Hubble (R < l_f; otherwise Theorem 6h-1/6i-11 apply) and
                            collapses within ~R/c: its interior expands by a fraction ~R/l_f < 1 before it is gone (TeV vacuum, R = 1e-18 m:
                            3e-16).
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.integrate import cumulative_trapezoid

G_SI = 6.67430e-11
C_SI = 2.99792458e8
HBAR_SI = 1.054571817e-34
L_P = 1.616255e-35
M_P_KG = 2.176434e-8
E_P_GEV = 1.220890e19
HBARC_GEV_M = 1.97327e-16
M_SUN = 1.98847e30
YEAR = 3.15576e7


# --------------------------------------------------------------------------- Part I


def compactness_identity(E_J, R_m) -> dict:
    """k = 2GE/(c^4 R) and E sqrt(rho_mean) with rho_mean = 3E/(4 pi R^3), against sqrt(3/32 pi) c^6/G^{3/2} k^{3/2}."""
    k = 2 * G_SI * E_J / (C_SI**4 * R_m)
    lhs = E_J * np.sqrt(3 * E_J / (4 * np.pi * R_m**3))
    rhs = np.sqrt(3 / (32 * np.pi)) * C_SI**6 / G_SI**1.5 * k**1.5
    return dict(k=float(k), lhs=float(lhs), rhs=float(rhs), invariant_k1=float(np.sqrt(3 / (32 * np.pi)) * C_SI**6 / G_SI**1.5))


@lru_cache(maxsize=1)
def misner_sharp_identity() -> dict:
    """Verify d_u m = r^2 e^{-2s}(G_uv r_u - G_uu r_v), d_v m = r^2 e^{-2s}(G_uv r_v - G_vv r_u), and the orthonormal-frame form
    e_n m = 4 pi r^2 (rho r_n + F r_t) with rho = e^{-2s}(T_uu + 2T_uv + T_vv), F = -e^{-2s}(T_vv - T_uu)."""
    u, v, th, ph = sp.symbols("u v theta phi")
    X = [u, v, th, ph]
    r, s = sp.Function("r")(u, v), sp.Function("s")(u, v)
    g = sp.zeros(4)
    g[0, 1] = g[1, 0] = -sp.exp(2 * s) / 2
    g[2, 2], g[3, 3] = r**2, r**2 * sp.sin(th) ** 2
    gi, n = g.inv(), 4
    Gam = [[[sp.simplify(sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(n)) / 2)
             for c in range(n)] for b in range(n)] for a in range(n)]

    def riem(a, b, c, d):
        return sp.diff(Gam[a][d][b], X[c]) - sp.diff(Gam[a][c][b], X[d]) + sum(
            Gam[a][c][k] * Gam[k][d][b] - Gam[a][d][k] * Gam[k][c][b] for k in range(n))

    Ric = sp.Matrix(n, n, lambda b, d: sp.simplify(sum(riem(a, b, a, d) for a in range(n))))
    Rs = sp.simplify(sum(gi[a, b] * Ric[a, b] for a in range(n) for b in range(n)))
    G = sp.simplify(Ric - Rs * g / 2)
    m = r / 2 * (1 + 4 * sp.exp(-2 * s) * sp.diff(r, u) * sp.diff(r, v))
    ok_u = sp.simplify(sp.diff(m, u) - r**2 * sp.exp(-2 * s) * (G[0, 1] * sp.diff(r, u) - G[0, 0] * sp.diff(r, v))) == 0
    ok_v = sp.simplify(sp.diff(m, v) - r**2 * sp.exp(-2 * s) * (G[0, 1] * sp.diff(r, v) - G[1, 1] * sp.diff(r, u))) == 0
    Tuu, Tuv, Tvv, ru, rv, S, R_ = sp.symbols("Tuu Tuv Tvv r_u r_v s r")
    mu = 8 * sp.pi * R_**2 * sp.exp(-2 * S) * (Tuv * ru - Tuu * rv)
    mv = 8 * sp.pi * R_**2 * sp.exp(-2 * S) * (Tuv * rv - Tvv * ru)
    rho = sp.exp(-2 * S) * (Tuu + 2 * Tuv + Tvv)
    F = -sp.exp(-2 * S) * (Tvv - Tuu)
    rt, rn = sp.exp(-S) * (ru + rv), sp.exp(-S) * (rv - ru)
    ok_frame = sp.simplify(sp.exp(-S) * (mv - mu) - 4 * sp.pi * R_**2 * (rho * rn + F * rt)) == 0
    return dict(d_u=bool(ok_u), d_v=bool(ok_v), orthonormal=bool(ok_frame))


def _profiles():
    """Time-symmetric radial energy-density shapes (arbitrary units; scaled by an amplitude until a marginal sphere appears)."""
    return {
        "uniform ball": lambda r: (r < 1.0).astype(float),
        "centrally peaked (Gaussian)": lambda r: np.exp(-(r / 0.3) ** 2),
        "hollow shell": lambda r: np.exp(-((r - 1.0) / 0.1) ** 2),
        "vacuum ball + negative dust shell": lambda r: (r < 1.0) - 0.6 * np.exp(-((r - 1.0) / 0.05) ** 2) / 0.05 * 0.1,
        "negative core + positive shell": lambda r: -0.5 * (r < 0.5) + 3 * np.exp(-((r - 1.0) / 0.1) ** 2),
    }


def which_density(n=20001, r_max=3.0) -> list:
    """For each profile, the smallest amplitude at which some sphere has 2m = r (time-symmetric data, m = int 4 pi r^2 rho dr, G = c = 1),
    and at that moment the peak positive density over rho_mean = 3 m(R) / (4 pi R^3), plus the positive energy needed over R/2."""
    r = np.linspace(1e-6, r_max, n)
    rows = []
    for name, f in _profiles().items():
        rho1 = f(r)
        m1 = cumulative_trapezoid(4 * np.pi * r**2 * rho1, r, initial=0.0)
        ratio = 2 * m1 / r
        k = int(np.argmax(ratio))
        A = 1.0 / ratio[k]                     # amplitude making the first marginal sphere at r[k]
        R = r[k]
        rho = A * rho1
        rho_mean = 3 * (A * m1[k]) / (4 * np.pi * R**3)
        inside = r <= R
        E_plus = np.trapezoid(4 * np.pi * r[inside] ** 2 * np.maximum(rho[inside], 0), r[inside])
        rows.append(dict(profile=name, R=float(R), peak_over_mean=float(rho[inside].max() / rho_mean),
                         positive_energy_over_R_half=float(E_plus / (R / 2))))
    return rows


def three_bounds(R_m, E_natural_J=0.0) -> dict:
    """Geometric (exact at k = 1), stored (>= Misner-Sharp, rho >= 0) and supplied (stored minus what nature provides)."""
    E_ms = R_m * C_SI**4 / (2 * G_SI)
    rho_mean = 3 * C_SI**4 / (8 * np.pi * G_SI * R_m**2)
    return dict(R_m=R_m, geometric_E_J=float(E_ms), geometric_rho_J_m3=float(rho_mean), stored_min_J=float(E_ms),
                supplied_min_J=float(max(E_ms - E_natural_J, 0.0)),
                note="supplied energy is bounded only by what the precursor lacks; the invariant bounds stored energy")


def hawking_lifetime_s(M_kg, g_eff=1.0) -> float:
    """t ~ 5120 pi G^2 M^3 / (hbar c^4) / g_eff (photons and gravitons only for g_eff = 1)."""
    return float(5120 * np.pi * G_SI**2 * M_kg**3 / (HBAR_SI * C_SI**4) / g_eff)


def precursor_window(E_sup_J=1e20, M_kg=None) -> dict:
    """The tolerance delta = E_sup / (M c^2) a natural precursor must sit within, for precursors that outlive the universe, and the
    near-Nariai decision time."""
    t_univ = 13.8e9 * YEAR
    M_min = (t_univ * HBAR_SI * C_SI**4 / (5120 * np.pi * G_SI**2)) ** (1 / 3)
    M = M_min if M_kg is None else M_kg
    r_N = 3 * G_SI * M / C_SI**2                      # M_N = r_N c^2 / 3G
    return dict(M_kg=float(M), M_min_survivor_kg=float(M_min), delta=float(E_sup_J / (M * C_SI**2)),
                nariai_decision_time_s=float(3.2 * r_N / C_SI), hawking_lifetime_s=hawking_lifetime_s(M))


def qei_negative_energy(R_m) -> dict:
    """Ford-Roman scaling for the negative energy a quantum field can hold in a region R over a time R/c: |E_-| ~ hbar c / R; the throat
    tipping law needs W_th c^4 / G with W_th = 1.85e-3 R."""
    E_neg = HBAR_SI * C_SI / R_m
    W_need = 1.85e-3 * R_m * C_SI**4 / G_SI
    return dict(R_m=R_m, E_neg_max_J=float(E_neg), E_needed_J=float(W_need), log10_shortfall=float(np.log10(W_need / E_neg)),
                shortfall_formula=float((R_m / L_P) ** 2 * 1.85e-3))


def penrose_form(A_m2, M_kg) -> dict:
    """Riemannian Penrose inequality M >= sqrt(A / 16 pi) c^2 / G (time-symmetric, rho >= 0): its compactness k_A = 2GM/(c^2 r_A) >= 1."""
    r_A = np.sqrt(A_m2 / (4 * np.pi))
    k_A = 2 * G_SI * M_kg / (C_SI**2 * r_A)
    return dict(r_A_m=float(r_A), k_A=float(k_A), satisfies=bool(k_A >= 1 - 1e-12))


# --------------------------------------------------------------------------- Group A


def strong_gravity_core(F, L_m=1.0) -> dict:
    """Core of radius L (Einstein frame) where F = (M_P,eff / M_P)^2.  Delta chi = sqrt(3/(16 pi)) |ln F| M_P; Dirichlet (capacitance)
    minimum of the gradient energy with chi fixed on the core and chi -> chi_inf at infinity: E = 2 pi L Delta chi^2 (G = c = hbar = 1).
    Compactness of that energy: 2E/L = (3/4) ln^2 F.  Gain in the invariant measured outside: F^{3/4}."""
    lnF = abs(np.log(F))
    dchi = np.sqrt(3 / (16 * np.pi)) * lnF
    k_core = 4 * np.pi * dchi**2                       # 2 E / L with E = 2 pi L dchi^2
    return dict(F=F, M_P_eff_GeV=float(np.sqrt(F) * E_P_GEV), delta_chi_over_MP=float(dchi), core_compactness=float(k_core),
                horizon_class=bool(k_core >= 1), invariant_gain=float(F**0.75),
                E_core_J=float(2 * np.pi * L_m / L_P * dchi**2 * M_P_KG * C_SI**2))


def strong_gravity_threshold() -> float:
    """F below which the core is horizon-class: (3/4) ln^2 F = 1."""
    return float(np.exp(-np.sqrt(4 / 3)))


def dirichlet_check(n=4000) -> dict:
    """Numerical check of the capacitance minimum: the harmonic profile chi = Delta L / r outside L has gradient energy 2 pi L Delta^2;
    a linear ramp over [L, 2L] (then constant) has more."""
    L, D = 1.0, 1.0
    r = np.linspace(L, 400 * L, 400000)
    harm = 0.5 * np.trapezoid(4 * np.pi * r**2 * (D * L / r**2) ** 2, r) + 0.5 * 4 * np.pi * L * D**2 * (L / r[-1])
    r2 = np.linspace(L, 2 * L, n)
    ramp = 0.5 * np.trapezoid(4 * np.pi * r2**2 * (D / L) ** 2, r2)
    return dict(harmonic=float(harm), formula=float(2 * np.pi * L * D**2), linear_ramp=float(ramp))


def pocket_expansion(vac_scale_GeV, R_m) -> dict:
    """Sub-Hubble pocket: l_f = sqrt(3/(8 pi G rho)) c; collapse time ~ R/c; interior expansion fraction ~ H R / c = R / l_f."""
    rho_J = vac_scale_GeV**4 / HBARC_GEV_M**3 * 1.602176634e-10
    l_f = C_SI * np.sqrt(3 / (8 * np.pi * G_SI * rho_J / C_SI**2))
    return dict(vac_scale_GeV=vac_scale_GeV, R_m=R_m, l_f_m=float(l_f), expansion_fraction=float(R_m / l_f),
                anti_trapped=bool(R_m >= l_f))
