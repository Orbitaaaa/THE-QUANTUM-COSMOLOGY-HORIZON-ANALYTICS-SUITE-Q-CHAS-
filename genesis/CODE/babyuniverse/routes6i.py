"""Round 6i, routes 4 and 6, and the ledger of all seven routes (state, trigger, completion).  G = c = hbar = 1 unless stated.

Route 4 (gauged ghost condensation; Cheng, Luty, Mukohyama & Thaler 2006, hep-th/0603010):
* ``nec_charge_identity``   In a ghost condensate rho + p = 2 X P'(X) M^4 and the shift-charge density is J^0 = 2 P'(X) phidot, so every NEC-
                            violating region (P' < 0) is charged: |rho + p| = |J^0| phidot.  [Theorem, one line]
* ``coulomb_ratio``         Gauged, the charge sources an electric field.  A uniform NEC-violating ball of radius R carries Coulomb energy
                            density / |rho + p| = (P''|delta| / 9) g^2 M^2 R^2: the violation is screened on scales R >~ 3 / (g M sqrt(P''|delta|)).
* ``gauged_window``         [Proposition 6i-8] Tipping a throat of radius r_N needs M^3 r_N >= K_lo (Round 6f window) and an unscreened
                            violation, r_N <~ c / (g M); CLMT need M <~ min(1e12 GeV, g^2 1e15 GeV).  Together M >= (K_lo^2 / (c^2 M_*))^(1/3)
                            ~ 1e17 GeV for c = 30: the window is empty, by >~ 5 orders of magnitude in M.  Gauging does not rescue the universal
                            condensate.
* ``vacuum_dependent``      [Proposition 6i-9, conditional] A condensate whose scale depends on the vacuum, M(chi): large M_in in the false
                            vacuum, M_out <= 100 GeV in ours (an overall factor F(chi) -> 0 instead makes the strong-coupling scale M F^2 vanish
                            in our vacuum).  The effective theory stays valid across the wall only if the wall is thicker than 1/M_out, so
                            r_N >= 1/M_out = 2e-18 m: M_N >= 9e8 kg with a ~3e10 GeV false vacuum, 30 orders of magnitude below the universal-
                            condensate floor (2.5e9 M_sun).  Assembling that background still needs ~c^5/G power (Proposition 6g-2, conditional).
Route 6 (classical transitions; Easther, Giblin, Hui & Lim 2009; Johnson & Yang 2010, arXiv:1005.3506):
* ``dtr_region_c``          Dray-'t Hooft-Redmount relation f_A f_D = f_B f_C at a spherical collision: a pocket of any Lambda_C is kinematically
                            allowed (Johnson-Yang: in planar collisions the new vacuum cannot exceed the parent's).
* ``pocket_completion``     [Proposition 6i-10] Where f_out > 0, beta_out^2 = f_out + Rdot^2 > 0, so beta_out keeps its sign.  A pocket wall
                            launched outward from a collision in a normal region has beta_out > 0 and can never reach the hidden inflating branch
                            (beta_out < 0); with positive tension it cannot grow into the parent either.  It collapses (Theorem 6h-1 in thin-wall
                            form).  In an anti-trapped (inflating) parent collisions can make vacua up to the parent's: the cosmological route.
Synthesis:
* ``bridge_invariant``      [Theorem 6i-11] Every hidden daughter contains a sphere with 2Gm/(c^2 r) >= 1 (the area radius has a minimum between
                            the wall, where beta_out < 0, and our asymptotic region).  From regular data its first appearance is a black-hole
                            formation event with E sqrt(rho_mean) = sqrt(3/(32 pi)) c^6/G^{3/2} = 2.3e65 J^{3/2} m^{-3/2}, whatever the route,
                            the vacuum scale or the NEC-violating ingredient; the best technology is 10^39.5 short.  Escapes: quantum
                            tunnelling (bounded by R10/R37), a pre-existing bridge (natural black holes; Round 6c needs horizon-scale false vacuum
                            inside them), a charged Cauchy horizon, an anti-trapped (cosmological) background.
* ``side_change_lemma``     Thin-shell form: beta_out changes sign only at m_s = 0 or where f_out <= 0 (3351 of 3351 random sign changes).
"""

from __future__ import annotations

import numpy as np

from .accessibility import WINDOW_ALPHA1, E_P_GEV, L_P, physical_scales, background_requirements

CLMT_M_MAX_GEV = 1e12          # M < min(1e12 GeV, g^2 1e15 GeV)
CLMT_M_STAR_GEV = 1e15
UNGAUGED_M_MAX_GEV = 100.0     # minimal ghost condensate (published abstract of hep-th/0603010)


# --------------------------------------------------------------------------- route 4


def nec_charge_identity(P1, X, phidot, M=1.0) -> dict:
    """rho + p and the shift-charge density for L = M^4 P(X) (X = phidot^2 / M^4 in the condensate frame)."""
    return dict(rho_plus_p=2 * X * P1 * M**4, J0=2 * P1 * phidot, ratio=2 * X * P1 * M**4 / (2 * P1 * phidot))


def coulomb_ratio(g, M, R, P2_delta=1.0) -> float:
    """Coulomb energy density of a uniform charged NEC-violating ball / |rho + p| = (P''|delta| / 9) g^2 M^2 R^2 (natural units)."""
    return float(P2_delta / 9 * g**2 * M**2 * R**2)


def gauged_window(c_screen=30.0, K_lo=WINDOW_ALPHA1[0]) -> dict:
    """Smallest condensate scale compatible with tipping (M^3 r_N >= K_lo), no screening (r_N <= c / (g M)) and the CLMT bound
    M <= g^2 M_*; compared with the CLMT ceiling 1e12 GeV."""
    M_star = CLMT_M_STAR_GEV / E_P_GEV
    M_min = (K_lo**2 / (c_screen**2 * M_star)) ** (1 / 3)
    return dict(M_min_GeV=float(M_min * E_P_GEV), M_max_GeV=CLMT_M_MAX_GEV, empty=bool(M_min * E_P_GEV > CLMT_M_MAX_GEV),
                log10_gap=float(np.log10(M_min * E_P_GEV / CLMT_M_MAX_GEV)), c_screen=c_screen)


def strong_coupling_scale_overall_factor(F, M=1.0) -> float:
    """L = F(chi) [M^4 P(X) - ...]: canonical normalisation multiplies the leading cubic coupling by F^{-1/2}; its dimensionless strength
    (E/M)^{1/4} (ghost-condensate scaling) becomes F^{-1/2} (E/M)^{1/4}, strong at E = M F^2."""
    return float(M * F**2)


def vacuum_dependent(M_out_GeV=UNGAUGED_M_MAX_GEV, c_wall=1.0) -> dict:
    """Condensate scale M(chi): M_out in our vacuum (observationally bounded), large inside.  Wall thickness >= c_wall / M_out and <= r_N
    gives r_N >= c_wall / M_out; the inside scale then follows from the tipping window."""
    r_N = c_wall * E_P_GEV / M_out_GeV                       # Planck lengths
    sc = physical_scales(r_N)
    M_in = (WINDOW_ALPHA1[0] / r_N) ** (1 / 3) * E_P_GEV
    return dict(r_N_m=sc["r_N_m"], M_N_kg=sc["M_N_kg"], vac_scale_GeV=sc["vac_scale_GeV"], M_in_GeV=float(M_in),
                assembly=background_requirements(sc["vac_scale_GeV"]))


# --------------------------------------------------------------------------- route 6


def dtr_region_c(r_c, M_A, M_B, Lam_C, M_D=0.0) -> dict:
    """Mass parameter of the pocket region C (SdS with Lambda_C) produced at radius r_c by a collision of two concentric walls in our
    vacuum: f_A f_D = f_B f_C with A outside, D inside, B between before the collision."""
    fA, fB, fD = 1 - 2 * M_A / r_c, 1 - 2 * M_B / r_c, 1 - 2 * M_D / r_c
    fC = fA * fD / fB
    M_C = r_c / 2 * (1 - Lam_C * r_c**2 / 3 - fC)
    return dict(f_C=float(fC), M_C=float(M_C), Lambda_C=Lam_C)


def pocket_completion(R_c, M, chi_C, k_w, Rdot=0.0) -> dict:
    """A pocket of vacuum C (dS, Hubble rate chi_C) with wall tension k_w / 4 pi at radius R_c, exterior mass M.  If the exterior is normal
    there (f_out > 0), beta_out > 0 for the outward orientation and keeps its sign while f_out > 0; the hidden inflating branch needs
    beta_out < 0.  Reports which branch R_c lies on at this mass."""
    from .stable_shells import Shell
    sh = Shell(chi=chi_C, M=M, k_w=k_w)
    f_out = float(sh.f_out(R_c))
    tps = sh.turning_points(1e-4 / chi_C, 50 / chi_C)
    branch = "bound" if (tps and R_c <= tps[0]) else ("inflating" if (tps and R_c >= tps[-1]) else "no barrier")
    return dict(f_out=f_out, normal=f_out > 0, beta_out_at_launch=float(np.sqrt(max(f_out, 0) + Rdot**2)) if f_out > 0 else None,
                branch=branch, beta_out_inflating=float(sh.beta_out(tps[-1] * 1.01)) if tps else None,
                completes=False if f_out > 0 else None)


def pocket_scan(chi_C=1.0, k_w=0.5, n=40) -> list:
    """Every launch point in the exterior static region at masses below the merger: bound branch or forbidden; never the hidden branch."""
    from .stable_shells import Shell
    rows = []
    for M in np.geomspace(1e-3, 0.3, n):
        sh = Shell(chi=chi_C, M=M, k_w=k_w)
        tps = sh.turning_points(1e-4, 50.0)
        if len(tps) < 2:
            continue
        R1, R2 = tps[0], tps[-1]
        rows.append(dict(M=float(M), R1=R1, R2=R2, f_out_R1=float(sh.f_out(R1)), beta_out_R1=float(sh.beta_out(R1)),
                         beta_out_R2=float(sh.beta_out(R2 * 1.001)), f_out_R2=float(sh.f_out(R2))))
    return rows


# --------------------------------------------------------------------------- the ledger

def ledger() -> list:
    """State / trigger / completion for each route of Round 6i, with the blocking point.  Numbers from ROADMAP_ROUND6I.md."""
    return [
        dict(route="1 Gauss-Bonnet curvature activation", state="needs horizon-class curvature (lab short by 19 orders)",
             trigger="activation drives F up: pro-collapse (Prop. 6i-1); scalarised Nariai throats collapse or lose predictivity",
             completion="not reached (kinetic degeneracy)", blocking="activation and daughter-tipping have opposite signs"),
        dict(route="2 Monopole core", state="natural if relic monopoles exist (none observed); stable below eta_c ~ 0.25-0.4 m_P",
             trigger="channel at the monopole's own mass ends disconnected from our universe (M << Q); connected channel needs a "
                     "magnetic black hole M > Q ~ 7 m_P, then B ~ 3e3-1.5e4", completion="hidden daughter after tunnelling",
             blocking="Planck-density assembly (M > Q) and e^{-B}"),
        dict(route="3 Breathing bubble", state="composite shell needs negative-energy dust (NEC violation) for the hidden case; every "
                     "well found crosses zero net surface energy at light speed", trigger="tunnelling, B = B_hat (l_f/l_P)^2, B_hat ~ 0.5-0.9 "
                     "(-> 0 only at merger points)", completion="hidden daughter after tunnelling", blocking="exotic dust; consistency of the "
                     "composite shell; e^{-B}"),
        dict(route="4 Gauged / vacuum-dependent ghost", state="gauged universal: window empty (Prop. 6i-8); vacuum-dependent M(chi): background "
                     ">= 9e8 kg (Prop. 6i-9, conditional)", trigger="condensate tipping (Rounds 6c-6f)", completion="hidden daughter (simulated "
                     "in Rounds 6e-6f)", blocking="assembly of the 1e9 kg background within r_N / c (~c^5/G power)"),
        dict(route="5 Resonant tunnelling", state="two-well shell potentials exist (charge + gas components)", trigger="resonance lowers the "
                     "exponent from B1 + B2 to max(B2, B1 - B2) at best", completion="as the underlying tunnelling", blocking="the barrier next "
                     "to the inflating branch is never removed"),
        dict(route="6 Classical transitions", state="pockets of higher vacuum allowed kinematically (planar: <= parent vacuum)",
             trigger="collision", completion="impossible in a normal parent (Prop. 6i-10, Theorem 6h-1); possible in an inflating parent",
             blocking="needs an anti-trapped parent: the cosmological route"),
        dict(route="7 Four-form membrane", state="flux pocket = Farhi-Guth wall (NEC saturated)", trigger="quantum creation, B(M) of earlier "
                     "rounds", completion="hidden daughter after tunnelling", blocking="Theorem 6h-1 classically; e^{-B} with B ~ (l_f/l_P)^2"),
    ]


# --------------------------------------------------------------------------- the synthesis: the bridge invariant (Theorem 6i-11)

G_SI = 6.67430e-11
C_SI = 2.99792458e8


def bridge_invariant() -> dict:
    """Theorem 6i-11.  (1) A hidden daughter has beta_out < 0 at its wall: r decreases outward there, and r -> infinity in our asymptotic
    region, so r has a minimum on every spacelike slice in between; at that sphere n.grad r = 0, grad r is timelike or null, and
    1 - 2Gm/(c^2 r) = (grad r)^2 <= 0.  (2) From regular data (2m/r < 1 everywhere) the first such sphere is marginal, 2Gm = c^2 r, a black-hole
    formation event, with E sqrt(rho_mean) = sqrt(3/(32 pi)) c^6 / G^{3/2} exactly (E = m c^2, rho_mean = 3 m c^2 / (4 pi r^3)): scale-free and
    independent of the mechanism.  Negative energy lowers m and cannot help."""
    val = np.sqrt(3 / (32 * np.pi)) * C_SI**6 / G_SI**1.5
    best = max(r["fom"] for r in __import__("babyuniverse.accessibility", fromlist=["x"]).accessibility_table() if r["kind"] == "technology")
    return dict(E_sqrt_rho_SI=float(val), best_technology=float(best), log10_gap=float(np.log10(val / best)))


def side_change_lemma(n=3000, seed=11) -> dict:
    """Thin-shell form: beta_out = (R (f_in - f_out) - m_s^2 / R) / (2 m_s) is continuous where m_s != 0 and beta_out^2 = f_out + Rdot^2, so
    on a classical trajectory beta_out changes sign only where m_s = 0 (crossed at the speed of light) or f_out <= 0 (inside a horizon).
    Checked on random shells: every sign change of beta_out along R within allowed regions sits at a zero of m_s or where f_out <= 0."""
    from .stable_shells import Shell
    rng = np.random.default_rng(seed)
    changes, explained = 0, 0
    for _ in range(n):
        sh = Shell(chi=1.0, M=10 ** rng.uniform(-3, 0), k_w=rng.choice([-1, 1]) * 10 ** rng.uniform(-1.5, 1),
                   m_d=rng.choice([-1, 1]) * 10 ** rng.uniform(-3, 0), k_1=rng.choice([0, 1]) * rng.uniform(-1, 1))
        R = np.geomspace(1e-3, 30.0, 6000)
        with np.errstate(all="ignore"):
            V, b, m, fo = sh.V(R), sh.beta_out(R), sh.ms(R), sh.f_out(R)
        allowed = V <= 0
        for i in np.where(allowed[:-1] & allowed[1:] & (np.sign(b[:-1]) != np.sign(b[1:])))[0]:
            changes += 1
            explained += int(m[i] * m[i + 1] <= 0 or min(fo[i], fo[i + 1]) <= 0)
    return dict(sign_changes=changes, explained=explained)
