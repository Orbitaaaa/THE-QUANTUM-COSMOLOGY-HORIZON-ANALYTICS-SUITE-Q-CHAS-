"""Round 6c: what the magnetic throat can and cannot do (P-R6f, Section III of the brief).

* ``rnds_f``, ``segment_Lmax``   Theorem 6c-1 (normal false vacuum is sub-horizon).  In a region of vacuum energy Lambda = x/Q^2 with radial
                                 magnetic field, d m = 4 pi r^2 (rho + rho_EM) dr along any curve, so f = 1 - 2m/r = f_C(r) =
                                 1 - C/r + Q^2/r^2 - Lambda r^2/3 exactly (one constant C per region).  On normal data |dr/dl| >= sqrt(f) along
                                 every spacelike radial curve, so a normal false-vacuum segment whose inner end has f = f_a has proper length
                                 L <= L_max = int_{r_a}^{r_c} dr / sqrt(f_C) (r_c: the local cosmological horizon), ~ 2 r_c sqrt(f_a)/W(r_c).
* ``growth_ratio_max``           The largest L_max / (2 l_H) over all depths, l_H = sqrt(3/Lambda): 0.78 (x = 0.3), 0.53 (0.5), 0.41 (1),
                                 for an inner end at 1.01 Q.  Normal data never contain a super-Hubble false-vacuum segment.
* ``qnec_c_needed``, ``species_window``   Corollary 6c-2 (no throat discount for growth).  Growth needs an anti-trapped false-vacuum
                                 region whose null extents exceed ~ l_H (verified with ``throat_evolution``); making one out of normal data
                                 needs the G1 budget at f = O(1), which QNEC-saturating matter supplies only for c >~ 12 pi sqrt(f) Q^2/eta,
                                 i.e. a species length sqrt(c) l_P above the pulse scale Q/eta: outside the effective theory.
* ``lifting_focusing``           Proposition 6c-3 (creating the vacuum energy afterwards costs the same).  Depositing rho over a longitudinal
                                 length L of the throat focuses the crossing rays by at least 4 pi G Q rho L = x L/(2Q): O(1) for L ~ l_H.
* ``corona_identity``            Proposition 6c-4 (field-restored phases).  At the decompactification point Lambda Q^2 = 1/4 the horizon
                                 radius is a = sqrt2 Q and rho_vac = B(a)^2/2 exactly; a phase restored by the throat's own field needs
                                 B^2/2 >= rho (thermodynamic stability of the restored phase), i.e. the AdS_2 condition.  Such coronae saturate
                                 but never exceed x = 1/4, and B ~ 1/r^2 ends them as soon as the region grows.
* ``hiding_growth``              Corollary 6c-5.  Growth needs the inner end at f >= f_* = O(1), i.e. outside the throat mouth of an extremal
                                 hole, where Proposition GSL-10 does not hide anything: the throat cannot both hide and grow a daughter.

Units G = c = hbar = 1.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def rnds_f(r, C, Q, Lam):
    return 1 - C / r + Q**2 / r**2 - Lam * r**2 / 3


def segment_Lmax(Q: float, x: float, r_a: float, f_a: float) -> dict:
    """Theorem 6c-1: maximal proper length of a normal false-vacuum segment starting at areal radius r_a with f = f_a (0 < f_a), and the
    local cosmological horizon r_c where f_C vanishes."""
    Lam = x / Q**2
    C = r_a * (1 + Q**2 / r_a**2 - Lam * r_a**2 / 3 - f_a)
    f = lambda r: rnds_f(r, C, Q, Lam)
    r, step = r_a, 1e-3 * Q
    while f(r) > 0:
        r += step
        step *= 1.2
    rc = brentq(f, r_a, r, xtol=1e-14 * Q)
    W = Lam * rc**2 + Q**2 / rc**2 - 1
    # substitution r = rc - s^2 removes the square-root endpoint singularity; f(rc - s^2) = s^2 [W/rc + O(s^2)], so the integrand is
    # evaluated as 2 / sqrt(f/s^2) with the ratio taken from its series where f itself underflows
    s_max = np.sqrt(rc - r_a)

    def integrand(s):
        if s < 1e-7 * np.sqrt(Q):
            return 2 / np.sqrt(W / rc)
        return 2 * s / np.sqrt(max(f(rc - s * s), 1e-300))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        val = quad(integrand, 0, s_max, limit=800, epsabs=0, epsrel=1e-9)[0]
    return dict(L_max=float(val), r_c=float(rc), C=float(C), approx=float(2 * rc * np.sqrt(f_a) / W), l_H=float(np.sqrt(3 / Lam)))


def growth_ratio_max(Q: float, x: float, r_a: float, fs=None) -> dict:
    """max over f_a of L_max/(2 l_H): below 1 means no normal super-Hubble false-vacuum segment starts at r_a."""
    fs = np.geomspace(1e-5, 0.95, 60) if fs is None else fs
    rows = [(fa, segment_Lmax(Q, x, r_a, fa)["L_max"]) for fa in fs]
    lH = np.sqrt(3 * Q**2 / x)
    ratios = np.array([L / (2 * lH) for _, L in rows])
    k = int(np.argmax(ratios))
    return dict(ratio_max=float(ratios[k]), f_at_max=float(rows[k][0]), l_H=float(lH))


def qnec_c_needed(Q: float, f: float, eta: float = 1.0, r: float = None) -> float:
    """Central charge for which a QNEC-saturating pulse of local scale Q/eta anti-traps spheres with expansion sqrt f (Prop. GSL-6):
    c > 12 pi sqrt(f) Q r / eta."""
    r = Q if r is None else r
    return float(12 * np.pi * np.sqrt(f) * Q * r / eta)


def species_window(Q: float, f: float, eta: float) -> dict:
    """Is the needed c compatible with the species bound?  The pulse scale Q/eta must exceed the species length sqrt(c) l_P."""
    c = qnec_c_needed(Q, f, eta)
    l_sp = np.sqrt(c)
    return dict(c_needed=c, species_length=float(l_sp), pulse_scale=float(Q / eta), inside_eft=bool(Q / eta > l_sp))


def best_eta_window(Q: float, f: float) -> dict:
    """(pulse scale)/(species length) = (Q/eta) / sqrt(12 pi sqrt f Q^2/eta) = 1/sqrt(12 pi sqrt f eta), largest at the smallest admissible
    eta (~1, a pulse no wider than the throat); at eta = 1 it exceeds 1 only for f < 1/(12 pi)^2 = 7.0e-4.  Cheap flips (small f) are inside
    the effective theory; the f = O(1) flips that growth needs are not."""
    ratio = 1 / np.sqrt(12 * np.pi * np.sqrt(f))
    return dict(ratio_at_eta1=float(ratio), inside_eft_at_eta1=bool(ratio > 1), f_limit=float(1 / (12 * np.pi) ** 2))


def lifting_focusing(x: float, L: float, Q: float) -> float:
    """Proposition 6c-3: minimal focusing 4 pi G Q rho L = x L / (2Q) of rays crossing a region of longitudinal length L whose vacuum
    energy rho = x/(8 pi Q^2) was deposited by NEC-satisfying fluxes."""
    return float(x * L / (2 * Q))


def corona_identity(Q: float = 1.0) -> dict:
    """Proposition 6c-4: at Lambda Q^2 = 1/4 the AdS_2/dS_2 branches meet at a^2 = 2 Q^2 and Q^2/a^4 = Lambda, i.e. 8 pi rho_vac = 8 pi rho_B(a)
    with rho_B = B^2/2 = Q^2/(8 pi a^4).  Returns both sides."""
    Lam = 1 / (4 * Q**2)
    a2 = 2 * Q**2
    return dict(a2=a2, rho_vac=float(Lam / (8 * np.pi)), rho_B=float(Q**2 / (8 * np.pi * a2**2)))


def corona_efolds(B_ratio: float) -> float:
    """A restored phase exists while B(r) = B(a) a^2/r^2 >= B_c: the radius can grow by sqrt(B(a)/B_c), i.e. ln(B(a)/B_c)/2 e-folds of r;
    near the decompactification point B(a)^2/2 ~ rho ~ B_c^2/2, so B(a)/B_c -> 1 and the e-folds vanish."""
    return float(0.5 * np.log(max(B_ratio, 1.0)))


def extremal_f_at(r_over_Q: float) -> float:
    return float((1 - 1 / r_over_Q) ** 2)


def hiding_growth(f_star: float) -> dict:
    """Corollary 6c-5: the inner end of a growable segment needs f >= f_star; in extremal RN, f = (1 - Q/r)^2 gives r >= Q/(1 - sqrt f_star).
    Proposition GSL-10 hides only regions below the mouth r = 2Q (f = 1/4) by ~2Q; growth needs r beyond this when f_star >= 1/4."""
    r_min = 1 / (1 - np.sqrt(f_star))
    return dict(r_min_over_Q=float(r_min), beyond_mouth=bool(r_min >= 2.0))
