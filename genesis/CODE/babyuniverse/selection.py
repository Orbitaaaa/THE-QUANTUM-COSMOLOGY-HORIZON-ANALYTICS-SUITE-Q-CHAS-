"""Problem A (Stage 2 tasks 8-11): which transition does the reduced quantum theory permit?

Theorem A1 (charge selection rule).  For the analytic constraint psi(u+ih) + psi(u-ih) = W_+(u) psi(u) (complex/Sorkin
angles, the convention selected by analyticity across the flip and by the full-theory superspace argument), the Casoratian
K[psi, psi^#] is a single ih-periodic analytic function Phi(q), q = exp(2 pi u/h).  A *geometric* solution in a region whose
net flip number is n = n_out - n_in (exterior flip counts +1, interior flip -1; a factor exp(-n pi u/h) times a slowly varying
wave) contributes only to the Laurent coefficient q^{-n}.  Requiring the solution to be geometric in the laboratory region
(n = 0) and in the daughter region (n_D) forces:
    n_D = 0 (double flip):  K = phi_0 everywhere, the usual flux balance |I|^2 - |R|^2 = (k_D/k_L)|D|^2 > 0;
    n_D = 1 (single flip):  K = 0 identically: total reflection, the daughter part is a zero-flux standing wave.
The squared constraint escapes the rule only by admitting non-geometric states (Sorkin angle mismatched by i pi).

Exact check (Theorem A2): the single-flip toy C(u) = -lambda u, i.e. psi(u+ih) + psi(u-ih) = -2 lambda u psi(u), is solved by
cylinder functions Z_nu(x) of order nu = -i u/h and argument x = i/(lambda h).  The unique solution bounded and oscillatory
in the laboratory region (u -> -inf) and bounded in the daughter region (u -> +inf) is the Hankel function H^(1)_nu(x).  Its
Casoratian vanishes identically, and its decay across the barrier follows the analytic angle (B_an), not arccos|C| (B_sq).
"""

from __future__ import annotations

import mpmath as mp

from .model import ThinWallModel
from . import daughter as D


def _toy(lam=1.0, h=0.1, dps=50):
    mp.mp.dps = dps
    lam, h = mp.mpf(lam), mp.mpf(h)
    x = 1j / (lam * h)
    nu = lambda u: -1j * u / h
    sols = {
        "J": lambda u: mp.besselj(nu(u), x),
        "Y": lambda u: mp.bessely(nu(u), x),
        "H1": lambda u: mp.hankel1(nu(u), x),
        "H2": lambda u: mp.hankel2(nu(u), x),
        "decaying": lambda u: mp.e ** (-1j * mp.pi * nu(u)) * mp.besselj(-nu(u), x),
    }
    return lam, h, x, sols


def toy_casoratian(f, u, h):
    g = lambda z: mp.conj(f(mp.conj(z)))
    return f(u) * g(u + 1j * h) - f(u + 1j * h) * g(u)


def toy_report(lam=1.0, h=0.1, us=(-3, -1.5, -0.5, 0, 0.5, 1.5, 3)) -> dict:
    """Residuals, |psi|, and the Laurent structure of K for each cylinder-function solution of the toy."""
    lam, h, x, sols = _toy(lam, h)
    out = {}
    for name, f in sols.items():
        u0 = mp.mpf("0.37")
        res = abs(f(u0 + 1j * h) + f(u0 - 1j * h) + 2 * lam * u0 * f(u0)) / abs(f(u0))
        rows = []
        for u in us:
            u = mp.mpf(u)
            K = toy_casoratian(f, u, h)
            rows.append(dict(u=float(u), abs_psi=float(abs(f(u))), K=complex(K),
                             K_times_q=float(abs(K * mp.e ** (2 * mp.pi * u / h))),
                             K_over_q=float(abs(K * mp.e ** (-2 * mp.pi * u / h)))))
        out[name] = dict(residual=float(res), rows=rows)
    return out


def toy_amplitude_check(lam=1.0, h=0.1, a=0.5) -> dict:
    """|H1(a/lam)| / |H1(-a/lam)| against exp(-int arccos(-lam u) du / h) = exp(-pi a/(lam h)) (analytic angle, B_an) and
    exp(-int arccos|lam u| du / h) (sector-blind, B_sq)."""
    lam_, h_, x, sols = _toy(lam, h)
    f = sols["H1"]
    ratio = abs(f(mp.mpf(a) / lam_)) / abs(f(-mp.mpf(a) / lam_))
    an = mp.pi * a / (lam_ * h_)
    sq = 2 * (a * mp.acos(a) - mp.sqrt(1 - a * a) + 1) / (lam_ * h_)
    return dict(measured_exponent=float(-mp.log(ratio)), analytic_exponent=float(an), sector_blind_exponent=float(sq))


def charge_signs(model: ThinWallModel) -> dict:
    """Signs of the conserved Casoratian charge (sign of the real-angle momentum, K ~ 2|a|^2 sinh(G P/R)) for the incident
    laboratory wave and for the physically expanding daughter, with the daughter's net flip number."""
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]
    rh = model.horizons_out()[0]
    Rl = 0.5 * (R1 + rh) if 0.95 * R1 <= rh else 0.95 * R1
    Rd = 1.05 * R2
    nD = (1 if model.A_out(Rd) < 0 else 0) - (1 if model.A_in(Rd) < 0 else 0)
    return dict(P_lab=float(D.lorentzian_momentum(model, Rl)), P_daughter=float(D.lorentzian_momentum(model, Rd)),
                n_daughter=nD, dRdt_daughter=float(D.group_velocity(model, Rd)))
