"""Reduced phase space of the thin wall and comparison of quantizations (QUANTIZATION.md).

Reduced phase space (constraints solved, clock = exterior Killing time t,
Hamiltonian = ADM mass M):

    P(R, Rdot) = (R/G) [ atanh(Rdot/A_out) - atanh(Rdot/A_in) ],
    (+,+) sector:  cosh(G P / R) = (f_in + f_out - kappa^2 R^2) / (2 sqrt(f_in f_out)),
    M = H(R, P) = (R/2G) [1 - (sqrt(f_in) cosh q - sqrt(f_in sinh^2 q + kappa^2 R^2))^2],
    q = G P / R   (eta = 0; for eta > 0 the same relation defines M implicitly).

In the conjugate pair (u, P_u) = (R^2/2, P/R) the Hamiltonian depends on the
momentum only through cosh(G P_u): the operator cosh(G P_u) -> (T_+ + T_-)/2,
T_pm psi(u) = psi(u -+ i hbar G), is a *shift by an imaginary amount in the area
variable* (Delta A = 8 pi Delta u = -+ 8 pi i hbar G).  The quantum constraint is
a second-order difference equation with imaginary step; see QUANTIZATION.md.

This module provides
  * the Lorentzian canonical momentum and the numerical Hamiltonian H(R, P)
    in any sign sector, with a check of Hamilton's equation dR/dt = dH/dP;
  * the (+,+) closed form and its flat-space limit;
  * "naive" effective-potential quantizations (which reproduce the classical
    orbit Rdot^2 + V = 0 but not the gravitational symplectic structure) and
    their WKB exponents, to document that they are inequivalent (dead end).
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import quad

from .model import ThinWallModel
from . import tunnelling


def canonical_momentum(model: ThinWallModel, R: float, Rdot: float) -> float:
    """Lorentzian P(R, Rdot) in a region where both f_in, f_out > 0 (and |Rdot| < |A|)."""
    Ain, Aout = float(model.A_in(R)), float(model.A_out(R))
    return R / model.G * (np.arctanh(Rdot / Aout) - np.arctanh(Rdot / Ain))


def hamiltonian_pp(model: ThinWallModel, R: float, P: float) -> float:
    """Closed form, (+,+) sector, eta = 0."""
    if model.eta != 0:
        raise NotImplementedError
    q = model.G * P / R
    fin = float(model.f_in(R))
    s = np.sqrt(fin) * np.cosh(q) - np.sqrt(fin * np.sinh(q) ** 2 + model.kappa**2 * R**2)
    return R / (2 * model.G) * (1 - s**2)


def hamiltonian_numeric(model: ThinWallModel, R: float, P: float, sector=None,
                        rdot_sign: int = +1, M0: float = None, tol: float = 1e-14) -> float:
    """H(R, P) by Newton iteration in M at fixed R, using the exact derivative
    dP/dM|_R = A_out / (Rdot f_out)  (Hamilton's equation, symbolic.py), and
    staying inside the static sign sector `sector` = (eps_in, eps_out) and on
    the orbit branch sign(Rdot) = rdot_sign.

    A global bracket search is not used: P(M) is defined only on the M-interval
    where (R, M) is classically allowed and in the given sector, and that
    interval can be narrow; the Newton step is damped to remain in it."""
    sector = sector or (int(np.sign(model.A_in(R))), int(np.sign(model.A_out(R))))
    Mv = model.M if M0 is None else M0

    def state(Mx):
        m = model.with_mass(Mx)
        V = float(m.V(R))
        if V >= 0 or m.f_in(R) <= 0 or m.f_out(R) <= 0:
            return None
        Ain, Aout = float(m.A_in(R)), float(m.A_out(R))
        if (int(np.sign(Ain)), int(np.sign(Aout))) != sector:
            return None
        Rd = rdot_sign * np.sqrt(-V)
        if abs(Rd) >= abs(Ain) or abs(Rd) >= abs(Aout):
            return None
        Pm = R / m.G * (np.arctanh(Rd / Aout) - np.arctanh(Rd / Ain))
        dP = Aout / (Rd * float(m.f_out(R)))
        return Pm, dP

    st = state(Mv)
    if st is None:
        raise ValueError("starting mass not in the requested sector/branch")
    for _ in range(200):
        Pm, dP = st
        step = (P - Pm) / dP
        lam = 1.0
        while True:
            Mn = Mv + lam * step
            sn = state(Mn) if Mn > 0 else None
            if sn is not None:
                break
            lam *= 0.5
            if lam < 1e-12:
                raise ValueError("Newton step cannot stay inside the sector")
        Mv, st = Mn, sn
        if abs(lam * step) < tol * max(1.0, abs(Mv)):
            return Mv
    raise ValueError("Newton iteration did not converge")


def hamilton_equation_check(model: ThinWallModel, R: float, sign: int = +1, h: float = 1e-6) -> dict:
    """Compare dR/dt_out = Rdot f_out / A_out with dH/dP (numerical)."""
    V = float(model.V(R))
    if V >= 0:
        raise ValueError("R not classically allowed")
    Rd = sign * np.sqrt(-V)
    P = canonical_momentum(model, R, Rd)
    dRdt = Rd * float(model.f_out(R)) / float(model.A_out(R))
    sec = (int(np.sign(model.A_in(R))), int(np.sign(model.A_out(R))))
    dP = h * max(1.0, abs(P))
    Hp = hamiltonian_numeric(model, R, P + dP, sector=sec, rdot_sign=sign)
    Hm = hamiltonian_numeric(model, R, P - dP, sector=sec, rdot_sign=sign)
    # second-order central difference; h chosen so truncation and round-off balance
    out = dict(R=R, P=P, dRdt=dRdt, dHdP=(Hp - Hm) / (2 * dP),
               H_at_P=hamiltonian_numeric(model, R, P, sector=sec, rdot_sign=sign), M=model.M,
               sector=sec)
    if model.eta == 0 and sec == (1, 1):
        out["H_closed"] = hamiltonian_pp(model, R, P)
    return out


# ------------------------------------------------------ inequivalent quantizations
def B_naive(model: ThinWallModel, mass="membrane") -> float:
    """WKB exponent of a point particle with Lagrangian L = m(R)(Rdot^2 - V)/2,
    which has the same classical orbits Rdot^2 + V = 0 but a different
    symplectic form:  B = 2 int m(R) sqrt(V) dR.
    mass = "membrane": m = 4 pi sigma R^2 ;  "unit": m = 1/G (Planck units)."""
    Ra, Rb = tunnelling.barrier(model)
    if mass == "membrane":
        m = lambda R: 4 * np.pi * model.sigma * R**2  # noqa: E731
    elif mass == "unit":
        m = lambda R: 1.0 / model.G  # noqa: E731
    else:
        raise ValueError(mass)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        val = quad(lambda R: m(R) * np.sqrt(max(model.V(R), 0.0)), Ra, Rb, limit=400)[0]
    return 2 * val


def area_shift_symbol(model: ThinWallModel, u: float, Pu: complex) -> complex:
    """Classical symbol of the (+,+) constraint in (u, P_u):
    C = 2 sqrt(f_in f_out) cosh(G P_u) - (f_in + f_out - 2 kappa^2 u).
    In the barrier P_u = i P_E/R and C = 0 reproduces the Euclidean angle
    relation cos(Theta_out - Theta_in) (symbolic.euclidean_angle_identity)."""
    R = np.sqrt(2 * u)
    fin, fout = float(model.f_in(R)), float(model.f_out(R))
    return 2 * np.sqrt(fin * fout) * np.cosh(model.G * Pu) - (fin + fout - 2 * model.kappa**2 * u)


def eikonal_check(model: ThinWallModel, n: int = 25) -> float:
    """max |C(u, i P_E/R)| over the barrier with the analytic-branch P_E:
    the difference-operator symbol vanishes on the Euclidean orbit."""
    Ra, Rb = tunnelling.barrier(model)
    Rs = np.linspace(Ra + 1e-6 * (Rb - Ra), Rb - 1e-6 * (Rb - Ra), n)
    errs = []
    for R in Rs:
        PE = float(tunnelling.euclidean_momentum(model, R))
        errs.append(abs(area_shift_symbol(model, R**2 / 2, 1j * PE / R)))
    return float(max(errs))
