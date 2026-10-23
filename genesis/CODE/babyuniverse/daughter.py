"""Roadmap Stage 2 (P-Q9, P-Q10, tasks 14-21): quantum sectors, the daughter edge, the daughter clock.

* ``orientation_is_configuration_data``  Prop. Q9.1: on any slice through the wall in the static domain the sign of
  R'/Lambda (the orientation eps) is boost invariant; the slice momentum G P_Lambda / R is the boosted Rdot.
* ``casoratian``                         Theorem Q9.3: Casoratian charges of the shift equation.  A wave with n flipped
  normals (Sorkin angle i pi n) contributes only to the Laurent coefficient q^{-n}, q = exp(2 pi u / h).
* ``regularized_momentum``               Theorem Q10.1: the interior-horizon "edge" of Theorem V-1 is a coordinate
  singularity of the static-frame momentum; P~ = P -+ (R/2G) ln|f_in| is continuous and the constraint is analytic in it.
* ``regularized_hamiltonian``            the true Hamiltonian H(u, P~), analytic across f_in = 0.
* ``clock_orientation``                  Theorem Q10.2: d t_inf / d tau = A_out / f_out > 0 on the parent, < 0 on every
  daughter (eps_out = -1, region III): the asymptotic clock folds back at the transition.
* ``flux_balance``                       Corollary Q10.3: bookkeeping of t_inf-fluxes (H = M as the generator).  In that clock the
  expanding daughter is an incoming wave, |R|^2 = |I|^2 + |D|^2.  Round 3 (F17) withdraws the physical reading
  ("Klein-type", "time-reversed daughter"): the conserved charge of the constraint is positive on the expanding
  double-flip daughter (selection.charge_signs), and every global time function increases along it.
* ``global_time_orientation``            Proposition A5 (round 3): along the parent and along the daughter, on both sides of
  the wall, both Kruskal null coordinates increase with the wall's proper time.  Only the Killing time t_inf folds back.

See ROADMAP_STAGE1_2_4.md, Parts C and D, and ROADMAP_ROUND3_PROBLEMS_ABC.md, Part A.
"""

from __future__ import annotations

import numpy as np

from .model import ThinWallModel
from . import contours as Cc


# ------------------------------------------------------------------ Q9.1
def orientation_is_configuration_data(A: float, Rdot: float, rapidities) -> dict:
    """(R'/Lambda, G P_Lambda / R) on a slice boosted by rapidity beta relative to the wall frame is
    (A cosh b + Rdot sinh b, Rdot cosh b + A sinh b).  For A^2 - Rdot^2 = f > 0 the sign of R'/Lambda equals sign(A)
    for every beta; the sign of the slice momentum does not (it is the time-orientation datum squared away by WDW)."""
    b = np.asarray(rapidities, dtype=float)
    Rp = A * np.cosh(b) + Rdot * np.sinh(b)
    Pl = Rdot * np.cosh(b) + A * np.sinh(b)
    return dict(f=A * A - Rdot * Rdot, sign_Rprime=np.sign(Rp), sign_PLambda=np.sign(Pl),
                invariant=bool(np.all(np.sign(Rp) == np.sign(A))),
                momentum_sign_changes=bool(len(set(np.sign(Pl))) > 1),
                interval=Rp**2 - Pl**2)


# ------------------------------------------------------------------ Q9.3
def wave(u, p: float, n: int, h: float):
    """Exact solution branch of psi(u+ih)+psi(u-ih) = W psi for constant W = 2(-1)^n cosh(p):
    psi = exp(i p u / h - n pi u / h)  (G = 1, hbar = h; Sorkin angle i pi n on top of the real rapidity p)."""
    u = np.asarray(u, dtype=complex)
    return np.exp(1j * p * u / h - n * np.pi * u / h)


def casoratian(psi1, psi2, u, h: float):
    """K[psi1, psi2](u) = psi1(u) psi2(u + ih) - psi1(u + ih) psi2(u)."""
    return psi1(u) * psi2(u + 1j * h) - psi1(u + 1j * h) * psi2(u)


def charge_of_wave(p: float, n: int, h: float, u=None) -> dict:
    """Casoratian of a branch with its sharp-conjugate psi^#(u) = conj(psi(conj u)); returns K(u) e^{2 pi n u / h}
    (constant iff the wave carries only the q^{-n} charge) and the ih-periodicity defect."""
    if u is None:
        u = np.linspace(0.1, 2.0, 7)
    ps = lambda z: wave(z, p, n, h)
    sharp = lambda z: np.conj(wave(np.conj(z), p, n, h))
    K = casoratian(ps, sharp, np.asarray(u, dtype=complex), h)
    Kp = casoratian(ps, sharp, np.asarray(u, dtype=complex) + 1j * h, h)
    scaled = K * np.exp(2 * np.pi * n * np.asarray(u) / h)
    return dict(n=n, scaled=scaled, scaled_spread=float(np.ptp(np.abs(scaled)) / np.max(np.abs(scaled))),
                periodicity_defect=float(np.max(np.abs(Kp - K)) / np.max(np.abs(K))),
                expected=float((-1) ** n * 2 * np.sinh(p)))


# ------------------------------------------------------------------ Q10.1
def _re_artanh(x):
    return 0.5 * np.log(np.abs((1 + x) / (1 - x)))


def lorentzian_momentum(model: ThinWallModel, R: float) -> float:
    """Real-angle momentum P = (R/G)[Re artanh(Rdot/A_out) - Re artanh(Rdot/A_in)] on a branch (Rdot >= 0)."""
    rd = np.sqrt(max(-float(model.V(R)), 0.0))
    return R / model.G * (_re_artanh(rd / model.A_out(R)) - _re_artanh(rd / model.A_in(R)))


def edge_sign(model: ThinWallModel) -> int:
    """s = +1 if the daughter crosses f_in = 0 with eps_in = +1 (P -> -inf), -1 if eps_in = -1 (P -> +inf)."""
    return 1 if model.A_in(1 / model.chi) > 0 else -1


def regularized_momentum(model: ThinWallModel, R: float) -> float:
    """P~ = P - s (R/2G) ln|f_in|: removes the static-frame divergence at the interior de Sitter horizon."""
    s = edge_sign(model)
    return lorentzian_momentum(model, R) - s * R / (2 * model.G) * np.log(abs(float(model.f_in(R))))


def regularized_constraint_residual(model: ThinWallModel, R: float) -> float:
    """Residual of the regularized (real-angle) constraint on the daughter:
    s = +1:  sqrt(f_out) (f_in e^{p} + e^{-p}) = eps_in eps_out (f_in + f_out - kappa^2 R^2),
    s = -1:  sqrt(f_out) (e^{p} + f_in e^{-p}) = eps_in eps_out (...),   p = G P~ / R.
    Analytic in f_in, hence meaningful on both sides of the interior horizon."""
    s = edge_sign(model)
    p = model.G * regularized_momentum(model, R) / R
    fi, fo = float(model.f_in(R)), float(model.f_out(R))
    sp = np.sign(model.A_in(R) if fi > 0 else model.A_in(1 / model.chi)) * np.sign(model.A_out(R))
    lhs = np.sqrt(fo) * ((fi * np.exp(p) + np.exp(-p)) if s > 0 else (np.exp(p) + fi * np.exp(-p)))
    rhs = sp * (fi + fo - model.kappa**2 * R**2)
    return float(lhs - rhs)


def regularized_hamiltonian(model: ThinWallModel, R: float, Pt: float, sign_product: int, root: int) -> float:
    """H(u, P~) = (R/2G)(1 - s^2), s = sqrt(f_out) the root of sigma s^2 - c1 s + sigma (f_in - kappa^2 R^2) = 0,
    c1 = f_in e^{p} + e^{-p} (edge sign +1) or e^{p} + f_in e^{-p} (edge sign -1).  Analytic across f_in = 0
    (the discriminant there is c1^2 + 4 kappa^2 R^2 > 0).  eta = 0."""
    s_e = edge_sign(model)
    p = model.G * Pt / R
    fi = float(model.f_in(R))
    c1 = fi * np.exp(p) + np.exp(-p) if s_e > 0 else np.exp(p) + fi * np.exp(-p)
    disc = c1**2 - 4 * (fi - model.kappa**2 * R**2)
    s = (c1 + root * np.sqrt(disc)) / (2 * sign_product)
    return float(R / (2 * model.G) * (1 - s * s))


def edge_report(model: ThinWallModel, deltas=(1e-2, 1e-4, 1e-6, -1e-6, -1e-4, -1e-2)) -> list:
    """P, P~, constraint residual, H(P~) - M and dt/dtau on both sides of R = 1/chi along the daughter."""
    out = []
    for d in deltas:
        R = (1 - d) / model.chi
        Pt = regularized_momentum(model, R)
        sp = int(np.sign(model.A_in(1 / model.chi)) * np.sign(model.A_out(R)))
        Hs = [regularized_hamiltonian(model, R, Pt, sp, r) for r in (1, -1)]
        out.append(dict(delta=d, R=R, P=lorentzian_momentum(model, R), P_reg=Pt,
                        residual=regularized_constraint_residual(model, R),
                        H_minus_M=min(abs(H - model.M) for H in Hs),
                        dt_dtau=float(model.A_out(R) / model.f_out(R))))
    return out


# ------------------------------------------------------------------ Q10.2 / Q10.3
def clock_orientation(model: ThinWallModel, n: int = 7) -> dict:
    """Sign of d t_inf / d tau = A_out / f_out along the parent branch (2GM, R1] and the daughter [R2, 3 R2]."""
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]
    rh = model.horizons_out()[0]
    Rp = np.linspace(rh + 0.05 * (R1 - rh), R1, n)
    Rd = np.linspace(R2, 3 * R2, n)
    dp = model.A_out(Rp) / model.f_out(Rp)
    dd = model.A_out(Rd) / model.f_out(Rd)
    return dict(parent=dp.tolist(), daughter=dd.tolist(), parent_future_directed=bool(np.all(dp > 0)),
                daughter_past_directed=bool(np.all(dd < 0)), flip_class=Cc.flip_class(model))


def group_velocity(model: ThinWallModel, R: float, expanding: bool = True) -> float:
    """dR/dt_inf = Rdot f_out / A_out for a wall moving with proper-time velocity +-sqrt(-V)."""
    rd = np.sqrt(max(-float(model.V(R)), 0.0)) * (1 if expanding else -1)
    return float(rd * model.f_out(R) / model.A_out(R))


def flux_balance(model: ThinWallModel, B: float) -> dict:
    """Stationary scattering at fixed M with the boundary clock.  Laboratory incident wave (expanding bounded branch)
    has v > 0; the physically expanding daughter has v < 0 (incoming in t_inf).  With transmission factor
    tau = exp(-B):  (i) physical daughter:  |R|^2 = |I|^2 + |D|^2, |D|^2 = tau |I|^2 (reflected > incident in t_inf-flux);
    (ii) t-outgoing daughter (contracting in its proper time): |I|^2 = |R|^2 + |T|^2.
    These are statements about the t_inf clock only (F17); the clock-independent conserved charge is selection.charge_signs."""
    R1, R2 = model.turning_points()[0], model.turning_points()[-1]
    v_lab = group_velocity(model, 0.5 * (R1 + model.horizons_out()[0]) if R1 > model.horizons_out()[0] else R1, True)
    v_dau = group_velocity(model, 1.5 * R2, True)
    tau = float(np.exp(-B))
    return dict(v_lab_incident=v_lab, v_daughter_expanding=v_dau, klein=bool(v_lab > 0 > v_dau),
                physical_daughter=dict(I2=1.0, D2=tau, R2=1.0 + tau),
                t_outgoing_daughter=dict(I2=1.0, T2=tau, R2=1.0 - tau))


def global_time_orientation(model: ThinWallModel, R: float, expanding: bool = True) -> dict:
    """Proposition A5.  At a static point of the wall (f_in, f_out > 0) with proper-time velocity Rdot = +-sqrt(-V):
    on each side, with A the side's A_in / A_out and f its metric function,
        du/dtau = (A - Rdot)/f,  dv/dtau = (A + Rdot)/f,  dt/dtau = A/f   (u, v, t the side's static null and Killing times).
    The wall sits in the static patch whose Killing field has orientation sign(A) (exterior: region I for A_out > 0,
    region III for A_out < 0, Prop. G1; interior: the de Sitter patch containing the centre for A_in > 0, the antipodal
    patch for A_in < 0).  Kruskal coordinates there are U = -sign(A) e^{-k u}, V = sign(A) e^{k v} (black-hole type) or
    U = sign(A) e^{k u}, V = -sign(A) e^{-k v} (cosmological type), so
        dU/dtau = k |U| |A - Rdot| / f > 0,  dV/dtau = k |V| |A + Rdot| / f > 0,
    because the constraint Rdot^2 = A^2 - f with f > 0 gives |Rdot| < |A|, hence sign(A -+ Rdot) = sign(A).
    Returns the signs computed from the junction data (no patch assignment is used beyond sign(A))."""
    f_in, f_out = float(model.f_in(R)), float(model.f_out(R))
    if f_in <= 0 or f_out <= 0:
        raise ValueError("static point required (f_in, f_out > 0)")
    rd = np.sqrt(max(-float(model.V(R)), 0.0)) * (1 if expanding else -1)
    out = {}
    for side, A, f in (("out", float(model.A_out(R)), f_out), ("in", float(model.A_in(R)), f_in)):
        s = np.sign(A)
        du, dv = (A - rd) / f, (A + rd) / f
        out[side] = dict(A=A, dt_dtau=A / f, du_dtau=du, dv_dtau=dv,
                         dU_sign=float(np.sign(s * du)), dV_sign=float(np.sign(s * dv)), Rdot_below_A=bool(abs(rd) < abs(A)))
    out["killing_clock_reversed"] = bool(out["out"]["dt_dtau"] < 0)
    out["global_future_directed"] = all(out[k]["dU_sign"] > 0 and out[k]["dV_sign"] > 0 for k in ("out", "in"))
    return out
