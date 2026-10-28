"""Round 6i, route 1: a curvature-activated field (scalar-Gauss-Bonnet coupling).  Can it tip a throat into a daughter or make the false vacuum?

Action (Doneva-Yazadjiev normalisation, G = c = 1):  S = (1/16 pi) int sqrt(-g) [R - 2 (grad phi)^2 - 2 V(phi) + lam^2 f(phi) GB],
GB = R^2 - 4 R_ab R^ab + R_abcd R^abcd.  Field equations  G_ab + Gamma_ab = 2 d_a phi d_b phi - g_ab (grad phi)^2 - g_ab V,
box phi = V'/2 - (lam^2/4) f' GB, with Gamma_ab linear in grad grad F, F = lam^2 f(phi) (``gb_stress_double_null``).

Results (ROADMAP_ROUND6I.md, section 1):

* ``gb_stress_double_null``  [Theorem, symbolic] Gamma_ab for ds^2 = -e^{2s} du dv + r^2 dOmega^2, verified by the off-shell identity
                             div Gamma = -(GB/2) grad F (exact, sympy).  Null component
                                 -Gamma_uu = 4 (2m/r^3) nabla_u nabla_u F + (16 e^{-2s}/r^2) (nabla_u nabla_u r)(F_u r_v + F_v r_u),
                             m the Misner-Sharp mass.
* ``throat_null_flux``       [Proposition 6i-1] On a vacuum-energy throat (r ~ r_N, 2m/r^3 = 1/r_N^2, nabla_u nabla_u r = O(phi^2)) the
                             tipping flux of Round 6c is N = 4 pi r^2 T_uu = r_N^2 phi_u^2 + 2 F_uu (quadratic order, affine u).  The
                             Gauss-Bonnet part is a total second derivative, so its weighted deficit is
                                 W_GB = -2 int F_uu e^{-u/tau} du = -(2/tau^2) int (F - F_0) e^{-u/tau} du - (boundary terms of the same sign
                                 when F grows),
                             negative (pro-collapse) whenever F increases.  In a vacuum-energy throat GB = C^2 + 8 Lambda^2/3 > 0, and a field
                             that the coupling destabilises (curvature activation: m^2 = -(lam^2/4) f'' GB < 0) rolls so as to increase F.
                             Hence curvature activation can only push the throat towards the black hole.  The opposite arrangement
                             (f'' < 0, so F decreases during the roll) needs a potential that is tachyonic at low curvature, i.e. a field
                             that cannot be stored in the false-vacuum ball or in our vacuum.
* ``ks_run``                 [Proposition-N 6i-2] Fully nonlinear homogeneous (Kantowski-Sachs) throat with the coupling, from exact Nariai
                             data: scalarisation (f = (1 - e^{-6 phi^2})/12, lam^2 = 0.02-2 r_N^2) ends in a black hole or in a degenerate
                             kinetic matrix (the homogeneous analogue of loss of hyperbolicity); never in a daughter.  The reversed
                             arrangement (f'' < 0, tachyonic potential with a positive-energy minimum) ends in a black hole or in the
                             degeneracy in all 144 scanned runs.  (Runs with an unbounded potential that reach large b are discarded: the
                             potential becomes negative before b doubles.)
* ``activation_report``      [Assessment] Activation needs lam^2 |GB| ~ 1/L^2 over the device size L.  Laboratory GB <~ 1e-43 m^-4 requires
                             lam ~ 1e21 m against the bound lam <~ 0.28 km (GW230529): short by ~19 orders of magnitude.  Neutron stars and
                             stellar black holes reach the activation curvature; laboratories do not.
* ``hysteresis_bound``       [Proposition 6i-3] A false vacuum lifted by curvature has Lambda_f <= lam^2 |Delta f| GB / 2.  With GB = c_G / l_curv^4
                             (c_G = 12 at a Schwarzschild horizon of radius l_curv) its Hubble length is l_f >= l_curv sqrt(6/(c_G |Delta f|))
                             (l_curv / lam), i.e. l_f >~ l_curv for |Delta f| <~ 1/2 and lam <~ l_curv (effective theory valid): it can fill a
                             daughter-sized region only around an object that is already horizon-class at the scale l_f.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

# --------------------------------------------------------------------------- symbolic: Gauss-Bonnet stress in double-null spherical symmetry


@lru_cache(maxsize=1)
def gb_stress_double_null() -> dict:
    """Gamma_ab (Doneva-Yazadjiev) for ds^2 = -e^{2s} du dv + r^2 dOmega^2 with F(u, v) = lam^2 f(phi), and the GB invariant."""
    u, v, th, ph = sp.symbols("u v theta phi")
    X = [u, v, th, ph]
    r, s, F = (sp.Function(n)(u, v) for n in ("r", "s", "F"))
    g = sp.zeros(4)
    g[0, 1] = g[1, 0] = -sp.exp(2 * s) / 2
    g[2, 2], g[3, 3] = r**2, r**2 * sp.sin(th) ** 2
    gi, n = g.inv(), 4
    Gam = [[[sp.simplify(sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(n)) / 2)
             for c in range(n)] for b in range(n)] for a in range(n)]

    def riem(a, b, c, d):
        e = sp.diff(Gam[a][d][b], X[c]) - sp.diff(Gam[a][c][b], X[d])
        return sp.simplify(e + sum(Gam[a][c][k] * Gam[k][d][b] - Gam[a][d][k] * Gam[k][c][b] for k in range(n)))

    R4 = [[[[riem(a, b, c, d) for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = sp.Matrix(n, n, lambda b, d: sp.simplify(sum(R4[a][b][a][d] for a in range(n))))
    Rs = sp.simplify(sum(gi[a, b] * Ric[a, b] for a in range(n) for b in range(n)))
    Rl = [[[[sp.simplify(sum(g[a, e] * R4[e][b][c][d] for e in range(n))) for d in range(n)] for c in range(n)] for b in range(n)]
          for a in range(n)]
    nz = [(a, a2) for a in range(n) for a2 in range(n) if gi[a, a2] != 0]
    riem2 = 0
    for (a, a2) in nz:
        for (b, b2) in nz:
            for (c, c2) in nz:
                for (d, d2) in nz:
                    riem2 += Rl[a][b][c][d] * gi[a, a2] * gi[b, b2] * gi[c, c2] * gi[d, d2] * Rl[a2][b2][c2][d2]
    Ricup = gi * Ric * gi
    GB = sp.simplify(sp.simplify(riem2) - 4 * sum(Ric[a, b] * Ricup[a, b] for a in range(n) for b in range(n)) + Rs**2)
    dF = [sp.diff(F, x) for x in X]
    H = sp.Matrix(n, n, lambda a, b: sp.simplify(sp.diff(F, X[a], X[b]) - sum(Gam[k][a][b] * dF[k] for k in range(n))))
    Hm = gi * H
    trH = sp.simplify(sum(gi[a, b] * H[a, b] for a in range(n) for b in range(n)))
    Gt = Ric - Rs * g / 2
    RH = sp.simplify(sum(Ricup[a, b] * H[a, b] for a in range(n) for b in range(n)))

    def gamma(m, q):
        e = -2 * Rs * H[m, q] - 4 * trH * Gt[m, q] - 4 * g[m, q] * RH
        e += 4 * sum(Ric[m, a] * Hm[a, q] + Ric[q, a] * Hm[a, m] for a in range(n))
        e += 4 * sum(R4[b][m][a][q] * Hm[a, b] for a in range(n) for b in range(n))
        return sp.simplify(e)

    Gm = sp.Matrix(n, n, lambda a, b: gamma(a, b))
    return dict(syms=(u, v, r, s, F), GB=GB, Gamma=Gm, ginv=gi, Christoffel=Gam, coords=X)


def verify_divergence_identity() -> bool:
    """div Gamma_{. nu} + (GB/2) d_nu F = 0 identically (off shell), nu = u, v."""
    d = gb_stress_double_null()
    X, Gam, gi, Gm = d["coords"], d["Christoffel"], d["ginv"], d["Gamma"]
    F = d["syms"][4]
    Gmix = gi * Gm
    ok = True
    for q in (0, 1):
        e = sum(sp.diff(Gmix[m, q], X[m]) for m in range(4))
        e += sum(Gam[m][m][k] * Gmix[k, q] for m in range(4) for k in range(4))
        e -= sum(Gam[k][m][q] * Gmix[m, k] for m in range(4) for k in range(4))
        ok &= sp.simplify(e + d["GB"] * sp.diff(F, X[q]) / 2) == 0
    return bool(ok)


def null_component_closed_form() -> bool:
    """-Gamma_uu = 4 (2m/r^3) nabla_u nabla_u F + (16 e^{-2s}/r^2)(nabla_u nabla_u r)(F_u r_v + F_v r_u), 2m/r^3 = (1 + 4 e^{-2s} r_u r_v)/r^2."""
    d = gb_stress_double_null()
    u, v, r, s, F = d["syms"]
    ddF = sp.diff(F, u, 2) - 2 * sp.diff(s, u) * sp.diff(F, u)
    ddr = sp.diff(r, u, 2) - 2 * sp.diff(s, u) * sp.diff(r, u)
    two_m_r3 = (1 + 4 * sp.exp(-2 * s) * sp.diff(r, u) * sp.diff(r, v)) / r**2
    closed = 4 * two_m_r3 * ddF + 16 * sp.exp(-2 * s) / r**2 * ddr * (sp.diff(F, u) * sp.diff(r, v) + sp.diff(F, v) * sp.diff(r, u))
    return bool(sp.simplify(-d["Gamma"][0, 0] - closed) == 0)


# --------------------------------------------------------------------------- Proposition 6i-1: the weighted Gauss-Bonnet flux in a throat


def throat_null_flux(phi_u, F_uu, r_N=1.0):
    """N = 4 pi r^2 T_uu on a vacuum-energy throat to quadratic order: r_N^2 phi_u^2 (canonical) + 2 F_uu (Gauss-Bonnet)."""
    return r_N**2 * np.asarray(phi_u) ** 2 + 2 * np.asarray(F_uu)


def weighted_gb_deficit(F, u, tau):
    """W_GB = -2 int F_uu e^{-u/tau} du computed two ways: directly (finite differences) and by the integrated-by-parts identity
    -(2/tau^2) int (F - F_0) e^{-u/tau} du - 2 [F_u e^{-u/tau}] - (2/tau) [(F - F_0) e^{-u/tau}] (upper limit; data start at rest)."""
    F = np.asarray(F, float)
    Fu = np.gradient(F, u)
    Fuu = np.gradient(Fu, u)
    w = np.exp(-u / tau)
    direct = -2 * np.trapezoid(Fuu * w, u)
    dF = F - F[0]
    parts = -(2 / tau**2) * np.trapezoid(dF * w, u) - 2 * Fu[-1] * w[-1] - (2 / tau) * dF[-1] * w[-1]
    return dict(direct=float(direct), by_parts=float(parts))


def nariai_gb_invariant(Lam=1.0) -> float:
    """GB of dS2 x S2 with both radii 1/sqrt(Lambda): 8 Lambda^2 (= C^2 + 8 Lambda^2 / 3 with C^2 = 16 Lambda^2 / 3)."""
    return 8.0 * Lam**2


def sds_gb_invariant(M, Lam, r):
    """GB of Schwarzschild-de Sitter: 48 M^2 / r^6 + 8 Lambda^2 / 3 (positive: activation in a vacuum throat raises F)."""
    return 48 * M**2 / r**6 + 8 * Lam**2 / 3


# --------------------------------------------------------------------------- Proposition-N 6i-2: the nonlinear homogeneous throat


@lru_cache(maxsize=1)
def ks_equations():
    """Kantowski-Sachs ds^2 = -N^2 dt^2 + a^2 dx^2 + b^2 dOmega^2 with homogeneous phi.  sqrt(-g) GB = d/dt[8 a' (1 + b'^2/N^2)/N] (checked),
    so L = N a b^2 (R + 2 phi'^2/N^2 - 2 V) - lam^2 f'(phi) phi' 8 a' (1 + b'^2/N^2)/N.  Returns lambdified accelerations (a'', b'', phi''),
    the Hamiltonian constraint and the determinant of the acceleration (kinetic) matrix, as functions of
    (a, b, phi, a', b', phi', lam^2, f', f'', f''', V, V')."""
    from sympy.calculus.euler import euler_equations
    t, lam2 = sp.symbols("t lam2")
    N, a, b, p = (sp.Function(n)(t) for n in ("N", "a", "b", "p"))
    f, V = sp.Function("f"), sp.Function("V")
    ad, bd, Nd, pd = (sp.diff(q, t) for q in (a, b, N, p))
    R = (2 / b**2 + 4 * sp.diff(b, t, 2) / (N**2 * b) + 2 * bd**2 / (N**2 * b**2) + 2 * sp.diff(a, t, 2) / (N**2 * a)
         + 4 * ad * bd / (N**2 * a * b) - 4 * Nd * bd / (N**3 * b) - 2 * Nd * ad / (N**3 * a))
    K = 8 * ad * (1 + bd**2 / N**2) / N
    L = N * a * b**2 * (R + 2 * pd**2 / N**2 - 2 * V(p)) - lam2 * sp.diff(f(p), p) * pd * K
    A2, B2, P2, A1, B1, P1, A0, B0, P0 = sp.symbols("A2 B2 P2 A1 B1 P1 A0 B0 P0")

    def rep(e):
        e = e.subs({sp.diff(a, t, 2): A2, sp.diff(b, t, 2): B2, sp.diff(p, t, 2): P2})
        return e.subs({ad: A1, bd: B1, pd: P1}).subs({a: A0, b: B0, p: P0})

    E = [sp.simplify(rep(e.lhs.subs(N, 1).doit())) for e in euler_equations(L, [N, a, b, p], t)]
    sol = sp.solve(E[1:], [A2, B2, P2], dict=True)[0]
    Mk = sp.Matrix([[sp.diff(e, q) for q in (A2, B2, P2)] for e in E[1:]])
    f1, f2, f3, V0, V1 = sp.symbols("f1 f2 f3 V0 V1")

    def clean(e):
        e = e.subs(sp.Derivative(f(P0), (P0, 3)), f3).subs(sp.Derivative(f(P0), (P0, 2)), f2).subs(sp.Derivative(f(P0), P0), f1)
        return e.subs(sp.Derivative(V(P0), P0), V1).subs(V(P0), V0)

    args = (A0, B0, P0, A1, B1, P1, lam2, f1, f2, f3, V0, V1)
    acc = [sp.lambdify(args, clean(sol[q]), "numpy") for q in (A2, B2, P2)]
    return dict(acc=acc, con=sp.lambdify(args, clean(E[0]), "numpy"), det=sp.lambdify(args, clean(sp.factor(Mk.det())), "numpy"),
                total_derivative_ok=_ks_gb_total_derivative())


def _ks_gb_total_derivative() -> bool:
    t, x, th, ph = sp.symbols("t x theta phi")
    N, a, b = (sp.Function(n)(t) for n in ("N", "a", "b"))
    X = [t, x, th, ph]
    g = sp.diag(-N**2, a**2, b**2, b**2 * sp.sin(th) ** 2)
    gi, n = g.inv(), 4
    Gam = [[[sp.simplify(sum(gi[A, d] * (sp.diff(g[d, B], X[c]) + sp.diff(g[d, c], X[B]) - sp.diff(g[B, c], X[d])) for d in range(n)) / 2)
             for c in range(n)] for B in range(n)] for A in range(n)]
    R4 = [[[[sp.simplify(sp.diff(Gam[A][d][B], X[c]) - sp.diff(Gam[A][c][B], X[d])
                         + sum(Gam[A][c][k] * Gam[k][d][B] - Gam[A][d][k] * Gam[k][c][B] for k in range(n)))
             for d in range(n)] for c in range(n)] for B in range(n)] for A in range(n)]
    Ric = sp.Matrix(n, n, lambda B, d: sp.simplify(sum(R4[A][B][A][d] for A in range(n))))
    Rs = sp.simplify(sum(gi[A, A] * Ric[A, A] for A in range(n)))
    Rmix = lambda A, B, c, d: gi[B, B] * R4[A][B][c][d]
    riem2 = sum(Rmix(A, B, c, d) * Rmix(c, d, A, B) for A in range(n) for B in range(n) for c in range(n) for d in range(n))
    Ricmix = gi * Ric
    GB = riem2 - 4 * sum(Ricmix[A, B] * Ricmix[B, A] for A in range(n) for B in range(n)) + Rs**2
    K = 8 * sp.diff(a, t) * (1 + sp.diff(b, t) ** 2 / N**2) / N
    return bool(sp.simplify(N * a * b**2 * GB - sp.diff(K, t)) == 0)


def coupling_dy(p):
    """Doneva-Yazadjiev scalarising coupling f = (1 - e^{-6 phi^2}) / 12 (f''(0) = 1): (f', f'', f''')."""
    e = np.exp(-6 * p * p)
    return p * e, (1 - 12 * p * p) * e, (-36 * p + 144 * p**3) * e


def coupling_reverse(p):
    """f = -phi^2 / 2 (f'' = -1): the coupling raises the mass in a GB > 0 throat (no activation)."""
    return -p, -1.0, 0.0


def coupling_reverse_saturating(ps):
    """f = -(ps^2/2)(1 - e^{-phi^2/ps^2}): reversed sign, with |lam^2 f'| bounded (stays inside the effective theory longer)."""
    def F(p):
        e = np.exp(-p * p / ps**2)
        return -p * e, -(1 - 2 * p * p / ps**2) * e, (6 * p / ps**2 - 4 * p**3 / ps**4) * e
    return F


def potential_flat(Lam=1.0):
    return lambda p: (Lam, 0.0)


def potential_hilltop(mu2, V_min=0.5, Lam=1.0):
    """V = Lam - mu2 phi^2/2 + k phi^4/4 with a positive-energy minimum V_min (tachyonic at phi = 0, mass^2 -mu2/2)."""
    k = mu2**2 / (4 * (Lam - V_min))
    return lambda p: (Lam - 0.5 * mu2 * p * p + 0.25 * k * p**4, -mu2 * p + k * p**3)


def ks_run(coupling, potential, lam2, phi0=0.0, phidot0=0.0, T=80.0, b_daughter=40.0, b_hole=0.3, det_floor=1e-3, path=False) -> dict:
    """Evolve from exact Nariai (Lambda = V(0) = 1: a = 1, a' = 0, b' = 0, b fixed by the constraint).  Outcomes: 'daughter' (b grows past
    b_daughter), 'black hole' (b falls below b_hole), 'kinetic degeneracy' (det of the acceleration matrix / (64 a b^4) < det_floor:
    the evolution stops being predictive), or 'undecided'."""
    eq = ks_equations()

    def args(y):
        A0, B0, P0, A1, B1, P1 = y
        f1, f2, f3 = coupling(P0)
        V0, V1 = potential(P0)
        return (A0, B0, P0, A1, B1, P1, lam2, f1, f2, f3, V0, V1)

    V0, _ = potential(phi0)
    y0 = [1.0, 1.0 / np.sqrt(phidot0**2 + V0), phi0, 0.0, 0.0, phidot0]

    def rhs(t, y):
        ar = args(y)
        return [y[3], y[4], y[5], eq["acc"][0](*ar), eq["acc"][1](*ar), eq["acc"][2](*ar)]

    def ev_d(t, y): return y[1] - b_daughter
    def ev_h(t, y): return y[1] - b_hole
    def ev_k(t, y): return eq["det"](*args(y)) / (64 * y[0] * y[1] ** 4) - det_floor
    for e in (ev_d, ev_h, ev_k):
        e.terminal = True
    s = solve_ivp(rhs, (0, T), y0, method="DOP853", rtol=1e-10, atol=1e-13, events=[ev_d, ev_h, ev_k])
    y = s.y[:, -1]
    out = ("daughter" if len(s.t_events[0]) else "black hole" if len(s.t_events[1]) else
           "kinetic degeneracy" if len(s.t_events[2]) else "undecided")
    con = max(abs(eq["con"](*args(s.y[:, i]))) / max(1.0, abs(s.y[0, i] * s.y[1, i] ** 2)) for i in range(s.y.shape[1]))
    Vmin = min(potential(q)[0] for q in s.y[2])
    res = dict(outcome=out, t_end=float(s.t[-1]), b_end=float(y[1]), phi_end=float(y[2]), H_b_end=float(y[4] / y[1]),
               min_V=float(Vmin), rel_constraint=float(con))
    if path:
        res["t"], res["b"], res["phi"] = s.t.tolist(), s.y[1].tolist(), s.y[2].tolist()
    return res


def scalarisation_scan(lam2s=(0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0), phi0=1e-3) -> list:
    """Curvature activation of a Nariai throat (DY coupling, flat potential Lambda = 1)."""
    return [dict(lam2=l, **ks_run(coupling_dy, potential_flat(), l, phi0, 0.0)) for l in lam2s]


def reverse_scan(lam2s=(1.0, 2.0, 5.0, 10.0), sats=(0.1, 0.2, 0.4, None), epss=(-0.2, -1.0, -3.0), V_mins=(0.9, 0.97, 0.99),
                 kick=1e-3) -> list:
    """Reversed sign (F decreases during the roll) with a hilltop potential mu2 = 4 lam2 (1 - eps): net mass^2 2 lam2 eps < 0 at the
    throat; a positive-energy minimum V_min.  sats = None: unsaturated f = -phi^2/2."""
    rows = []
    for l in lam2s:
        for ps in sats:
            for e in epss:
                for vm in V_mins:
                    cpl = coupling_reverse if ps is None else coupling_reverse_saturating(ps)
                    r = ks_run(cpl, potential_hilltop(4 * l * (1 - e), vm), l, 0.0, kick)
                    rows.append(dict(lam2=l, phi_sat=ps, eps=e, V_min=vm, **r))
    return rows


# --------------------------------------------------------------------------- Assessment: where does curvature activate the field?

L_GB_BOUND_M = 280.0      # sqrt(alpha_GB) <~ 0.28 km (GW230529; Sanger et al. 2024, arXiv:2406.03568)
G_SI = 6.674e-11
C_SI = 2.99792458e8


def environment_gb(kind: str) -> float:
    """|GB| in m^-4 (orders of magnitude): Kretschmann 48 (GM/c^2)^2/r^6 outside matter, (8 pi G rho / c^2)^2 inside."""
    rs = lambda M: G_SI * M / C_SI**2
    if kind == "earth surface":
        return 48 * rs(5.972e24) ** 2 / 6.371e6**6
    if kind == "dense lab mass":                    # tungsten, rho = 1.9e4 kg/m^3, interior
        return (8 * np.pi * G_SI * 1.93e4 / C_SI**2) ** 2
    if kind == "neutron star surface":
        return 48 * rs(1.4 * 1.98847e30) ** 2 / 1.2e4**6
    if kind == "10 Msun black hole horizon":
        M = 10 * 1.98847e30
        return 48 * rs(M) ** 2 / (2 * rs(M)) ** 6
    raise ValueError(kind)


def activation_report(L_device=1.0) -> list:
    """Coupling length lam needed for lam^2 |GB| = 1/L^2 (the activation condition over the region size L), against the bound."""
    rows = []
    for kind, L in (("earth surface", L_device), ("dense lab mass", L_device), ("neutron star surface", 1.2e4),
                    ("10 Msun black hole horizon", 2 * G_SI * 10 * 1.98847e30 / C_SI**2)):
        gb = environment_gb(kind)
        lam_req = 1 / (L * np.sqrt(gb))
        rows.append(dict(environment=kind, GB_m4=gb, region_m=L, lam_required_m=lam_req,
                         log10_shortfall=float(np.log10(lam_req / L_GB_BOUND_M))))
    return rows


def hysteresis_bound(lam_over_lcurv: float = 1.0, delta_f: float = 0.5, c_G: float = 12.0) -> dict:
    """Proposition 6i-3.  V_eff = V - lam^2 f GB / 2 (same units as Lambda): curvature can raise V by at most lam^2 |Delta f| GB / 2.
    With GB = c_G / l_curv^4 (c_G = 12 at a Schwarzschild horizon of radius l_curv) the lifted vacuum has Lambda_f <= lam^2 |Delta f| c_G
    / (2 l_curv^4) and Hubble length l_f = sqrt(3 / Lambda_f) >= l_curv sqrt(6 / (c_G |Delta f|)) / (lam / l_curv)."""
    lf_over_lcurv = np.sqrt(6.0 / (c_G * abs(delta_f))) / lam_over_lcurv
    return dict(l_f_over_l_curv_min=float(lf_over_lcurv), horizon_class=bool(lf_over_lcurv > 0.3))
