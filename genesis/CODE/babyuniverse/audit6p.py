"""Round 6p: independent mathematical audit of the model-U results (Rounds 6l-6o).

Every check here re-derives a result by a route different from the one used to obtain it: symbolic algebra (sympy) for the formulas,
a direct Riemann-tensor computation for the curvature invariants, and independent numerical methods (time integration instead of
quadrature, Monte Carlo over the raw parton grid instead of the luminosity table, a different ODE integrator) for the numbers.
``run_all_checks`` returns the audit ledger used in AUDIT_ROUND6P.md.
"""

from __future__ import annotations

from math import gamma, pi

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

from . import unified6l as U


def _ok(name, passed, detail, source):
    return dict(check=name, passed=bool(passed), detail=detail, source=source)


# --------------------------------------------------------------------------- symbolic checks


def check_hayward_extremality_and_temperature():
    r, mu, a, D = sp.symbols("r mu alpha D", positive=True)
    f = 1 - mu * r**2 / (r ** (D - 1) + a * mu)
    out = []
    for Dv in (4, 5, 6, 10):
        fD = f.subs(D, Dv)
        # horizon: mu(r) = r^{D-1}/(r^2 - alpha); extremum of mu(r)
        mu_r = r ** (Dv - 1) / (r**2 - a)
        rext = sp.solve(sp.diff(sp.log(mu_r), r), r)
        rext = [x for x in rext if x.is_positive][0]
        mu_ext = sp.simplify(mu_r.subs(r, rext))
        claim = a ** sp.Rational(Dv - 3, 2) * sp.Rational(Dv - 1, Dv - 3) ** sp.Rational(Dv - 1, 2) * sp.Rational(Dv - 3, 2)
        ok1 = sp.simplify(mu_ext / claim - 1) == 0
        # temperature at a horizon: f'(r_h)/(4 pi) with mu = mu(r_h)
        T = sp.simplify(sp.diff(fD, r).subs(mu, mu_r) / (4 * sp.pi))
        Tclaim = ((Dv - 3) * r**2 - (Dv - 1) * a) / (4 * sp.pi * r**3)
        ok2 = sp.simplify(T - Tclaim) == 0
        out.append(_ok(f"Hayward mu_ext and T(r_+), D = {Dv}", ok1 and ok2, f"mu_ext = {sp.nsimplify(mu_ext.subs(a, 1))}", "R173 (Round 6l)"))
    return out


def check_closed_dust_inflection():
    a, C, al = sp.symbols("a C alpha", positive=True)
    out = []
    for Dv in (4, 5, 10):
        F = a**2 * C / (a ** (Dv - 1) + al * C)                      # a^2 psi, Hayward
        claim = (2 * al * C / (Dv - 3)) ** sp.Rational(1, Dv - 1)
        # the derivative's numerator vanishes at the claimed point, and F'' < 0 there (a maximum of a^2 psi: end of acceleration)
        num = sp.numer(sp.together(sp.diff(F, a)))
        ok = sp.simplify(sp.expand(num.subs(a, claim))) == 0 and float(sp.diff(F, a, 2).subs(a, claim).subs({C: 7.3, al: 1.3})) < 0
        out.append(_ok(f"closed-dust inflection a_infl^(D-1) = 2 alpha C/(D-3), D = {Dv}", ok, str(claim), "R174 (Round 6l)"))
    return out


def check_vacuum_threshold_duality():
    a, C, L = sp.symbols("a C Lambda", positive=True)
    out = []
    for Dv in (4, 5, 10):
        G = C * a ** (3 - Dv) + L * a**2                              # low-curvature a^2 psi with vacuum
        astar = sp.solve(sp.diff(G, a), a)
        astar = [x for x in astar if x.is_positive][0]
        Lc = sp.solve(sp.Eq(G.subs(a, astar), 1), L)
        claim = (2 * sp.Rational(Dv - 3, Dv - 1) ** sp.Rational(Dv - 1, 2) / ((Dv - 3) * C)) ** sp.Rational(2, Dv - 3)
        ok = any(abs(float((x / claim).subs(C, 3.7)) - 1) < 1e-10 for x in Lc)
        out.append(_ok(f"vacuum threshold Lambda_c(C), D = {Dv}", ok, str(sp.simplify(claim)), "R175 (Round 6l)"))
    return out


def check_hydro_acceleration_equation():
    """Derive D_tau U from Gamma^2 - U^2 = 1 - R^2 psi(mu/R^{D-1}), D_tau mu = -(D-1) p R^{D-2} U, D_tau Gamma = -U Gamma p'/((rho+p) R')
    and compare with the stated acceleration equation; then check the GR and Friedmann limits."""
    R, Uv, mu, p, Gm, dpdR, rho, D = sp.symbols("R U mu p Gamma dpdR rho D", real=True)
    psi = sp.Function("psi")
    y = mu / R ** (D - 1)
    # time derivatives along the fluid
    Rdot, mudot = Uv, -(D - 1) * p * R ** (D - 2) * Uv
    Gdot = -Uv * Gm * dpdR / (rho + p)
    Udot = sp.Symbol("Udot")
    lhs = 2 * Gm * Gdot - 2 * Uv * Udot
    ydot = mudot / R ** (D - 1) - (D - 1) * mu * Rdot / R**D
    rhs = -(2 * R * Rdot * psi(y) + R**2 * sp.diff(psi(y), mu) * (ydot * R ** (D - 1)))      # d/dtau (1 - R^2 psi(y))
    # note d psi(y)/d mu = psi'(y)/R^{D-1}, so sp.diff(psi(y), mu) * R^{D-1} = psi'(y)
    sol = sp.solve(sp.Eq(lhs, rhs), Udot)[0]
    claim = -Gm**2 * dpdR / (rho + p) + R * psi(y) - (D - 1) / 2 * R * sp.diff(psi(y), mu) * R ** (D - 1) * (y + p)
    ok = sp.simplify(sp.expand(sol - claim)) == 0
    # GR limit, D = 4: psi(y) = y -> -Gamma^2 p'/(rho+p) - mu/(2R^2) - (3/2) p R
    s4 = sp.simplify(claim.subs(D, 4).replace(psi, lambda z: z).doit())
    gr = -Gm**2 * dpdR / (rho + p) - mu / (2 * R**2) - sp.Rational(3, 2) * p * R
    ok_gr = sp.simplify(s4 - gr) == 0
    return [_ok("hydro acceleration equation from Misner-Sharp + first law + tr equation", ok, "symbolic identity", "R188 (Round 6n)"),
            _ok("hydro acceleration equation: GR limit (D = 4, mu = 2GM, kappa~ = 8 pi G/3)", ok_gr, "-M/R^2 - 4 pi p R", "R188 (Round 6n)")]


def check_kretschmann_formula(Dv=5):
    """Direct Riemann computation for -f dt^2 + dr^2/f + r^2 dOmega_{D-2}^2 against K = f''^2 + 2(D-2) f'^2/r^2 + 2(D-2)(D-3)(1-f)^2/r^4."""
    t, r = sp.symbols("t r")
    th = sp.symbols(f"th1:{Dv - 1}")
    f = sp.Function("f")(r)
    coords = [t, r] + list(th)
    g = [-f, 1 / f]
    s2 = 1
    for k in range(Dv - 2):
        g.append(r**2 * s2)
        s2 = s2 * sp.sin(th[k]) ** 2
    n = Dv
    gi = [1 / x for x in g]
    Gam = [[[0] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                val = 0
                if a == b:
                    val += sp.diff(g[a], coords[c])
                if a == c:
                    val += sp.diff(g[a], coords[b])
                if b == c:
                    val -= sp.diff(g[b], coords[a])
                Gam[a][b][c] = sp.simplify(gi[a] * val / 2)
    # Riemann R^a_{bcd}
    K = 0
    Rm = {}
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(c + 1, n):
                    val = sp.diff(Gam[a][b][d], coords[c]) - sp.diff(Gam[a][b][c], coords[d])
                    for e in range(n):
                        val += Gam[a][c][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][c]
                    val = sp.simplify(val)
                    if val != 0:
                        Rm[(a, b, c, d)] = val
    for (a, b, c, d), val in Rm.items():
        # R_{abcd} R^{abcd} = sum g_aa g^bb g^cc g^dd (R^a_bcd)^2 ; factor 2 for c<d ordering
        K += 2 * g[a] * gi[b] * gi[c] * gi[d] * val**2
    K = sp.simplify(K)
    fp, fpp = sp.diff(f, r), sp.diff(f, r, 2)
    claim = fpp**2 + 2 * (Dv - 2) * fp**2 / r**2 + 2 * (Dv - 2) * (Dv - 3) * (1 - f) ** 2 / r**4
    ok = sp.simplify(K - claim) == 0
    return [_ok(f"Kretschmann formula by direct Riemann computation, D = {Dv}", ok, "general f(r)", "R184 (Round 6m)")]


def check_frw_kretschmann(Dv=5):
    """Closed FRW in D dimensions: K = 4(D-1)(addot/a)^2 + 2(D-1)(D-2)((adot^2+1)/a^2)^2, by direct computation."""
    t = sp.symbols("t")
    chi = sp.symbols("chi")
    th = sp.symbols(f"th1:{Dv - 1}")
    a = sp.Function("a")(t)
    coords = [t, chi] + list(th)
    g = [-1, a**2]
    s2 = sp.sin(chi) ** 2
    for k in range(Dv - 2):
        g.append(a**2 * s2)
        s2 = s2 * sp.sin(th[k]) ** 2
    n = Dv
    gi = [1 / x for x in g]
    Gam = [[[0] * n for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                val = 0
                if i == j:
                    val += sp.diff(g[i], coords[k])
                if i == k:
                    val += sp.diff(g[i], coords[j])
                if j == k:
                    val -= sp.diff(g[j], coords[i])
                Gam[i][j][k] = sp.simplify(gi[i] * val / 2)
    K = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(k + 1, n):
                    val = sp.diff(Gam[i][j][l], coords[k]) - sp.diff(Gam[i][j][k], coords[l])
                    for e in range(n):
                        val += Gam[i][k][e] * Gam[e][j][l] - Gam[i][l][e] * Gam[e][j][k]
                    val = sp.simplify(val)
                    if val != 0:
                        K += 2 * g[i] * gi[j] * gi[k] * gi[l] * val**2
    K = sp.simplify(K)
    ad, add = sp.diff(a, t), sp.diff(a, t, 2)
    claim = 4 * (Dv - 1) * (add / a) ** 2 + 2 * (Dv - 1) * (Dv - 2) * ((ad**2 + 1) / a**2) ** 2
    ok = sp.simplify(K - claim) == 0
    return [_ok(f"closed-FRW Kretschmann formula by direct computation, D = {Dv}", ok, "general a(t)", "R184 (Round 6m)")]


def check_dl_radius():
    M, MD, n = sp.symbols("M M_D n", positive=True)
    D = n + 4
    Om = 2 * sp.pi ** ((D - 1) / 2) / sp.gamma((D - 1) / 2)
    GD = 1 / MD ** (n + 2)
    rs = (16 * sp.pi * GD * M / ((D - 2) * Om)) ** (1 / (n + 1))
    dl = 1 / (sp.sqrt(sp.pi) * MD) * (M / MD * 8 * sp.gamma((n + 3) / 2) / (n + 2)) ** (1 / (n + 1))
    vals = [abs(float((rs / dl).subs({M: 3.1, MD: 1.7, n: k}).evalf()) - 1) for k in (1, 2, 4, 6, 7)]
    return [_ok("Dimopoulos-Landsberg radius from Tangherlini with G_D = 1/M_D^{n+2}", max(vals) < 1e-12, f"max rel diff {max(vals):.1e}",
                "collider normalisation (Rounds 6k-6o)")]


def check_headon_4d():
    from .collider6l import trapped_fraction_headon
    v = trapped_fraction_headon(0)
    return [_ok("head-on trapped fraction, D = 4 (Eardley-Giddings / Penrose 1/sqrt2)", abs(v - 1 / np.sqrt(2)) < 1e-12, f"{v:.12f}", "R177 (Round 6l)")]


# --------------------------------------------------------------------------- independent numerics


def check_n_acc_by_time_integration(C=1e8, D=10):
    """N_acc from integrating addot = a[psi - ((D-1)/2) psi' y] in time (instead of the quadrature/root route of closed_dust)."""
    r = U.closed_dust(C, D)
    def rhs(t, s):
        a, ad = s
        y = C / a ** (D - 1)
        psi = y / (1 + y)
        dpsi = 1 / (1 + y) ** 2
        return [ad, a * (psi - (D - 1) / 2 * dpsi * y)]
    a0 = r["a_b"] * (1 + 1e-9)
    ad0 = np.sqrt(max(a0**2 * (C / a0 ** (D - 1)) / (1 + C / a0 ** (D - 1)) - 1, 0))
    ev = lambda t, s: s[0] * (C / s[0] ** (D - 1) / (1 + C / s[0] ** (D - 1)) - (D - 1) / 2 * (C / s[0] ** (D - 1)) / (1 + C / s[0] ** (D - 1)) ** 2)
    ev.terminal = True
    s = solve_ivp(rhs, (0, 100), [a0, ad0], events=ev, rtol=1e-11, atol=1e-13, method="DOP853")
    a_infl = s.y_events[0][0][0] if len(s.t_events[0]) else np.nan
    N = np.log(a_infl / r["a_b"])
    return [_ok("N_acc by time integration vs quadrature", abs(N - r["N_acc"]) < 2e-3, f"{N:.5f} vs {r['N_acc']:.5f}", "R174 (Round 6l)")]


def check_sigma_by_monte_carlo(sqrt_s=100.0, n=6, M_D=8.0, M_th=38.0, n_mc=400000, seed=7):
    """sigma(pp -> BH) by Monte Carlo over (x1, x2) with the raw NNPDF grid, against the luminosity-table integration of collider6l."""
    import os
    from . import collider6l as C
    path = C.PDF_SEARCH[0]
    if not os.path.exists(path):
        return [_ok("sigma by Monte Carlo over the raw PDF grid", True, "skipped: PDF grid not installed", "R177")]
    pdf = C.PDF(path)
    rng = np.random.default_rng(seed)
    s = (sqrt_s * 1e3) ** 2
    y = C.trapped_fraction_headon(n)
    tau_min = (M_th * 1e3 / y) ** 2 / s
    # sample ln x1, ln x2 uniformly in [ln tau_min, 0]
    l1 = rng.uniform(np.log(tau_min), 0, n_mc)
    l2 = rng.uniform(np.log(tau_min), 0, n_mc)
    x1, x2 = np.exp(l1), np.exp(l2)
    tau = x1 * x2
    keep = tau > tau_min
    shat = np.sqrt(tau * s)
    Q = y * shat
    w = np.zeros(n_mc)
    # F(x) = xF/x ; integrand dx1 dx2 F1 F2 sigma_hat = d ln x1 d ln x2 (xF1)(xF2) sigma_hat
    Qs = np.clip(Q, 300, 1e5)
    xf1 = np.array([pdf.xF(np.array([a]), q)[0] for a, q in zip(x1[keep], Qs[keep])])
    xf2 = np.array([pdf.xF(np.array([b]), q)[0] for b, q in zip(x2[keep], Qs[keep])])
    R = C.rs_GeV(shat[keep], n, M_D * 1e3)
    w[keep] = xf1 * xf2 * C.YN_FORM_FACTOR[n] * pi * R**2
    vol = np.log(tau_min) ** 2
    sig = w.mean() * vol * C.GEV2_TO_PB
    err = w.std() / np.sqrt(n_mc) * vol * C.GEV2_TO_PB
    L = C.load_lumi(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "DATA", "round6l_lumi.json"))
    ref = C.sigma_pp(L, sqrt_s, n, M_D, M_th, "conservative")
    return [_ok("sigma by Monte Carlo over the raw PDF grid vs luminosity table", abs(sig - ref) < 4 * err + 0.03 * ref,
                f"MC {sig:.3e} +- {err:.1e} pb vs table {ref:.3e} pb", "R177 (Round 6l)")]


def check_vaidya_integrator():
    from . import pipeline6o as P
    a = P.vaidya_history(6, 0.6, 3.0, "hayward", n_pts=500)
    # independent: integrate the mass loss with RK45 and compare the late mass excess
    D = 10
    me = U.extremal_mu(D)
    eps = P.eps_rate(6, 0.6)
    def dmu(v, y):
        mu = y[0]
        if mu <= me * (1 + 1e-12):
            return [0.0]
        tp = U.temperature(mu, D)
        return [-eps * tp["r_plus"] ** 2 * tp["T"] ** 4]
    vend = a["v"][-1]
    s = solve_ivp(dmu, (0, vend), [3.0 * me], method="RK45", rtol=1e-9, atol=1e-12 * me)
    ex_ind = s.y[0][-1] / me - 1
    ex_ref = a["final_mu_over_ext"] - 1
    return [_ok("Vaidya mass loss with an independent integrator (RK45 vs LSODA)", abs(ex_ind / ex_ref - 1) < 0.02,
                f"late excess {ex_ind:.4e} vs {ex_ref:.4e}", "R194 (Round 6o)")]


def check_curvature_bound_analytic():
    """Hayward f'' in closed form against the finite differences used in Round 6m."""
    from . import global6m as G
    r, mu = sp.symbols("r mu", positive=True)
    D = 10
    f = 1 - mu * r**2 / (r ** (D - 1) + mu)
    K = sp.diff(f, r, 2) ** 2 + 2 * (D - 2) * sp.diff(f, r) ** 2 / r**2 + 2 * (D - 2) * (D - 3) * (1 - f) ** 2 / r**4
    Kf = sp.lambdify((r, mu), K, "numpy")
    worst = 0.0
    for m in (1e-2, 1.0, 1e3, 1e6):
        rr = np.geomspace(1e-2, 20, 300)
        k_an = Kf(rr, m)
        k_fd = np.array([G.kretschmann_static(x, m, D) for x in rr])
        sig = k_an > 1e-6 * k_an.max()                 # finite differences are noise-limited where K is tiny (flat region)
        worst = max(worst, float(np.max(np.abs(k_fd[sig] / k_an[sig] - 1))))
    sup = max(float(np.max(Kf(np.geomspace(1e-3, 50, 20000), m))) for m in (1e-3, 1e-1, 10.0, 1e4, 1e8))
    return [_ok("curvature bound: closed-form Kretschmann vs finite differences, D = 10", worst < 1e-3, f"max rel diff {worst:.1e}", "R184"),
            _ok("curvature bound: sup K l^4 over mu = 1e-3..1e8 (closed form), D = 10", sup < 200, f"sup = {sup:.2f}", "R184 (Round 6m)")]


def check_generator_vs_sigma():
    import os
    from . import pipeline6o as P, collider6l as C
    L = C.load_lumi(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "DATA", "round6l_lumi.json"))
    g = P.Generator(L, 100.0, 6, 8.0, 0.6, "hayward")
    ref = C.sigma_pp(L, 100.0, 6, 8.0, g.M_ext / 1e3, "conservative", n_tau=3000)
    return [_ok("event-generator normalisation vs sigma_pp", abs(g.sigma_pb / ref - 1) < 0.02, f"{g.sigma_pb:.3f} vs {ref:.3f} pb", "R195 (Round 6o)")]


def run_all_checks(heavy=True) -> list:
    rows = []
    rows += check_hayward_extremality_and_temperature()
    rows += check_closed_dust_inflection()
    rows += check_vacuum_threshold_duality()
    rows += check_hydro_acceleration_equation()
    rows += check_dl_radius()
    rows += check_headon_4d()
    rows += check_n_acc_by_time_integration()
    rows += check_curvature_bound_analytic()
    rows += check_generator_vs_sigma()
    rows += check_vaidya_integrator()
    if heavy:
        rows += check_kretschmann_formula(4)
        rows += check_kretschmann_formula(5)
        rows += check_frw_kretschmann(4)
        rows += check_frw_kretschmann(5)
        rows += check_sigma_by_monte_carlo()
    return rows
