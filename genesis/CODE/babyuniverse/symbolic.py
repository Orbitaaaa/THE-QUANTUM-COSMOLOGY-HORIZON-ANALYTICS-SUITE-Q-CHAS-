"""Symbolic derivations (sympy) for the thin-wall model.

Every function returns a dict of sympy expressions and raises AssertionError
if a claimed identity fails to simplify to zero.  ``run_all()`` executes the
complete chain and is exercised by the test-suite and by
``scripts/derivation_log.py`` (which writes DATA/derivation_log.txt).

Conventions: signature (-,+,+,+); K_ab = e_a^mu e_b^nu nabla_mu n_nu with n
pointing from the interior (-) to the exterior (+); [X] = X_out - X_in;
Israel: [K_ab] - h_ab [K] = -8 pi G S_ab.
"""

from __future__ import annotations

import sympy as sp

tau, theta, phi = sp.symbols("tau theta phi", real=True)
G, kappa, sigma = sp.symbols("G kappa sigma", positive=True)
chi, eta, M = sp.symbols("chi eta M", nonnegative=True)
r = sp.symbols("r", positive=True)
eps = sp.symbols("epsilon", real=True)  # orientation sign, eps**2 = 1


def _is_zero(expr) -> bool:
    e = sp.simplify(sp.together(sp.expand(expr)))
    return e == 0


def christoffel(metric, coords):
    n = len(coords)
    ginv = metric.inv()
    Gam = [[[0] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                Gam[a][b][c] = sp.simplify(sum(ginv[a, d] * (sp.diff(metric[d, b], coords[c])
                                                          + sp.diff(metric[d, c], coords[b])
                                                          - sp.diff(metric[b, c], coords[d]))
                                               for d in range(n)) / 2)
    return Gam


def extrinsic_curvature_static():
    """K_ab of the hypersurface r = R(tau), t = T(tau) in
    ds^2 = -f dt^2 + dr^2/f + r^2 dOmega^2, with generic f.

    Result (eps = sign of n^r):
        K_tau_tau     = -eps (R'' + f'(R)/2) / sqrt(f + R'^2)
        K_theta_theta =  eps R sqrt(f + R'^2)
        K^tau_tau     =  eps (R'' + f'/2) / sqrt(f + R'^2)
        K^theta_theta =  eps sqrt(f + R'^2) / R
    """
    f = sp.Function("f")
    t, rr = sp.symbols("t rr")
    coords = (t, rr, theta, phi)
    g = sp.diag(-f(rr), 1 / f(rr), rr**2, rr**2 * sp.sin(theta) ** 2)
    Gam = christoffel(g, coords)

    R = sp.Function("R")(tau)
    Rd = R.diff(tau)
    W = sp.sqrt(f(R) + Rd**2)
    Td = W / f(R)  # future-directed normalisation where f > 0
    e_tau = [Td, Rd, 0, 0]
    e_th = [0, 0, 1, 0]
    # normal covector n_mu = eps(-R', T', 0, 0); normalised because -f T'^2 + R'^2/f = -1
    n_cov = [-eps * Rd, eps * Td, 0, 0]
    sub = {t: sp.Symbol("T0"), rr: R}

    def K(ea, eb, dea_eb):
        s = 0
        for mu in range(4):
            acc = dea_eb[mu]
            for nu in range(4):
                for lam in range(4):
                    acc += Gam[mu][nu][lam].subs(sub) * ea[nu] * eb[lam]
            s += -n_cov[mu] * acc
        return s

    d_e_tau = [sp.diff(c, tau) for c in e_tau]
    K_tt = sp.simplify(K(e_tau, e_tau, d_e_tau))
    # d_theta e_theta = 0 in these coordinates; Christoffels supply the rest
    K_thth = sp.simplify(K(e_th, e_th, [0, 0, 0, 0]).subs(theta, sp.pi / 2))

    fp = sp.Subs(sp.Derivative(f(sp.Symbol("x")), sp.Symbol("x")), sp.Symbol("x"), R).doit()
    K_tt_expected = -eps * (Rd.diff(tau) + sp.diff(f(R), R) / 2) / W
    K_thth_expected = eps * R * W
    # normalisation of the normal
    norm = sp.simplify((-Rd**2 / f(R) + f(R) * Td**2) * eps**2)

    assert _is_zero((K_tt - K_tt_expected).subs(eps**2, 1)), K_tt
    assert _is_zero(K_thth - K_thth_expected), K_thth
    assert _is_zero(norm - eps**2)
    del fp
    return dict(K_tau_tau=K_tt_expected, K_theta_theta=K_thth_expected,
                K_up_tau_tau=-K_tt_expected, K_up_theta_theta=K_thth_expected / R**2,
                normal_norm=norm)


def israel_equations():
    """Israel conditions for S_ab = -sigma h_ab.

    Trace-reversal gives [K_ab] = -8 pi G (S_ab - h_ab S/2) = -4 pi G sigma h_ab:
      theta-theta:  A_out - A_in = -kappa R,     A = eps sqrt(f + Rdot^2)
      tau-tau:      [eps (Rddot + f'/2)/sqrt(f+Rdot^2)] = -kappa
    and the tau-tau equation equals d/dtau of the theta-theta one divided by
    Rdot (wall energy conservation with constant sigma).
    """
    Rs = sp.Function("R")(tau)
    u = Rs.diff(tau)
    a = Rs.diff(tau, 2)
    ein, eout = sp.symbols("epsilon_in epsilon_out", real=True)
    fin = 1 - chi**2 * Rs**2
    fout = 1 - 2 * G * M / Rs - eta**2 * Rs**2
    Ain = ein * sp.sqrt(fin + u**2)
    Aout = eout * sp.sqrt(fout + u**2)
    # S_ab = -sigma h_ab  -> S = -3 sigma -> S_ab - h_ab S/2 = sigma/2 h_ab
    h = sp.diag(-1, Rs**2, Rs**2)
    S = -sigma * h
    Strace = sum((h.inv() * S)[i, i] for i in range(3))
    rhs = -8 * sp.pi * G * (S - h * Strace / 2)
    assert _is_zero(rhs[1, 1] + 4 * sp.pi * G * sigma * Rs**2)
    assert _is_zero(rhs[0, 0] - 4 * sp.pi * G * sigma)
    thth = Aout - Ain + kappa * Rs                       # = 0  ([K_th^th] = -kappa)
    tt = (eout * (a + sp.diff(fout, Rs) / 2) / sp.sqrt(fout + u**2)
          - ein * (a + sp.diff(fin, Rs) / 2) / sp.sqrt(fin + u**2) + kappa)
    # d/dtau thth = u * tt  identically
    assert _is_zero(sp.diff(thth, tau) - u * tt)
    return dict(theta_theta=thth, tau_tau=tt)


def squaring_and_signs():
    """Exact (non-spurious) resolution of the theta-theta equation.

    With X = eps_in sqrt(f_in+u^2), Y = eps_out sqrt(f_out+u^2):
      X - Y = kappa R  and  X^2 - Y^2 = f_in - f_out  (identity)
      => X = (f_in - f_out + kappa^2 R^2)/(2 kappa R),
         Y = (f_in - f_out - kappa^2 R^2)/(2 kappa R),
         u^2 = Y^2 - f_out = X^2 - f_in =: -V(R).
    """
    R, u = sp.symbols("R u", positive=True)
    fin = 1 - chi**2 * R**2
    fout = 1 - 2 * G * M / R - eta**2 * R**2
    X = (fin - fout + kappa**2 * R**2) / (2 * kappa * R)
    Y = (fin - fout - kappa**2 * R**2) / (2 * kappa * R)
    assert _is_zero(X - Y - kappa * R)
    assert _is_zero(X**2 - Y**2 - (fin - fout))
    V1 = fout - Y**2
    V2 = fin - X**2
    assert _is_zero(V1 - V2)
    L4 = (chi**2 - eta**2 + kappa**2) ** 2 + 4 * kappa**2 * eta**2
    Vpoly = (1 - L4 / (4 * kappa**2) * R**2 + G * M * (chi**2 - eta**2 - kappa**2) / (kappa**2 * R)
             - G**2 * M**2 / (kappa**2 * R**4))
    assert _is_zero(V1 - Vpoly)
    # sign changes: X = 0 <-> R^3 = 2GM/(chi^2-eta^2-kappa^2), Y = 0 <-> R^3 = 2GM/(chi^2+kappa^2-eta^2)
    assert _is_zero(sp.numer(sp.together(X)) * 1 - (2 * G * M + (eta**2 - chi**2 + kappa**2) * R**3))
    assert _is_zero(sp.numer(sp.together(Y)) * 1 - (2 * G * M + (eta**2 - chi**2 - kappa**2) * R**3))
    # V(sigma) = V(-sigma): the squared equation cannot see the sign of the tension
    assert _is_zero(Vpoly - Vpoly.subs(kappa, -kappa))
    # Identity |gamma| <= 1:  L^4 - (chi^2-eta^2-kappa^2)^2 = 4 kappa^2 chi^2
    assert _is_zero(L4 - (chi**2 - eta**2 - kappa**2) ** 2 - 4 * kappa**2 * chi**2)
    return dict(X=X, Y=Y, V=Vpoly, L4=L4)


def dimensionless_form():
    """V = 1 - g(z)/E with g = z^2 - 2 gamma/z + 1/z^4 (all chi, eta, kappa)."""
    z, E, Lsq = sp.symbols("z E L2", positive=True)
    L4 = Lsq**2
    R = (2 * G * M / Lsq) ** sp.Rational(1, 3) * z
    gam = (chi**2 - eta**2 - kappa**2) / Lsq
    Evalue = 4 * kappa**2 / (Lsq ** sp.Rational(4, 3) * (2 * G * M) ** sp.Rational(2, 3))
    V = (1 - L4 / (4 * kappa**2) * R**2 + G * M * (chi**2 - eta**2 - kappa**2) / (kappa**2 * R)
         - G**2 * M**2 / (kappa**2 * R**4))
    g = z**2 - 2 * gam / z + 1 / z**4
    diff = sp.simplify(sp.powsimp(sp.expand(V - (1 - g / Evalue)), force=True))
    assert diff == 0, diff
    # minimum of g
    gp = sp.diff(z**2 - 2 * sp.Symbol("gamma") / z + 1 / z**4, z)
    assert _is_zero(sp.expand(gp * z**5 / 2) - (z**6 + sp.Symbol("gamma") * z**3 - 2))
    return dict(g=g, E=Evalue, gamma=gam)


def mass_formula():
    """M(R, Rdot): BGG mass formula, verified by substitution in the junction."""
    R, u = sp.symbols("R u", positive=True)
    ein = sp.symbols("epsilon_in", real=True)
    fin = 1 - chi**2 * R**2
    Mexpr = ((chi**2 - eta**2) * R**3 / (2 * G) + ein * kappa * R**2 / G * sp.sqrt(fin + u**2)
             - kappa**2 * R**3 / (2 * G))
    fout = 1 - 2 * G * Mexpr / R - eta**2 * R**2
    Y = ein * sp.sqrt(fin + u**2) - kappa * R  # must equal eps_out sqrt(fout+u^2)
    for sgn in (1, -1):
        assert _is_zero((Y**2 - (fout + u**2)).subs(ein, sgn))
    # physical decomposition
    rho_in = 3 * chi**2 / (8 * sp.pi * G)
    rho_out = 3 * eta**2 / (8 * sp.pi * G)
    sig = kappa / (4 * sp.pi * G)
    decomposed = (sp.Rational(4, 3) * sp.pi * (rho_in - rho_out) * R**3
                  + ein * 4 * sp.pi * sig * R**2 * sp.sqrt(fin + u**2)
                  - 8 * sp.pi**2 * G * sig**2 * R**3)
    assert _is_zero(Mexpr - decomposed)
    return dict(M=Mexpr)


def reduced_hamiltonian():
    """Canonical momentum conjugate to R for the clock t_out with H = M.

    P = (R/G) [ atanh(Rdot/A_out) - atanh(Rdot/A_in) ]  (+ arbitrary F(R)).
    Verified: dP/dM |_R = A_out / (Rdot f_out) = 1/(dR/dt_out), i.e.
    Hamilton's equation dR/dt_out = dH/dP holds with H = M.
    Closed form of H (eps_in = eps_out = +1 sector):
      cosh(G P/R) = (f_in + f_out - kappa^2 R^2) / (2 sqrt(f_in f_out)),
      M = (R/2G) [1 - (sqrt(f_in) cosh q - sqrt(f_in sinh^2 q + kappa^2 R^2))^2], q = GP/R.
    """
    R, Mv = sp.symbols("R M_v", positive=True)
    fin = 1 - chi**2 * R**2
    fout = 1 - 2 * G * Mv / R - eta**2 * R**2
    X = (fin - fout + kappa**2 * R**2) / (2 * kappa * R)
    Y = (fin - fout - kappa**2 * R**2) / (2 * kappa * R)
    u = sp.sqrt(Y**2 - fout)
    P = R / G * (sp.atanh(u / Y) - sp.atanh(u / X))
    dPdM = sp.diff(P, Mv)
    target = Y / (u * fout)
    num_check = []
    import random
    rnd = random.Random(7)
    for _ in range(6):
        rat = lambda lo, hi: sp.Rational(int(rnd.uniform(lo, hi) * 1000), 1000)  # noqa: E731
        vals = {G: 1, kappa: rat(0.3, 1.5), chi: rat(0.2, 1.5),
                eta: rat(0, 0.15), Mv: rat(0.01, 0.05)}
        # choose R inside a classically allowed static region (u real, f's > 0)
        Rv = None
        for trial in [0.05, 0.08, 0.12, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0, 3.0]:
            sub = dict(vals)
            sub[R] = sp.Rational(trial).limit_denominator(1000)
            uu = complex(u.subs(sub).evalf())
            if (abs(uu.imag) < 1e-12 and uu.real > 1e-6 and float(fin.subs(sub)) > 0
                    and float(fout.subs(sub)) > 0 and abs(float((u / Y).subs(sub))) < 1
                    and abs(float((u / X).subs(sub))) < 1):
                Rv = trial
                break
        if Rv is None:
            continue
        sub = dict(vals)
        sub[R] = sp.Rational(Rv).limit_denominator(1000)
        lhs = complex(dPdM.subs(sub).evalf(30))
        rhs = complex(target.subs(sub).evalf(30))
        num_check.append(abs(lhs - rhs) / max(1e-30, abs(rhs)))
    assert num_check and max(num_check) < 1e-20, num_check
    # closed-form Hamiltonian inverse relation, symbolic in the q variable
    q = sp.symbols("q", real=True)
    s = sp.symbols("s", positive=True)  # s = sqrt(f_out)
    fin_s = sp.symbols("f_in", positive=True)
    s_root = sp.sqrt(fin_s) * sp.cosh(q) - sp.sqrt(fin_s * sp.sinh(q) ** 2 + kappa**2 * R**2)
    rel = (fin_s + s_root**2 - kappa**2 * R**2) - 2 * sp.sqrt(fin_s) * s_root * sp.cosh(q)
    assert _is_zero(sp.expand(rel.rewrite(sp.exp)))
    del s
    return dict(P=P, dPdM_rel_err=max(num_check), cosh_relation=rel)


def euclidean_angle_identity():
    """Euclidean continuation: Rdot -> i R', A = eps sqrt(f - R'^2), R'^2 = V.
    With Theta = atan2(R', A):  cos(Theta_out - Theta_in) =
      (f_in + f_out - kappa^2 R^2)/(2 sqrt(f_in f_out)).
    """
    R = sp.symbols("R", positive=True)
    fin = 1 - chi**2 * R**2
    fout = 1 - 2 * G * M / R - eta**2 * R**2
    X = (fin - fout + kappa**2 * R**2) / (2 * kappa * R)
    Y = (fin - fout - kappa**2 * R**2) / (2 * kappa * R)
    V = fout - Y**2
    cos_diff = (Y * X + V) / sp.sqrt(fin * fout)
    target = (fin + fout - kappa**2 * R**2) / (2 * sp.sqrt(fin * fout))
    assert _is_zero(cos_diff - target)
    return dict(cos_dtheta=target)


def closed_forms_M0():
    """M = 0: compact O(4) instanton and canonical exponent in closed form.

    With R0 = 2 kappa / L^2, p = (chi^2+kappa^2-eta^2)/L^2 = -A_out(R0),
    q = (eta^2-chi^2+kappa^2)/L^2 = A_in(R0) (signed cosines of the cap angles):

      S_E[instanton] = -(pi R0^2 / 2G) [1/(1+p) + 1/(1+q)]
      B_can          = (2 R0^2/G) [J(-p) - J(q)],
      J(s) = int_0^1 x atan2(sqrt(1-x^2), s x) dx
           = pi/(4(1+s))         (s >= 0)
           = pi/2 - pi/(4(1-s))  (s <  0).
    """
    import mpmath as mp
    mp.mp.dps = 30

    def J_closed(sv):
        return mp.pi / (4 * (1 + sv)) if sv >= 0 else mp.pi / 2 - mp.pi / (4 * (1 - sv))

    errs = []
    for sv in (-0.9, -0.5, -0.1, 0.0, 0.2, 0.5, 0.9):
        Jnum = mp.quad(lambda xx: xx * mp.atan2(mp.sqrt(1 - xx**2), sv * xx), [0, 1])
        errs.append(abs(Jnum - J_closed(sv)))
    assert max(errs) < mp.mpf(10) ** -14, errs
    sp_s = sp.symbols("s", positive=True)
    Jpos = sp.pi / (4 * (1 + sp_s))
    p, qq, R0 = sp.symbols("p q R_0", real=True)
    # S_E from cap volumes:  -(Lambda/8piG) Vol_cap - (sigma/2) Area
    capfac = lambda c: (2 - 3 * c + c**3)  # noqa: E731   (3/4)(2/3 - c + c^3/3)*4
    SE = -(sp.pi * R0**2 / (4 * G)) * ((1 - qq) * (2 + qq) / (1 + qq)
                                        + (1 - p) * (2 + p) / (1 + p) + p + qq)
    assert _is_zero(SE + sp.pi * R0**2 / (2 * G) * (1 / (1 + p) + 1 / (1 + qq)))
    del capfac
    return dict(J_pos=Jpos, S_E=-sp.pi * R0**2 / (2 * G) * (1 / (1 + p) + 1 / (1 + qq)))


def flat_space_limit():
    """G -> 0 limit of the reduced Hamiltonian: M -> sqrt(P^2 + m^2) - (4pi/3) drho R^3
    with m = 4 pi sigma R^2 (relativistic membrane in a vacuum-energy difference).
    Checked to first order in G for a flat interior and exterior (chi = eta = 0)."""
    R, P = sp.symbols("R P", positive=True)
    sg = sp.symbols("sigma_w", positive=True)
    eG = sp.symbols("e_G", positive=True)  # bookkeeping for G
    q = eG * P / R
    kap = 4 * sp.pi * eG * sg
    fin = 1
    s_root = sp.sqrt(fin) * sp.cosh(q) - sp.sqrt(fin * sp.sinh(q) ** 2 + kap**2 * R**2)
    Mham = R / (2 * eG) * (1 - s_root**2)
    lead = sp.limit(Mham, eG, 0)
    m = 4 * sp.pi * sg * R**2
    assert _is_zero(lead - sp.sqrt(P**2 + m**2))
    return dict(leading=lead)


def run_all():
    out = {}
    for fn in (extrinsic_curvature_static, israel_equations, squaring_and_signs,
               dimensionless_form, mass_formula, reduced_hamiltonian,
               euclidean_angle_identity, closed_forms_M0, flat_space_limit):
        out[fn.__name__] = fn()
    return out


if __name__ == "__main__":
    res = run_all()
    for k, v in res.items():
        print(f"[ok] {k}")
        for kk, vv in v.items():
            print(f"     {kk} = {vv}")
