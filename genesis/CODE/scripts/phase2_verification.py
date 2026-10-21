"""Phase II: independent high-precision verification.

Everything here is re-implemented from the formulas in the documents with
mpmath (40 digits), *without importing* babyuniverse.tunnelling:

  (a) B by the arccos formula  B = (2/G) int R arccos C dR,
  (b) the covariant Smarr pieces (interior volume via Stokes, wall area,
      exterior Killing time, exterior excised volume for eta > 0),
  (c) the Smarr identity residual,
  (d) the regime-B critical-limit angle formula versus direct quadrature,
  (e) the M = 0 closed form,
across extreme parameters (small M, gamma -> +-1, de Sitter exteriors)."""
import mpmath as mp

from _common import dump

mp.mp.dps = 40


class M:
    def __init__(s, chi, kap, Mass, eta=0, G=1):
        s.chi, s.kap, s.M, s.eta, s.G = map(mp.mpf, (chi, kap, Mass, eta, G))
        s.L4 = (s.chi**2 - s.eta**2 + s.kap**2) ** 2 + 4 * s.kap**2 * s.eta**2
        s.L2 = mp.sqrt(s.L4)

    def fin(s, R):
        return 1 - s.chi**2 * R**2

    def fout(s, R):
        return 1 - 2 * s.G * s.M / R - s.eta**2 * R**2

    def Ain(s, R):
        return (s.fin(R) - s.fout(R) + s.kap**2 * R**2) / (2 * s.kap * R)

    def Aout(s, R):
        return (s.fin(R) - s.fout(R) - s.kap**2 * R**2) / (2 * s.kap * R)

    def V(s, R):
        return s.fout(R) - s.Aout(R) ** 2

    def dV(s, R):
        return mp.diff(s.V, R)

    def turning(s):
        if s.M == 0:
            return mp.mpf(0), 2 * s.kap / s.L2
        GM = s.G * s.M
        coeffs = [-s.L4 / 4, 0, s.kap**2, GM * (s.chi**2 - s.eta**2 - s.kap**2), 0, 0, -GM**2]
        rts = mp.polyroots(coeffs, maxsteps=200, extraprec=200)
        pos = sorted(mp.re(r) for r in rts if abs(mp.im(r)) < mp.mpf(10) ** -25 and mp.re(r) > 0)
        pos = [mp.findroot(s.V, r) for r in pos]
        return pos[0], pos[-1]


def barrier_integral(m, fun):
    """int fun(R) dR / sqrt(V) over the barrier, R = Ra + (Rb-Ra) sin^2(s/2)."""
    Ra, Rb = m.turning()
    span = Rb - Ra

    def g(s):
        R = Ra + span * mp.sin(s / 2) ** 2
        V = m.V(R)
        if V <= 0:
            return mp.mpf(0)
        return fun(R) * span * mp.sin(s) / 2 / mp.sqrt(V)
    return mp.quad(g, [0, mp.pi / 2, mp.pi])


def B_arccos(m):
    Ra, Rb = m.turning()

    def integrand(R):
        C = (m.fin(R) + m.fout(R) - m.kap**2 * R**2) / (2 * mp.sqrt(m.fin(R) * m.fout(R)))
        C = max(min(C, mp.mpf(1)), mp.mpf(-1))
        return R * mp.acos(C)
    return 2 * mp.quad(integrand, [Ra, (Ra + Rb) / 2, Rb]) / m.G


def smarr(m):
    Ra, Rb = m.turning()
    G = m.G
    sig = m.kap / (4 * mp.pi * G)
    B = B_arccos(m)
    dtW = 2 * barrier_integral(m, lambda R: m.Aout(R) / m.fout(R))
    AW = 2 * 4 * mp.pi * barrier_integral(m, lambda R: R**2)
    rh = 1 / m.chi
    Din = barrier_integral(m, lambda R: m.Ain(R) / m.fin(R))
    W = barrier_integral(m, lambda R: (R**3 - rh**3) / 3 * m.Ain(R) / m.fin(R))
    ax = Din + (mp.pi / m.chi if m.Ain(Rb) < 0 else 0) - (mp.pi / m.chi if (Ra > 0 and m.Ain(Ra) < 0) else 0)
    Vol = 2 * 4 * mp.pi * (W + rh**3 / 3 * ax)
    Wout = 2 * 4 * mp.pi * barrier_integral(m, lambda R: R**3 / 3 * m.Aout(R) / m.fout(R)) if m.eta > 0 else 0
    # endpoint (Morse) terms
    def n(R):
        return (1 if m.Aout(R) < 0 else 0) - (1 if m.Ain(R) < 0 else 0)
    end = n(Rb) * mp.pi * Rb**2 / G - (n(Ra) * mp.pi * Ra**2 / G if Ra > 0 else 0)
    if m.Ain(Rb) < 0:
        end += mp.pi / (G * m.chi**2)
    if Ra > 0 and m.Ain(Ra) < 0:
        end -= mp.pi / (G * m.chi**2)
    if m.M == 0 and m.eta > 0 and m.Aout(Rb) > 0:
        end += mp.pi / (G * m.eta**2)
    cov = (-(3 * m.chi**2) / (8 * mp.pi * G) * Vol - sig / 2 * AW - m.M / 2 * dtW
           + (3 * m.eta**2) / (8 * mp.pi * G) * Wout)
    return dict(B=B, Bt=B - end, cov=cov, residual=B - end - cov, dtW=dtW)


def model_from_gamma(gam, E, kap=1):
    gam, E, kap = mp.mpf(gam), mp.mpf(E), mp.mpf(kap)
    chi = kap * mp.sqrt((1 + gam) / (1 - gam))
    L2 = chi**2 + kap**2
    twoGM = (4 * kap**2 / (E * L2 ** (mp.mpf(4) / 3))) ** mp.mpf(1.5)
    return M(chi, kap, twoGM / 2)


rows = []
cases = [("gamma=-0.5, E=4.8", model_from_gamma(-0.5, 4.8)),
         ("gamma=0.6, E=1.2", model_from_gamma(0.6, 1.2)),
         ("gamma=-0.95, E=50 g1 (small M)", model_from_gamma(-0.95, 50 * 3.9)),
         ("gamma=0.95, E=1.5 g1", model_from_gamma(0.95, 1.5 * 0.1)),
         ("gamma=0, E=1e4 (M/Mcr~1e-6)", model_from_gamma(0.0, 1e4)),
         ("SdS chi=1 kap=1 M=0.1 eta=0.4", M(1, 1, 0.1, 0.4)),
         ("SdS chi=2 kap=1 M=0.12 eta=0.8", M(2, 1, 0.12, 0.8)),
         ("M=0 chi=1 kap=0.5", M(1, 0.5, 0)),
         ("M=0 chi=1 kap=1.5 eta=0.5", M(1, 1.5, 0, 0.5))]
for name, m in cases:
    r = smarr(m)
    rows.append(dict(case=name, B=mp.nstr(r["B"], 25), residual=mp.nstr(r["residual"], 5),
                     rel=mp.nstr(abs(r["residual"]) / max(1, abs(r["B"])), 5)))
    print(f"{name:34s} B={mp.nstr(r['B'], 18):>24s}  Smarr residual={mp.nstr(r['residual'], 3)}")

# closed form M = 0 at high precision
m0 = []
for chi, kap, eta in ((1, 0.5, 0), (1, 1.5, 0), (1, 0.5, 0.3), (1, 1.5, 0.5), (1, 3.56e-3, 0)):
    m = M(chi, kap, 0, eta)
    R0 = 2 * m.kap / m.L2
    p, q = -m.Aout(R0), m.Ain(R0)
    Bc = mp.pi * R0**2 / 2 * (p + q + 2 * max(p, 0) * max(q, 0)) / ((1 + abs(p)) * (1 + abs(q)))
    Bn = B_arccos(m)
    m0.append(dict(chi=chi, kappa=kap, eta=eta, B_closed=mp.nstr(Bc, 25), rel_diff=mp.nstr(abs(Bn / Bc - 1), 3)))
    print(f"M=0 closed form chi={chi} kappa={kap} eta={eta}: rel diff {mp.nstr(abs(Bn / Bc - 1), 3)}")

# regime-B critical-limit angle formula vs direct quadrature close to M_cr
crit = []
for gam in (-0.8, 0.0, 0.6, 0.9):
    gam_m = mp.mpf(gam)
    zs3 = (-gam_m + mp.sqrt(gam_m**2 + 8)) / 2
    zs = mp.cbrt(zs3)
    gmin = zs**2 - 2 * gam_m / zs + zs**-4
    m = model_from_gamma(gam, gmin * (1 + mp.mpf(10) ** -8))
    Dout = barrier_integral(m, lambda R: m.Aout(R) / m.fout(R))
    th = abs(Dout / (4 * m.G * m.M)) / mp.pi
    t = zs3
    F = 3 * mp.sqrt(3) * (t + 1) ** 2 / (4 * (t + 2) * mp.sqrt(t**2 + 2))
    crit.append(dict(gamma=gam, quadrature=mp.nstr(th, 15), closed_form=mp.nstr(F, 15)))
    print(f"gamma={gam}: |theta_D|/pi near M_cr: quadrature {mp.nstr(th, 12)}, closed form {mp.nstr(F, 12)}")
dump("phase2_high_precision.json", dict(smarr=rows, M0=m0, theta_D_critical=crit,
                                        dps=mp.mp.dps))
