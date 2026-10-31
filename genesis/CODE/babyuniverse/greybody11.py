"""Round 11 (G11-12): brane greybody factors for model U and the Hawking emission envelope used by the collider generator.

The emitted Standard Model fields live on the brane, which sees the induced 4D metric ds^2 = -f dt^2 + dr^2/f + r^2 dOmega_2 with the
model-U f(r) = 1 - r^2 psi(mu / r^{D-1}) (alpha = 1, lengths in units of l).  Since g_tt g_rr = -1, the standard 4D master equations
hold for any f:  d^2 psi/dr*^2 + (omega^2 - V) psi = 0 with
    spin 0 (minimal):  V = f [ l(l+1)/r^2 + f'/r ],                    l >= 0
    spin 1/2:          V = W^2 + dW/dr*,  W = k sqrt(f)/r,  k = j + 1/2 >= 1   (the partner potential has the same transmission)
    spin 1:            V = f l(l+1)/r^2,                                l >= 1
Purely ingoing data at the outer horizon; transmission Gamma = omega / |A|^2 with the incident amplitude read off adiabatically
(|A|^2 = |sqrt(k) psi + i psi'/sqrt(k)|^2/4, as in greybody.py).

Checks (tests/test_round11.py): Gamma_0(s=0) -> 4 r_+^2 omega^2 at low frequency (sigma_abs -> A_H); sum_j (2j+1) Gamma_j ->
(omega b_c)^2 at high frequency with b_c the photon-sphere impact parameter (geometric optics).

Per degree of freedom the number spectrum is dN/dt domega = (1/2pi) sum_j (2j+1) Gamma_j(omega) / (e^{omega/T} -+ 1).  The generator
of Round 6o used x^3/(e^x -+ 1) energies with blackbody species weights; ``emission_tables`` replaces both with the greybody result at a
given mu/mu_ext.  Bulk graviton emission is not included (it adds missing energy; recorded as unresolved).
"""

from __future__ import annotations

import types
from math import pi

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

from . import unified6l as U

SPIN_OF = {"quark": 0.5, "top": 0.5, "charged_lepton": 0.5, "neutrino": 0.5, "gluon": 1, "W": 1, "Z": 1, "photon": 1, "higgs": 0}


class Metric:
    """Brane metric function of model U (kind 'hayward' / 'dymnikova') or Tangherlini (kind 'gr'), D = n + 4, alpha = 1."""

    def __init__(self, D=10, mu=None, mu_over_ext=3.0, kind="hayward"):
        self.D, self.kind = D, kind
        if mu is None:
            mu = mu_over_ext * (U.extremal_mu(D, 1.0, kind) if kind != "gr" else 1.0)
        self.mu = float(mu)
        if kind == "gr":
            self.r_plus = self.mu ** (1.0 / (D - 3)); self.r_minus = None
            self.T = (D - 3) / (4 * pi * self.r_plus)
        else:
            t = U.temperature(self.mu, D, 1.0, kind)
            self.r_plus, self.r_minus, self.T = t["r_plus"], t["r_minus"], t["T"]
        rr = minimize_scalar(lambda r: -self.r_over_sqrtf_inv(r), bounds=(self.r_plus * 1.0001, self.r_plus * 20), method="bounded")
        self.r_c = float(rr.x)
        self.b_c = float(self.r_c / np.sqrt(self.f(self.r_c)))

    def f(self, r):
        x = self.mu / np.asarray(r, float) ** (self.D - 1)
        return 1.0 - np.asarray(r, float) ** 2 * U.h_inv(x, 1.0, self.kind)

    def fp(self, r, eps=1e-6):
        return (self.f(r * (1 + eps)) - self.f(r * (1 - eps))) / (2 * eps * r)

    def r_over_sqrtf_inv(self, r):
        return np.sqrt(max(self.f(r), 0.0)) / r


def potential(m: Metric, r, s, l):
    f = m.f(r)
    if s == 0:
        return f * (l * (l + 1) / r**2 + m.fp(r) / r)
    if s == 1:
        return f * l * (l + 1) / r**2
    k = l                                                   # for s = 1/2 the argument l is k = j + 1/2 = 1, 2, ...
    sf = np.sqrt(max(f, 0.0))
    W = k * sf / r
    dW = k * (sf * m.fp(r) / (2 * r) - f * sf / r**2)        # f dW/dr
    return W * W + dW


def transmission(m: Metric, s, l, omega, r_match_fac=60.0):
    """Gamma for spin s (0, 0.5, 1); for s = 0.5, l is k = j + 1/2."""
    if s == 1 and l < 1:
        return 0.0
    rp = m.r_plus
    Vscale = (l * (l + 1) + 2) / rp**2
    f0 = min(1e-6, 1e-3 * omega**2 / Vscale)
    r0 = brentq(lambda r: m.f(r) - f0, rp * (1 + 1e-14), m.r_c)
    r1 = max(r_match_fac * rp, 40.0 * (l + 1) / omega, m.r_c * 3)

    def rhs(r, y):
        ps = y[0] + 1j * y[1]; ch = y[2] + 1j * y[3]
        fr = m.f(r)
        dps = ch / fr
        dch = (potential(m, r, s, l) - omega**2) * ps / fr
        return [dps.real, dps.imag, dch.real, dch.imag]

    k0 = np.sqrt(max(omega**2 - potential(m, r0, s, l), 1e-300))
    y0 = [1.0, 0.0, 0.0, -k0]                                 # psi = e^{-i k r*}
    sol = solve_ivp(rhs, (r0, r1), y0, method="DOP853", rtol=1e-9, atol=1e-12)
    ps = sol.y[0, -1] + 1j * sol.y[1, -1]; ch = sol.y[2, -1] + 1j * sol.y[3, -1]
    V1 = potential(m, r1, s, l)
    if omega**2 <= V1:
        return 0.0
    kk = np.sqrt(omega**2 - V1)
    A2 = abs(np.sqrt(kk) * ps + 1j * ch / np.sqrt(kk)) ** 2 / 4.0
    flux = k0                                                 # transmitted flux |psi|^2 k0 at the horizon side
    return float(min(flux / A2, 1.0))


def sum_partial(m: Metric, s, omega, tol=1e-5, l_max=200):
    """sum_j (2j+1) Gamma_j(omega) over all partial waves."""
    tot = 0.0
    l0 = {0: 0, 0.5: 1, 1: 1}[s]
    for l in range(l0, l_max):
        deg = 2 * l if s == 0.5 else 2 * l + 1               # 2j+1 = 2k for spin 1/2
        g = transmission(m, s, l, omega)
        tot += deg * g
        if l > omega * m.b_c + 2 and deg * g < tol * max(tot, 1e-300):
            break
    return tot


def spectrum(m: Metric, s, x_grid):
    """Per-dof number spectrum dN/(dt dx) on x = omega/T, and the geometric-optics blackbody reference (omega b_c)^2."""
    stat = -1.0 if s in (0, 1) else 1.0
    T = m.T
    S = np.array([sum_partial(m, s, x * T) for x in x_grid])
    occ = 1.0 / (np.exp(x_grid) + stat)
    dN = S * occ * T / (2 * pi)
    geo = (x_grid * T * m.b_c) ** 2 * occ * T / (2 * pi)
    return dict(x=x_grid.tolist(), sum_partial=S.tolist(), dN=dN.tolist(), dN_geo=geo.tolist())


def _spin_job(args):
    D, q, kind, s, x_grid = args
    m = Metric(D, mu_over_ext=q, kind=kind)
    sp = spectrum(m, s, np.asarray(x_grid))
    x = np.array(sp["x"]); dN = np.array(sp["dN"]); geo = np.array(sp["dN_geo"])
    N = np.trapezoid(dN, x); Ng = np.trapezoid(geo, x)
    E = np.trapezoid(dN * x, x); Eg = np.trapezoid(geo * x, x)
    sp.update(N=float(N), N_geo=float(Ng), E=float(E), E_geo=float(Eg), mean_x=float(E / N), mean_x_geo=float(Eg / Ng),
              number_ratio=float(N / Ng), energy_ratio=float(E / Eg))
    return q, s, sp


def emission_tables(D=10, mu_over_ext=(1.05, 1.5, 3.0, 10.0), kind="hayward", x_grid=None, workers=1) -> dict:
    """Greybody spectra for spins 0, 1/2, 1 at several mu/mu_ext, with number and energy fluxes relative to the geometric blackbody."""
    x_grid = np.concatenate([np.geomspace(0.05, 1.0, 8), np.linspace(1.4, 16.0, 18)]) if x_grid is None else np.asarray(x_grid)
    out = {}
    for q in mu_over_ext:
        m = Metric(D, mu_over_ext=q, kind=kind)
        out["%g" % q] = dict(r_plus=m.r_plus, r_minus=m.r_minus, T=m.T, b_c=m.b_c, T_rplus=m.T * m.r_plus, spins={})
    jobs = [(D, q, kind, s, x_grid.tolist()) for q in mu_over_ext for s in (0, 0.5, 1)]
    if workers > 1:
        import multiprocessing as mp
        with mp.Pool(workers) as pool:
            res = pool.map(_spin_job, jobs)
    else:
        res = [_spin_job(j) for j in jobs]
    for q, s, sp in res:
        out["%g" % q]["spins"][str(s)] = sp
    return out


def make_emission(tables: dict, mode="greybody"):
    """Return a patch for pipeline6o.Generator: mode 'greybody' (spin- and mu-dependent spectra and species weights from ``tables``) or
    'geometric' (blackbody species weights, energies from the number spectrum x^2/(e^x -+ 1) instead of x^3)."""
    from . import pipeline6o as P
    keys = sorted(tables, key=float)
    qs = np.array([float(k) for k in keys])
    cdfs, rates = {}, {}
    for k in keys:
        for s, sp in tables[k]["spins"].items():
            x = np.array(sp["x"]); dN = np.maximum(np.array(sp["dN"]), 0)
            c = np.concatenate([[0.0], np.cumsum(0.5 * (dN[1:] + dN[:-1]) * np.diff(x))])
            cdfs[(k, s)] = (x, c / c[-1])
            rates[(k, s)] = sp["N"] / tables[k]["T"]           # per-dof number rate in units of T
    xg = np.linspace(1e-4, 45.0, 6000)
    geo_cdf = {}
    for st, sgn in (("B", -1.0), ("F", 1.0)):
        pdf = xg**2 / (np.exp(xg) + sgn)
        c = np.concatenate([[0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(xg))])
        geo_cdf[st] = c / c[-1]

    def patch(g):
        def _emit(self, T):
            if mode == "geometric":
                w = np.array([d * (1.0 if st == "B" else 0.75) * (np.exp(-mm / T) if mm > 0 else 1.0) for (_, d, st, mm, _, _) in P.SPECIES])
            else:
                mu = getattr(self, "_mu_now", 3.0 * self.me)
                q = mu / self.me
                k = keys[int(np.argmin(np.abs(np.log(qs) - np.log(max(q, 1.0001)))))]
                w = np.array([d * rates[(k, str(SPIN_OF[nm]))] * (np.exp(-mm / T) if mm > 0 else 1.0)
                              for (nm, d, st, mm, _, _) in P.SPECIES])
            for _ in range(50):
                i = self.rng.choice(len(P.SPECIES), p=w / w.sum())
                name, d, st, mm, qq, vis = P.SPECIES[i]
                for _ in range(20):
                    if mode == "geometric":
                        x = float(np.interp(self.rng.random(), geo_cdf[st], xg))
                    else:
                        xx, cc = cdfs[(k, str(SPIN_OF[name]))]
                        x = float(np.interp(self.rng.random(), cc, xx))
                    E = x * T
                    if E > mm:
                        return name, E, vis, qq
            return "photon", 3 * T, True, 0

        orig_T = g.T_of_M

        def T_of_M(self, M):
            self._mu_now = M / (self.A * self.MDg)
            return orig_T(M)

        g._emit = types.MethodType(_emit, g)
        g.T_of_M = types.MethodType(T_of_M, g)
        return g

    return patch
