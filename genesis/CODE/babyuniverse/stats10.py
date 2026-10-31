"""Round 10C: statistical analysis.

* ``Model``  Binned likelihood L = prod_i Pois(n_i | nu_i) x prod_k N(theta_k), nu_i = mu s_i(theta) + sum_p b_pi(theta), with nuisances:
  a log-normal normalisation per background process (sigma_p), the jet/MET energy scale (template morphing, +-delta), the luminosity,
  the signal normalisation (PDFs at high x, shower) and one Barlow-Beeston-lite parameter per bin for the Monte-Carlo statistics.
* ``q0`` / ``qmu``  Profile-likelihood test statistics; ``Z_asimov`` and ``CLs_asimov`` with the asymptotic formulae of Cowan, Cranmer,
  Gross & Vitells (EPJC 71, 1554 (2011)); ``toys`` validates them with pseudo-experiments.
* ``llr_events_for_5sigma``  Exact distribution of the summed per-event log-likelihood ratio for N events (FFT convolution of the
  single-event distribution, so tail probabilities of 3e-7 need no toys): the number of events for a median 5 sigma separation of two
  simple hypotheses (e.g. model U against a GR hole at the same threshold), with backgrounds mixed in.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm


def met_templates(met, w, edges, scale=0.02):
    """Nominal and MET-scale-shifted histograms (the energy-scale nuisance)."""
    h0 = np.histogram(met, edges, weights=w)[0]
    hu = np.histogram(met * (1 + scale), edges, weights=w)[0]
    hd = np.histogram(met * (1 - scale), edges, weights=w)[0]
    return h0, hu, hd


class Model:
    def __init__(self, sig, bkgs, sig_unc=0.3, lumi_unc=0.02, bkg_unc=None, mc_rel=None):
        """sig: (h0, hu, hd); bkgs: {name: (h0, hu, hd)}; bkg_unc: {name: sigma}; mc_rel: relative MC-stat error of the total background."""
        self.s0, self.su, self.sd = [np.asarray(x, float) for x in sig]
        self.names = list(bkgs)
        self.b0 = np.array([bkgs[k][0] for k in self.names], float)
        self.bu = np.array([bkgs[k][1] for k in self.names], float)
        self.bd = np.array([bkgs[k][2] for k in self.names], float)
        self.bunc = np.array([(bkg_unc or {}).get(k, 0.1) for k in self.names])
        self.sig_unc, self.lumi_unc = sig_unc, lumi_unc
        self.nb = len(self.s0)
        self.mc_rel = np.zeros(self.nb) if mc_rel is None else np.minimum(np.asarray(mc_rel, float), 1.0)
        self.npar = len(self.names) + 3 + self.nb        # process norms, JES, lumi, signal norm, per-bin gammas
        self.aux = np.zeros(self.npar)                    # global observables of the Gaussian constraints (randomised in toys)

    @staticmethod
    def _morph(h0, hu, hd, a):
        return h0 + np.where(a >= 0, a * (hu - h0), -a * (hd - h0))

    def nu(self, mu, th):
        k = len(self.names)
        norms, jes, lum, sn, gam = th[:k], th[k], th[k + 1], th[k + 2], th[k + 3:]
        # log-normal normalisation nuisances (HistFactory "normsys" with exponential interpolation): yields stay positive and a
        # signal cannot be switched off at a bounded cost
        b = sum(np.maximum(self._morph(self.b0[p], self.bu[p], self.bd[p], jes), 0.0) * (1 + self.bunc[p]) ** norms[p] for p in range(k))
        b = b * (1 + self.mc_rel) ** gam
        s = np.maximum(self._morph(self.s0, self.su, self.sd, jes), 0.0) * (1 + self.sig_unc) ** sn
        return np.maximum((mu * s + b) * (1 + self.lumi_unc) ** lum, 1e-12)

    def nll(self, mu, th, n):
        """Negative log-likelihood relative to the saturated model (Poisson deviance / 2): O(1) at the optimum even with 1e7 events
        per bin, which keeps the minimiser precise."""
        v = self.nu(mu, th)
        n = np.asarray(n, float)
        logterm = np.where(n > 0, n * np.log(np.where(n > 0, n, 1.0) / v), 0.0)
        return float(np.sum(v - n + logterm) + 0.5 * np.sum((th - self.aux) ** 2))

    def nll_grad(self, mu, th, n):
        """NLL and its analytic gradient with respect to (mu, theta)."""
        k = len(self.names)
        norms, jes, lum, sn, gam = th[:k], th[k], th[k + 1], th[k + 2], th[k + 3:]
        Lf = (1 + self.lumi_unc) ** lum
        Bp = (1 + self.bunc) ** norms
        G = (1 + self.mc_rel) ** gam
        bm = np.array([np.maximum(self._morph(self.b0[p], self.bu[p], self.bd[p], jes), 0.0) for p in range(k)])
        dbm = np.array([(self.bu[p] - self.b0[p]) if jes >= 0 else (self.b0[p] - self.bd[p]) for p in range(k)])
        sm = np.maximum(self._morph(self.s0, self.su, self.sd, jes), 0.0)
        dsm = (self.su - self.s0) if jes >= 0 else (self.s0 - self.sd)
        S_ = (1 + self.sig_unc) ** sn
        bk = (bm * Bp[:, None]).sum(0)
        nu = np.maximum((mu * sm * S_ + bk * G) * Lf, 1e-12)
        n = np.asarray(n, float)
        logterm = np.where(n > 0, n * np.log(np.where(n > 0, n, 1.0) / nu), 0.0)
        f = float(np.sum(nu - n + logterm) + 0.5 * np.sum((th - self.aux) ** 2))
        w = 1.0 - n / nu
        g_mu = float(np.sum(w * Lf * sm * S_))
        g = np.zeros(self.npar)
        g[:k] = [np.sum(w * Lf * bm[p] * Bp[p] * np.log(1 + self.bunc[p]) * G) for p in range(k)]
        g[k] = np.sum(w * Lf * (mu * dsm * S_ + (dbm * Bp[:, None]).sum(0) * G))
        g[k + 1] = np.sum(w * nu * np.log(1 + self.lumi_unc))
        g[k + 2] = np.sum(w * Lf * mu * sm * S_ * np.log(1 + self.sig_unc))
        g[k + 3:] = w * Lf * bk * G * np.log(1 + self.mc_rel)
        g += th - self.aux
        return f, g_mu, g

    def nu_jac(self, mu, th):
        """Expected yields and the Jacobian d nu / d(mu, theta)."""
        k = len(self.names)
        norms, jes, lum, sn, gam = th[:k], th[k], th[k + 1], th[k + 2], th[k + 3:]
        Lf = (1 + self.lumi_unc) ** lum
        Bp = (1 + self.bunc) ** norms
        G = (1 + self.mc_rel) ** gam
        a = jes
        bm = np.maximum(self.b0 + (a * (self.bu - self.b0) if a >= 0 else -a * (self.bd - self.b0)), 0.0)
        dbm = (self.bu - self.b0) if a >= 0 else (self.b0 - self.bd)
        sm = np.maximum(self.s0 + (a * (self.su - self.s0) if a >= 0 else -a * (self.sd - self.s0)), 0.0)
        dsm = (self.su - self.s0) if a >= 0 else (self.s0 - self.sd)
        S_ = (1 + self.sig_unc) ** sn
        bk = (bm * Bp[:, None]).sum(0)
        nu = np.maximum((mu * sm * S_ + bk * G) * Lf, 1e-12)
        J = np.zeros((self.nb, 1 + self.npar))
        J[:, 0] = Lf * sm * S_
        J[:, 1:k + 1] = (Lf * bm * Bp[:, None] * np.log(1 + self.bunc)[:, None] * G).T
        J[:, k + 1] = Lf * (mu * dsm * S_ + (dbm * Bp[:, None]).sum(0) * G)
        J[:, k + 2] = nu * np.log(1 + self.lumi_unc)
        J[:, k + 3] = Lf * mu * sm * S_ * np.log(1 + self.sig_unc)
        J[:, k + 4:] = np.diag(Lf * bk * G * np.log(1 + self.mc_rel))
        return nu, J

    def fit(self, n, mu=None, n_iter=60):
        """Damped Gauss-Newton (Fisher scoring) on the analytic Jacobian, with bounds mu >= 0, |theta| <= 5; several starts for mu."""
        n = np.asarray(n, float)
        free = mu is None
        best = None
        for m0 in ([0.0, 1.0, 3.0] if free else [mu]):
            x = np.concatenate([[m0], np.zeros(self.npar)])
            f = self.nll(x[0], x[1:], n)
            lam = 1e-3
            for _ in range(n_iter):
                nu, J = self.nu_jac(x[0], x[1:])
                w = 1.0 - n / nu
                g = J.T @ w; g[1:] += x[1:] - self.aux
                H = J.T @ (J * (n / nu**2)[:, None]); H[1:, 1:] += np.eye(self.npar)
                idx = slice(0, None) if free else slice(1, None)
                Hs, gs = H[idx, idx], g[idx]
                improved = False
                for _t in range(20):
                    try:
                        step = np.linalg.solve(Hs + lam * np.diag(np.diag(Hs) + 1e-12), gs)
                    except np.linalg.LinAlgError:
                        lam *= 10; continue
                    xn = x.copy(); xn[idx] = x[idx] - step
                    xn[0] = max(xn[0], 0.0); xn[1:] = np.clip(xn[1:], -5, 5)
                    fn = self.nll(xn[0], xn[1:], n)
                    if fn <= f + 1e-12:
                        improved = (f - fn) > 1e-10
                        x, f = xn, fn; lam = max(lam / 10, 1e-9)
                        break
                    lam *= 10
                if not improved:
                    break
            if best is None or f < best[0]:
                best = (f, x.copy())
        return (best[0], best[1]) if free else (best[0], best[1][1:])

    def fit_lbfgs(self, n, mu=None):
        """Minimise over nuisances (and mu if mu is None, mu >= 0) with analytic gradients; several starting points, best kept."""
        best = None
        starts = [None] if mu is not None else [0.0, 1.0, 3.0]
        for m0 in starts:
            if mu is None:
                x0 = np.zeros(self.npar + 1); x0[0] = m0

                def f(x):
                    v, gm, g = self.nll_grad(x[0], x[1:], n)
                    return v, np.concatenate([[gm], g])
                bnds = [(0.0, None)] + [(-5, 5)] * self.npar
            else:
                x0 = np.zeros(self.npar)

                def f(x):
                    v, gm, g = self.nll_grad(mu, x, n)
                    return v, g
                bnds = [(-5, 5)] * self.npar
            r = minimize(f, x0, jac=True, method="L-BFGS-B", bounds=bnds, options=dict(maxiter=3000))
            if best is None or r.fun < best.fun:
                best = r
        return best.fun, best.x

    def q0(self, n):
        l0, _ = self.fit(n, 0.0)
        lh, xh = self.fit(n, None)
        return max(2 * (l0 - lh), 0.0) if xh[0] > 0 else 0.0

    def qmu(self, n, mu=1.0):
        lm, _ = self.fit(n, mu)
        lh, xh = self.fit(n, None)
        if xh[0] > mu:
            return 0.0
        return max(2 * (lm - lh), 0.0)

    def asimov(self, mu):
        return self.nu(mu, np.zeros(self.npar))

    def Z_asimov(self):
        return float(np.sqrt(self.q0(self.asimov(1.0))))

    def CLs_asimov(self, mu=1.0):
        qA = self.qmu(self.asimov(0.0), mu)
        sA = np.sqrt(max(qA, 1e-12))
        # median expected CLs under background only: CLs = (1 - Phi(sqrt(q_mu) )) / Phi(sqrt(q_A) - sqrt(q_mu)) at q_mu = q_A
        return float((1 - norm.cdf(sA)) / max(norm.cdf(0.0), 1e-300))

    def toys(self, mu_true=1.0, n_toys=200, seed=1, stat="q0"):
        """Pseudo-experiments: Poisson main measurements and Gaussian-fluctuated global observables of every constraint."""
        rng = np.random.default_rng(seed)
        nu = self.asimov(mu_true)
        out = []
        for _ in range(n_toys):
            n = rng.poisson(nu)
            self.aux = rng.normal(0.0, 1.0, self.npar)
            out.append(self.q0(n) if stat == "q0" else self.qmu(n, 1.0))
        self.aux = np.zeros(self.npar)
        return np.array(out)


def scale_model(m: Model, f_sig=1.0, f_lumi=1.0):
    """Copy with the signal scaled by f_sig and all yields by f_lumi (luminosity scan)."""
    import copy
    c = copy.deepcopy(m)
    c.s0, c.su, c.sd = c.s0 * f_sig * f_lumi, c.su * f_sig * f_lumi, c.sd * f_sig * f_lumi
    c.b0, c.bu, c.bd = c.b0 * f_lumi, c.bu * f_lumi, c.bd * f_lumi
    return c


def lumi_for_Z(m: Model, L_ref, Z_target=5.0, lo=1e-6, hi=1e3):
    """Luminosity (same units as L_ref) at which the median Asimov significance reaches Z_target (bisection in log L)."""
    f = lambda L: scale_model(m, 1.0, L / L_ref).Z_asimov() - Z_target
    if f(hi * L_ref) < 0:
        return np.inf
    a, b = np.log(lo * L_ref), np.log(hi * L_ref)
    for _ in range(40):
        c = 0.5 * (a + b)
        a, b = (c, b) if f(np.exp(c)) < 0 else (a, c)
    return float(np.exp(b))


def binned_cls_limit(m: Model, n_obs, cl=0.95):
    """Observed CLs upper limit on mu (asymptotic, q~_mu with the Asimov background-only dataset for sigma)."""
    nA = m.asimov(0.0)

    def cls(mu):
        q = m.qmu(n_obs, mu); qA = m.qmu(nA, mu)
        sq, sA = np.sqrt(max(q, 0)), np.sqrt(max(qA, 1e-12))
        pmu = 1 - norm.cdf(sq)
        pb = norm.cdf(sA - sq)
        return pmu / max(pb, 1e-300)

    lo, hi = 0.0, 1.0
    while cls(hi) > 1 - cl:
        hi *= 2
        if hi > 1e8:
            return np.inf
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if cls(mid) > 1 - cl else (lo, mid)
    return hi


# --------------------------------------------------------------------------------------------------------------- simple vs simple

def llr_events_for_5sigma(pA, pB, alpha=2.87e-7, n_max=5000, grid=4096):
    """pA, pB: binned per-event probabilities (same binning) of hypotheses A (true) and B (null).  Per-event LLR t = ln(pA/pB).
    For N events the LLR under B is the N-fold convolution (FFT on a fine grid); returns the smallest N with median_A(T) beyond the
    (1 - alpha) quantile of T under B, and the exact-distribution p-value at that N."""
    pA = np.asarray(pA, float); pB = np.asarray(pB, float)
    eps = 1e-12
    pA = (pA + eps) / (pA + eps).sum(); pB = (pB + eps) / (pB + eps).sum()
    t = np.log(pA / pB)
    kl = float(np.sum(pA * t))
    lo_t, hi_t = t.min(), t.max()

    def dist(N, p):
        # distribution of the sum of N iid t's with probabilities p, on a grid
        span = N * (hi_t - lo_t) + 1e-9
        dx = span / (grid - 1)
        idx = np.round((t - lo_t) / dx * 1.0).astype(int)
        # single-event pmf on a grid with resolution dx/N? use characteristic function directly on the summed grid
        g = np.zeros(grid); xs = N * lo_t + dx * np.arange(grid)
        single = np.zeros(grid)
        np.add.at(single, np.clip(np.round((t - lo_t) / dx).astype(int), 0, grid - 1), p)
        F = np.fft.rfft(single, 2 * grid)
        tot = np.fft.irfft(F**N, 2 * grid)[:grid]
        tot = np.maximum(tot, 0); tot /= tot.sum()
        return xs, tot

    N = max(1, int(25.0 / max(kl, 1e-9) * 0.05))
    best = None
    for N in sorted(set(np.unique(np.geomspace(1, n_max, 200).astype(int)))):
        xs, dB = dist(N, pB)
        xA, dA = dist(N, pA)
        medA = xA[np.searchsorted(np.cumsum(dA), 0.5)]
        pval = float(dB[xs >= medA].sum())
        if pval < alpha:
            best = (N, pval)
            break
    return dict(kl_per_event=kl, events_5sigma=best[0] if best else None, p_value=best[1] if best else None,
                gaussian_estimate=float(25.0 / max(2 * kl, 1e-12)))
