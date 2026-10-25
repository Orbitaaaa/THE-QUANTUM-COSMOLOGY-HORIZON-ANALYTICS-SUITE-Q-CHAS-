"""Round 6f (item 3): the ghost condensate beyond its IR fluid limit, and non-spherical perturbations of the tipping geometry.

(1) Nonlinear 1+1 evolution of the full condensate EFT along the expanding direction of the Nariai band (dS_2 x S^2, flat dS_2 slicing,
    a = e^{t/r_N}, units r_N = 1).  phi = M^2 (t + pi), X = (1 + pi_t)^2 - a^{-2} pi_x^2, P(X) = (X - 1)^2 / 2 (P'' = 1), plus the leading
    higher-derivative term -(alpha/2M^2)(lap pi)^2 of the ghost condensate and the condensate's Jeans coupling J (from mixing with gravity):
        Q = 2 a P'(X) (1 + pi_t)                (charge per comoving length, units M^4)
        Q_t = d_x[2 a^{-1} P'(X) pi_x] - 2 a^{-1} J pi_xx - (alpha/M^2) a^{-3} pi_xxxx - 2 a sigma,
    and pi_t from Q by the cubic 2 a ((1 + pi_t)^2 - a^{-2} pi_x^2 - 1)(1 + pi_t) = Q.  sigma is the false-vacuum tilt source (sigma = 0:
    the NEC phase dilutes as e^{-t}; sigma > 0: it is regenerated, steady eps = -sigma).  Its exact linearisation is
    pi_tt + (1 + A'/A) pi_t + [eps k^2 + (alpha/2M^2) k^4 - J k^2]/A pi = 0 with k = k_comoving / a and A = P' + 2 X P'' = 2 + 3 eps; the
    mode equation of ``ghost_condensate.mode_growth`` (Round 6d) uses A = 2 + eps and no A'/A, an O(eps) difference (< 1% for |eps| < 3e-3).
    Pseudo-spectral in x (periodic, comoving), RK4 in t.  Dimensionless inputs: mu = M r_N, alpha, eps0, J = 4 pi alpha A (M/m_P)^2.

(2) Non-spherical perturbations: the l-mode master equation of a massless test field on the simulated spherical background,
        psi_uv = (r_uv / r) psi - (e^{2s}/4) [l(l+1)/r^2 + kappa 6 m / r^3 ... ] psi,
    with kappa = 0 (scalar proxy) or the Regge-Wheeler-type term (kappa = 1, odd-parity gravitational proxy, with the Misner-Sharp mass m).
    Integrated with the same diamond scheme on the stored grid; returns the largest amplification of each l.
"""

from __future__ import annotations

import numpy as np

from .throat_evolution import njit


def _y_from_Q(Q, px2a, a, y0):
    """Solve 2 a (y^2 - px2a - 1) y = Q for y near y0 (Newton, vectorised); px2a = a^{-2} pi_x^2."""
    y = y0.copy()
    for _ in range(30):
        f = 2 * a * (y * y - px2a - 1) * y - Q
        fp = 2 * a * (3 * y * y - px2a - 1)
        dy = f / fp
        y -= dy
        if np.max(np.abs(dy)) < 1e-15:
            break
    return y


def evolve_1p1(mu, alpha, eps0, J=0.0, sigma=0.0, L=1.0, N=512, t_end=10.0, seed_amp=1e-8, seed=0, dt=None, n_out=50, k_seed=None):
    """Returns times, the rms and max of the perturbation pi - <pi>, the background eps(t), and a breakdown flag (cubic without root, i.e.
    X crosses zero: a caustic)."""
    rng = np.random.default_rng(seed)
    k = 2 * np.pi * np.fft.rfftfreq(N, d=L / N)
    C = alpha / mu**2
    # initial perturbation: band-limited random field (all comoving modes up to N/4)
    phase = rng.uniform(0, 2 * np.pi, len(k))
    amp = np.where((k > 0) & (k < (k.max() / 2 if k_seed is None else k_seed)), 1.0, 0.0)
    pi = np.fft.irfft(amp * np.exp(1j * phase), n=N)
    pi *= seed_amp / max(np.std(pi), 1e-300)
    a = 1.0
    A = 2.0 + 3.0 * eps0
    y = np.sqrt(1.0 + eps0) * np.ones(N)                 # background X - 1 = eps0 exactly
    X1 = eps0
    Q = 2 * a * X1 * y + 0 * pi
    t = 0.0
    if dt is None:
        dt = 0.5 / (np.sqrt(C / A) * k.max() ** 2 + 1.0)

    def deriv(pi, Q, t, y_guess):
        a = np.exp(t)
        pk = np.fft.rfft(pi)
        px = np.fft.irfft(1j * k * pk, n=N)
        pxx = np.fft.irfft(-(k**2) * pk, n=N)
        pxxxx = np.fft.irfft(k**4 * pk, n=N)
        px2a = px * px / (a * a)
        if np.any(2 * a * (3 * y_guess**2 - px2a - 1) <= 0):
            return None
        y = _y_from_Q(Q, px2a, a, y_guess)
        X1 = y * y - px2a - 1
        flux = 2 * X1 * px / a
        dflux = np.fft.irfft(1j * k * np.fft.rfft(flux), n=N)
        dQ = dflux - 2 * J * pxx / a - C * pxxxx / a**3 - 2 * a * sigma
        return y - 1, dQ, y
    out = dict(t=[], rms=[], max=[], eps=[], broke=False)
    n_steps = int(np.ceil(t_end / dt))
    every = max(1, n_steps // n_out)
    for step in range(n_steps + 1):
        if step % every == 0:
            d = pi - pi.mean()
            out["t"].append(t)
            out["rms"].append(float(np.std(d)))
            out["max"].append(float(np.max(np.abs(d))))
            out["eps"].append(float(np.mean(Q) / (2 * np.exp(t))))
        k1 = deriv(pi, Q, t, y)
        if k1 is None:
            out["broke"] = True
            break
        k2 = deriv(pi + 0.5 * dt * k1[0], Q + 0.5 * dt * k1[1], t + 0.5 * dt, k1[2])
        k3 = deriv(pi + 0.5 * dt * k2[0], Q + 0.5 * dt * k2[1], t + 0.5 * dt, k2[2]) if k2 is not None else None
        k4 = deriv(pi + dt * k3[0], Q + dt * k3[1], t + dt, k3[2]) if k3 is not None else None
        if k4 is None or not np.all(np.isfinite(k4[0])):
            out["broke"] = True
            break
        pi = pi + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        Q = Q + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        y = k4[2]
        t += dt
    for key in ("t", "rms", "max", "eps"):
        out[key] = np.array(out[key])
    out["amplification"] = float(np.max(out["rms"]) / out["rms"][0]) if len(out["rms"]) else np.nan
    return out


def _background(eps0, sigma, t):
    """Exact homogeneous background of the 1+1 EFT: w = eps sqrt(1 + eps) obeys d(a w)/dt = -a sigma, so w = -sigma + (w0 + sigma) e^{-t};
    returns eps(t), the kinetic coefficient A = 2(1 + eps) + eps = 2 + 3 eps (from P' + 2 X P'') and dA/dt."""
    w0 = eps0 * np.sqrt(1 + eps0)
    w = -sigma + (w0 + sigma) * np.exp(-t)
    wd = -(w0 + sigma) * np.exp(-t)
    eps = np.asarray(w, float).copy()
    for _ in range(40):
        eps = eps - (eps * np.sqrt(1 + eps) - w) / (np.sqrt(1 + eps) + eps / (2 * np.sqrt(1 + eps)))
    epsd = wd / (np.sqrt(1 + eps) + eps / (2 * np.sqrt(1 + eps)))
    return eps, 2 + 3 * eps, 3 * epsd


def linear_amplification(mu, alpha, eps0, J=0.0, sigma=0.0, L=1.0, N=512, t_end=10.0):
    """The same set of comoving modes evolved with the linear equation (eps(t) as in the nonlinear background): max over modes."""
    from scipy.integrate import solve_ivp
    k = 2 * np.pi * np.fft.rfftfreq(N, d=L / N)[1: N // 4]
    C = alpha / mu**2
    best = 1.0
    for kk in k[:: max(1, len(k) // 60)]:
        def f(t, z, kk=kk):
            eps, A, Ad = _background(eps0, sigma, t)
            q2 = (kk * np.exp(-t)) ** 2
            return [z[1], -(1 + Ad / A) * z[1] - (eps * q2 + 0.5 * C * q2 * q2 - J * q2) / A * z[0]]
        s = solve_ivp(f, (0, t_end), [1.0, 0.0], rtol=1e-8, atol=1e-12, max_step=0.05)
        best = max(best, float(np.max(np.abs(s.y[0]))))
    return best


# ------------------------------------------------------------------ non-spherical test-field perturbations on the simulated background
@njit(cache=True)
def _lmode(R, S, h, ell, kappa, psi_a, psi_b, Uv):
    """Diamond integration of psi_uv = (r_uv/r) psi - (e^{2s}/4)(l(l+1)/r^2 - kappa 6 M/r^3) psi on a stored background (r_uv and the
    Misner-Sharp m from finite differences of R, S).  M = m - (4 pi/3) U r^3 is the mass parameter without the local vacuum energy U
    (Uv, per cell; the Regge-Wheeler term of Schwarzschild-de Sitter carries M, not m).  psi_a on u = 0, psi_b on v = 0."""
    M, N = R.shape
    P = np.full((M, N), np.nan)
    for j in range(N):
        P[0, j] = psi_a[j]
    for i in range(1, M):
        P[i, 0] = psi_b[i]
        for j in range(1, N):
            rN, rE, rW, rS = R[i, j], R[i - 1, j], R[i, j - 1], R[i - 1, j - 1]
            if not (rN == rN and rE == rE and rW == rW and rS == rS and P[i - 1, j] == P[i - 1, j] and P[i, j - 1] == P[i, j - 1]):
                break
            rc = 0.25 * (rN + rE + rW + rS)
            sc = 0.25 * (S[i, j] + S[i - 1, j] + S[i, j - 1] + S[i - 1, j - 1])
            ruv = (rN - rE - rW + rS) / (h * h)
            ru = (rN + rW - rE - rS) / (2 * h)
            rv = (rN + rE - rW - rS) / (2 * h)
            a2 = np.exp(2.0 * sc)
            m = 0.5 * rc * (1.0 + 4.0 * ru * rv / a2) - 4.0 * np.pi / 3.0 * Uv[i, j] * rc * rc * rc
            V = ruv / rc - 0.25 * a2 * (ell * (ell + 1.0) / (rc * rc) - kappa * 6.0 * m / (rc * rc * rc))
            # psi_N - psi_E - psi_W + psi_S = h^2 V psi_c, psi_c = average (implicit in psi_N)
            Pc3 = P[i - 1, j] + P[i, j - 1] + P[i - 1, j - 1]
            num = P[i - 1, j] + P[i, j - 1] - P[i - 1, j - 1] + h * h * V * 0.25 * Pc3
            den = 1.0 - h * h * V * 0.25
            P[i, j] = num / den
    return P


def lmode_amplification(res, ells=(1, 2, 3, 5, 10, 20), kappa=0.0, width=2.0, v0=2.0, pot=None, r_cut=None, power=1):
    """Largest |psi| / r^power relative to its largest value on the data, for a Gaussian l-pulse placed on the outgoing data ray near
    the corner (and zero on the ingoing ray), for each l.  power = 1: the scalar amplitude phi = psi / r; power = 3: the odd-parity
    curvature (perturbed Weyl scalar ~ psi / r^3) of the Regge-Wheeler proxy.  pot: the scalar potential (its U(phi) is removed from the
    Regge-Wheeler mass); r_cut: only cells with r > r_cut count (None: the whole domain, including the crushing region r -> 0)."""
    R, S, h = res["r"], res["s"], res["h"]
    Uv = np.zeros_like(R) if pot is None else np.nan_to_num(pot.U(res["phi"]))
    v = np.arange(R.shape[1]) * h
    pa = np.exp(-((v - v0) / width) ** 2) * R[0]          # psi = r phi with a unit phi pulse
    pb = np.zeros(R.shape[0])
    pb[0] = pa[0]
    out = {}
    for ell in ells:
        P = _lmode(np.ascontiguousarray(R), np.ascontiguousarray(S), h, float(ell), float(kappa), pa, pb, np.ascontiguousarray(Uv))
        F = np.abs(P) / np.abs(R) ** power / np.max(np.abs(pa) / R[0] ** power)
        if r_cut is not None:
            F = np.where(R > r_cut, F, np.nan)
        out[int(ell)] = float(np.nanmax(F))                            # amplification of phi = psi / r
    return out


def linear_rms_1p1(mu, alpha, eps0, J=0.0, sigma=0.0, L=1.0, N=512, t_end=10.0, seed_amp=1e-8, seed=0, k_seed=None, n_out=50,
                   dt_max=0.02, wdt=0.2):
    """The same initial field as ``evolve_1p1`` evolved mode by mode with the linear equation (all modes at once, RK4 with
    omega_max dt <= wdt); rms(t) by Parseval."""
    rng = np.random.default_rng(seed)
    k = 2 * np.pi * np.fft.rfftfreq(N, d=L / N)
    phase = rng.uniform(0, 2 * np.pi, len(k))
    amp = np.where((k > 0) & (k < (k.max() / 2 if k_seed is None else k_seed)), 1.0, 0.0)
    pi0 = np.fft.irfft(amp * np.exp(1j * phase), n=N)
    c0 = np.fft.rfft(pi0 * seed_amp / max(np.std(pi0), 1e-300))
    keep = np.abs(c0) > 0
    k, w = k[keep], 2 * np.abs(c0[keep]) ** 2
    C = alpha / mu**2

    def w2(t):
        eps, A, _ = _background(eps0, sigma, t)
        q2 = (k * np.exp(-t)) ** 2
        return (eps * q2 + 0.5 * C * q2 * q2 - J * q2) / A

    def f(t, z, y):
        _, A, Ad = _background(eps0, sigma, t)
        return y, -(1 + Ad / A) * y - w2(t) * z

    ts = np.linspace(0, t_end, n_out + 1)
    z, y, t = np.ones_like(k), np.zeros_like(k), 0.0
    pw = np.zeros(len(ts))
    pw[0] = np.sum(w)
    for n in range(1, len(ts)):
        while t < ts[n] - 1e-14:
            dt = min(dt_max, wdt / np.sqrt(max(np.max(np.abs(w2(t))), 1e-300)), ts[n] - t)
            k1z, k1y = f(t, z, y)
            k2z, k2y = f(t + dt / 2, z + dt / 2 * k1z, y + dt / 2 * k1y)
            k3z, k3y = f(t + dt / 2, z + dt / 2 * k2z, y + dt / 2 * k2y)
            k4z, k4y = f(t + dt, z + dt * k3z, y + dt * k3y)
            z = z + dt / 6 * (k1z + 2 * k2z + 2 * k3z + k4z)
            y = y + dt / 6 * (k1y + 2 * k2y + 2 * k3y + k4y)
            t += dt
        pw[n] = np.sum(w * z**2)
    rms = np.sqrt(pw) / N
    return dict(t=ts, rms=rms, amplification=float(rms.max() / rms[0]))


def growth_exponent(mu, alpha, eps0, J=0.0, sigma=0.0, t_end=10.0, n_k=400, n_t=4000):
    """Exponential growth only: for each comoving k, G_k = int max(0, -1/2 + sqrt(1/4 - omega_k^2(t))) dt, the WKB growth exponent of
    pi_tt + (1 + A'/A) pi_t + omega^2 pi = 0 with omega^2 = [eps(t) q^2 + (C/2) q^4 - J q^2] / A, q = k e^{-t}, C = alpha / mu^2,
    A = 2 + 3 eps (the exact linearisation of ``evolve_1p1``).  Stable modes
    (omega^2 > 0) contribute nothing; their bounded algebraic growth (at most omega_0^{1/4}, absent for vacuum-normalised seeds) is not
    counted.  Returns max_k G_k and the k that attains it."""
    C = alpha / mu**2
    eps_max = max(abs(eps0), abs(sigma)) + abs(J)
    qb = np.sqrt(2 * eps_max / C) if eps_max > 0 else 1.0                # edge of the unstable band in physical q
    k = np.geomspace(1e-3 * qb, 10 * qb * np.exp(t_end), n_k)
    t = np.linspace(0, t_end, n_t + 1)
    eps_t, A, Ad = _background(eps0, sigma, t)
    g = 0.5 * (1 + Ad / A)                                               # half the friction
    q2 = (k[None, :] * np.exp(-t)[:, None]) ** 2
    w2 = (eps_t[:, None] * q2 + 0.5 * C * q2 * q2 - J * q2) / A[:, None]
    lam = np.maximum(0.0, -g[:, None] + np.sqrt(np.maximum(g[:, None] ** 2 - w2, 0.0)))
    G = np.trapezoid(lam, t, axis=0)
    i = int(np.argmax(G))
    return dict(G=float(G[i]), k=float(k[i]), valid=bool(np.all(A > 0)))


def growth_exponent_history(M3rN, r_N, alpha, x, nu, rr, decay_after=True, t_end=9.6, n_k=60, n_l=60, n_t=2000, A=2.0):
    """The full linear stability measure on a measured history: modes along the expanding dS_2 direction (comoving k, redshifting) and
    along the S^2 (L = sqrt(l(l+1)), physical q = L / r(t), not redshifted while r ~ r_N), both at once:
        omega^2 = [eps(t) q^2 + (C/2) q^4 - J q^2] / A,   q^2 = k^2 e^{-2t} + L^2 / rho(t)^2,   units r_N = 1,
    with eps = nu r_N^2 / (2 mu^4) from the measured charge nu(x) = n r_N^2 at x = t / r_N, rho = r / r_N, mu = (M^3 r_N r_N^2)^{1/3},
    C = alpha / mu^2 and the Jeans term J = 8 pi alpha (mu / r_N)^2.  After the measured history ends the charge decays as e^{-(x - x_m)}
    (decay_after) or stays (self-generated), and rho grows as e^{(x - x_m)/sqrt(3)} in the self-generated case (Round 6e convention).
    Returns max over (k, L) of G = int max(0, -1/2 + sqrt(1/4 - omega^2)) dt.  (A = 2 + 3 eps + O(eps^2) differs from 2 by < 1% here.)"""
    mu = (M3rN * r_N**2) ** (1 / 3)
    C = alpha / mu**2
    J = 8 * np.pi * alpha * (mu / r_N) ** 2
    xm = float(np.max(x))
    t = np.linspace(0, t_end, n_t + 1)
    if decay_after:
        nut = np.where(t <= xm, np.interp(t, x, nu), nu[-1] * np.exp(-(t - xm)))
        rho = np.where(t <= xm, np.interp(t, x, rr), rr[-1])
    else:
        nut = np.where(t <= xm, np.interp(t, x, nu), nu[-1])
        rho = np.where(t <= xm, np.interp(t, x, rr), rr[-1] * np.exp((t - xm) / np.sqrt(3)))
    eps = nut * r_N**2 / (2 * mu**4)
    qb = np.sqrt(2 * (np.max(np.abs(eps)) + J) / C)
    ks = np.concatenate([[0.0], np.geomspace(1e-2 * qb, 10 * qb * np.exp(t_end), n_k)])
    Ls = np.concatenate([[0.0], np.geomspace(1.0, 10 * qb * np.max(rho), n_l)])
    best, arg = 0.0, (0.0, 0.0)
    for L in Ls:
        q2 = (ks[None, :] * np.exp(-t)[:, None]) ** 2 + (L / rho)[:, None] ** 2
        w2 = (eps[:, None] * q2 + 0.5 * C * q2 * q2 - J * q2) / A
        lam = np.maximum(0.0, -0.5 + np.sqrt(np.maximum(0.25 - w2, 0.0)))
        G = np.trapezoid(lam, t, axis=0)
        i = int(np.argmax(G))
        if G[i] > best:
            best, arg = float(G[i]), (float(ks[i]), float(L))
    return dict(G=best, k=arg[0], L=arg[1], eps_max=float(np.max(np.abs(eps))))
