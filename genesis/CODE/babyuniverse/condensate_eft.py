"""Round 6f (item 3): the ghost condensate beyond its IR fluid limit, and non-spherical perturbations of the tipping geometry.

(1) Nonlinear 1+1 evolution of the full condensate EFT along the expanding direction of the Nariai band (dS_2 x S^2, flat dS_2 slicing,
    a = e^{t/r_N}, units r_N = 1).  phi = M^2 (t + pi), X = (1 + pi_t)^2 - a^{-2} pi_x^2, P(X) = (X - 1)^2 / 2 (P'' = 1), plus the leading
    higher-derivative term -(alpha/2M^2)(lap pi)^2 of the ghost condensate and the condensate's Jeans coupling J (from mixing with gravity):
        Q = 2 a P'(X) (1 + pi_t)                (charge per comoving length, units M^4)
        Q_t = d_x[2 a^{-1} P'(X) pi_x] - 2 a^{-1} J pi_xx - (alpha/M^2) a^{-3} pi_xxxx - 2 a sigma,
    and pi_t from Q by the cubic 2 a ((1 + pi_t)^2 - a^{-2} pi_x^2 - 1)(1 + pi_t) = Q.  sigma is the false-vacuum tilt source (sigma = 0:
    the NEC phase dilutes as e^{-t}; sigma > 0: it is regenerated, steady eps = -sigma).  The linearisation is exactly the mode equation of
    ``ghost_condensate.mode_growth``: pi_tt + pi_t + [eps k^2 + (alpha/2M^2) k^4 - J k^2]/A pi = 0 with k = k_comoving / a and A = 2 + eps.
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
    A = 2.0 + eps0
    y = np.ones(N)
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


def linear_amplification(mu, alpha, eps0, J=0.0, sigma=0.0, L=1.0, N=512, t_end=10.0):
    """The same set of comoving modes evolved with the linear equation (eps(t) as in the nonlinear background): max over modes."""
    from scipy.integrate import solve_ivp
    k = 2 * np.pi * np.fft.rfftfreq(N, d=L / N)[1: N // 4]
    C = alpha / mu**2
    A = 2.0 + eps0
    eps_t = (lambda t: eps0 * np.exp(-t)) if sigma == 0 else (lambda t: -sigma + (eps0 + sigma) * np.exp(-t))
    best = 1.0
    for kk in k[:: max(1, len(k) // 60)]:
        f = lambda t, z: [z[1], -z[1] - (eps_t(t) * (kk * np.exp(-t)) ** 2 + 0.5 * C * (kk * np.exp(-t)) ** 4 - J * (kk * np.exp(-t)) ** 2) / A * z[0]]
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


def lmode_amplification(res, ells=(1, 2, 3, 5, 10, 20), kappa=0.0, width=2.0, v0=2.0, pot=None, r_cut=None):
    """Largest |phi| = |psi/r| relative to the initial pulse amplitude, for a unit Gaussian l-pulse placed on the outgoing data ray near
    the corner (and zero on the ingoing ray), for each l.  pot: the scalar potential (its U(phi) is removed from the Regge-Wheeler mass);
    r_cut: only cells with r > r_cut count (r_cut = None: the whole domain, including the crushing region r -> 0 at the singularities)."""
    R, S, h = res["r"], res["s"], res["h"]
    Uv = np.zeros_like(R) if pot is None else np.nan_to_num(pot.U(res["phi"]))
    v = np.arange(R.shape[1]) * h
    pa = np.exp(-((v - v0) / width) ** 2) * R[0]          # psi = r phi with a unit phi pulse
    pb = np.zeros(R.shape[0])
    pb[0] = pa[0]
    out = {}
    for ell in ells:
        P = _lmode(np.ascontiguousarray(R), np.ascontiguousarray(S), h, float(ell), float(kappa), pa, pb, np.ascontiguousarray(Uv))
        F = np.abs(P / R)
        if r_cut is not None:
            F = np.where(R > r_cut, F, np.nan)
        out[int(ell)] = float(np.nanmax(F))                            # amplification of phi = psi / r
    return out


def linear_rms_1p1(mu, alpha, eps0, J=0.0, sigma=0.0, L=1.0, N=512, t_end=10.0, seed_amp=1e-8, seed=0, k_seed=None, n_out=50):
    """The same initial field as ``evolve_1p1`` evolved mode by mode with the linear equation; rms(t) by Parseval."""
    from scipy.integrate import solve_ivp
    rng = np.random.default_rng(seed)
    k = 2 * np.pi * np.fft.rfftfreq(N, d=L / N)
    phase = rng.uniform(0, 2 * np.pi, len(k))
    amp = np.where((k > 0) & (k < (k.max() / 2 if k_seed is None else k_seed)), 1.0, 0.0)
    pi0 = np.fft.irfft(amp * np.exp(1j * phase), n=N)
    c0 = np.fft.rfft(pi0 * seed_amp / max(np.std(pi0), 1e-300))
    C = alpha / mu**2
    A = 2.0 + eps0
    eps_t = (lambda t: eps0 * np.exp(-t)) if sigma == 0 else (lambda t: -sigma + (eps0 + sigma) * np.exp(-t))
    ts = np.linspace(0, t_end, n_out + 1)
    pw = np.zeros(len(ts))
    for kk, ck in zip(k, c0):
        if abs(ck) == 0:
            continue
        f = lambda t, z: [z[1], -z[1] - (eps_t(t) * (kk * np.exp(-t)) ** 2 + 0.5 * C * (kk * np.exp(-t)) ** 4 - J * (kk * np.exp(-t)) ** 2) / A * z[0]]
        sol = solve_ivp(f, (0, t_end), [1.0, 0.0], t_eval=ts, rtol=1e-9, atol=1e-12, max_step=0.02)
        pw += 2 * np.abs(ck) ** 2 * sol.y[0] ** 2
    rms = np.sqrt(pw) / N
    return dict(t=ts, rms=rms, amplification=float(rms.max() / rms[0]))
