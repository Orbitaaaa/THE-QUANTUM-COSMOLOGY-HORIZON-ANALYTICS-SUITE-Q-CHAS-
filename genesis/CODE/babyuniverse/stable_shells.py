"""Round 6i, routes 2, 3, 5, 7: long-lived thin-shell objects as daughter precursors (monopole core, breathing bubble, resonant tunnelling,
four-form membrane).  G = c = hbar = 1.

One model covers all four.  A spherical shell of radius R(tau) separates a de Sitter interior f_in = 1 - chi^2 R^2 from a
Reissner-Nordstrom-de Sitter exterior f_out = 1 - 2M/R + Q^2/R^2 - eta^2 R^2.  The shell's total surface energy is
m_s(R) = 4 pi R^2 sigma(R) = k_w R^2 + m_d + k_1 R + k_3 / R, made of a domain wall (k_w = 4 pi sigma_w, p = -sigma), dust (rest mass m_d,
p = 0), and optional string-like (k_1) and radiation-like (k_3) gases.  Israel junction: beta_in - beta_out = m_s / R, beta = eps sqrt(f + Rdot^2),
hence Rdot^2 + V(R) = 0 with
    V = f_out - beta_out^2,  beta_out = (R (f_in - f_out) - m_s^2 / R) / (2 m_s),  beta_in = beta_out + m_s / R.
A hidden daughter is a shell on the inflating branch with beta_out < 0 (beyond the exterior Einstein-Rosen bridge); beta_out > 0 at large R
is a bubble that grows into the parent.

Semiclassical exponent between two turning points (the standard thin-shell WKB / Fischler-Morgan-Polchinski form, as in tunnelling.py):
    B = 2 int R [Theta_out - Theta_in] dR,  Theta = atan2(sqrt V, beta)  (analytic prescription);
Bohr-Sommerfeld phase of a classically allowed well S = 2 int R [atanh(Rdot/beta_out) - atanh(Rdot/beta_in)] dR (static regions only).
Assessment: the canonical momentum of the constant-tension shell is used unchanged for shells with dust or gas components; for dust it agrees
with the reduced Hamiltonian of Kraus and Wilczek, and it is the convention of the breathing-bubble literature.

Results (ROADMAP_ROUND6I.md, sections 2, 3, 5, 7):
* ``breathing_random_scan`` Proposition-N 6i-4 (breathing bubble, Guendelman & Sakai 2008).  Bounded wells need k_w m_d < 0.  The hidden-daughter
                           case (beta_out < 0 on the large-R branch) needs k_w > 0, hence m_d < 0: negative-energy dust, which violates the NEC
                           (sigma + p = sigma_d < 0); k_w < 0 gives a bubble that eats the parent (V is even in m_s, beta_out odd).  In 6000
                           random shells all 713 wells contain the radius where m_s = 0, which the shell crosses at the speed of light: the
                           change of side of the bridge happens there.
* ``breathing_barrier_vs_mass`` B = B_hat (l_f/l_P)^2 with B_hat ~ 0.5-0.9 along typical families; B_hat -> 0 only at merger points.
* ``monopole_shell``       Proposition-N 6i-5 (monopole core, thin shell after Sakai, Nakao, Ishihara & Kobayashi 2006): static stable cores
                           exist below eta_c = 0.25-0.26 (lam = 1, e = 0.5), 0.32-0.33 (e = 1), 0.39-0.40 (lam = 0.1) m_P (literature: O(0.3)
                           m_P).  At the monopole's own mass B ~ 650/lam (eta-independent; ~1.5 x the flat thin-wall CdL action 421/lam of the
                           same wall), but M << Q there: the exterior has no bridge and the beta_out < 0 branch is disconnected from our
                           asymptotic region (forbidden by Gauss's law).  The connected channel needs M > Q ~ 7 m_P (a magnetic black hole)
                           and has B ~ 3e3-3e4.  Negative dust on the shell lowers B by <= 20% before the stable core disappears.
* ``resonance_scan``       Proposition 6i-6 (resonant tunnelling): two-well shell potentials exist (charge + gas components; 13 clean cases in
                           20000 samples).  On resonance the exponent is max(B_2, B_1 - B_2) >= B_2: the barrier next to the inflating
                           branch is never removed.
* ``membrane_branches``    Proposition 6i-7 (four-form membrane): a Brown-Teitelboim membrane is a constant-tension wall (NEC saturated) and
                           the branches are those of Farhi-Guth: the hidden inflating branch exists for M < M_cr but classically starts at the
                           white-hole singularity (Theorem 6h-1, Theorem 6i-11).  B(M -> 0) = 1.03 x the flat thin-wall CdL action: the "no
                           free lunch" bound (R10, R37) holds.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar


@dataclass(frozen=True)
class Shell:
    chi: float                 # interior Hubble rate
    M: float                   # exterior (ADM) mass
    k_w: float = 0.0           # 4 pi sigma_wall
    m_d: float = 0.0           # dust rest mass on the shell
    Q2: float = 0.0            # exterior charge^2 (magnetic monopole: 4 pi / e^2)
    eta: float = 0.0           # exterior Hubble rate
    k_1: float = 0.0           # string-gas term (m_s ~ R)
    k_3: float = 0.0           # radiation-gas term (m_s ~ 1/R)

    def ms(self, R):
        R = np.asarray(R, float)
        return self.k_w * R**2 + self.m_d + self.k_1 * R + self.k_3 / R

    def f_in(self, R):
        return 1 - self.chi**2 * np.asarray(R, float) ** 2

    def f_out(self, R):
        R = np.asarray(R, float)
        return 1 - 2 * self.M / R + self.Q2 / R**2 - self.eta**2 * R**2

    def beta_out(self, R):
        R = np.asarray(R, float)
        m = self.ms(R)
        return (R * (self.f_in(R) - self.f_out(R)) - m**2 / R) / (2 * m)

    def beta_in(self, R):
        return self.beta_out(R) + self.ms(R) / np.asarray(R, float)

    def V(self, R):
        return self.f_out(R) - self.beta_out(R) ** 2

    # ------------------------------------------------------------------ structure of V
    def turning_points(self, R_lo=1e-3, R_hi=50.0, n=20000):
        """Sign changes of V on a log grid, refined by brentq; excludes zeros of m_s (where V -> -inf)."""
        R = np.geomspace(R_lo, R_hi, n)
        with np.errstate(all="ignore"):
            v = self.V(R)
            m = self.ms(R)
        idx = np.where(np.isfinite(v[:-1]) & np.isfinite(v[1:]) & (v[:-1] * v[1:] < 0) & (m[:-1] * m[1:] > 0))[0]
        return [brentq(lambda x: float(self.V(x)), R[i], R[i + 1], xtol=1e-14, rtol=1e-13) for i in idx]

    def regions(self, R_lo=1e-3, R_hi=50.0):
        """Alternating classically allowed (V < 0) and forbidden intervals between turning points."""
        tp = [R_lo] + self.turning_points(R_lo, R_hi) + [R_hi]
        regs = []
        for a, b in zip(tp[:-1], tp[1:]):
            mid = np.sqrt(a * b)
            regs.append(dict(R=(a, b), allowed=bool(self.V(mid) < 0), beta_out=float(self.beta_out(mid)),
                             beta_in=float(self.beta_in(mid)), ms=float(self.ms(mid))))
        return regs

    # ------------------------------------------------------------------ semiclassical integrals
    def B(self, Ra, Rb):
        """Euclidean exponent across the forbidden interval (Ra, Rb): 2 int R (Theta_out - Theta_in) dR."""
        def integrand(R):
            y = np.sqrt(max(float(self.V(R)), 0.0))
            return R * (np.arctan2(y, float(self.beta_out(R))) - np.arctan2(y, float(self.beta_in(R))))
        return 2 * quad(integrand, Ra, Rb, limit=400, epsabs=0, epsrel=1e-10)[0]

    def phase(self, Ra, Rb):
        """Bohr-Sommerfeld phase over a classically allowed interval inside static regions: 2 int R [atanh(Rdot/b_out) - atanh(Rdot/b_in)]."""
        def integrand(R):
            rd = np.sqrt(max(-float(self.V(R)), 0.0))
            bo, bi = float(self.beta_out(R)), float(self.beta_in(R))
            return R * (np.arctanh(np.clip(rd / bo, -1 + 1e-15, 1 - 1e-15)) - np.arctanh(np.clip(rd / bi, -1 + 1e-15, 1 - 1e-15)))
        return abs(2 * quad(integrand, Ra, Rb, limit=400, epsabs=0, epsrel=1e-10)[0])

    def well_frequency(self, R_min):
        """Small-oscillation proper angular frequency about a minimum of V: omega^2 = V''(R_min) / 2."""
        h = 1e-4 * R_min
        d2 = (float(self.V(R_min + h)) - 2 * float(self.V(R_min)) + float(self.V(R_min - h))) / h**2
        return float(np.sqrt(max(d2, 0.0) / 2))


# --------------------------------------------------------------------------- route 3: the breathing bubble (Guendelman & Sakai 2008)


def classify_breathing(sh: Shell, R_lo=1e-3, R_hi=50.0) -> dict:
    """A bounded allowed well, separated by barriers from the collapse region and from the large-R branch.  Flags: whether the well stays in
    static regions (f_in, f_out > 0) and whether it contains the radius where the total surface energy m_s vanishes (there V -> -inf and the
    shell moves at the speed of light).  Outcome of the large-R branch: hidden daughter (beta_out < 0) or parent-eating bubble (beta_out > 0)."""
    regs = sh.regions(R_lo, R_hi)
    out = dict(breathing=False)
    for i in range(1, len(regs) - 2):
        w, bar, inf = regs[i], regs[i + 1], regs[-1]
        if not (w["allowed"] and not regs[i - 1]["allowed"] and not bar["allowed"] and inf["allowed"] and inf is regs[i + 2]):
            continue
        a, b = w["R"]
        Rs = np.geomspace(a, b, 400)[1:-1]
        m = sh.ms(Rs)
        crosses = bool(m.min() < 0 < m.max())
        static = bool(np.all(sh.f_in(Rs) > 0) and np.all(sh.f_out(Rs) > 0))
        out = dict(breathing=True, crosses_zero_energy=crosses, static=static, well=(a, b), barrier=bar["R"], inflating_from=inf["R"][0],
                   beta_out_large_R=float(sh.beta_out(R_hi)),
                   outcome="hidden daughter" if sh.beta_out(R_hi) < 0 else "parent-eating bubble",
                   B=float(sh.B(*bar["R"])), E_vac_outer=float(sh.chi**2 * b**3 / 2), E_wall_outer=float(sh.k_w * b**2), E_dust=float(sh.m_d))
        break
    return out


def breathing_family(chi=1.0, k_ws=(0.5, 1.0, 2.0, -0.5, -1.0), m_ds=None, Ms=None) -> list:
    """Scan (k_w, m_d, M) for breathing bubbles; record which signs give stable oscillation and which outcome follows."""
    m_ds = np.linspace(-0.6, 0.6, 25) if m_ds is None else m_ds
    Ms = np.linspace(0.005, 0.6, 60) if Ms is None else Ms
    rows = []
    for kw in k_ws:
        for md in m_ds:
            if md == 0:
                continue
            for M in Ms:
                c = classify_breathing(Shell(chi=chi, M=M, k_w=kw, m_d=md))
                if c.get("breathing"):
                    rows.append(dict(k_w=kw, m_d=float(md), M=float(M), **c))
    return rows


def breathing_random_scan(n=6000, seed=1, chi=1.0) -> dict:
    """Random two-component shells (k_w, m_d of either sign, log-uniform magnitudes; M log-uniform): every bounded well found, and whether
    it contains the zero of m_s."""
    rng = np.random.default_rng(seed)
    wells, crossing, signs = 0, 0, set()
    for _ in range(n):
        kw = rng.choice([-1, 1]) * 10 ** rng.uniform(-1.5, 1)
        md = rng.choice([-1, 1]) * 10 ** rng.uniform(-2.5, 0)
        M = 10 ** rng.uniform(-3, 0)
        sh = Shell(chi=chi, M=M, k_w=kw, m_d=md)
        regs = sh.regions(1e-3, 30.0)
        for j in range(1, len(regs) - 1):
            r = regs[j]
            if r["allowed"] and not regs[j - 1]["allowed"] and not regs[j + 1]["allowed"]:
                m = sh.ms(np.geomspace(*r["R"], 300))
                wells += 1
                crossing += int(m.min() < 0 < m.max())
                signs.add((int(np.sign(kw)), int(np.sign(md))))
    return dict(samples=n, wells=wells, wells_crossing_zero_energy=crossing, sign_pairs=sorted(signs))


def breathing_barrier_vs_mass(k_w, m_d, chi=1.0, n=60):
    """B(M) along a breathing family at fixed (k_w, m_d): find the mass interval with a well, and B approaching the merger point."""
    Ms = np.linspace(1e-3, 1.0, 400)
    ok = [M for M in Ms if classify_breathing(Shell(chi=chi, M=M, k_w=k_w, m_d=m_d)).get("breathing")]
    if not ok:
        return []
    lo, hi = min(ok), max(ok)
    rows = []
    for M in np.linspace(lo, hi, n):
        c = classify_breathing(Shell(chi=chi, M=M, k_w=k_w, m_d=m_d))
        if c.get("breathing"):
            rows.append(dict(M=float(M), B=c["B"], well=c["well"], barrier=c["barrier"], outcome=c["outcome"]))
    return rows


# --------------------------------------------------------------------------- route 2: the monopole core (Sakai, Nakao, Ishihara & Kobayashi 2006)


def monopole_shell(eta, lam=1.0, e=0.5, M=None, m_d=0.0) -> Shell:
    """Thin-shell monopole: interior false vacuum rho = lam eta^4 / 4 (chi^2 = 8 pi rho / 3), wall tension (sqrt2/3) sqrt(lam) eta^3
    (Higgs kink from 0 to eta), magnetic RN exterior Q^2 = 4 pi / e^2 (unit 't Hooft-Polyakov charge 4 pi / e, Heaviside-Lorentz).
    m_d < 0 models an NEC-violating admixture on the shell."""
    chi = np.sqrt(8 * np.pi / 3 * lam * eta**4 / 4)
    k_w = 4 * np.pi * np.sqrt(2) / 3 * np.sqrt(lam) * eta**3
    return Shell(chi=chi, M=0.0 if M is None else M, k_w=k_w, m_d=m_d, Q2=4 * np.pi / e**2)


def monopole_static(eta, lam=1.0, e=0.5, m_d=0.0):
    """The ground state: the smallest M for which V has a local minimum with V_min = 0 (a static shell in static regions).
    Returns None if no static stable core exists (topological inflation regime)."""
    sh0 = monopole_shell(eta, lam, e, 0.0, m_d)
    R_hub = 1 / sh0.chi
    R_lo, R_hi = 1e-3 * min(np.sqrt(sh0.Q2), R_hub), 3 * R_hub

    def vmin(M):
        sh = replace(sh0, M=M)
        Rs = np.geomspace(R_lo, min(R_hi, 0.999 / sh.chi), 4000)
        with np.errstate(all="ignore"):
            v = sh.V(Rs)
            ok = (sh.f_out(Rs) > 0) & (sh.f_in(Rs) > 0) & (sh.ms(Rs) > 0) & np.isfinite(v)
        if ok.sum() < 5:
            return None
        idx = np.where(ok)[0]
        vv = v[idx]
        # interior local minima of V restricted to static region
        loc = [k for k in range(1, len(vv) - 1) if vv[k] < vv[k - 1] and vv[k] < vv[k + 1] and idx[k + 1] - idx[k - 1] == 2]
        if not loc:
            return None
        k = min(loc, key=lambda k: vv[k])
        return float(vv[k]), float(Rs[idx[k]])

    # flat-space scale of the monopole mass ~ 4 pi eta / e; search M on a log grid for the static (V_min = 0) state
    Ms = np.geomspace(1e-4 * 4 * np.pi * eta / e, 50 * 4 * np.pi * eta / e, 300)
    prev = None
    for M in Ms:
        r = vmin(M)
        if r is None:
            prev = None
            continue
        if prev is not None and prev[1][0] > 0 >= r[0]:
            Ms_ = brentq(lambda m: (vmin(m) or (np.nan,))[0], prev[0], M, xtol=1e-12 * M)
            v, R = vmin(Ms_)
            return dict(M=float(Ms_), R=R, chi_R=float(R * sh0.chi), V_min=v)
        prev = (M, r)
    return None


def monopole_critical_eta(lam=1.0, e=0.5, m_d_over_M=0.0, eta_grid=None) -> dict:
    """Largest eta with a static stable core (thin-shell model).  m_d_over_M: negative dust on the shell as a fraction of the flat-space
    monopole mass 4 pi eta / e (the NEC-field admixture)."""
    eta_grid = np.geomspace(1e-3, 3.0, 60) if eta_grid is None else eta_grid
    last = None
    for eta in eta_grid:
        md = m_d_over_M * 4 * np.pi * eta / e
        st = monopole_static(eta, lam, e, md)
        if st is None:
            return dict(eta_c_between=(float(last) if last is not None else None, float(eta)), lam=lam, e=e, m_d_over_M=m_d_over_M)
        last = eta
    return dict(eta_c_between=(float(last), None), lam=lam, e=e, m_d_over_M=m_d_over_M)


def monopole_tunnelling(eta, lam=1.0, e=0.5, m_d=0.0) -> dict:
    """B from the static core (its own mass) to the inflating branch, when that branch exists at the same M (analytic prescription;
    eps_out flips inside the barrier, so the Euclidean region is a pseudomanifold as in Farhi-Guth-Guven, R4)."""
    st = monopole_static(eta, lam, e, m_d)
    if st is None:
        return dict(static=False)
    sh = replace(monopole_shell(eta, lam, e, m_d=m_d), M=st["M"])
    R_hub = 1 / sh.chi
    tps = [x for x in sh.turning_points(1e-4 * min(np.sqrt(sh.Q2), R_hub), 50 * R_hub, n=60000) if x > st["R"] * 1.01]
    if not tps:
        return dict(static=True, M=st["M"], inflating_branch=False)
    Rb = tps[0]
    Ra = max(x for x in sh.turning_points(1e-4 * min(np.sqrt(sh.Q2), R_hub), 50 * R_hub, n=60000) if x < Rb)
    Ra = max(Ra, st["R"])
    return dict(static=True, M=st["M"], R_static=st["R"], R_inflate=float(Rb), B=float(sh.B(Ra, Rb)),
                beta_out_inflating=float(sh.beta_out(Rb * 1.01)))


# --------------------------------------------------------------------------- route 5: resonant tunnelling through a two-well shell potential


def double_barrier_rate(B1, B2, omega1=1.0, omega2=1.0, resonant=True) -> float:
    """WKB decay rate of a level in well 1 through barrier 1, well 2 and barrier 2 (T = e^{-B}).  Off resonance: omega T1 T2 / 4.
    On resonance (levels aligned): min(omega2 T2 / 2, omega1 T1 / T2) (coherent sharing vs incoherent hopping)."""
    T1, T2 = np.exp(-B1), np.exp(-B2)
    if not resonant:
        return float(omega1 * T1 * T2 / 4)
    return float(min(omega2 * T2 / 2, omega1 * T1 / T2))


def resonance_exponent(B1, B2) -> dict:
    """Exponent of the decay rate: B1 + B2 off resonance, max(B2, B1 - B2) on resonance."""
    return dict(off=float(B1 + B2), on=float(max(B2, B1 - B2)), gain=float(B1 + B2 - max(B2, B1 - B2)))


def find_two_well(n=20000, seed=7, chi=1.0, R_hi=30.0) -> list:
    """Random wall + dust + charge + string + radiation shells whose potential ends in well, barrier, well, barrier, inflating branch with
    m_s > 0 throughout the last five regions (no zero-energy crossing) and a hidden-daughter large-R branch."""
    rng = np.random.default_rng(seed)
    found = []
    for _ in range(n):
        sh = Shell(chi=chi, M=10 ** rng.uniform(-3, 0.3), k_w=rng.choice([-1, 1]) * 10 ** rng.uniform(-1.5, 1),
                   m_d=rng.choice([-1, 1, 0]) * 10 ** rng.uniform(-3, 0), Q2=rng.choice([0, 1]) * 10 ** rng.uniform(-4, -0.5),
                   k_1=rng.choice([0, 1]) * rng.uniform(-1, 1), k_3=rng.choice([0, 1]) * rng.choice([-1, 1]) * 10 ** rng.uniform(-4, -1))
        try:
            regs = sh.regions(1e-3, R_hi)
        except Exception:
            continue
        pat = "".join("A" if r["allowed"] else "F" for r in regs)
        if not pat.endswith("FAFAFA"):
            continue
        last = regs[-5:]
        Rs = np.geomspace(last[0]["R"][0], R_hi, 2000)
        if np.all(sh.ms(Rs) > 0) and sh.beta_out(R_hi) < 0:
            found.append(dict(shell=sh, regions=last))
    return found


def resonance_scan(sh: Shell) -> dict:
    """Barriers, well phases and the resonance exponents for a two-well shell."""
    w1, b1, w2, b2, inf = sh.regions(1e-3, 30.0)[-5:]
    B1, B2 = sh.B(*b1["R"]), sh.B(*b2["R"])
    return dict(B1=float(B1), B2=float(B2), phase1=sh.phase(*w1["R"]), phase2=sh.phase(*w2["R"]),
                outcome="hidden daughter" if sh.beta_out(30.0) < 0 else "parent-eating bubble", **resonance_exponent(B1, B2))


def resonant_masses(sh: Shell, L_over_lP: float, dM: float = 0.02, n: int = 400) -> dict:
    """Bohr-Sommerfeld levels of the two wells versus M (phases scale as (l_f/l_P)^2): masses at which a level of well 1 is degenerate with a
    level of well 2 within the level spacing of well 2, i.e. where resonant transmission occurs."""
    Ms = np.linspace(sh.M * (1 - dM), sh.M * (1 + dM), n)
    n1, n2 = [], []
    for M in Ms:
        s2 = replace(sh, M=M)
        regs = s2.regions(1e-3, 30.0)[-5:]
        if "".join("A" if r["allowed"] else "F" for r in regs) != "AFAFA":
            n1.append(np.nan); n2.append(np.nan)
            continue
        n1.append(L_over_lP**2 * s2.phase(*regs[0]["R"]) / (2 * np.pi) - 0.5)
        n2.append(L_over_lP**2 * s2.phase(*regs[2]["R"]) / (2 * np.pi) - 0.5)
    n1, n2 = np.array(n1), np.array(n2)
    f1 = n1 - np.round(n1)        # well-1 level when f1 crosses 0; resonance if well 2 also has a level there
    lev = [i for i in range(n - 1) if np.isfinite(f1[i]) and np.isfinite(f1[i + 1]) and f1[i] * f1[i + 1] < 0 and abs(f1[i] - f1[i + 1]) < 0.5]
    res = [float(Ms[i]) for i in lev if abs(n2[i] - np.round(n2[i])) < 0.05]
    return dict(levels_well1=len(lev), resonances=res, dn1_dM=float(np.nanmean(np.gradient(n1, Ms))), dn2_dM=float(np.nanmean(np.gradient(n2, Ms))))


# --------------------------------------------------------------------------- route 7: a Brown-Teitelboim membrane


def membrane_shell(rho_out, flux_E, q, tension, M) -> Shell:
    """Interior flux E + q (upward), Delta rho = q (E + q/2) (four-form energy density E^2 / 2); constant tension: a Farhi-Guth wall."""
    rho_in = rho_out + q * (flux_E + q / 2)
    return Shell(chi=np.sqrt(8 * np.pi / 3 * rho_in), M=M, k_w=4 * np.pi * tension, eta=np.sqrt(8 * np.pi / 3 * max(rho_out, 0.0)))


def membrane_branches(rho_out=0.0, flux_E=0.05, q=0.05, tension=1e-3, n=80) -> dict:
    """Branch structure versus M: the bound branch (R < R_1) and the hidden inflating branch (R > R_2, beta_out < 0), the critical mass
    where they merge, and B(M) for quantum creation.  NEC on the membrane: sigma + p = 0 (saturated), so Theorem 6h-1 applies classically."""
    sh0 = membrane_shell(rho_out, flux_E, q, tension, 1e-6)
    rows = []
    Mcr = None
    for M in np.geomspace(1e-5, 2.0, n) / sh0.chi:
        sh = replace(sh0, M=M)
        tp = sh.turning_points(1e-4 / sh.chi, 30 / sh.chi)
        if len(tp) >= 2:
            rows.append(dict(M=float(M), R1=tp[0], R2=tp[-1], B=float(sh.B(tp[0], tp[-1])),
                             hidden=bool(sh.beta_out(tp[-1] * 1.001) < 0)))
        elif rows and Mcr is None:
            Mcr = float(M)
    return dict(rows=rows, M_cr_upper=Mcr, chi=float(sh0.chi), drho=float(q * (flux_E + q / 2)), nec_saturated=True)


def cdl_flat(sigma, drho) -> float:
    """Flat-space thin-wall Coleman-De Luccia exponent 27 pi^2 sigma^4 / (2 drho^3) (gravity neglected)."""
    return float(27 * np.pi**2 * sigma**4 / (2 * drho**3))


def monopole_cdl(lam=1.0) -> float:
    """The same for the monopole wall: sigma = (sqrt2/3) sqrt(lam) eta^3, drho = lam eta^4 / 4  ->  (6912 pi^2 / 162) / lam ~ 421 / lam."""
    eta = 1.0
    return cdl_flat(np.sqrt(2) / 3 * np.sqrt(lam) * eta**3, lam * eta**4 / 4)
