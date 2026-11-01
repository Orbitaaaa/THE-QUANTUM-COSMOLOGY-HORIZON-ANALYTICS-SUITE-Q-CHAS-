"""Round 13: robust finite-daughter formation.

1. Causal domain (``in_domain``, ``certify_case``).  In the Oppenheimer-Snyder spacetime (closed FRW patch chi <= chi_s glued to the static
   model-U exterior), with Sigma the moment of time symmetry, an interior event (eta, chi) lies in the domain of dependence D(Sigma) iff
   eta + chi < eta_E + chi_s, where E is the post-bounce crossing of r = r_- by the surface (Theorem 13-1, ROADMAP_ROUND13.md).
   ``certify_case`` evaluates every event of the cycle against this condition, with a tolerance and refinement study.

2. Realistic matter (``run_cycle``).  The Round 6n Lagrangian code (pressure, gradients, shocks, rezoning; time = proper time of the
   pressure-free surface) is extended with:
   - E: the surface's post-bounce crossing of r_-(mu_total) (the exterior is static model U by Birkhoff);
   - the Cauchy-horizon cone, an ingoing radial null ray from E traced through the fluid in the rest-mass coordinate.  A zone of areal
     thickness dR has proper thickness dR / Gamma; the ray covers proper distance e^Phi dt;
   - cycle events of tracked fluid elements (bounce, maximum, second bounce), each compared with the cone's arrival there.

3. Collider formation (``collider_content``, ``open_frw_no_bounce``, ``anisotropy_amplification``).
"""

from __future__ import annotations

from math import gamma, pi

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from . import daughter12 as Dg
from . import hydro6n as H
from . import unified6l as U


# --------------------------------------------------------------------------------------------------------------- 1. causal domain

def in_domain(eta, chi, eta_E, chi_s):
    """Theorem 13-1: interior event (eta, chi) is in D(Sigma) iff eta + chi < eta_E + chi_s (and it is to the future of Sigma)."""
    return eta + chi < eta_E + chi_s


def certify_case(D=10, mu_over_ext=1.05, chi_s=1.4, kind="dymnikova", rtols=(1e-8, 1e-10, 1e-12)) -> dict:
    """Every event of the cycle checked against Theorem 13-1, with the numerical error estimated from a tolerance scan and the
    independent quadrature.  Returns, for each event, eta - eta_b, the largest chi at which it lies in D, and the margin."""
    runs = []
    for rt in rtols:
        r = Dg.cdd(D, mu_over_ext, chi_s, kind)
        ev = Dg.evolve_ball(r["C"], D, kind, rtol=rt)
        sol = ev["sol"]
        aE_t = brentq(lambda t: sol.sol(t)[0] * np.sin(chi_s) - r["r_minus"], ev["tau_b"], ev["tau_max"], xtol=1e-14)
        etaE = float(sol.sol(aE_t)[2])
        runs.append(dict(rtol=rt, eta_b=ev["eta_b"], eta_E=etaE, eta_max=ev["eta_max"], eta_b2=ev["eta_b2"]))
    q = Dg.cdd_fast(D, mu_over_ext, chi_s, kind)
    best = runs[-1]
    budget = best["eta_E"] - best["eta_b"] + chi_s
    events = {"bounce": 0.0, "maximum": best["eta_max"] - best["eta_b"], "second bounce": best["eta_b2"] - best["eta_b"]}
    err_ode = max(abs((x["eta_E"] - x["eta_b"] + chi_s) - budget) for x in runs) + max(abs((x["eta_b2"] - x["eta_b"]) - events["second bounce"]) for x in runs)
    err_quad = abs(q["cone_budget"] - budget) + abs(2 * q["half_period"] - events["second bounce"])
    err = max(err_ode, err_quad, 1e-12)
    rows = {}
    for name, de in events.items():
        chi_max = budget - de                      # the event (simultaneous in eta) lies in D for chi < chi_max
        cin = float(min(max(chi_max, 0.0), chi_s))
        rows[name] = dict(eta_minus_eta_b=de, chi_in_D=cin, centre_margin=float(chi_max), certified_at_centre=bool(chi_max > 10 * err),
                          volume_fraction_in_D=float(Dg._I(cin, D) / Dg._I(chi_s, D)))
    return dict(D=D, kind=kind, mu_over_ext=mu_over_ext, chi_s=chi_s, cone_budget=float(budget), numerical_error=float(err),
                runs=runs, events=rows, untrapped_initially=q["untrapped_initially"], R0_over_rplus=q["R0_over_rplus"])


def full_cycle_boundary(D=10, kind="dymnikova", mus=None, n_chi=60) -> list:
    """For each mu/mu_ext, the interval of chi_s for which the second bounce is in D at the centre (budget > 2 x half period)."""
    mus = mus if mus is not None else [1.002, 1.005, 1.01, 1.02, 1.03, 1.05, 1.07, 1.1, 1.13, 1.16, 1.2, 1.25, 1.3]
    out = []
    chis = np.linspace(0.8, 1.569, n_chi)
    for m in mus:
        ok = []
        best = -np.inf
        for c in chis:
            r = Dg.cdd_fast(D, m, c, kind)
            if not r["valid"] or not r["untrapped_initially"]:
                continue
            margin = r["cone_budget"] - 2 * r["half_period"]
            best = max(best, margin)
            if margin > 0:
                ok.append(float(c))
        out.append(dict(mu_over_ext=m, chi_min=min(ok) if ok else None, chi_max=max(ok) if ok else None, best_margin=float(best)))
    return out


# --------------------------------------------------------------------------------------------------------------- 2. realistic matter

def _rest_mass_coord(star):
    return np.concatenate([[0.0], np.cumsum(star.db)])


def run_cycle(star, track_x=(0.025, 0.05, 0.1, 0.2, 0.3, 0.6), t_max=None, max_steps=200000, rezone=True, frac=0.1, wall=900.0) -> dict:
    """Evolve a (Conservative)Star, locate E and trace the Cauchy-horizon cone, and record cycle events of tracked fluid elements
    (labelled by their rest-mass fraction).  An event is 'in D' if it happens before the cone reaches that element."""
    import time
    t0 = time.time()
    D, kind = star.D, star.kind
    mu_tot = float(star.mu[-1])
    hz = U.static_horizons(mu_tot, D, 1.0, kind)
    if hz is None:
        return dict(error="below extremality: no black hole")
    rm, rp = hz
    star.N0 = getattr(star, "N0", star.N)
    M0 = float(np.sum(star.db))
    m_init = _rest_mass_coord(star)
    R_init = star.R.copy()
    # tracked fluid elements labelled by their initial areal-radius fraction x = R(0)/R_s (rest-mass coordinate is Lagrangian)
    m_of = {x: float(np.interp(x * R_init[-1], R_init, m_init)) for x in track_x}
    track_fracs = tuple(track_x)
    track = {f: dict(R=[], U=[], t=[], events=[]) for f in track_fracs}
    surf_bounced = False
    tE = None
    t_plus = None                              # surface crosses r_+ inward: every later event is inside the event horizon
    cone_m = None                              # rest-mass coordinate of the cone (decreasing to 0)
    cone_t_at = {}                             # arrival time of the cone at each tracked element
    t_cone_centre = None
    status = "ok"
    k = 0
    prevU = {f: None for f in track_fracs}
    prev_surface = (star.R[-1], star.Ub[-1], star.t)
    hist = []
    t_max = t_max or 200.0
    s = star.state(star.R, star.Ub, star.mu) if isinstance(star, H.ConservativeStar) else star.state(star.R, star.Ub, star.mu, star.eps)
    while star.t < t_max and k < max_steps and time.time() - t0 < wall:
        dt, s = star.step()
        k += 1
        if not np.all(np.isfinite(star.R)) or np.min(np.diff(star.R)) <= 0:
            status = "breakdown (shell crossing)"
            break
        if rezone and isinstance(star, H.ConservativeStar):
            if H.merge_thin_zones(star, frac):
                s = star.state(star.R, star.Ub, star.mu)
        if np.min(s["G2"][1:]) < 1e-10:
            status = "comoving slicing fails (Gamma -> 0: neck/throat)"
            break
        mcoord = _rest_mass_coord(star)
        # tracked elements: interpolate R and U at fixed rest-mass fraction
        for f in track_fracs:
            Rf = float(np.interp(m_of[f], mcoord, star.R)); Uf = float(np.interp(m_of[f], mcoord, star.Ub))
            tr = track[f]
            if prevU[f] is not None and tr.get("collapsing"):
                names = [e[0] for e in tr["events"]]
                if prevU[f] < 0 <= Uf:
                    tr["events"].append(("bounce" if "bounce" not in names else "second bounce", star.t, Rf))
                elif prevU[f] > 0 >= Uf and "bounce" in names:
                    tr["events"].append(("maximum", star.t, Rf))
            if Uf < -1e-6:
                tr["collapsing"] = True
            prevU[f] = Uf
            if k % 10 == 0:
                tr["R"].append(Rf); tr["U"].append(Uf); tr["t"].append(star.t)
        # surface: bounce, then the crossing of r_- outward = E (on the exterior Cauchy horizon v = infinity)
        Rs, Us = float(star.R[-1]), float(star.Ub[-1])
        if t_plus is None and Rs < rp and Us < 0:
            t_plus = float(star.t)
        if not surf_bounced and Us > 0:
            surf_bounced = True
        if surf_bounced and tE is None and Rs >= rm and prev_surface[0] < rm:
            w = (rm - prev_surface[0]) / max(Rs - prev_surface[0], 1e-300)
            tE = prev_surface[2] + w * (star.t - prev_surface[2])
            cone_m = M0
        prev_surface = (Rs, Us, star.t)
        # advance the cone inward: proper distance e^Phi dt through zones of proper thickness dR/Gamma
        if cone_m is not None and cone_m > 0:
            Gz = 0.5 * (s["Gam"][1:] + s["Gam"][:-1])
            ell = np.diff(star.R) / np.maximum(Gz, 1e-12)
            lapz = np.exp(s["Phi_z"])
            remaining = dt
            i = int(np.searchsorted(mcoord, cone_m) - 1)
            i = min(max(i, 0), star.N - 1)
            while remaining > 0 and cone_m > 0:
                lo, hi = mcoord[i], mcoord[i + 1]
                frac_in = (cone_m - lo) / max(hi - lo, 1e-300)
                t_need = frac_in * ell[i] / lapz[i]
                if t_need <= remaining:
                    remaining -= t_need
                    cone_m = lo
                    i -= 1
                    if i < 0:
                        cone_m = 0.0
                        break
                else:
                    cone_m = lo + (hi - lo) * (frac_in - remaining * lapz[i] / ell[i])
                    remaining = 0.0
            for f in track_fracs:
                if f not in cone_t_at and cone_m <= m_of[f]:
                    cone_t_at[f] = star.t
            if cone_m <= 0 and t_cone_centre is None:
                t_cone_centre = star.t
        if k % 50 == 0:
            d = star.diagnostics(s)
            hist.append(dict(t=d["t"], R_surface=d["R_surface"], U_surface=d["U_surface"], rho_c=d["rho_c"], n_trapped=d["n_trapped"],
                             n_anti=d["n_anti"], psi_max=d["psi_max"], energy_error=d["constraint"], N=star.N))
        if t_cone_centre is not None and star.t > t_cone_centre + 1.0:
            break
    # verdict per tracked element
    results = {}
    for f, tr in track.items():
        tc = cone_t_at.get(f)
        # in the daughter region: before the cone reaches the element (determined) AND after the surface has crossed r_+ (hidden;
        # sufficient: inside the matter the event horizon reaches the surface at t_plus, so later events lie inside it)
        # significant cycle only: a maximum must re-expand the element by >= AMP over its bounce radius, a second bounce must
        # recontract it by >= AMP below the maximum (filters shock/pressure wiggles; the exact full-cycle case re-expands by 11 %)
        AMP = 0.02
        seq, b1 = [], None
        for name, te, Re in tr["events"]:
            if name == "bounce" and b1 is None:
                b1 = (te, Re); seq.append(("bounce", te, Re))
            elif name == "maximum" and b1 is not None and not any(q[0] == "maximum" for q in seq) and Re >= (1 + AMP) * b1[1]:
                seq.append(("maximum", te, Re))
            elif name == "second bounce" and any(q[0] == "maximum" for q in seq) and not any(q[0] == "second bounce" for q in seq):
                Rm = [q[2] for q in seq if q[0] == "maximum"][0]
                if Re <= Rm / (1 + AMP):
                    seq.append(("second bounce", te, Re))
        evs = [dict(event=e[0], t=float(e[1]), R=float(e[2]), determined=bool(tc is None or e[1] < tc),
                    hidden=bool(t_plus is not None and e[1] >= t_plus),
                    in_D=bool((tc is None or e[1] < tc) and t_plus is not None and e[1] >= t_plus)) for e in seq]
        names = [e["event"] for e in evs if e["in_D"]]
        results["%.3f" % f] = dict(mass_fraction=m_of[f] / M0, cone_arrival=tc, events=evs, max_in_D=("maximum" in names), second_bounce_in_D=("second bounce" in names),
                                   t=tr["t"], R=tr["R"], U=tr["U"])
    return dict(D=D, kind=kind, mu_over_ext=mu_tot / U.extremal_mu(D, 1.0, kind), r_minus=float(rm), r_plus=float(rp), status=status,
                R0_over_rplus=float(R_init[-1] / rp), t_plus=t_plus, black_hole_formed=bool(t_plus is not None),
                steps=k, t_end=float(star.t), t_E=tE, t_cone_centre=t_cone_centre, elements=results, history=hist,
                energy_error=float(max([h["energy_error"] for h in hist] + [0.0])))


def make_d10_star(mu_over_ext=1.05, chi_s=1.4, profile="uniform", eps0=0.0, kind="dymnikova", N=120, ripple=0.0, seed=1, modes=3):
    """D = 10 ball at rest with surface at R^2 psi = sin^2 chi_s (f(R_s) = cos^2 chi_s), optional pressure (eps0) and a smooth, resolved
    density perturbation of relative amplitude `ripple`: rho -> rho (1 + ripple sum_k a_k cos(k pi x + phi_k)), k = 1..modes, random
    a_k (normalised to max |sum| = 1) and phases (inhomogeneity on top of the profile; total mass renormalised)."""
    st = H.make_conservative(10, N, mu_over_ext, profile, eps0, np.sin(chi_s) ** 2, kind)
    if ripple > 0:
        rng = np.random.default_rng(seed)
        D = st.D
        xz = 0.5 * (st.R[1:] + st.R[:-1]) / st.R[-1]
        pert = sum(rng.uniform(0.3, 1.0) * np.cos(k * np.pi * xz + rng.uniform(0, 2 * np.pi)) for k in range(1, modes + 1))
        fac = 1 + ripple * pert / np.max(np.abs(pert))
        rho = (st.mu[1:] - st.mu[:-1]) / ((D - 1) * (st.R[1:] ** (D - 1) - st.R[:-1] ** (D - 1)) / (D - 1)) * fac
        mu = np.concatenate([[0.0], np.cumsum(rho * (st.R[1:] ** (D - 1) - st.R[:-1] ** (D - 1)))])
        mu *= st.mu[-1] / mu[-1]
        st.mu = mu
        s = st.state(st.R, st.Ub, st.mu)
        if np.any(s["G2"][1:] <= 0):
            raise ValueError("rippled data trapped")
        Gz = 0.5 * (s["Gam"][1:] + s["Gam"][:-1])
        V = (st.R[1:] ** (D - 1) - st.R[:-1] ** (D - 1)) / (D - 1)
        rho_new = (st.mu[1:] - st.mu[:-1]) / ((D - 1) * V)
        x = 0.5 * (st.R[1:] + st.R[:-1]) / st.R[-1]
        eps = (eps0 * (1 - x**2)) if np.isscalar(eps0) else np.asarray(eps0)
        n = rho_new / (1 + eps)
        st.db = n * V / Gz
        st._mu0 = float(st.mu[-1])
    st.N0 = st.N
    return st


# --------------------------------------------------------------------------------------------------------------- 3. collider formation

def open_frw_no_bounce(C=50.0, D=10, kind="hayward") -> dict:
    """Unbound matter: an open (k = -1) FRW interior has adot^2 = a^2 psi + 1 > 0, so it never bounces; with psi <= 1/alpha it reaches
    a -> 0 in finite proper time (adot -> 1) while rho ~ C a^{-(D-1)} diverges: a matter caustic with bounded psi."""
    def rhs(t, y):
        a = y[0]
        psi = float(H.psi_and_derivative(C / a ** (D - 1), 1.0, kind)[0])
        return [-np.sqrt(a * a * psi + 1.0)]
    stop = lambda t, y: y[0] - 1e-6
    stop.terminal = True
    a0 = 5.0
    sol = solve_ivp(rhs, (0, 100), [a0], events=[stop], rtol=1e-10, atol=1e-14, max_step=0.01)
    a = sol.y[0]
    psi_end = float(H.psi_and_derivative(C / a[-1] ** (D - 1), 1.0, kind)[0])
    return dict(C=C, D=D, kind=kind, reaches_zero=bool(sol.status == 1), t_crunch=float(sol.t[-1]), adot_end=float(np.sqrt(a[-1] ** 2 * psi_end + 1)),
                psi_end_alpha=psi_end, rho_end=float(C / a[-1] ** (D - 1)))


def unbound_lemma(Gamma=2.0, D=10, mu_over_ext=1.2, kind="dymnikova", n=2000) -> dict:
    """A shell with Gamma > 1 (unbound; massless content is the limit Gamma -> infinity) has U^2 = Gamma^2 - f > 0 at every r in the core,
    where f <= 1: no turning point, hence no bounce of that matter in the regular core."""
    me = U.extremal_mu(D, 1.0, kind)
    rm, rp = U.static_horizons(mu_over_ext * me, D, 1.0, kind)
    r = np.linspace(1e-4, rm, n)
    f = U.static_f(r, mu_over_ext * me, D, 1.0, kind)
    return dict(Gamma=Gamma, min_U2_in_core=float(np.min(Gamma**2 - f)), f_max_core=float(np.max(f)), no_turning_point=bool(np.all(Gamma**2 - f > 0)))


def collider_content(sqrt_s_TeV=100.0, n=6, M_D_TeV=8.0, c_l=0.70, kind="hayward", n_events=4000, seed=3) -> dict:
    """Holes from the Round 6o generator: distribution of mu/mu_ext and of the dimensionless spin; the infalling content is two partons
    (massless, Gamma -> infinity relative to the hole) plus gravitational radiation, so it is unbound (unbound_lemma)."""
    from . import collider6l as C
    from . import pipeline6o as P
    import os
    L = C.load_lumi(os.path.join(os.path.dirname(__file__), "..", "..", "DATA", "round6l_lumi.json"))
    g = P.Generator(L, sqrt_s_TeV, n, M_D_TeV, c_l, kind, seed=seed)
    M, a = [], []
    for _ in range(n_events):
        e = g.event()
        if e["horizon"]:
            M.append(e["M"] / g.M_ext); a.append(e["a_star"])
    M = np.array(M); a = np.array(a)
    q = np.percentile(M, [5, 25, 50, 75, 95])
    return dict(sqrt_s_TeV=sqrt_s_TeV, n=n, M_D_TeV=M_D_TeV, c_l=c_l, M_ext_TeV=g.M_ext / 1e3, n_holes=int(len(M)),
                mu_over_ext_percentiles=q.tolist(), frac_below_1p1=float(np.mean(M < 1.1)), frac_below_1p2=float(np.mean(M < 1.2)),
                a_star_median=float(np.median(a)), trapped_fraction_headon=float(C.trapped_fraction_headon(n)),
                hist_mu=np.histogram(M, bins=np.linspace(1, 2, 41))[0].tolist())


def anisotropy_amplification(D=10, mu_over_ext=1.05, chi_s=1.4, kind="dymnikova") -> dict:
    """Effective estimate: shear energy scales like a stiff component, sigma^2 ~ a^{-2(D-1)} against dust rho ~ a^{-(D-1)}, so a small
    anisotropy ratio eps at maximum expansion is amplified by (a_max / a_b)^{D-1} at the bounce (GR scaling, assumed to carry over)."""
    r = Dg.cdd_fast(D, mu_over_ext, chi_s, kind)
    amp = (r["a_max"] / r["a_b"]) ** (D - 1)
    return dict(D=D, mu_over_ext=mu_over_ext, chi_s=chi_s, kind=kind, a_max_over_a_b=r["a_max"] / r["a_b"], amplification=float(amp),
                note="effective-fluid scaling; the tower's anisotropic (Bianchi) sector is not known")
