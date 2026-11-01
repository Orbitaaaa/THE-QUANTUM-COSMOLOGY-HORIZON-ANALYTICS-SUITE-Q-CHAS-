"""Round 13 (GENESIS-13): robust finite-daughter formation.
1. causal-domain certification of every cycle event (Theorem 13-1) and the full-cycle region;
2. realistic matter in D = 10 (pressure, gradients, smooth inhomogeneities) with the Cauchy-horizon cone traced through the fluid;
3. collider formation (content, kinematics, anisotropy).
Writes DATA/round13.json and FIGURES/round13.png.   Usage: python round13.py [--figure-only]"""
import json
import sys
import time
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import robust13 as R
from babyuniverse import daughter12 as Dg

warnings.simplefilter("ignore")
CONFIGS = [("dymnikova", 1.05, 1.4), ("dymnikova", 1.1, 1.3), ("hayward", 1.2, 1.0), ("hayward", 1.05, 1.4)]
PROFILES = ["uniform", "cored", "parabolic", "smooth", "gaussian"]
# central density / volume-averaged density in D = 10 for each profile shape, and the shape's central curvature s''(0)/s(0)
SHAPE_CONC = {"uniform": 1.00, "cored": 6.91, "parabolic": 3.79, "smooth": 17.88, "gaussian": 10.31}
SHAPE_CURV = {"uniform": 0.0, "cored": 0.0, "parabolic": -1.8, "smooth": -4.0, "gaussian": -6.0}
EPS = [0.0, 0.1, 0.5, 2.0]
RIPPLES = [(0.0, 0), (0.05, 1), (0.05, 2)]


def classify(run):
    if "error" in run:
        return "error"
    if not run.get("black_hole_formed"):
        return "no black hole (pressure-supported)"
    e = run["elements"]["0.025"]
    names = [x["event"] for x in e["events"] if x["in_D"]]
    if "second bounce" in names:
        return "full cycle"
    if "maximum" in names:
        return "maximum + recollapse"
    if "bounce" in names:
        return "expands, ends"
    return "no bounce in D (%s)" % run["status"]


def largest_x(run, event):
    xs = [float(f) for f, e in run["elements"].items() if any(x["event"] == event and x["in_D"] for x in e["events"])]
    return max(xs) if xs else 0.0


def survey_job(spec):
    kind, m, chi, prof, eps, rip, seed = spec
    row = dict(kind=kind, mu_over_ext=m, chi_s=chi, profile=prof, eps0=eps, ripple=rip, seed=seed)
    try:
        st = R.make_d10_star(m, chi, prof, eps, kind, N=100, ripple=rip, seed=seed)
        x = 0.5 * (st.R[1:] + st.R[:-1]) / st.R[-1]
        rho = (st.mu[1:] - st.mu[:-1]) / (st.R[1:] ** 9 - st.R[:-1] ** 9)
        r = R.run_cycle(st, t_max=15, wall=120)
        p_over_rho = float(np.max((1.0 / 9) * eps * (1 - x**2) / (1 + eps * (1 - x**2))))
        row.update(status=r["status"], verdict=classify(r), black_hole=r.get("black_hole_formed"), t_plus=r.get("t_plus"),
                   R0_over_rplus=r.get("R0_over_rplus"), t_E=r.get("t_E"),
                   t_cone_centre=r.get("t_cone_centre"), central_concentration=float(SHAPE_CONC[prof]),
                   p_over_rho_central=p_over_rho, x_max_in_D=largest_x(r, "maximum"), x_second_bounce_in_D=largest_x(r, "second bounce"),
                   energy_error=r.get("energy_error"))
    except Exception as e:
        row.update(status="error: %s" % e, verdict="error")
    return row


if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round13.json"))
else:
    out = {}
    # ---- 1. causal-domain certification
    out["certified"] = [R.certify_case(D, m, c, k) for (D, m, c, k) in
                        [(10, 1.05, 1.4, "dymnikova"), (10, 1.01, 1.5, "dymnikova"), (10, 1.1, 1.4, "dymnikova"), (10, 1.2, 1.0, "hayward"),
                         (10, 1.05, 1.4, "hayward"), (5, 20.0, 0.5, "hayward")]]
    for c in out["certified"]:
        print("certify", c["D"], c["kind"], c["mu_over_ext"], c["chi_s"], {k: (round(v["chi_in_D"], 4), v["certified_at_centre"]) for k, v in c["events"].items()},
              "err %.1e" % c["numerical_error"], flush=True)
    out["full_cycle_boundary"] = R.full_cycle_boundary(10, "dymnikova")
    out["full_cycle_boundary_D9"] = R.full_cycle_boundary(9, "dymnikova", mus=[1.005, 1.02, 1.05, 1.1])
    print("boundary done", flush=True)
    # ---- 2a. hydro validation against the exact dust solution, and resolution
    val = []
    for N in (60, 120, 200):
        st = R.make_d10_star(1.05, 1.4, "uniform", 0.0, "dymnikova", N=N)
        r = R.run_cycle(st, t_max=15)
        val.append(dict(N=N, t_E=r["t_E"], t_cone_centre=r["t_cone_centre"], centre=r["elements"]["0.025"]["events"][:3]))
    ev = Dg.evolve_ball(Dg.cdd(10, 1.05, 1.4, "dymnikova")["C"], 10, "dymnikova")
    out["hydro_validation"] = dict(runs=val, analytic=dict(tau_bounce=ev["tau_b"], tau_max=ev["tau_max"], tau_b2=ev["tau_b2"]))
    res = []
    for N in (60, 120, 200):
        for cfg in [("dymnikova", 1.05, 1.4, "cored", 0.1, 0.0, 0), ("hayward", 1.2, 1.4, "parabolic", 0.1, 0.0, 0), ("dymnikova", 1.05, 1.4, "uniform", 0.5, 0.05, 1)]:
            st = R.make_d10_star(cfg[1], cfg[2], cfg[3], cfg[4], cfg[0], N=N, ripple=cfg[5], seed=cfg[6])
            r = R.run_cycle(st, t_max=15)
            res.append(dict(N=N, config=list(cfg), status=r["status"], verdict=classify(r), t_cone_centre=r["t_cone_centre"],
                            centre=r["elements"]["0.025"]["events"][:3]))
    out["resolution"] = res
    print("validation + resolution done", flush=True)
    # ---- 2b. survey (parallel)
    import multiprocessing as mp
    specs = [(kind, m, chi, prof, eps, rip, seed) for (kind, m, chi) in CONFIGS for prof in PROFILES for eps in EPS for (rip, seed) in RIPPLES]
    t0 = time.time()
    with mp.Pool(4) as pool:
        survey = pool.map(survey_job, specs, chunksize=2)
    print("survey %d runs in %.0fs" % (len(survey), time.time() - t0), flush=True)
    out["survey"] = survey
    # example histories for the figure
    ex = {}
    for lab, cfg in [("uniform + pressure (full cycle)", ("dymnikova", 1.05, 1.4, "uniform", 0.5, 0.05, 1)),
                     ("cored + pressure (maximum in D)", ("dymnikova", 1.05, 1.4, "cored", 0.1, 0.0, 0)),
                     ("gaussian (ends while expanding)", ("dymnikova", 1.05, 1.4, "gaussian", 0.1, 0.0, 0))]:
        st = R.make_d10_star(cfg[1], cfg[2], cfg[3], cfg[4], cfg[0], N=120, ripple=cfg[5], seed=cfg[6])
        r = R.run_cycle(st, t_max=15)
        e = r["elements"]["0.025"]
        ex[lab] = dict(t=e["t"], R=e["R"], cone=e["cone_arrival"], events=e["events"], t_E=r["t_E"])
    out["examples"] = ex
    # ---- 3. collider formation
    out["collider_content"] = R.collider_content(n_events=6000)
    out["unbound_lemma"] = [R.unbound_lemma(G, 10, m, k) for G in (1.01, 2.0, 100.0) for m, k in ((1.05, "dymnikova"), (1.2, "hayward"))]
    out["open_frw"] = [R.open_frw_no_bounce(C, 10, k) for C in (5.0, 50.0) for k in ("hayward", "dymnikova")]
    out["anisotropy"] = [R.anisotropy_amplification(10, m, c, k) for (k, m, c) in CONFIGS]
    dump("round13.json", out)
    print("wrote DATA/round13.json", flush=True)

# ------------------------------------------------------------------ figure
INK, MUTED, C1, C2, C3, C4 = "#0b0b0b", "#6b6a66", "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6"
plt.rcParams.update({"font.size": 9})
fig, ax = plt.subplots(2, 3, figsize=(18, 10.5))
# (a) full-cycle window, D = 10 Dymnikova
a0 = ax[0, 0]
fb = out["full_cycle_boundary"]
for row in fb:
    if row["chi_min"] is not None:
        a0.plot([row["chi_min"], row["chi_max"]], [row["mu_over_ext"]] * 2, color=C3, lw=5, solid_capstyle="butt")
a0.set_xlim(0.8, 1.58); a0.set_ylim(1.0, 1.32)
a0.set_xlabel("chi_s"); a0.set_ylabel("mu / mu_ext")
cert = out["certified"][0]
a0.set_title("(a) D = 10 Dymnikova: initial data whose daughter centre completes the full cycle inside D(Sigma)\n"
             "(green bars; certified case: margin %.4f vs numerical error %.0e)" % (cert["events"]["second bounce"]["centre_margin"], cert["numerical_error"]), fontsize=8.5)
# (b) hydro validation
a1 = ax[0, 1]
for lab, e in out["examples"].items():
    l, = a1.plot(e["t"], e["R"], lw=1.6, label=lab)
    if e["cone"]:
        a1.axvline(e["cone"], color=l.get_color(), ls="--", lw=1)
a1.set_xlabel("proper time of the surface [l]"); a1.set_ylabel("areal radius of the central fluid element (x = 0.025)")
a1.set_title("(b) D = 10 fluid runs (pressure, gradients, inhomogeneity): central element vs the cone (dashed)", fontsize=8.5)
a1.legend(fontsize=7)
# (c, d) survey heatmaps
cat = {"full cycle": 3, "maximum + recollapse": 2, "expands, ends": 1}
from matplotlib.colors import ListedColormap
cm = ListedColormap(["#d9d7cf", C1, C4, C3])
for axx, cfgs, title in ((ax[0, 2], CONFIGS[:2], "(c) survey, Dymnikova"), (ax[1, 0], CONFIGS[2:], "(d) survey, Hayward")):
    rows, labels = [], []
    for (kind, m, chi) in cfgs:
        for prof in PROFILES:
            vals = []
            for eps in EPS:
                sub = [s for s in out["survey"] if s["kind"] == kind and s["mu_over_ext"] == m and s["chi_s"] == chi and s["profile"] == prof and s["eps0"] == eps
                   and s.get("status") == "ok"]
                vals.append(np.mean([cat.get(s["verdict"], 0) for s in sub]) if sub else np.nan)
            rows.append(vals); labels.append("%s mu=%.2f chi=%.1f %s" % (kind[:4], m, chi, prof))
    axx.imshow(np.array(rows), cmap=cm, vmin=-0.5, vmax=3.5, aspect="auto")
    axx.set_yticks(range(len(labels))); axx.set_yticklabels(labels, fontsize=6.5)
    axx.set_xticks(range(len(EPS))); axx.set_xticklabels(["eps0=%g" % e for e in EPS])
    axx.set_title(title + " (D = 10), completed runs: verdict at the centre (mean over\n3 perturbation seeds). Green: full cycle; purple: max + recollapse; blue: expands, ends", fontsize=8)
# (e) collider holes
a4 = ax[1, 1]
cc = out["collider_content"]
edges = np.linspace(1, 2, 41)
a4.bar(0.5 * (edges[1:] + edges[:-1]), cc["hist_mu"], width=np.diff(edges), color=C1, alpha=0.7, label="holes at 100 TeV (Round 6o generator)")
a4.axvspan(1.0, 1.1, color=C3, alpha=0.15, label="full-cycle window for BOUND matter (D = 10 Dymnikova)")
a4.set_xlabel("mu / mu_ext of produced holes"); a4.set_ylabel("events")
a4.set_title("(e) collider holes are near threshold (median %.3f, %.0f%% below 1.1), but their content is two massless\n"
             "partons: unbound (Gamma -> infinity), so no bounce in the core and no massive daughter (Lemma 13-3)" % (cc["mu_over_ext_percentiles"][2], 100 * cc["frac_below_1p1"]), fontsize=8)
a4.legend(fontsize=7)
# (f) bound vs unbound interiors
a5 = ax[1, 2]
from scipy.integrate import solve_ivp
from babyuniverse import hydro6n as H


def frw(C, k, a0=3.0, T=12.0):
    def rhs(t, y):
        a, ad = y
        ps, dps = H.psi_and_derivative(C / a**9, 1.0, "dymnikova")
        return [ad, a * (ps - 4.5 * dps * C / a**9)]
    ad0 = -np.sqrt(max(a0 * a0 * float(H.psi_and_derivative(C / a0**9, 1.0, "dymnikova")[0]) - k, 0))
    stop = lambda t, y: y[0] - 0.02
    stop.terminal = True
    return solve_ivp(rhs, (0, T), [a0, ad0], max_step=0.01, events=[stop])


s1 = frw(4.4, 1.0, a0=1.11); s2 = frw(4.4, -1.0, a0=1.11)
a5.plot(s1.t, s1.y[0], color=C3, lw=1.8, label="bound (closed, k = +1): bounces at a_b >= l and cycles")
a5.plot(s2.t, s2.y[0], color=C2, lw=1.8, label="unbound (open, k = -1): no bounce, crunches to a matter caustic")
a5.axhline(1.0, color=MUTED, lw=0.8, ls=":")
a5.set_xlabel("proper time [l]"); a5.set_ylabel("interior scale factor a [l]"); a5.set_ylim(0, 1.5)
a5.set_title("(f) bound vs unbound interior, same tower and density (D = 10)", fontsize=8.5)
a5.legend(fontsize=7)
fig.tight_layout()
fig.savefig(f"{FIG}/round13.png", dpi=130)
print("wrote FIGURES/round13.png")
