"""Round 12: a finite daughter that forms, expands and recollapses inside the Cauchy development of a model-U black hole, plus the
five loopholes.  Writes DATA/round12.json, FIGURES/round12.png and FIGURES/round12_causal.png.
Usage: python round12.py [--figure-only]"""
import json
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import daughter12 as Dg

warnings.simplefilter("ignore")
CASES = {"A: D=5 Hayward, mu=20 mu_ext, chi_s=0.5 (forms, inflates, ends while expanding)": (5, 20.0, 0.5, "hayward"),
         "B: D=10 Hayward, mu=1.2 mu_ext, chi_s=1.0 (maximum and recollapse inside D)": (10, 1.2, 1.0, "hayward"),
         "C: D=10 Dymnikova, mu=1.05 mu_ext, chi_s=1.4 (full cycle inside D)": (10, 1.05, 1.4, "dymnikova")}
L_COLLIDER_M = 0.70 / 8.0 * 1.97327e-19           # l = c_l / M_D for c_l = 0.70, M_D(DL) = 8 TeV, in metres

if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round12.json"))
else:
    out = {"cases": {}}
    for lab, args in CASES.items():
        r = Dg.cdd(*args)
        r["quadrature_check"] = Dg.cdd_fast(*args)
        r["collider_units"] = dict(l_m=L_COLLIDER_M, a_b_m=r["a_b"] * L_COLLIDER_M, a_max_m=r["a_max"] * L_COLLIDER_M,
                                   tau_daughter_center_s=r["tau_daughter_center"] * L_COLLIDER_M / 2.99792458e8 if np.isfinite(r["tau_daughter_center"]) else None)
        out["cases"][lab] = r
        print(lab, "cone %.4f max %.4f b2 %.4f" % (r["eta_cone_minus_eta_b"], r["eta_max_minus_eta_b"], r["eta_b2_minus_eta_b"]), flush=True)
    # timing law and expansion scaling (Proposition-N 12-1/12-2)
    rows = []
    for D in (5, 6, 7, 8, 9, 10):
        for kind in ("hayward", "dymnikova"):
            for m in (1.05, 1.5, 3.0, 10.0, 100.0):
                for chi in (0.5, 1.0, 1.4):
                    try:
                        q = Dg.cdd_fast(D, m, chi, kind)
                    except Exception:
                        continue
                    if not q["valid"]:
                        continue
                    rows.append(dict(D=D, kind=kind, mu_over_ext=m, chi_s=chi, untrapped=q["untrapped_initially"], budget=q["cone_budget"],
                                     half=q["half_period"], max_inside=q["max_inside_D"], full_inside=q["full_cycle_inside_D"],
                                     chi_reaches_max=q["chi_reaches_max"], chi_full=q["chi_full_cycle"], a_max_over_a_b=q["a_max"] / q["a_b"]))
    out["scan"] = rows
    print("scan", len(rows), flush=True)
    large = [Dg.cdd(5, m, 0.5) for m in (10.0, 100.0, 1e4)]
    out["scaling_D5"] = [dict(mu_over_ext=x["mu_over_ext"], C=x["C"], a_cone=x["a_cone"], a_cone_over_C_root=x["a_cone_over_C_root"],
                              efolds_accelerating=x["efolds_accelerating"], eta_cone_minus_eta_b=x["eta_cone_minus_eta_b"],
                              log10_V=x["log10_V_max_over_V_bounce"]) for x in large]
    out["maps"] = {}
    for D, kind in ((10, "dymnikova"), (10, "hayward"), (8, "hayward"), (5, "hayward")):
        out["maps"]["D%d_%s" % (D, kind)] = Dg.cycle_map(D, kind)
        print("map", D, kind, out["maps"]["D%d_%s" % (D, kind)]["frac_max_inside"], out["maps"]["D%d_%s" % (D, kind)]["frac_full_inside"], flush=True)
    out["lemma_core_monotone"] = [Dg.core_monotone(D, m, k) for D in (5, 10) for m in (1.05, 20.0) for k in ("hayward", "dymnikova")]
    # loophole 1 and 2
    out["L1_compact_trapped_region"] = [Dg.compact_trapped_region(D=5, mu0_over_ext=5.0, T_life=T) for T in (5.0, 20.0, 60.0)]
    out["L2_attractor"] = Dg.extremality_attractor(alpha_decay=1e-3)
    out["L2_attractor_our_universe"] = dict(alpha_decay_per_l=2.2e-18 * L_COLLIDER_M / 2.99792458e8,
                                            note="H_0 in units of the collider l; v* = ((1-q) c / alpha)^{1/q}")
    q = out["L2_attractor"]["kappa_minus_power_q"]; c = out["L2_attractor"]["Lambda_coefficient"]
    out["L2_attractor_our_universe"]["v_star_in_l"] = float(((1 - q) * c / out["L2_attractor_our_universe"]["alpha_decay_per_l"]) ** (1 / q))
    out["L2_attractor_our_universe"]["age_of_universe_in_l"] = float(4.35e17 * 2.99792458e8 / L_COLLIDER_M)
    out["L2_finite_time_extremality"] = [Dg.finite_time_extremality(q_law=qq, v_max=2000.0) for qq in (0.5, 0.9, 2.0)]
    dump("round12.json", out)
    print("wrote DATA/round12.json", flush=True)

# ------------------------------------------------------------------ figure
INK, MUTED, C1, C2, C3, C4 = "#0b0b0b", "#6b6a66", "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6"
plt.rcParams.update({"font.size": 9})
fig, ax = plt.subplots(2, 3, figsize=(18, 10.5))
cols = [C1, C2, C3]
# (a) a(eta) with the cone
a0 = ax[0, 0]
for (lab, r), col in zip(out["cases"].items(), cols):
    s = r["series"]; eta = np.array(s["eta"]); a = np.array(s["a"])
    a0.plot(eta, a / r["a_b"], color=col, lw=1.8, label=lab.split(" (")[0])
    a0.axvline(r["eta_cone_minus_eta_b"], color=col, ls="--", lw=1)
    a0.plot([r["eta_max_minus_eta_b"]], [r["a_max"] / r["a_b"]], "o", color=col)
a0.set_yscale("log"); a0.set_xlim(0, 3.0); a0.set_ylim(0.9, 60)
a0.set_xlabel("conformal time since the bounce, eta - eta_b"); a0.set_ylabel("a / a_b (daughter scale factor)")
a0.set_title("(a) daughter scale factor; dashed: the Cauchy-horizon cone reaches the centre (end of D)\n"
             "cone arrives at eta - eta_b ~ pi/2 (the core's de Sitter horizon); dots: maximum expansion", fontsize=8.5)
a0.legend(fontsize=7, loc="upper right")
# (b) exact (chi, eta) diagram of the daughter region for case C
a1 = ax[0, 1]
rC = list(out["cases"].values())[2]
chi_s = rC["chi_s"]; etaE = rC["eta_E"] - rC["eta_b"]; etac = rC["eta_cone_minus_eta_b"]
chis = np.linspace(0, chi_s, 200)
cone = etaE + (chi_s - chis)
a1.fill_between(chis, 0, np.minimum(cone, 2.2), color=C3, alpha=0.25, label="daughter region D (in the Cauchy development)")
a1.plot(chis, cone, color=C2, lw=2, label="Cauchy-horizon cone from E")
for y, lab in ((0.0, "bounce"), (rC["eta_max_minus_eta_b"], "maximum expansion"), (rC["eta_b2_minus_eta_b"], "second bounce")):
    a1.axhline(y, color=INK, lw=0.8, ls=":"); a1.text(chi_s * 0.98, y + 0.02, lab, ha="right", fontsize=8)
a1.plot([chi_s], [etaE], "s", color=C2); a1.text(chi_s - 0.02, etaE - 0.08, "E: surface crosses r_- (v = inf)", ha="right", fontsize=7.5)
a1.set_xlim(0, chi_s); a1.set_ylim(-0.05, 2.0)
a1.set_xlabel("comoving radius chi inside the ball (centre at 0, surface at chi_s)"); a1.set_ylabel("eta - eta_b")
a1.set_title("(b) case C, exact causal diagram of the dust interior: the centre completes\nbounce -> maximum -> recollapse -> second bounce before the cone arrives", fontsize=8.5)
a1.legend(fontsize=7, loc="upper left")


def regime_map(axx, key, title):
    m = out["maps"][key]
    mus = np.array(m["mu_over_ext"]); chis_ = np.array(m["chi_s"])
    b = np.array(m["cone_budget"], float); h = np.array(m["half_period"], float); u = np.array(m["untrapped"], bool)
    reg = np.where(b > 2 * h, 2, np.where(b > h, 1, 0)).astype(float)
    reg[~u] = np.nan
    from matplotlib.colors import ListedColormap
    cm = ListedColormap(["#d9d7cf", C1, C3])
    axx.pcolormesh(chis_, mus, reg, cmap=cm, vmin=-0.5, vmax=2.5, shading="nearest")
    axx.set_yscale("log"); axx.set_xlabel("chi_s (compactness of the initial ball)"); axx.set_ylabel("mu / mu_ext")
    axx.set_title(title + "\ngrey: expands, then ends while expanding; blue: maximum + recollapse inside D;\ngreen: full cycle (second bounce) inside D; white: initially trapped (excluded)", fontsize=8)


regime_map(ax[0, 2], "D10_dymnikova", "(c) D = 10, Dymnikova tower")
regime_map(ax[1, 0], "D10_hayward", "(d) D = 10, Hayward tower")
# (e) loophole 1
a4 = ax[1, 1]
L1 = out["L1_compact_trapped_region"][1]
v = np.array(L1["v"]); a4.plot(v, L1["r_plus"], color=INK, lw=1.5, label="outer apparent horizon")
a4.plot(v, L1["r_minus"], color=C2, lw=1.5, label="inner apparent horizon")
for ray in L1["rays"][::3]:
    a4.plot(ray["v"], ray["r"], color=C1, lw=0.5, alpha=0.6)
a4.plot(L1["body"]["v"], L1["body"]["r"], color=C3, lw=1.6, label="bounced matter (test body)")
a4.set_ylim(0, 12); a4.set_xlim(0, v[-1]); a4.set_xlabel("advanced time v [l]"); a4.set_ylabel("r [l]")
a4.set_title("(e) loophole 1: compact trapped region (D = 5): no event horizon, no Cauchy horizon;\n"
             "blueshift exponent int kappa_- dv = %.0f (RSET ~ e^%.0f); bounced matter returns to the parent" %
             (L1["blueshift_exponent"], L1["blueshift_exponent"]), fontsize=8)
a4.legend(fontsize=7, loc="upper right")
# (f) loophole 2
a5 = ax[1, 2]
L2 = out["L2_attractor"]
vv = np.array(L2["v"]); km = np.array(L2["kappa_minus"]); Lam = np.array(L2["Lambda"])
a5.loglog(vv[1:], km[1:], color=C2, label="kappa_-(v): -> 0 like v^-%.2f (attractor)" % L2["kappa_minus_power_q"])
a5.loglog(vv[1:], Lam[1:], color=C4, label="blueshift exponent int kappa_- dv ~ %.1f v^0.5 (unbounded)" % L2["Lambda_coefficient"])
a5.set_xlabel("advanced time v [l] (D = 10 evaporating hole)"); a5.legend(fontsize=7)
fte = out["L2_finite_time_extremality"]
a5.set_title("(f) loophole 2: evaporation drives kappa_- -> 0, but too slowly to stop mass inflation;\n"
             "finite-time extremality (L ~ (mu - mu_ext)^q, q < 1) would bound it: exponent %.0f (q = 0.5)" % fte[0]["Lambda_total"], fontsize=8)
fig.tight_layout()
fig.savefig(f"{FIG}/round12.png", dpi=130)
print("wrote FIGURES/round12.png")

# ------------------------------------------------------------------ schematic causal diagram
fig, ax = plt.subplots(figsize=(9, 8))
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-0.5, 2.5); ax.set_ylim(-1.6, 2.0)
i0, ip = (2.2, -0.6), (1.4, 0.2)
ax.plot([0, 0], [-1.45, 1.9], color=INK, lw=1.0); ax.text(-0.08, -0.5, "centre of the ball (r = 0)", rotation=90, fontsize=8, ha="right")
ax.plot([0, i0[0]], [-1.45, i0[1]], color=MUTED, lw=0.8, ls=":"); ax.text(1.25, -1.15, "initial slice (ball at rest, untrapped)", fontsize=7.5, color=MUTED, rotation=20)
ax.plot([i0[0], ip[0]], [i0[1], ip[1]], color=INK, lw=1.2); ax.text(1.85, -0.15, "I+ (parent)", fontsize=8)
ax.plot([0, ip[0]], [ip[1] - ip[0], ip[1]], color=INK, lw=1.2); ax.text(0.95, -0.35, "event horizon", rotation=45, fontsize=7.5)
E = (0.55, ip[1] + (ip[0] - 0.55))
ax.plot([ip[0], E[0]], [ip[1], E[1]], color=C2, lw=2.4)
ax.text(1.02, 0.82, "Cauchy horizon N1 (v = infinity)", rotation=-45, color=C2, fontsize=8)
apex = (0.0, E[1] + E[0])
ax.plot([E[0], apex[0]], [E[1], apex[1]], color=C2, lw=1.8, ls="--")
ax.text(0.62, 1.42, "N1 inside the ball:\ningoing cone from E", color=C2, fontsize=7.5)
yb = 0.55
ts = np.linspace(0, 1, 100)
fall_x = 0.8 + (0.25 - 0.8) * ts; fall_y = -1.25 + (yb + 1.25) * ts
rise_x = 0.25 + (E[0] - 0.25) * np.sin(ts * np.pi / 2); rise_y = yb + (E[1] - yb) * ts
ax.plot(fall_x, fall_y, color=C1, lw=2); ax.plot(rise_x, rise_y, color=C1, lw=2)
ax.plot([E[0], E[0] + 0.12], [E[1], E[1] + 0.35], color=C1, lw=1.2, ls=":")
ax.text(0.62, -0.95, "ball surface", color=C1, fontsize=8)
ax.plot(*E, "s", color=C2, ms=7); ax.text(E[0] + 0.06, E[1] - 0.06, "E", color=C2, fontsize=9)
poly_x = np.concatenate([[0.0, 0.0], [apex[0], E[0]], rise_x[::-1]])
poly_y = np.concatenate([[yb, apex[1]], [apex[1], E[1]], rise_y[::-1]])
ax.fill(poly_x, poly_y, color=C3, alpha=0.28, lw=0)
rC = list(out["cases"].values())[2]
fmax, fb2 = rC["eta_max_minus_eta_b"] / rC["eta_cone_minus_eta_b"], rC["eta_b2_minus_eta_b"] / rC["eta_cone_minus_eta_b"]
for frac, lab in ((0.0, "bounce, a = a_b ~ l"), (fmax, "maximum expansion"), (fb2, "second bounce (case C)")):
    y = yb + (apex[1] - yb) * frac
    xr = min(E[0] - (y - E[1]) if y > E[1] else np.interp(y, rise_y, rise_x), 0.6)
    ax.plot([0, xr], [y, y], color=INK, lw=0.7, ls=":"); ax.text(0.03, y - 0.07 if frac > 0.9 else y + 0.02, lab, fontsize=7)
ax.text(0.03, 0.85, "DAUGHTER (region D)", fontsize=8.5, color="#0d6b49", weight="bold")
ax.text(0.78, 0.05, "black-hole interior\n(hidden from the parent)", fontsize=7.5, color=MUTED)
ax.text(1.35, -0.75, "parent universe", fontsize=8.5)
ax.set_title("Round 12: the daughter forms, expands, reaches its maximum and recollapses\n"
             "inside the Cauchy development (schematic conformal diagram; timings from case C)", fontsize=9.5)
fig.tight_layout()
fig.savefig(f"{FIG}/round12_causal.png", dpi=140)
print("wrote FIGURES/round12_causal.png")
