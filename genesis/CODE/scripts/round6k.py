"""Round 6k: the two doors.  Y1 (black-hole interiors that become expanding daughters) in definite theories, Y2 (black holes at collider
energies), the conditional machine Y1 + Y2, and natural bridges within reach (a primordial black hole as Planet Nine).  Writes
DATA/round6k.json and FIGURES/round6k.png; `--figure-only` redraws the figure."""
import json
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from _common import FIG, DATA, dump
from babyuniverse import doors6k as D

warnings.simplefilter("ignore")


def ecsk_path(l_d, A=1.0, S=1e-6, T=40.0):
    V = lambda a: A / a - S / a**4 - 1 + (0 if l_d is None else a**2 / l_d**2)
    if l_d is None:
        a0, ad0 = 0.999999, -1e-9
    else:
        a0, ad0 = A, -np.sqrt(V(A))
    rhs = lambda t, y: [y[1], -A / (2 * y[0] ** 2) + 2 * S / y[0] ** 5 + (0 if l_d is None else y[0] / l_d**2)]
    big = lambda t, y: y[0] - 30.0
    big.terminal = True
    s = solve_ivp(rhs, (0, T), [a0, ad0], rtol=1e-10, atol=1e-13, max_step=0.01, events=[big])
    k = max(1, len(s.t) // 1500)
    return dict(t=s.t[::k].tolist(), a=s.y[0][::k].tolist())


if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round6k.json"))
else:
    out = {}
    b = D.ecsk_bounce()
    out["ecsk_bounce"] = b
    out["ecsk_min_hole"] = [D.ecsk_min_hole(b["rho_b_kg_m3"]), D.ecsk_min_hole(1e45)]
    out["ecsk_criterion"] = D.ecsk_inflation_criterion()
    out["ecsk_interior"] = [D.ecsk_interior(1.0, 1e-6, ld) for ld in (None, 5.0, 3.0, 2.5, 2.0, 1.0)]
    out["ecsk_paths"] = {"none": ecsk_path(None), "2.0": ecsk_path(2.0), "3.0": ecsk_path(3.0)}
    out["fmm_D"] = {str(n): [D.fmm_D(n, x) for x in np.geomspace(0.5, 100, 60)] for n in (0, 4, 6)}
    out["machine"] = [D.add_daughter_machine(n, md) for n in (4, 5, 6, 7) for md in (6.0, 8.0, 10.0, 12.0, 15.0, 20.0)]
    out["machine_safe"] = [D.add_daughter_machine(n, 8.0, safety=3.0) for n in (4, 6)]
    out["reach"] = [D.collider_reach(n, md) for n in (4, 5, 6, 7) for md in (6.0, 8.0, 10.0, 12.0, 15.0, 20.0)]
    out["regular_core"] = D.regular_core_temperature(6, 1.0, 250)
    out["pbh"] = [D.pbh_partner(m) for m in (0.5, 1.0, 3.3, 5.0, 10.0, 20.0)]
    out["y1_theories"] = D.y1_theories()
    dump("round6k.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
ax = ax.ravel()
# (a) ECSK interior
for key, col, lab in (("none", C1, "no vacuum energy: bounce, then cycles"), ("3.0", C3, "vacuum l_d = 3.0 A: bounce, recollapse"),
                      ("2.0", C2, "vacuum l_d = 2.0 A: bounce, then inflates")):
    p = out["ecsk_paths"][key]
    ax[0].semilogy(p["t"], p["a"], color=col, label=lab)
ax[0].set_xlabel("proper time (turnaround units)"); ax[0].set_ylabel("scale factor of the collapsed interior")
b = out["ecsk_bounce"]
ax[0].set_title(f"(a) Einstein-Cartan interior of a collapsing star: torsion bounce at {b['rho_b_kg_m3']:.1e} kg/m^3\n"
                f"(any hole above {out['ecsk_min_hole'][0]['M_min_kg']:.1e} kg); inflates if l_d < 2.6 A", fontsize=8)
ax[0].legend(frameon=False, fontsize=7)
# (b) D-dimensional limiting-curvature junction
for n, col in (("0", MUTED), ("4", C1), ("6", C2)):
    rows = [r for r in out["fmm_D"][n] if r["valid"]]
    ax[1].loglog([r["r_s"] for r in rows], [r["r0_over_l"] for r in rows], color=col, label=f"D = {4 + int(n)}")
ax[1].axvline(1.0, color=INK, lw=0.8, ls="--")
ax[1].axhline(1.0, color=INK, lw=0.6)
ax[1].set_xlabel("horizon radius r_s / limiting length l"); ax[1].set_ylabel("de Sitter region at the junction, r_0 / l")
ax[1].set_title("(b) limiting curvature: Schwarzschild(-Tangherlini) interior joins de Sitter at\nr_0 = (r_s^{n+1} l^2)^{1/(n+3)}; exists iff r_s > l, starts super-Hubble", fontsize=8)
ax[1].legend(frameon=False, fontsize=7)
# (c) the conditional machine: daughter threshold against collider reach
for n, col in ((4, C1), (5, C3), (6, C2), (7, C4)):
    rows = [m for m in out["machine"] if m["n"] == n]
    ax[2].plot([m["M_D_TeV"] for m in rows], [m["M_Y_TeV"] for m in rows], "o-", color=col, label=f"n = {n}: daughter threshold R_S = 1/M_D")
ax[2].axhline(np.sqrt(0.5) * 100, color=INK, ls="--", lw=1)
ax[2].text(6.2, 73, "100 TeV collider reach (tau = 0.5 scaling)", fontsize=7)
ax[2].axhspan(0, 11.4, color="#f2a07e", alpha=0.4)
ax[2].text(6.2, 5, "excluded by CMS (black holes < 8.4-11.4 TeV)", fontsize=7)
ax[2].set_xlabel("fundamental gravity scale M_D (TeV)"); ax[2].set_ylabel("black-hole mass (TeV)")
ax[2].set_title("(c) Y1 + Y2: a black hole holds a limiting-curvature daughter once R_S > 1/M_D\n(M_BH ~ 4-5 M_D): within a 100 TeV collider's reach for M_D <~ 15 TeV", fontsize=8)
ax[2].legend(frameon=False, fontsize=7, loc="upper left")
# (d) regular-core class: evaporation stops at a remnant
rc = out["regular_core"]
mu = np.array([r["mu"] for r in rc["rows"]]); T = np.array([r["T"] for r in rc["rows"]]); rh = np.array([r["r_h"] for r in rc["rows"]])
ax[3].semilogx(mu, T, color=C2, label="static de Sitter core (Hayward-type, D = 10)")
ax[3].semilogx(mu, 7 / (4 * np.pi * mu ** (1 / 7)), color=C1, ls="--", label="general relativity, T = (n+1)/(4 pi r_s)")
ax[3].axvline(rc["mu_remnant"], color=MUTED, ls=":")
ax[3].set_xlabel("mass parameter mu = r_s^{n+1} (units of l)"); ax[3].set_ylabel("Hawking temperature (1/l)")
ax[3].set_title("(d) a collider test of the interior class: a static core makes T fall to zero\nat a remnant; GR (and the expanding FMM class outside) keep rising", fontsize=8)
ax[3].legend(frameon=False, fontsize=7)
# (e) natural bridges within reach: primordial black holes and their partner vacuum
pb = out["pbh"]
ax[4].loglog([p["M_earth"] for p in pb], [p["partner_vac_GeV"] for p in pb], "o-", color=C1, label="Nariai partner vacuum of the hole")
ax[4].axhline(pb[0]["ew_vac_GeV"], color=C2, ls="--", label="electroweak vacuum energy V(0)^{1/4} = 104 GeV")
ax[4].axvspan(5, 10, color=C3, alpha=0.2, label="Planet Nine (5-10 Earth masses, ~500 AU)")
ax[4].set_xlabel("hole mass (Earth masses)"); ax[4].set_ylabel("vacuum scale (GeV)")
ax[4].set_title("(e) a primordial black hole as Planet Nine would be a natural bridge ~8 years away\nat 0.001c, and the Nariai partner of the electroweak scale (to a factor ~2)", fontsize=8)
ax[4].legend(frameon=False, fontsize=7)
# (f) premises
prem = [("P1  gravity strong at 6-20 TeV (n >= 4)", "testable: 100 TeV collider", 1),
        ("P2  black-hole interiors -> expanding daughters", "ECSK, limiting curvature: yes; LQG, mimetic: other", 1),
        ("P3  interior forms before evaporation", "computed: lifetime / interior time ~ 7", 2),
        ("P4  threshold reachable", "M_Y = 24-68 TeV < 71 TeV reach", 2),
        ("P5  semiclassical control", "needs R_S >= 3 l: 1e3-1e4 M_D (10-100 PeV)", 0),
        ("Natural bridge nearby (Planet Nine PBH)", "unconfirmed; probes could check in ~10 years", 1)]
cm = matplotlib.colors.ListedColormap(["#f2a07e", "#f6dc8a", "#9fd8bd"])
ax[5].imshow(np.array([p[2] for p in prem])[:, None], cmap=cm, vmin=0, vmax=2, aspect="auto")
ax[5].set_yticks(range(len(prem))); ax[5].set_yticklabels([p[0] for p in prem], fontsize=7)
ax[5].set_xticks([]); ax[5].grid(False)
for i, p in enumerate(prem):
    ax[5].text(0, i, p[1], ha="center", va="center", fontsize=7, color=INK)
ax[5].set_title("(f) the conditional machine: green computed and satisfied, yellow open\n(testable or theory-dependent), red not met at 100 TeV", fontsize=8)
fig.tight_layout()
fig.savefig(f"{FIG}/round6k.png", dpi=150, facecolor="#fcfcfb")
print("done")
