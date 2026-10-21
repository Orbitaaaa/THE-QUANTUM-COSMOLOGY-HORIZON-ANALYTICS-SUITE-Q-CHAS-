"""Phase IV: Einstein-scalar dynamics of false-vacuum bubbles (double-null evolution).

1. validation of the evolution code;
2. laboratory (parent) bubbles: collapse, horizon formation, thin-wall comparison;
3. Farhi-Guth in field theory: no anti-trapped region from regular laboratory data (scan);
4. daughter configurations: bounce branch (eps_out = -1) and white-hole branch, inflation;
5. thin -> thick wall transition;
6. stability: radial robustness (non-linear) and aspherical growth (probe Jacobi equation).
Writes DATA/phase4_dynamics.json and FIGURES/phase4_*.png.
"""
import os
import time
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from _common import FIG, dump
from babyuniverse import doublenull as D, stability as S
from babyuniverse.thickwall import QuarticPotential
from babyuniverse.model import ThinWallModel

warnings.simplefilter("ignore")
G, CHI = 1.0, 1.0
RHO_F = 3 * CHI**2 / (8 * np.pi * G)
out = {}
T0 = time.time()


def potential(w, kappa):
    q = QuarticPotential.from_width_tension(w, kappa / (4 * np.pi * G), RHO_F)
    return q, D.Potential.from_quartic(q)


def thin_wall_curve(model, R_start, sign, tau_max, n=800):
    """R(tau) from rest at a turning point: Rddot = -V'(R)/2, moving in direction sign."""
    eps = 1e-9 * R_start
    R0 = R_start + sign * eps
    v0 = sign * np.sqrt(max(-model.V(R0), 0.0))

    def rhs(t, y):
        return [y[1], -0.5 * float(model.dV(y[0]))]

    def hit(t, y):
        return y[0] - 1e-3
    hit.terminal = True
    sol = solve_ivp(rhs, (0, tau_max), [R0, v0], t_eval=np.linspace(0, tau_max, n), rtol=1e-10, atol=1e-12,
                    events=hit)
    return sol.t, sol.y[0]


def leading_stretch(mask):
    """The first contiguous run of True values."""
    mask = np.asarray(mask, dtype=bool)
    out = np.zeros_like(mask)
    idx = np.nonzero(mask)[0]
    if len(idx) == 0:
        return out
    k0 = idx[0]
    k1 = k0
    while k1 + 1 < len(mask) and mask[k1 + 1]:
        k1 += 1
    out[k0:k1 + 1] = True
    return out


def bondi_mass(res):
    mb = res["m_boundary"]
    return float(np.nanmedian(mb[1:max(3, len(mb) // 50)]))


# ------------------------------------------------------------------ 1. validation
val = {}
pot = D.Potential(lam=50.0, v=0.1, rho_f=RHO_F)
errs = []
for h in (0.02, 0.01, 0.005, 0.0025):
    k = int(round(0.04 / h))
    res = D.evolve(lambda r: -pot.v * np.ones_like(r), pot, h, v_max=1.6, ks=k, kv=k)
    U, V = np.meshgrid(res["u"], res["v"], indexing="ij")
    rex, sex = D.de_sitter_exact(CHI, U, V)
    ok = np.isfinite(res["r"])
    off = ok & (rex > 0.15)
    errs.append(dict(h=h, r_max=float(np.nanmax(np.abs(res["r"] - rex)[ok])),
                     s_max=float(np.nanmax(np.abs(res["s"] - sex)[ok])),
                     s_offaxis=float(np.nanmax(np.abs(res["s"] - sex)[off]))))
val["de_sitter"] = errs
pot0 = D.Potential(lam=0.0, v=1.0, rho_f=0.0)
prof = lambda r: 0.01 * np.exp(-((r - 0.5) / 0.1) ** 2)
res = D.evolve(prof, pot0, 0.005, v_max=3.0, G=1e-12)
U, V = np.meshgrid(res["u"], res["v"], indexing="ij")
F = lambda x: (x / 2) * prof(x / 2)
with np.errstate(invalid="ignore", divide="ignore"):
    ex = (F(V) - F(U)) / ((V - U) / 2)
val["flat_pulse_max_error"] = float(np.nanmax(np.abs(res["phi"] - ex)[(V - U) > 0.05]))
sols = {}
for h in (0.004, 0.002, 0.001, 0.0005):
    k = int(round(0.016 / h))
    sols[h] = D.evolve(lambda r: 0.15 * np.exp(-((r - 0.5) / 0.1) ** 2), pot0, h, v_max=3.2, ks=k, kv=k)
hs = sorted(sols, reverse=True)
mask = sols[hs[0]]["r"] > 0.05
for hh in hs:
    rv = np.gradient(sols[hh]["r"], 0.016, axis=1)
    mask &= np.isfinite(sols[hh]["r"]) & (rv > 0.1)
val["self_convergence"] = {}
for key in ("r", "phi", "s"):
    d = [float(np.nanmax(np.abs(sols[hs[n]][key] - sols[hs[n + 1]][key])[mask])) for n in range(3)]
    val["self_convergence"][key] = dict(diffs=d, ratios=[d[0] / d[1], d[1] / d[2]])
coll = []
for A in (0.03, 0.04, 0.05, 0.15):
    rr = {}
    for h in (0.002, 0.001):
        res = D.evolve(lambda r, A=A: A * np.exp(-((r - 0.5) / 0.1) ** 2), pot0, h, v_max=3.2, ks=10, kv=10)
        rah = res["r_ah"][np.isfinite(res["r_ah"])]
        rr[h] = dict(trapped=bool(np.nanmax(res["n_trapped"]) > 0), r_ah=float(rah.max()) if len(rah) else None,
                     M_initial=bondi_mass(res))
    coll.append(dict(A=A, **{f"h={h}": v for h, v in rr.items()}))
val["collapse_threshold"] = coll
out["validation"] = val
print("validation done", round(time.time() - T0), "s")
for e in errs:
    print("  dS", e)
print("  flat pulse", val["flat_pulse_max_error"])
print("  self-convergence", {k: np.round(v["ratios"], 2).tolist() for k, v in val["self_convergence"].items()})
for c in coll:
    print("  collapse", c)

# ------------------------------------------------------------------ 2. laboratory parent
par = []
store = {}
for w in (0.04, 0.02, 0.01, 0.005):
    q, potb = potential(w, 1.0)
    h = w / 24
    res = D.evolve(lambda r: q.v * np.tanh((r - 0.4) / w), potb, h, v_max=2.4, ks=max(1, int(0.004 / h)),
                   kv=max(1, int(0.004 / h)))
    tr = D.wall_trajectory(res)
    M = bondi_mass(res)
    tm = ThinWallModel(CHI, 1.0, M)
    R1 = tm.turning_points()[0]
    i_ah = np.nonzero(np.isfinite(res["r_ah"]))[0]
    early = (tr["r"] > 0.95 * np.nanmax(tr["r"]))
    israel = (tr["A_in"] - tr["A_out"]) / (1.0 * tr["r"])
    rec = dict(w=w, M=M, R_max=float(np.nanmax(tr["r"])), R1_thin=float(R1),
               rel_dev_Rmax=float(np.nanmax(tr["r"]) / R1 - 1), two_GM=2 * G * M,
               r_AH_first=float(res["r_ah"][i_ah[0]]) if len(i_ah) else None,
               u_AH_first=float(res["u_ray"][i_ah[0]]) if len(i_ah) else None,
               max_antitrapped=float(np.nanmax(res["n_antitrapped"])),
               israel_ratio_early=float(np.nanmedian(israel[early])),
               eps_out_early=float(np.sign(np.nanmedian(tr["A_out"][early]))),
               rays_done=res["n_done"], rays=res["n_rays"])
    # thin-wall R(tau) comparison during collapse, aligned at the maximum, until R = 1.3 * 2GM
    k0 = int(np.nanargmax(tr["r"]))
    tt = tr["tau"][k0:] - tr["tau"][k0]
    Rt = tr["r"][k0:]
    sel = Rt > 1.3 * 2 * G * M
    tth, Rth = thin_wall_curve(tm, R1, -1, tt[sel][-1] if sel.any() else 1.0)
    rec["max_rel_dev_R_tau"] = float(np.max(np.abs(np.interp(tt[sel], tth, Rth) - Rt[sel]) / Rt[sel])) if sel.any() else None
    par.append(rec)
    store[w] = (res, tr, tm, tt, Rt, tth, Rth)
    print("  parent", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in rec.items()})
out["parent"] = par

# ------------------------------------------------------------------ 3. Farhi-Guth scan
fg = []
for kappa in (0.5, 1.0, 2.0):
    for rw in (0.15, 0.25, 0.4, 0.55, 0.7, 0.85):
        w = 0.02
        q, potb = potential(w, kappa)
        h = w / 16
        try:
            res = D.evolve(lambda r: q.v * np.tanh((r - rw) / w), potb, h, v_max=3.0, ks=8, kv=8)
        except Exception as e:                                   # pragma: no cover
            fg.append(dict(kappa=kappa, r_w=rw, error=str(e)))
            continue
        M = bondi_mass(res)
        init_antitrapped = bool(np.nanmax(res["n_antitrapped"][1:4]) > 0)
        trapped = bool(np.nanmax(res["n_trapped"]) > 0)
        fg.append(dict(kappa=kappa, r_w=rw, M=M, two_GM_over_rw=2 * G * M / rw,
                       laboratory=not init_antitrapped and 2 * G * M < rw,
                       max_antitrapped=float(np.nanmax(res["n_antitrapped"])), black_hole=trapped,
                       rays_done=res["n_done"], rays=res["n_rays"]))
        print("  FG", fg[-1])
out["farhi_guth_scan"] = fg
lab = [x for x in fg if x.get("laboratory")]
out["farhi_guth_summary"] = dict(n_laboratory=len(lab), n_with_antitrapped=sum(x["max_antitrapped"] > 0 for x in lab),
                                 n_black_hole=sum(x["black_hole"] for x in lab))
print("FG summary", out["farhi_guth_summary"], round(time.time() - T0), "s")

# ------------------------------------------------------------------ 4. daughters
dau = []
dstore = {}
for w in (0.02, 0.01, 0.005):
    q, potb = potential(w, 1.0)
    M0 = 0.1545
    tm0 = ThinWallModel(CHI, 1.0, M0)
    R2 = tm0.turning_points()[-1]
    Ain = np.sqrt(1 - CHI**2 * R2**2)
    wv = 2 * w * Ain
    h = wv / 16
    vw = 2 * R2
    r0, s0, p0 = D.initial_cone_alpha1(lambda v: q.v * np.tanh((v - vw) / wv), h, v_max=vw + 0.3)
    res = D.evolve_cone(r0, s0, p0, potb, h, ks=16, kv=16)
    tr = D.wall_trajectory(res)
    good = leading_stretch(np.isfinite(tr["A_out"]) & np.isfinite(tr["A_in"]) & tr["resolved"])
    for key in ("u", "v", "r", "tau", "A_out", "A_in", "Rdot", "m_out"):
        tr[key] = tr[key][good]
    M = float(np.nanmedian(tr["m_out"][: max(3, len(tr["m_out"]) // 5)]))
    tm = ThinWallModel(CHI, 1.0, M)
    R2m = tm.turning_points()[-1]
    sel = (tr["r"] > 1.05 * R2m) & (tr["r"] < 1.8)
    vel_dev = (np.abs(tr["Rdot"][sel]) - np.sqrt(np.clip(-tm.V(tr["r"][sel]), 0, None))) / np.abs(tr["Rdot"][sel])
    tm_thin = ThinWallModel(CHI, 4 * np.pi * G * q.sigma_thin(), M)
    vel_dev_thin = (np.abs(tr["Rdot"][sel]) - np.sqrt(np.clip(-tm_thin.V(tr["r"][sel]), 0, None))) / np.abs(tr["Rdot"][sel])
    israel = (tr["A_in"] - tr["A_out"]) / tr["r"]
    rec = dict(w=w, M=M, R2_thin=float(R2m), R_start=float(tr["r"][0]), R_end=float(tr["r"][-1]),
               tau_end=float(tr["tau"][-1]), eps_in=float(np.sign(np.median(tr["A_in"]))),
               eps_out=float(np.sign(np.median(tr["A_out"]))), eps_out_all_negative=bool(np.all(tr["A_out"] < 0)),
               israel_ratio_median=float(np.median(israel)), israel_ratio_max_dev=float(np.max(np.abs(israel - 1))),
               median_signed_rel_dev_Rdot_sigma0=float(np.median(vel_dev)) if sel.any() else None,
               median_signed_rel_dev_Rdot_sigma_thin=float(np.median(vel_dev_thin)) if sel.any() else None,
               israel_ratio_p10_p90=[float(np.percentile(israel, 10)), float(np.percentile(israel, 90))],
               max_antitrapped=float(np.nanmax(res["n_antitrapped"])), max_trapped=float(np.nanmax(res["n_trapped"])))
    dau.append(rec)
    dstore[w] = (res, tr, tm)
    print("  daughter (bounce)", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in rec.items()},
          round(time.time() - T0), "s")
out["daughter_bounce"] = dau

# white-hole (regime C) daughter from anti-trapped data on two null rays
w = 0.02
q, potb = potential(w, 1.0)
h = w / 16
# null-cone width 8w: Rdot + A_in = r_v^in w_v / w = 4 at R = 1.7, so Rdot + A_out > 0 (exterior
# anti-trapped: white-hole region IV) and the outgoing ray stays regular
du, dv = D.anti_trapped_data(lambda v: q.v * np.tanh((v - 1.0) / (8 * w)), CHI, 1.2, h, v_max=40.0, u_max=1.0)
resC = D.evolve_rect(du, dv, potb, h, ks=8, kv=8, r_inf=1e3)
trC = D.wall_trajectory(resC)
goodC = leading_stretch(np.isfinite(trC["A_out"]) & np.isfinite(trC["A_in"]) & trC["resolved"])
MC = float(np.nanmedian(trC["m_out"][goodC][:20]))
tmC = ThinWallModel(CHI, 1.0, MC)
rC = trC["r"][goodC]
selC = rC > 0
devC = (np.abs(trC["Rdot"][goodC]) - np.sqrt(np.clip(-tmC.V(rC), 0, None))) / np.abs(trC["Rdot"][goodC])
out["daughter_white_hole"] = dict(M=MC, M_critical=float(tmC.M_critical),
                                  regime=("C" if MC > tmC.M_critical else "A/B"), R_start=float(rC[0]), R_end=float(rC[-1]),
                                  two_GM=2 * G * MC, A_out_sign_start=float(np.sign(trC["A_out"][goodC][0])),
                                  A_in_sign_start=float(np.sign(trC["A_in"][goodC][0])),
                                  israel_ratio_median=float(np.median((trC["A_in"] - trC["A_out"])[goodC] / rC)),
                                  median_signed_rel_dev_Rdot_sigma0=float(np.median(devC)),
                                  max_antitrapped=float(np.nanmax(resC["n_antitrapped"])),
                                  end_kinds=np.bincount(np.nan_to_num(resC["end_kind"]).astype(int)).tolist())
print("  daughter (white hole)", out["daughter_white_hole"])

# ------------------------------------------------------------------ 5. thick walls
thick = []
for w in (0.005, 0.01, 0.02, 0.04, 0.08, 0.16):
    q, potb = potential(w, 1.0)
    h = min(w / 16, 0.002)
    res = D.evolve(lambda r: q.v * np.tanh((r - 0.4) / w), potb, h, v_max=3.0, ks=max(1, int(0.008 / h)),
                   kv=max(1, int(0.008 / h)))
    tr = D.wall_trajectory(res)
    M = bondi_mass(res)
    tm = ThinWallModel(CHI, 1.0, M)
    tps = tm.turning_points()
    i_ah = np.nonzero(np.isfinite(res["r_ah"]))[0]
    thick.append(dict(w=w, w_over_R=w / 0.4, M=M, R_max=float(np.nanmax(tr["r"])) if len(tr["r"]) else None,
                      R1_thin=float(tps[0]) if len(tps) else None,
                      black_hole=bool(np.nanmax(res["n_trapped"]) > 0),
                      r_AH_over_2GM=float(res["r_ah"][i_ah[0]] / (2 * G * M)) if len(i_ah) else None,
                      max_antitrapped=float(np.nanmax(res["n_antitrapped"]))))
    print("  thick", thick[-1])
out["thick_wall"] = thick

# ------------------------------------------------------------------ 6. stability
stab = {}
w = 0.01
q, potb = potential(w, 1.0)
h = w / 16
rob = []
base = lambda r: q.v * np.tanh((r - 0.4) / w)
variants = {"tanh": base,
            "erf profile": lambda r: q.v * __import__("scipy.special", fromlist=["erf"]).erf((r - 0.4) / w),
            "plus outgoing pulse": lambda r: base(r) + 0.01 * q.v * np.exp(-((r - 0.8) / 0.05) ** 2),
            "plus interior pulse": lambda r: base(r) + 0.05 * q.v * np.exp(-((r - 0.2) / 0.04) ** 2)}
for name, f in variants.items():
    res = D.evolve(f, potb, h, v_max=2.4, ks=8, kv=8)
    tr = D.wall_trajectory(res)
    i_ah = np.nonzero(np.isfinite(res["r_ah"]))[0]
    rob.append(dict(variant=name, M=bondi_mass(res), R_max=float(np.nanmax(tr["r"])),
                    black_hole=bool(np.nanmax(res["n_trapped"]) > 0),
                    r_AH=float(res["r_ah"][i_ah[0]]) if len(i_ah) else None,
                    max_antitrapped=float(np.nanmax(res["n_antitrapped"]))))
stab["parent_robustness"] = rob
robd = []
for Rdot0 in (-0.05, 0.0, 0.05):
    M0 = 0.1545
    R2 = ThinWallModel(CHI, 1.0, M0).turning_points()[-1]
    Ain = np.sqrt(1 - CHI**2 * R2**2 + Rdot0**2)
    wv = 2 * w * (Rdot0 + Ain)
    hh = wv / 16
    r0, s0, p0 = D.initial_cone_alpha1(lambda v: q.v * np.tanh((v - 2 * R2) / wv), hh, v_max=2 * R2 + 0.2)
    res = D.evolve_cone(r0, s0, p0, potb, hh, ks=16, kv=16)
    tr = D.wall_trajectory(res)
    g = leading_stretch(np.isfinite(tr["A_out"]) & tr["resolved"])
    robd.append(dict(Rdot0=Rdot0, R_min=float(np.min(tr["r"][g])), R_end=float(tr["r"][g][-1]),
                     eps_out_negative=bool(np.all(tr["A_out"][g] < 0)),
                     inflating=bool(tr["r"][g][-1] > 1.5 * tr["r"][g][0]),
                     max_antitrapped=float(np.nanmax(res["n_antitrapped"]))))
stab["daughter_robustness"] = robd
traj = S.cmc_trajectory(1.0 / (4 * np.pi * G), RHO_F, par[1]["R_max"], 0.02)
stab["probe_zero_modes"] = S.zero_mode_check(traj)
stab["probe_amplification_to_2GM"] = S.amplification(traj, R_stop=par[1]["two_GM"])
stab["probe_amplification_to_half_2GM"] = S.amplification(traj, R_stop=0.5 * par[1]["two_GM"])
out["stability"] = stab
print("  stability", stab)

dump("phase4_dynamics.json", out)

# ------------------------------------------------------------------ figures
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for w, (res, tr, tm, tt, Rt, tth, Rth) in store.items():
    keep = leading_stretch(Rt > 0.5 * 2 * G * par[-1]["M"])
    ax[0].plot(tt[keep], Rt[keep], lw=1.2, label=f"field theory, w = {w}")
tm_last = store[0.005][2]
tth, Rth = thin_wall_curve(tm_last, tm_last.turning_points()[0], -1, 0.3)
keep = Rth > 0.5 * 2 * G * par[-1]["M"]
ax[0].plot(tth[keep], Rth[keep], "k--", lw=1.0, label="thin wall (Israel), w -> 0")
M = par[-1]["M"]
ax[0].axhline(2 * G * M, color="grey", ls=":", label="2GM")
ax[0].set_xlabel("wall proper time from maximum")
ax[0].set_ylabel("wall radius R")
ax[0].set_title("parent: laboratory bubble collapses to a black hole", fontsize=9)
ax[0].legend(fontsize=7)
for w, (res, tr, tm) in dstore.items():
    kk = slice(0, max(1, len(tr["r"]) - 5))
    ax[1].plot(tr["tau"][kk], tr["r"][kk], lw=1.2, label=f"field theory, w = {w}")
tm = dstore[0.005][2]
tth, Rth = thin_wall_curve(tm, tm.turning_points()[-1], +1, float(dstore[0.005][1]["tau"][-1]))
ax[1].plot(tth, Rth, "k--", lw=1.0, label="thin wall bounce branch")
ax[1].set_xlabel("wall proper time from bounce")
ax[1].set_ylabel("R")
ax[1].set_title("daughter: bounce branch inflates (eps_out = -1)", fontsize=9)
ax[1].legend(fontsize=7)
res, tr, tm = dstore[0.01]
kk = slice(0, max(1, len(tr["r"]) - 5))
ax[2].plot(tr["r"][kk], tr["A_in"][kk], label="A_in (field)")
ax[2].plot(tr["r"][kk], tr["A_out"][kk], label="A_out (field)")
ax[2].set_ylim(-2.5, 1.0)
Rg = np.linspace(tr["r"].min(), tr["r"].max(), 200)
ax[2].plot(Rg, tm.A_in(Rg), "k--", lw=0.8, label="thin wall A_in, A_out")
ax[2].plot(Rg, tm.A_out(Rg), "k--", lw=0.8)
ax[2].axhline(0, color="k", lw=0.5)
ax[2].set_xlabel("R")
ax[2].set_title("daughter: signs of the extrinsic curvature (Israel)", fontsize=9)
ax[2].legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "phase4_bubbles.png"), dpi=150)
plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for a, (res, title) in zip(ax, ((store[0.01][0], "parent (regular centre): no anti-trapped region"),
                                (dstore[0.01][0], "daughter (bounce branch): anti-trapped region forms"))):
    r = res["r"]
    hu, hv = res["h"] * res["ks"], res["h"] * res["kv"]
    ru = np.gradient(r, hu, axis=0)
    rv = np.gradient(r, hv, axis=1)
    cls = np.full(r.shape, np.nan)
    cls[(ru < 0) & (rv > 0)] = 0
    cls[(ru > 0) & (rv < 0)] = 1
    cls[(ru < 0) & (rv < 0)] = 2
    cls[(ru > 0) & (rv > 0)] = 3
    Uu, Vv = np.meshgrid(res["u"], res["v"], indexing="ij")
    T, X = (Uu + Vv) / 2, (Vv - Uu) / 2
    a.pcolormesh(X, T, cls, cmap=matplotlib.colors.ListedColormap(["#dfe8f5", "#f5e6c8", "#e8a3a3", "#9fd39f"]),
                 vmin=-0.5, vmax=3.5, shading="auto")
    ph = res["phi"]
    a.contour(X, T, ph, levels=[0.0], colors="k", linewidths=1.2)
    a.set_xlabel("(v - u)/2")
    a.set_ylabel("(u + v)/2")
    a.set_title(title, fontsize=9)
handles = [matplotlib.patches.Patch(color=c, label=l) for c, l in
           zip(["#dfe8f5", "#f5e6c8", "#e8a3a3", "#9fd39f"],
               ["normal (r_u<0<r_v)", "normal, reversed (r_v<0<r_u)", "trapped", "anti-trapped"])]
ax[1].legend(handles=handles + [matplotlib.lines.Line2D([], [], color="k", label="wall (phi = 0)")], fontsize=7,
             loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "phase4_causal_regions.png"), dpi=150)
plt.close(fig)
print("done", round(time.time() - T0), "s")
