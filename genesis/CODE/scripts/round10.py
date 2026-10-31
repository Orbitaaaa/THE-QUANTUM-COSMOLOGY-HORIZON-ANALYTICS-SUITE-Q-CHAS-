"""Round 10 (10A-10D and the four strengthening items) on the samples of round10a.py.  Writes DATA/round10.json and FIGURES/round10.png."""
import json
import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import FIG, DATA, dump
from babyuniverse import analysis10 as A
from babyuniverse import lhcdata10 as LD
from babyuniverse import pipeline10 as P10
from babyuniverse import stats10 as S

warnings.simplefilter("ignore")
EV = os.path.join(DATA, "round10_events")
L_AVAIL = {100.0: 30000.0, 85.0: 20000.0}          # fb^-1


def lhc_reinterpretation(K13):
    """Full likelihood over the 13 exclusive ATLAS MET bins (post-fit background +- error per bin, uncorrelated) for each 13 TeV model-U
    sample; observed and expected CLs mu95.  Signal: our generator -> Pythia shower -> ATLAS-like detector -> exact selection."""
    a = LD.em_arrays()
    edges = list(a["lo"]) + [1e9]
    bk = {("bin%02d" % i): (np.eye(13)[i] * a["bkg"][i],) * 3 for i in range(13)}
    unc = {("bin%02d" % i): a["err"][i] / a["bkg"][i] for i in range(13)}
    rows = []
    for name, smp in A.load(EV, "l13_*").items():
        ev = smp["ev"]
        if len(ev) == 0:
            continue
        m = P10.atlas_monojet(ev)
        w = ev["weight_pb"][m] * 1e3 * LD.LUMI_FB
        st = S.met_templates(ev["met"][m], w, edges, 0.02)
        model = S.Model(st, bk, sig_unc=0.3, lumi_unc=0.017, bkg_unc=unc)
        mu_obs = S.binned_cls_limit(model, a["obs"])
        mu_exp = S.binned_cls_limit(model, model.asimov(0.0))
        _, n_, md, cl, trap = name.split("_")
        rows.append(dict(sample=name, n=int(n_[1:]), M_D_TeV=float(md[2:]), c_l=float(cl), trapping=trap, M_ext_TeV=float(smp["M_ext_TeV"]),
                         sigma_pb=float(smp["sigma_pb"]), acc=float(m.mean()), S_in_IM12=float(w[ev["met"][m] > 1200].sum()),
                         mu95_obs=float(mu_obs), mu95_exp=float(mu_exp)))
    edges_out = {}
    for key in sorted({(r["n"], r["M_D_TeV"], r["trapping"]) for r in rows}):
        rr = sorted([r for r in rows if (r["n"], r["M_D_TeV"], r["trapping"]) == key], key=lambda r: r["M_ext_TeV"])
        edges_out["n%d_MD%g_%s" % key] = dict(
            M_ext_excluded_below_TeV=A.excl_from_rows([dict(M_ext_TeV=r["M_ext_TeV"], mu95_expected=r["mu95_obs"]) for r in rr]),
            expected_TeV=A.excl_from_rows([dict(M_ext_TeV=r["M_ext_TeV"], mu95_expected=r["mu95_exp"]) for r in rr]))
    return rows, edges_out


def cms_acceptance():
    """Model-U acceptance for the CMS black-hole search signature (N objects >= N_min with pT > 70 GeV-ish, S_T > S_T,min):
    S_T and multiplicity of our 13 TeV samples (objects: jets pT > 30 within |eta| < 2.8 plus leptons/photons)."""
    out = []
    for name, smp in A.load(EV, "l13_n6_MD4_*_opt").items():
        ev = smp["ev"]
        for Nmin in (2, 3, 5):
            for ST in (2000.0, 4000.0, 6000.0):
                out.append(dict(sample=name, M_ext_TeV=float(smp["M_ext_TeV"]), N_min=Nmin, ST_min=ST,
                                acceptance=float(np.mean((ev["n_jets"] + ev["n_e"] + ev["n_mu"] + ev["n_gamma"] >= Nmin) & (ev["ST"] > ST)))))
    return out


def detector_studies(sig):
    """10B: remnant time of flight against the 25 ns bunch spacing, ionisation, timing-layer beta and mass resolution, trigger acceptance."""
    b = np.concatenate([s["ev"]["beta_rem"] for s in sig.values()])
    b = b[(b > 0) & (b < 1)]
    c = 0.299792458                                   # m/ns
    out = {"beta_quantiles": np.quantile(b, [0.1, 0.5, 0.9]).tolist()}
    for R, lab in ((1.5, "timing layer (R = 1.5 m)"), (4.0, "calorimeter (R ~ 4 m)"), (10.0, "muon system (R ~ 10 m)")):
        dt = R / c * (1 / b - 1)
        out[lab] = dict(median_delay_ns=float(np.median(dt)), frac_beyond_25ns=float(np.mean(dt > 25.0)), frac_beyond_75ns=float(np.mean(dt > 75.0)))
    # timing-layer beta resolution and mass: sigma_t = 30 ps, R = 1.5 m; p resolution (FCC-hh tracker, B = 4 T) ~ 2e-5 pT[GeV] (+) 0.5 %
    st = 0.030
    bb = np.median(b)
    t = 1.5 / (c * bb)
    sig_b = bb * st / t
    gam = 1 / np.sqrt(1 - bb**2)
    dlnbg = gam**2 * sig_b / bb                       # d ln(beta gamma) = gamma^2 d beta / beta
    out["timing"] = dict(sigma_t_ps=30, median_beta=float(bb), sigma_beta=float(sig_b), mass_resolution_from_timing=float(dlnbg))
    out["dEdx_over_MIP_median"] = float(np.median(1 / b**2))
    # MET trigger acceptance: offline-selected signal with an L1 MET (10 % resolution) above 1 TeV
    rng = np.random.default_rng(1)
    acc = []
    for s in sig.values():
        ev = s["ev"]; m = P10.fcc_monojet(ev)
        l1 = ev["met"][m] * (1 + 0.1 * rng.normal(size=m.sum()))
        acc.append(float(np.mean(l1 > 1000.0 * 0.8)))
    out["L1_MET_trigger_acceptance"] = float(np.mean(acc))
    return out


def timing_demonstrator(sigma_plane_ps=35.0, L_m=1.0, n_mu=5000, seed=3):
    """10D: two LGAD timing planes L apart.  Cosmic muons (beta ~ 1) calibrate the time-of-flight resolution; a beam test with protons
    and pions at 0.5-3 GeV/c provides slow-particle proxies (beta 0.47-0.95)."""
    rng = np.random.default_rng(seed)
    c = 0.299792458
    sig = np.sqrt(2) * sigma_plane_ps * 1e-3                              # ns
    tof_mu = L_m / c + rng.normal(0, sig, n_mu)
    beta_mu = L_m / (c * tof_mu)
    rows = []
    for p in (0.5, 1.0, 2.0, 3.0):
        bp = p / np.hypot(p, 0.938); bpi = p / np.hypot(p, 0.1396)
        dt = L_m / c * (1 / bp - 1 / bpi)
        rows.append(dict(p_GeV=p, beta_proton=float(bp), beta_pion=float(bpi), dTOF_ps=float(dt * 1e3), separation_sigma=float(dt / sig)))
    return dict(sigma_tof_ps=float(sig * 1e3), beta_resolution_muons=float(np.std(beta_mu)),
                beta_0p4_separation_sigma=float((L_m / c * (1 / 0.4 - 1)) / sig), beam_test=rows,
                note="hardware is not built here; this is the specification and expected performance")


if "--figure-only" in sys.argv:
    out = json.load(open(f"{DATA}/round10.json"))
else:
    out = {}
    v = A.validation_13(EV)
    out["validation_13"] = v
    K = v["K"]
    # FCC-hh reach, two systematics scenarios
    for rs in (100.0, 85.0):
        for sc in ("lhc_like", "conservative"):
            rows = A.fcc_reach(EV, rs, L_AVAIL[rs], K=K, scenario=sc)
            out["reach_%g_%s" % (rs, sc)] = dict(rows=rows, discovery_TeV=A.reach_from_Z(rows), discovery_TeV_from_L5=A.reach_from_rows(rows, L_AVAIL[rs]),
                                                  exclusion_TeV=A.excl_from_rows(rows))
        for var in ("Udis", "Uspin", "Udym"):
            rows = A.fcc_reach(EV, rs, L_AVAIL[rs], K=K, scenario="conservative", pattern="s%g_" + var + "_*")
            out["variant_%g_%s" % (rs, var)] = rows
    # toys at one point (validate asymptotics)
    tem, mcrel = A.fcc_templates(EV, 100.0, K=K)
    smp = A.load(EV, "s100_U_0.70")["s100_U_0.70"]
    st, eff = A.signal_template(smp)
    model = S.scale_model(S.Model(st, tem, sig_unc=0.3, bkg_unc=A.BKG_UNC["conservative"], mc_rel=mcrel), 1.0, 30000.0)
    q1 = model.toys(1.0, 300, seed=5); q0 = model.toys(0.0, 500, seed=6)
    out["toys"] = dict(L_fb=30000.0, Z_asimov=model.Z_asimov(), toy_median_Z=float(np.sqrt(np.median(q1))),
                       bkg_only_frac_q0_gt_1=float(np.mean(q0 > 1)), asymptotic_frac_q0_gt_1=0.1587,
                       bkg_only_frac_q0_gt_4=float(np.mean(q0 > 4)), asymptotic_frac=0.02275)
    # alternatives and the exact separation (shape only, with backgrounds mixed in at the analysis S/B)
    bk = np.concatenate([A.stack(A.load(EV, "b100_%s_*" % p), "b100_%s" % p) for p in ("znunu", "wlnu", "ttbar", "qcd")])
    mb = P10.fcc_monojet(bk)
    bpmf = A.pmf(bk, mb)
    sig = A.load(EV, "s100_*")
    alt = A.load(EV, "a100_*")
    led = A.stack(alt, "a100_LED"); gr62 = alt.get("a100_GR_0.62", {}).get("ev"); gr70 = alt.get("a100_GR_0.70", {}).get("ev")
    dm = {k: A.stack(alt, k) for k in ("a100_DM_20_8", "a100_DM_40_15", "a100_DM_12_5")}
    U62 = sig["s100_U_0.62"]["ev"]; U70 = sig["s100_U_0.70"]["ev"]
    sep = {}
    for sb in (1.0, 0.5):
        tag = "S/B=%s" % ("inf" if sb == 1.0 else "1")
        sep[tag] = {}
        pairs = [("U(0.62) vs GR same threshold", U62, gr62), ("U(0.70) vs GR same threshold", U70, gr70),
                 ("U(0.70) vs ADD graviton emission", U70, led)] + [("U(0.70) vs DM %s" % k[8:], U70, dm[k]) for k in dm] + [
                 ("U-spin(0.70) vs GR", sig["s100_Uspin_0.70"]["ev"], gr70), ("U-discharge(0.70) vs GR", sig["s100_Udis_0.70"]["ev"], gr70),
                 ("U(0.70) vs U-spin(0.70)", U70, sig["s100_Uspin_0.70"]["ev"]), ("Hayward vs Dymnikova", U70, sig["s100_Udym_0.82"]["ev"])]
        for lab, ea, eb in pairs:
            if ea is None or eb is None:
                continue
            r = A.separation(ea, eb, bpmf if sb < 1 else None, sb if sb < 1 else None)
            sep[tag][lab] = r
    out["separation"] = sep
    # LHC reinterpretation (item 4)
    rows, edges = lhc_reinterpretation(K)
    out["lhc_reinterpretation"] = dict(rows=rows, edges=edges)
    out["cms_acceptance"] = cms_acceptance()
    out["onset_rate_precision_100"] = A.onset_rate_precision(EV, 100.0)
    out["onset_shape_fit_100"] = {nm: {k: v for k, v in A.onset_precision(EV, 100.0, nm, L_AVAIL[100.0], K=K, n_scan=60).items()
                                       if k not in ("scan", "q")} for nm in ("s100_U_0.66", "s100_U_0.69")}
    out["detector"] = detector_studies({k: v for k, v in sig.items() if "_U_" in k})
    out["timing_demonstrator"] = timing_demonstrator()
    dump("round10.json", out)

# ------------------------------------------------------------------ figure
C1, C2, C3, C4, INK, MUTED, GRIDC = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 1.8, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
ax = ax.ravel()
v = out["validation_13"]
x = np.array(v["edges"][:-1], float)
ax[0].step(x, np.array(v["atlas_bkg"]), where="post", color=INK, label="ATLAS post-fit SM (139 fb^-1)")
ax[0].step(x, np.array(v["pred_LO"]) * v["K"], where="post", color=C2, label="this pipeline: Pythia LO x K = %.2f, mu = 34" % v["K"])
ax[0].errorbar(x + 25, v["atlas_obs"], yerr=np.sqrt(v["atlas_obs"]), fmt="o", ms=3, color=C1, label="ATLAS data")
ax[0].set_yscale("log"); ax[0].set_xlabel("MET bin lower edge (GeV)"); ax[0].set_ylabel("events per bin")
ax[0].set_title("(a) 13 TeV validation of the full chain (generator, shower, detector, pile-up, exact ATLAS\nselection): shape residual %.0f %% rms above 300 GeV" % (100 * v["rms_shape_residual"]), fontsize=8)
ax[0].legend(frameon=False, fontsize=7)
for rs, col in ((100.0, C2), (85.0, C1)):
    for sc, ls in (("conservative", "-"), ("lhc_like", "--")):
        key = "reach_%g_%s" % (rs, sc)
        if key not in out:
            continue
        rr = out[key]["rows"]
        ax[1].semilogy([r["M_ext_TeV"] for r in rr], [r["L5_fb"] / 1e3 for r in rr], "o" + ls, color=col, ms=3,
                       label="%g TeV, %s syst.: 5 sigma to %.0f TeV" % (rs, sc.replace("_", "-"), out[key]["discovery_TeV"]))
    ax[1].axhline(L_AVAIL[rs] / 1e3, color=col, ls=":", lw=1)
ax[1].set_ylim(1e-6, 1e3); ax[1].set_xlabel("model-U threshold M_ext (TeV)"); ax[1].set_ylabel("luminosity for median 5 sigma (ab^-1)")
ax[1].set_title("(b) profile-likelihood discovery luminosity with nuisance parameters (process norms, MET scale,\nluminosity, signal norm 30 %%, MC stat.); toys: Z = %.2f vs asymptotic %.2f" % (out["toys"]["toy_median_Z"], out["toys"]["Z_asimov"]), fontsize=8)
ax[1].legend(frameon=False, fontsize=6.5)
for var, col, lab in (("Udis", C3, "colour-discharge partons"), ("Uspin", C4, "spin-down quanta"), ("Udym", C1, "Dymnikova tower")):
    rr = out.get("variant_100_%s" % var, [])
    ax[2].semilogy([r["M_ext_TeV"] for r in rr], [r["L5_fb"] / 1e3 for r in rr], "s", color=col, ms=6, label=lab)
rr = out["reach_100_conservative"]["rows"]
ax[2].semilogy([r["M_ext_TeV"] for r in rr], [r["L5_fb"] / 1e3 for r in rr], "o-", color=INK, ms=3, label="baseline model U")
ax[2].set_xlabel("M_ext (TeV)"); ax[2].set_ylabel("5 sigma luminosity (ab^-1), 100 TeV")
ax[2].set_title("(c) signal-model variants (P-R9-b discharge, P-R6n-b rotation, tower): the discovery\nchannel is robust - extra visible energy raises MET acceptance", fontsize=8)
ax[2].legend(frameon=False, fontsize=7)
sp = out["separation"]["S/B=1"]
labs = list(sp)
vals = [sp[k]["events_5sigma"] or np.nan for k in labs]
ax[3].barh(range(len(labs)), vals, color=[C2 if "GR" in k else (C4 if "DM" in k or "ADD" in k else C1) for k in labs])
ax[3].set_yticks(range(len(labs))); ax[3].set_yticklabels([l.replace(" same threshold", " (same M_th)").replace("graviton emission", "KK emission") for l in labs], fontsize=6.3); ax[3].set_xscale("log")
ax[3].set_xlabel("selected events for a median 5 sigma separation (exact LLR, S/B = 1)")
ax[3].set_title("(d) what the signature establishes (S/B = 1): GR at the same threshold ~13 events; DM / ADD\n240-560 events (shape) plus the onset-rate test; the two towers ~170 events", fontsize=8)
rows = out["lhc_reinterpretation"]["rows"]
for (n_, md, trap), col, mk in (((6, 4.0, "opt"), C2, "o"), ((6, 4.0, "con"), C1, "s"), ((6, 2.0, "opt"), C3, "^"), ((4, 4.0, "opt"), C4, "v")):
    rr = sorted([r for r in rows if r["n"] == n_ and r["M_D_TeV"] == md and r["trapping"] == trap], key=lambda r: r["M_ext_TeV"])
    ax[4].semilogy([r["M_ext_TeV"] for r in rr], [r["mu95_obs"] for r in rr], mk + "-", color=col, ms=4,
                   label="n=%d, M_D(DL)=%g, %s trapping" % (n_, md, "full" if trap == "opt" else "head-on"))
ax[4].axhline(1, color=INK, lw=0.8)
ax[4].set_ylim(1e-3, 1e3); ax[4].set_xlabel("M_ext (TeV)"); ax[4].set_ylabel("observed 95 % CL mu (excluded below 1)")
ax[4].set_title("(e) item 4: full-likelihood reinterpretation of the ATLAS 139 fb^-1 monojet bins with the exact\nmodel-U signal (CMS high-multiplicity limits do not apply: acceptance <~ few %)", fontsize=8)
ax[4].legend(frameon=False, fontsize=6.5)
d = out["detector"]; td = out["timing_demonstrator"]
txt = ["10B  remnant beta (10/50/90 %%): %.2f / %.2f / %.2f" % tuple(d["beta_quantiles"])]
for lab in ("timing layer (R = 1.5 m)", "calorimeter (R ~ 4 m)", "muon system (R ~ 10 m)"):
    txt.append("  %s: median +%.0f ns; %.0f %% > 25 ns" % (lab, d[lab]["median_delay_ns"], 100 * d[lab]["frac_beyond_25ns"]))
txt.append("  timing 30 ps: sigma_beta = %.4f; mass from timing %.1f %%" % (d["timing"]["sigma_beta"], 100 * d["timing"]["mass_resolution_from_timing"]))
txt.append("  dE/dx ~ %.0f x MIP; L1 MET trigger acceptance %.0f %%" % (d["dEdx_over_MIP_median"], 100 * d["L1_MET_trigger_acceptance"]))
txt.append("10D  timing demonstrator: 2 LGAD planes, 1 m: sigma_TOF = %.0f ps" % td["sigma_tof_ps"])
txt.append("  cosmic muons: sigma_beta = %.3f; beta = 0.4 separated at %.0f sigma" % (td["beta_resolution_muons"], td["beta_0p4_separation_sigma"]))
for r in td["beam_test"]:
    txt.append("  beam test p = %.1f GeV: p/pi dTOF = %.0f ps (%.0f sigma)" % (r["p_GeV"], r["dTOF_ps"], r["separation_sigma"]))
ax[5].axis("off")
ax[5].text(-0.05, 1.0, "\n".join(txt), va="top", family="monospace", fontsize=6.3, transform=ax[5].transAxes)
ax[5].set_title("(f) 10B detector and remnant studies; 10D timing-detector demonstrator", fontsize=8)
fig.subplots_adjust(left=0.09, right=0.98, top=0.93, bottom=0.06, wspace=0.42, hspace=0.38)
fig.savefig(f"{FIG}/round10.png", dpi=150, facecolor="#fcfcfb")
print("done")
