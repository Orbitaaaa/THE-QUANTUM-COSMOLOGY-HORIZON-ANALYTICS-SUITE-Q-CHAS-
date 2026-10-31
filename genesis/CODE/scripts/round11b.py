"""Round 11 collider package (G11-09 ... G11-15): greybody-emission signal variants, the uncertainty budget of the FCC-hh model-U
discovery and exclusion reach, remnant phenomenology, alternative-hypothesis separation under emission variations, and the
reproducibility manifest.  Writes DATA/round11_collider.json, DATA/round11_manifest.json, ROUND11/collider/analysis_config.json and
FIGURES/round11_collider.png.
Usage: python round11b.py [--workers 3] [--figure-only]"""
import json
import os
import sys
import time
import warnings
import multiprocessing as mp

import numpy as np

from _common import DATA, FIG
from babyuniverse import analysis10 as A
from babyuniverse import collider11 as C11
from babyuniverse import collider6l as C6
from babyuniverse import greybody11 as G
from babyuniverse import pipeline10 as P10
from babyuniverse import pipeline6o as P6
import round10a as R

warnings.simplefilter("ignore")
EV = R.EV
ROOT = os.path.dirname(DATA)
GB = json.load(open(f"{DATA}/round11_greybody.json")) if os.path.exists(f"{DATA}/round11_greybody.json") else None


def det(seed):
    return P10.Detector("fcchh", mu_pu=1000.0, pu_lib=R.pu_lib(100.0), eps_chs=0.95, r_pu=0.3, seed=seed)


def gen_job(spec):
    mode, cl, nev, seed = spec
    name = "s100%s_U_%.2f" % (mode, cl)
    out = os.path.join(EV, name + ".npz")
    if os.path.exists(out):
        return name, "cached"
    t = time.time()
    L = C6.load_lumi(os.path.join(DATA, "round6l_lumi.json"))
    patch = G.make_emission(GB["hayward"], mode="greybody" if mode == "gb" else "geometric")
    ev, sig, info = P10.signal(L, 100.0, 6, 8.0, cl, "hayward", nev, det(seed), seed=seed, emission=patch)
    np.savez_compressed(out, ev=ev, sigma_pb=sig, M_ext_TeV=info["M_ext_TeV"], M_th_TeV=info["M_th_TeV"])
    return name, "%.0f s" % (time.time() - t)


def signal_pdf_ratio(cls):
    L23 = C6.load_lumi(os.path.join(DATA, "round6l_lumi.json"))
    L31 = C6.load_lumi(os.path.join(DATA, "round6l_lumi_nnpdf31.json"))
    out = {}
    for cl in cls:
        g23 = P6.Generator(L23, 100.0, 6, 8.0, cl, "hayward"); g31 = P6.Generator(L31, 100.0, 6, 8.0, cl, "hayward")
        out["%.2f" % cl] = g31.sigma_pb / g23.sigma_pb if g23.sigma_pb > 0 else 1.0
    return out


def budget(K):
    base = C11.Config(K=K)
    nom = C11.reach(EV, base)
    res = dict(config_nominal=base.to_json(), nominal=nom, shifts_discovery={}, shifts_exclusion={}, variants={}, groups={})
    groups = {"background theory": [], "signal theory": [], "detector": []}

    def add(label, cfg, group):
        r = C11.reach(EV, cfg)
        res["variants"][label] = dict(config=cfg.to_json(), discovery_TeV=r["discovery_TeV"], exclusion_TeV=r["exclusion_TeV"], rows=r["rows"])
        res["shifts_discovery"][label] = r["discovery_TeV"] - nom["discovery_TeV"]
        res["shifts_exclusion"][label] = r["exclusion_TeV"] - nom["exclusion_TeV"]
        groups[group].append(label)
        print("%-34s disc %.2f (%+.2f)  excl %.2f (%+.2f)" % (label, r["discovery_TeV"], res["shifts_discovery"][label], r["exclusion_TeV"],
                                                            res["shifts_exclusion"][label]), flush=True)

    # background theory: scales and PDFs as smoothed ratios on V+jets; K-factor; EW Sudakov estimate
    res["smoothed_ratios"] = {}
    for var in ("scaleUp", "scaleDn", "pdfCTEQ6L1", "pdfNNPDF31"):
        ratios = {}
        for p in ("znunu", "wlnu"):
            sr = C11.smooth_ratio(EV, "b100%s_%s" % (var, p), "b100_%s" % p, base)
            ratios[p] = (sr["a"], sr["b"]); res["smoothed_ratios"]["%s_%s" % (var, p)] = sr
        add("bkg %s (smoothed)" % var, C11.Config(K=K, bkg_ratio=ratios), "background theory")
        raw = C11.reach(EV, C11.Config(K=K, bkg_prefix={p: "b100%s_%s" % (var, p) for p in ("znunu", "wlnu")}))
        res["variants"]["bkg %s (raw samples, cross-check)" % var] = dict(discovery_TeV=raw["discovery_TeV"], exclusion_TeV=raw["exclusion_TeV"])
    add("bkg K x1.2", C11.Config(K=K, bkg_scale=1.2), "background theory")
    add("bkg K x0.8", C11.Config(K=K, bkg_scale=0.8), "background theory")
    add("bkg EW Sudakov (LL estimate)", C11.Config(K=K, ew_a=0.0059), "background theory")
    add("bkg norm lhc_like priors", C11.Config(K=K, scenario="lhc_like"), "background theory")
    # signal theory
    add("signal norm x0.7", C11.Config(K=K, sig_scale={"all": 0.7}), "signal theory")
    add("signal norm x1.3", C11.Config(K=K, sig_scale={"all": 1.3}), "signal theory")
    pr = signal_pdf_ratio(C11.CL_COMMON); res["signal_pdf_ratio_nnpdf31"] = pr
    add("signal PDF NNPDF3.1 lumi", C11.Config(K=K, sig_scale=pr), "signal theory")
    add("signal PDF NNPDF3.1 mirrored", C11.Config(K=K, sig_scale={k: 1.0 / v for k, v in pr.items()}), "signal theory")
    if os.path.exists(os.path.join(EV, "s100gb_U_0.70.npz")):
        add("emission: greybody (model U)", C11.Config(K=K, sig_prefix="s100gb_U_"), "signal theory")
        add("emission: geometric x^2", C11.Config(K=K, sig_prefix="s100geo_U_"), "signal theory")
    # detector
    for dv in ("tail10", "puWorse"):
        add("detector %s" % dv, C11.Config(K=K, bkg_prefix={p: "b100%s_%s" % (dv, p) for p in ("znunu", "wlnu")},
                                            sig_prefix="s100%s_U_" % dv), "detector")
    add("detector JES 5%", C11.Config(K=K, jes=0.05), "detector")
    res["groups"] = groups
    # combination: within a group, two-sided pairs give the larger |shift| per side
    comb = {}
    for g, labs in groups.items():
        comb[g] = dict(discovery=C11.combine({l: res["shifts_discovery"][l] for l in labs}),
                       exclusion=C11.combine({l: res["shifts_exclusion"][l] for l in labs}))
    tot_d = C11.combine({l: v for l, v in res["shifts_discovery"].items() if "lhc_like" not in l and "EW" not in l})
    tot_e = C11.combine({l: v for l, v in res["shifts_exclusion"].items() if "lhc_like" not in l and "EW" not in l})
    res["combined"] = dict(groups=comb, total_discovery=tot_d, total_exclusion=tot_e,
                           rule="quadrature sum of one-at-a-time shifts, positive and negative separately; the EW Sudakov and lhc_like-prior rows "
                                "are reported but excluded from the total (one-sided improvements a real analysis would earn with NLO EW and "
                                "control regions)")
    return res


def separations():
    """U vs alternatives with the nominal and the greybody / geometric emission (shape-only, S/B = 1 with backgrounds mixed in)."""
    bk = np.concatenate([A.stack(A.load(EV, "b100_%s_*" % p), "b100_%s" % p) for p in C11.PROCS])
    bpmf = A.pmf(bk, P10.fcc_monojet(bk))
    alt = A.load(EV, "a100_*")
    gr70 = alt["a100_GR_0.70"]["ev"]; led = A.stack(alt, "a100_LED")
    out = {}
    for mode in ("", "gb", "geo"):
        nm = "s100%s_U_0.70" % mode
        s = A.load(EV, nm).get(nm)
        if s is None:
            continue
        out[mode or "nominal"] = {"vs GR same threshold": A.separation(s["ev"], gr70, bpmf, 0.5),
                                  "vs ADD graviton emission": A.separation(s["ev"], led, bpmf, 0.5)}
    return out


def figure(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9})
    fig, ax = plt.subplots(2, 3, figsize=(18, 10))
    b = out["budget"]
    # (a) tornado, discovery
    labs = list(b["shifts_discovery"]); vals = np.array([b["shifts_discovery"][l] for l in labs])
    o = np.argsort(np.abs(vals)); labs = [labs[i] for i in o]; vals = vals[o]
    col = ["#c0392b" if v < 0 else "#2471a3" for v in vals]
    ax[0, 0].barh(range(len(labs)), vals, color=col)
    ax[0, 0].set_yticks(range(len(labs))); ax[0, 0].set_yticklabels(labs, fontsize=7)
    ax[0, 0].axvline(0, color="k", lw=0.8)
    ax[0, 0].set_xlabel("shift of the 5 sigma discovery reach in M_ext [TeV]")
    tot = b["combined"]["total_discovery"]
    ax[0, 0].set_title("(a) G11-10/11 uncertainty budget, FCC-hh 100 TeV, 30 ab^-1\nnominal %.1f TeV  +%.1f / -%.1f TeV (total)"
                       % (b["nominal"]["discovery_TeV"], tot["up"], tot["down"]))
    # (b) Z vs M for nominal and envelope
    rows = b["nominal"]["rows"]
    M = [r["M_ext_TeV"] for r in rows]; Z = [max(r["Z_at_L"], 1e-3) for r in rows]
    ax[0, 1].semilogy(M, Z, "k-o", label="nominal")
    for lab, st in (("bkg scaleUp (smoothed)", "--"), ("signal norm x0.7", ":"), ("detector puWorse", "-."), ("emission: greybody (model U)", "-")):
        if lab in b["variants"]:
            rr = b["variants"][lab]["rows"]
            ax[0, 1].semilogy([r["M_ext_TeV"] for r in rr], [max(r["Z_at_L"], 1e-3) for r in rr], st, label=lab)
    ax[0, 1].axhline(5, color="grey", lw=0.8); ax[0, 1].set_ylim(1e-2, 1e5)
    ax[0, 1].set_xlabel("M_ext [TeV]"); ax[0, 1].set_ylabel("median Asimov Z at 30 ab^-1"); ax[0, 1].legend(fontsize=7)
    ax[0, 1].set_title("(b) discovery significance, nominal and leading variations")
    # (c) greybody: per-spin number flux relative to geometric blackbody vs mu/mu_ext
    if out.get("greybody"):
        gb = out["greybody"]["hayward"]
        qs = sorted(gb, key=float)
        for s, c in (("0", "#1b9e77"), ("0.5", "#d95f02"), ("1", "#7570b3")):
            ax[0, 2].semilogx([float(q) for q in qs], [gb[q]["spins"][s]["number_ratio"] for q in qs], "o-", color=c, label="spin %s, number" % s)
            ax[0, 2].semilogx([float(q) for q in qs], [gb[q]["spins"][s]["energy_ratio"] for q in qs], "s--", color=c, label="spin %s, energy" % s)
        if "gr" in out["greybody"]:
            g1 = out["greybody"]["gr"]["1"]["spins"]
            for s, c in (("0", "#1b9e77"), ("0.5", "#d95f02"), ("1", "#7570b3")):
                ax[0, 2].axhline(g1[s]["number_ratio"], color=c, lw=0.6, ls=":")
        ax[0, 2].set_xlabel("mu / mu_ext (Hayward, D = 10)"); ax[0, 2].set_ylabel("flux / geometric-optics blackbody")
        ax[0, 2].legend(fontsize=7); ax[0, 2].set_title("(c) G11-12 brane greybody factors of model U\n(dotted: Tangherlini n = 6)")
    # (d) greybody spectra at q = 1.05 and 3
    if out.get("greybody"):
        for q, ls in (("1.05", "-"), ("3", "--")):
            for s, c in (("0", "#1b9e77"), ("0.5", "#d95f02"), ("1", "#7570b3")):
                sp = gb[q]["spins"][s]; x = np.array(sp["x"]); dN = np.array(sp["dN"])
                ax[1, 0].plot(x, dN / np.trapezoid(dN, x), ls, color=c, label="q=%s spin %s" % (q, s))
        xx = np.linspace(0.01, 16, 300)
        for st, sg in (("B", -1), ("F", 1)):
            y = xx**3 / (np.exp(xx) + sg); ax[1, 0].plot(xx, y / np.trapezoid(y, xx), color="grey", lw=0.8, ls=":")
        ax[1, 0].set_xlabel("x = omega / T"); ax[1, 0].set_ylabel("normalised number spectrum")
        ax[1, 0].set_title("(d) per-quantum energy spectra (grey: Round 6o x^3 blackbody)"); ax[1, 0].legend(fontsize=6, ncol=2)
    # (e) remnant velocity and stopping
    rp = out["remnants"]
    ax[1, 1].axis("off")
    txt = ["G11-14 remnant phenomenology (s100_U_0.70, M_rem ~ %.1f TeV)" % rp["M_TeV"], "",
           "beta percentiles 5/25/50/75/95: " + ", ".join("%.2f" % v for v in rp["beta_percentiles"].values()),
           "fraction beta < 0.5: %.2f    beta < 0.9: %.2f" % (rp["frac_beta_lt_0p5"], rp["frac_beta_lt_0p9"]),
           "if singly charged:",
           "  stops in %.0f g/cm^2 Fe-equivalent: %.3f" % (rp["assumptions"]["depth_g_cm2"], rp["frac_stopping"]),
           "  median CSDA range: %.2e g/cm^2" % rp["median_range_g_cm2"],
           "  median dE/dx / MIP: %.1f" % rp["median_dEdx_over_MIP"],
           "  median TOF delay at %.1f m: %.2f ns (frac > 1 ns: %.2f)" % (rp["assumptions"]["R_timing_m"], rp["median_delay_ns"], rp["frac_delay_gt_1ns"]),
           "  TOF mass resolution (30 ps): %.4f" % rp["median_tof_mass_resolution"], "",
           "G11-15 separation from alternatives (S/B = 1, shape only):"]
    for mode, d in out.get("separations", {}).items():
        for k, v in d.items():
            nn = v.get("events_for_5sigma", v.get("n5", None))
            txt.append("  %-8s %-26s events for 5 sigma: %s" % (mode, k, "%.0f" % nn if isinstance(nn, (int, float)) else str(nn)))
    ax[1, 1].text(0, 1, "\n".join(txt), va="top", family="monospace", fontsize=7.5)
    ax[1, 1].set_title("(e) remnants and alternatives")
    # (f) group budget
    comb = b["combined"]["groups"]
    gs = list(comb)
    ax[1, 2].bar(np.arange(len(gs)) - 0.2, [comb[g]["discovery"]["up"] for g in gs], 0.4, label="+ (discovery)", color="#2471a3")
    ax[1, 2].bar(np.arange(len(gs)) + 0.2, [-comb[g]["discovery"]["down"] for g in gs], 0.4, label="- (discovery)", color="#c0392b")
    ax[1, 2].set_xticks(range(len(gs))); ax[1, 2].set_xticklabels(gs); ax[1, 2].axhline(0, color="k", lw=0.8)
    ax[1, 2].set_ylabel("reach shift [TeV]"); ax[1, 2].legend(fontsize=8)
    ax[1, 2].set_title("(f) grouped budget (quadrature); exclusion nominal %.1f TeV +%.1f/-%.1f" %
                       (b["nominal"]["exclusion_TeV"], b["combined"]["total_exclusion"]["up"], b["combined"]["total_exclusion"]["down"]))
    fig.tight_layout()
    fig.savefig(f"{FIG}/round11_collider.png", dpi=130)


if __name__ == "__main__":
    if "--figure-only" in sys.argv:
        out = json.load(open(f"{DATA}/round11_collider.json"))
        figure(out); sys.exit()
    nw = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 3
    R.pu_lib(100.0)
    specs = [(mode, cl, 1500, 1450 + int(round(cl * 100)) + (0 if mode == "gb" else 7)) for mode in ("gb", "geo") for cl in C11.CL_COMMON]
    with mp.Pool(nw) as pool:
        for name, st in pool.imap_unordered(gen_job, specs):
            print(name, st, flush=True)
    K = json.load(open(f"{DATA}/round10.json"))["validation_13"]["K"]
    out = dict(K13=K, greybody=GB)
    out["budget"] = budget(K)
    out["remnants"] = C11.remnant_phenomenology(A.load(EV, "s100_U_0.70")["s100_U_0.70"])
    print("remnants", out["remnants"], flush=True)
    out["separations"] = separations()
    print("separations done", flush=True)
    C11.write_json(out, f"{DATA}/round11_collider.json")
    os.makedirs(os.path.join(ROOT, "ROUND11", "collider"), exist_ok=True)
    C11.write_json(dict(nominal=out["budget"]["config_nominal"], variants={k: v.get("config") for k, v in out["budget"]["variants"].items()},
                        luminosity_fb=30000.0, sqrt_s_TeV=100.0, K13=K,
                        entry_point="python CODE/scripts/round11b.py  (after round10a.py, round11a.py, round11_greybody.py)"),
                   os.path.join(ROOT, "ROUND11", "collider", "analysis_config.json"))
    C11.write_json(C11.manifest(EV, extra=dict(generation_scripts=["CODE/scripts/round10a.py", "CODE/scripts/round11a.py", "CODE/scripts/round11b.py"],
                                               seeds_note="seeds are fixed in the spec lists of each script; sample names encode process, slice and variant")),
                   f"{DATA}/round11_manifest.json")
    figure(out)
    print("done", flush=True)
