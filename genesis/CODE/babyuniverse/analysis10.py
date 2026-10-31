"""Round 10 analysis on the cached samples of scripts/round10a.py: 13 TeV validation against ATLAS, FCC-hh profile-likelihood reach,
alternative hypotheses, exact U-vs-GR separation, and the LHC reinterpretation (item 4)."""

from __future__ import annotations

import glob
import os

import numpy as np

from . import lhcdata10 as LD
from . import pipeline10 as P10
from . import stats10 as S


def load(ev_dir, pattern):
    out = {}
    for f in sorted(glob.glob(os.path.join(ev_dir, pattern + ".npz"))):
        d = np.load(f, allow_pickle=True)
        out[os.path.basename(f)[:-4]] = {k: d[k] for k in d.files}
    return out


def stack(samples, prefix):
    evs = [v["ev"] for k, v in samples.items() if k.startswith(prefix)]
    return np.concatenate(evs) if evs else None


# --------------------------------------------------------------------------------------------------------------- validation

def validation_13(ev_dir, procs=("znunu", "wlnu", "ttbar", "qcd")) -> dict:
    """ATLAS monojet selection on our 13 TeV backgrounds (LO Pythia, mu = 34) against the published post-fit prediction."""
    a = LD.em_arrays()
    edges = list(a["lo"]) + [1e9]
    per = {}
    for p in procs:
        ev = stack(load(ev_dir, "b13_%s_*" % p), "b13_%s" % p)
        if ev is None:
            continue
        n, e = P10.binned(ev, P10.atlas_monojet(ev), edges, LD.LUMI_FB)
        per[p] = (n, e)
    tot = sum(v[0] for v in per.values()); err = np.sqrt(sum(v[1] ** 2 for v in per.values()))
    ratio = a["bkg"] / np.maximum(tot, 1e-9)
    sel = a["lo"] >= 300                                     # fit the normalisation above the trigger-dominated lowest bins
    w = 1.0 / ((err[sel] / np.maximum(tot[sel], 1e-9)) ** 2 + (a["err"][sel] / a["bkg"][sel]) ** 2)
    K = float(np.exp(np.sum(w * np.log(ratio[sel])) / np.sum(w)))
    resid = ratio / K - 1.0
    frac = {p: float(per[p][0][sel].sum() / tot[sel].sum()) for p in per}
    return dict(edges=edges, pred_LO=tot.tolist(), mc_err=err.tolist(), atlas_bkg=a["bkg"].tolist(), atlas_err=a["err"].tolist(),
                atlas_obs=a["obs"].tolist(), K=K, shape_residual=resid.tolist(), rms_shape_residual=float(np.sqrt(np.mean(resid[sel] ** 2))),
                composition=frac, per_process={p: per[p][0].tolist() for p in per})


# --------------------------------------------------------------------------------------------------------------- FCC-hh

BKG_UNC = {"lhc_like": {"znunu": 0.05, "wlnu": 0.07, "ttbar": 0.2, "qcd": 1.0},
           "conservative": {"znunu": 0.15, "wlnu": 0.2, "ttbar": 0.3, "qcd": 1.0}}


def fcc_templates(ev_dir, rs, edges=P10.FCC_EDGES, L_fb=1.0, K=1.0, lead_min=500.0, jes=0.02):
    tem = {}; mcrel = None; var = 0.0; tot = 0.0
    for p in ("znunu", "wlnu", "ttbar", "qcd"):
        ev = stack(load(ev_dir, "b%g_%s_*" % (rs, p)), "b%g_%s" % (rs, p))
        if ev is None:
            continue
        m = P10.fcc_monojet(ev, lead_min)
        w = ev["weight_pb"][m] * 1e3 * L_fb * K
        tem[p] = S.met_templates(ev["met"][m], w, edges, jes)
        var = var + np.histogram(ev["met"][m], edges, weights=w**2)[0]
        tot = tot + tem[p][0]
    mcrel = np.sqrt(var) / np.maximum(tot, 1e-12)
    return tem, mcrel


def signal_template(sample, edges=P10.FCC_EDGES, L_fb=1.0, lead_min=500.0, jes=0.02):
    ev = sample["ev"]
    m = P10.fcc_monojet(ev, lead_min)
    w = ev["weight_pb"][m] * 1e3 * L_fb
    return S.met_templates(ev["met"][m], w, edges, jes), float(m.mean())


def fcc_reach(ev_dir, rs, L_avail_fb, K=1.0, scenario="conservative", sig_unc=0.3, pattern="s%g_U_*") -> list:
    tem, mcrel = fcc_templates(ev_dir, rs, L_fb=1.0, K=K)
    sigs = load(ev_dir, pattern % rs)
    rows = []
    for name, smp in sigs.items():
        st, eff = signal_template(smp, L_fb=1.0)
        model = S.Model(st, tem, sig_unc=sig_unc, bkg_unc=BKG_UNC[scenario], mc_rel=mcrel)
        mL = S.scale_model(model, 1.0, L_avail_fb)
        Z = mL.Z_asimov()
        L5 = S.lumi_for_Z(model, 1.0, 5.0, lo=1e-8, hi=1e7)
        mu95 = S.scale_model(model, 1.0, L_avail_fb)
        lim = S.binned_cls_limit(mu95, mu95.asimov(0.0))
        rows.append(dict(sample=name, M_ext_TeV=float(smp["M_ext_TeV"]), sigma_pb=float(smp["sigma_pb"]), eff=eff,
                         Z_at_L=float(Z), L5_fb=float(L5), mu95_expected=float(lim)))
    return sorted(rows, key=lambda r: r["M_ext_TeV"])


def reach_from_rows(rows, L_avail_fb, key="L5_fb"):
    pts = sorted((r["M_ext_TeV"], r[key]) for r in rows)
    best = 0.0
    for (m0, L0), (m1, L1) in zip(pts, pts[1:] + [(None, np.inf)]):
        if L0 > L_avail_fb:
            break
        best = m0
        if m1 is not None and L1 > L_avail_fb and np.isfinite(L1):
            best = m0 + (m1 - m0) * (np.log(L_avail_fb) - np.log(L0)) / (np.log(L1) - np.log(L0))
            break
    return float(best)


def reach_from_Z(rows, Z_target=5.0):
    """Largest M_ext with median Z >= Z_target at the available luminosity (log-linear interpolation of Z between grid points)."""
    pts = sorted((r["M_ext_TeV"], max(r["Z_at_L"], 1e-6)) for r in rows)
    best = 0.0
    for (m0, z0), (m1, z1) in zip(pts, pts[1:] + [(None, 0.0)]):
        if z0 < Z_target:
            break
        best = m0
        if m1 is not None and z1 < Z_target:
            best = m0 + (m1 - m0) * (np.log(z0) - np.log(Z_target)) / (np.log(z0) - np.log(z1))
            break
    return float(best)


def excl_from_rows(rows):
    pts = sorted((r["M_ext_TeV"], r["mu95_expected"]) for r in rows)
    best = 0.0
    for (m0, u0), (m1, u1) in zip(pts, pts[1:] + [(None, np.inf)]):
        if u0 > 1:
            break
        best = m0
        if m1 is not None and u1 > 1 and np.isfinite(u1):
            best = m0 + (m1 - m0) * (np.log(1.0) - np.log(u0)) / (np.log(u1) - np.log(u0))
            break
    return float(best)


# --------------------------------------------------------------------------------------------------------------- discrimination

def features(ev, mask):
    """Per-event discriminating features after selection: MET bin, jet multiplicity (0-4+), MET/S_T bin."""
    met = ev["met"][mask]; nj = np.minimum(ev["n_jets"][mask], 4).astype(int)
    r = ev["met"][mask] / np.maximum(ev["ST"][mask], 1.0)
    mb = np.clip(np.digitize(met, [2000, 3000, 4500, 7000, 11000]), 0, 5)
    rb = np.clip(np.digitize(r, [0.5, 0.8, 0.95, 1.1]), 0, 4)
    return mb * 25 + nj * 5 + rb


def pmf(ev, mask, w=None, nbin=150):
    f = features(ev, mask)
    ww = (ev["weight_pb"][mask] if w is None else w)
    h = np.bincount(f, weights=ww, minlength=nbin)[:nbin]
    return h / max(h.sum(), 1e-300)


def separation(evA, evB, bkg_pmf=None, sb=None, lead_min=500.0):
    """Events (signal + background mixture at signal fraction sb) for a median 5 sigma separation of A from B (exact LLR)."""
    mA = P10.fcc_monojet(evA, lead_min); mB = P10.fcc_monojet(evB, lead_min)
    pA, pB = pmf(evA, mA), pmf(evB, mB)
    if bkg_pmf is not None and sb is not None:
        pA = sb * pA + (1 - sb) * bkg_pmf; pB = sb * pB + (1 - sb) * bkg_pmf
    return S.llr_events_for_5sigma(pA, pB)


def onset_precision(ev_dir, rs, true_name, L_fb, K=1.0, scenario="conservative", sig_unc=0.3, n_scan=200):
    """P-R9-c: precision on M_ext (hence c_l, given n and M_D from the KK monojet) from rate and MET shape.  Asimov data from the true
    point; the signal templates (nominal and MET-scale shifted) are interpolated log-linearly in M_ext between grid points; for each
    hypothesis, -2 ln L is profiled over all nuisances (signal normalisation 30 % included); 1 sigma where -2 Delta ln L = 1."""
    tem, mcrel = fcc_templates(ev_dir, rs, L_fb=L_fb, K=K)
    sigs = {k: v for k, v in load(ev_dir, "s%g_U_*" % rs).items() if len(v["ev"])}
    grid = sorted((float(v["M_ext_TeV"]), signal_template(v, L_fb=L_fb)[0]) for v in sigs.values())
    Ms = np.array([g[0] for g in grid])
    T = np.array([[np.maximum(h, 1e-9) for h in g[1]] for g in grid])        # (n_grid, 3, n_bins)

    def templ(M):
        j = int(np.clip(np.searchsorted(Ms, M) - 1, 0, len(Ms) - 2))
        t = (M - Ms[j]) / (Ms[j + 1] - Ms[j])
        return tuple(np.exp((1 - t) * np.log(T[j, k]) + t * np.log(T[j + 1, k])) for k in range(3))

    M0 = float(sigs[true_name]["M_ext_TeV"])
    data = S.Model(templ(M0), tem, sig_unc=sig_unc, bkg_unc=BKG_UNC[scenario], mc_rel=mcrel).asimov(1.0)
    scan = np.linspace(M0 * 0.85, M0 * 1.15, n_scan)
    q = np.array([S.Model(templ(M), tem, sig_unc=sig_unc, bkg_unc=BKG_UNC[scenario], mc_rel=mcrel).fit(data, 1.0)[0] for M in scan])
    q = 2 * (q - q.min())
    inside = scan[q <= 1.0]
    lo, hi = float(inside.min()), float(inside.max())
    return dict(true_M_ext_TeV=M0, scan=scan.tolist(), q=q.tolist(), lo_TeV=lo, hi_TeV=hi, rel_precision=float((hi - lo) / 2 / M0))


def onset_rate_precision(ev_dir, rs, sig_norm_unc=0.3):
    """Normalisation-limited precision on M_ext: delta M / M = sigma_norm / |d ln(sigma x eff) / d ln M| (shape information ignored, so
    robust against ISR-shape modelling).  Returned for every grid point (central differences in ln M)."""
    sigs = {k: v for k, v in load(ev_dir, "s%g_U_*" % rs).items() if len(v["ev"])}
    pts = sorted((float(v["M_ext_TeV"]), float(v["sigma_pb"]) * float(np.mean(P10.fcc_monojet(v["ev"])))) for v in sigs.values())
    M = np.array([p[0] for p in pts]); y = np.log(np.maximum([p[1] for p in pts], 1e-300))
    slope = np.gradient(y, np.log(M))
    return [dict(M_ext_TeV=float(m), dlnsig_dlnM=float(s_), rel_precision=float(sig_norm_unc / abs(s_)), c_l_rel_precision=float(sig_norm_unc / abs(s_) / 7.0))
            for m, s_ in zip(M, slope)]
