"""Round 11 (G11-09 ... G11-15): uncertainty budget for the FCC-hh model-U reach, remnant phenomenology, and the reproducibility manifest.

Every reach number is computed by ``reach`` from one ``Config`` (selection, binning, nuisances, luminosity, sample prefixes, weights); the
budget repeats it with one ingredient changed at a time and records the shift.  Background generator variations (scales, PDF sets) are
applied as smoothed per-process ratios to the higher-statistics nominal templates (fit of log ratio linear in log MET), so that MC noise
of the variation samples does not enter the budget; the unsmoothed shift is stored as a cross-check.
"""

from __future__ import annotations

import dataclasses
import hashlib
import inspect
import json
import os
from math import pi

import numpy as np

from . import analysis10 as A
from . import pipeline10 as P10
from . import stats10 as S

PROCS = ("znunu", "wlnu", "ttbar", "qcd")
CL_COMMON = (0.58, 0.62, 0.66, 0.68, 0.69, 0.70, 0.71, 0.72)
M_W = 80.4


@dataclasses.dataclass
class Config:
    rs: float = 100.0
    L_fb: float = 30000.0
    K: float = 1.0
    scenario: str = "conservative"
    sig_unc: float = 0.3
    jes: float = 0.02
    lead_min: float = 500.0
    edges: tuple = tuple(P10.FCC_EDGES)
    bkg_prefix: dict = dataclasses.field(default_factory=dict)      # process -> sample prefix (default "b100_<proc>")
    bkg_ratio: dict = dataclasses.field(default_factory=dict)       # process -> (a, b): weight x exp(a + b ln(MET/1 TeV))
    bkg_scale: float = 1.0
    ew_a: float = 0.0                                               # LL Sudakov weight 1 - ew_a ln^2(MET^2/M_W^2) on V+jets
    sig_prefix: str = "s100_U_"
    sig_scale: dict = dataclasses.field(default_factory=dict)       # c_l string -> factor (or "all")
    cls: tuple = CL_COMMON

    def to_json(self):
        d = dataclasses.asdict(self)
        d["edges"] = list(self.edges)
        d["bkg_ratio"] = {k: list(v) for k, v in self.bkg_ratio.items()}
        return d


def ew_weight(met, a):
    if a <= 0:
        return np.ones_like(met)
    L2 = np.log(np.maximum(met, M_W) ** 2 / M_W**2) ** 2
    return np.clip(1.0 - a * L2, 0.1, 1.0)


def _bkg_events(ev_dir, prefix):
    return A.stack(A.load(ev_dir, prefix + "_*"), prefix)


def templates(ev_dir, cfg: Config):
    tem = {}; var = 0.0; tot = 0.0
    for p in PROCS:
        pre = cfg.bkg_prefix.get(p, "b%g_%s" % (cfg.rs, p))
        ev = _bkg_events(ev_dir, pre)
        if ev is None:
            continue
        m = P10.fcc_monojet(ev, cfg.lead_min)
        w = ev["weight_pb"][m] * 1e3 * cfg.K * cfg.bkg_scale
        met = ev["met"][m]
        if p in cfg.bkg_ratio:
            a, b = cfg.bkg_ratio[p]
            w = w * np.exp(a + b * np.log(met / 1000.0))
        if p in ("znunu", "wlnu"):
            w = w * ew_weight(met, cfg.ew_a)
        tem[p] = S.met_templates(met, w, list(cfg.edges), cfg.jes)
        var = var + np.histogram(met, list(cfg.edges), weights=w**2)[0]
        tot = tot + tem[p][0]
    return tem, np.sqrt(var) / np.maximum(tot, 1e-12)


def reach(ev_dir, cfg: Config) -> dict:
    tem, mcrel = templates(ev_dir, cfg)
    rows = []
    for cl in cfg.cls:
        name = "%s%.2f" % (cfg.sig_prefix, cl)
        smp = A.load(ev_dir, name).get(name)
        if smp is None or len(smp["ev"]) == 0:
            continue
        ev = smp["ev"]
        m = P10.fcc_monojet(ev, cfg.lead_min)
        f = cfg.sig_scale.get("%.2f" % cl, cfg.sig_scale.get("all", 1.0))
        w = ev["weight_pb"][m] * 1e3 * f
        st = S.met_templates(ev["met"][m], w, list(cfg.edges), cfg.jes)
        model = S.Model(st, tem, sig_unc=cfg.sig_unc, bkg_unc=A.BKG_UNC[cfg.scenario], mc_rel=mcrel)
        mL = S.scale_model(model, 1.0, cfg.L_fb)
        lim = S.binned_cls_limit(mL, mL.asimov(0.0))
        rows.append(dict(sample=name, M_ext_TeV=float(smp["M_ext_TeV"]), sigma_pb=float(smp["sigma_pb"]) * f, eff=float(m.mean()),
                         Z_at_L=float(mL.Z_asimov()), mu95_expected=float(lim)))
    rows.sort(key=lambda r: r["M_ext_TeV"])
    return dict(rows=rows, discovery_TeV=A.reach_from_Z(rows), exclusion_TeV=A.excl_from_rows(rows))


def smooth_ratio(ev_dir, var_prefix, nom_prefix, cfg: Config):
    """Fit of the variation/nominal ratio of selected events: ln r = a + b ln(MET / 1 TeV), weighted by the MC statistical errors."""
    out = []
    for pre in (var_prefix, nom_prefix):
        ev = _bkg_events(ev_dir, pre)
        m = P10.fcc_monojet(ev, cfg.lead_min)
        w = ev["weight_pb"][m]
        h = np.histogram(ev["met"][m], list(cfg.edges), weights=w)[0]
        e = np.sqrt(np.histogram(ev["met"][m], list(cfg.edges), weights=w**2)[0])
        out.append((h, e))
    (hv, ev_), (hn, en) = out
    ok = (hv > 0) & (hn > 0)
    r = hv[ok] / hn[ok]
    sr = r * np.sqrt((ev_[ok] / hv[ok]) ** 2 + (en[ok] / hn[ok]) ** 2)
    c = np.array([np.sqrt(a * b) if np.isfinite(b) and b < 1e8 else a * 1.3 for a, b in zip(cfg.edges[:-1], cfg.edges[1:])])[ok]
    x = np.log(c / 1000.0); y = np.log(r); wt = (r / np.maximum(sr, 1e-12)) ** 2
    Xm = np.vstack([np.ones_like(x), x]).T
    W = np.diag(wt)
    cov = np.linalg.inv(Xm.T @ W @ Xm)
    a, b = cov @ Xm.T @ W @ y
    tot = np.log(hv.sum() / hn.sum())
    return dict(a=float(a), b=float(b), sa=float(np.sqrt(cov[0, 0])), sb=float(np.sqrt(cov[1, 1])), total_log_ratio=float(tot),
                bin_ratio=r.tolist(), bin_err=sr.tolist())


def combine(shifts: dict) -> dict:
    up = np.sqrt(sum(max(v, 0) ** 2 for v in shifts.values()))
    dn = np.sqrt(sum(min(v, 0) ** 2 for v in shifts.values()))
    return dict(up=float(up), down=float(dn))


# --------------------------------------------------------------------------------------------------------------- remnants (G11-14)

def bethe_dEdx(beta_gamma, Z_over_A=0.4656, I_eV=286.0):
    """Mean energy loss of a singly charged heavy particle [MeV cm^2/g] (Bethe, no density or shell corrections; iron by default).
    Below beta gamma ~ 0.05 the formula is not reliable (it is clipped there)."""
    bg = np.maximum(np.asarray(beta_gamma, float), 0.05)
    b2 = bg**2 / (1 + bg**2)
    Wmax = 2 * 0.511e6 * bg**2                                     # eV, heavy-particle limit
    L = 0.5 * np.log(2 * 0.511e6 * bg**2 * Wmax / I_eV**2) - b2
    return 0.307075 * Z_over_A / b2 * np.maximum(L, 0.1)


def range_g_cm2(M_GeV, beta, Z_over_A=0.4656, I_eV=286.0, n=400):
    """CSDA range [g/cm^2] of a singly charged particle of mass M from velocity beta to beta gamma = 0.05."""
    beta = np.asarray(beta, float)
    out = np.zeros_like(beta)
    for i, b in np.ndenumerate(beta):
        bg0 = b / np.sqrt(1 - b * b)
        if bg0 <= 0.05:
            continue
        bg = np.linspace(0.05, bg0, n)
        g = np.sqrt(1 + bg**2)
        dE_dbg = M_GeV * 1e3 * bg / g                              # MeV per unit beta gamma
        out[i] = np.trapezoid(dE_dbg / bethe_dEdx(bg, Z_over_A, I_eV), bg)
    return out


def remnant_phenomenology(sample: dict, depth_g_cm2=1450.0, R_timing_m=1.5, sigma_t_ps=30.0) -> dict:
    """Velocity spectrum of model-U remnants; if charged: stopping fraction in a calorimeter of the given iron-equivalent depth, ionisation
    relative to a MIP, time-of-flight delay at a timing layer, and the TOF contribution to the mass resolution."""
    ev = sample["ev"]
    beta = ev["beta_rem"][ev["beta_rem"] > 0]
    M = float(sample["M_ext_TeV"]) * 1e3
    q = np.percentile(beta, [5, 25, 50, 75, 95])
    R = range_g_cm2(M, beta[:600])
    c = 0.299792458
    delay = R_timing_m / c * (1 / beta - 1)                        # ns
    sig_beta = beta**2 * c * sigma_t_ps * 1e-3 / R_timing_m
    gam = 1 / np.sqrt(1 - beta**2)
    dEdx_rel = bethe_dEdx(beta * gam) / 1.451                      # iron MIP minimum ~1.451 MeV cm^2/g
    return dict(M_TeV=M / 1e3, n=int(len(beta)), beta_percentiles=dict(zip(("p5", "p25", "p50", "p75", "p95"), q.tolist())),
                frac_beta_lt_0p5=float(np.mean(beta < 0.5)), frac_beta_lt_0p9=float(np.mean(beta < 0.9)),
                frac_stopping=float(np.mean(R < depth_g_cm2)), median_range_g_cm2=float(np.median(R)),
                median_delay_ns=float(np.median(delay)), frac_delay_gt_1ns=float(np.mean(delay > 1.0)),
                median_dEdx_over_MIP=float(np.median(dEdx_rel)),
                median_tof_mass_resolution=float(np.median(gam**2 * sig_beta / beta)),
                assumptions=dict(depth_g_cm2=depth_g_cm2, R_timing_m=R_timing_m, sigma_t_ps=sigma_t_ps,
                                 note="iron-equivalent calorimeter depth and timing radius are assumptions, not FCC-hh CDR values"))


# --------------------------------------------------------------------------------------------------------------- manifest (G11-09)

def sha256(path, block=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def manifest(ev_dir, pdf_dir="/usr/local/share/Pythia8/pdfdata", extra=None) -> dict:
    import platform
    import importlib.metadata as md
    vers = {}
    for p in ("numpy", "scipy", "numba", "sympy", "matplotlib", "pythia8mc"):
        try:
            vers[p] = md.version(p)
        except Exception:
            vers[p] = "unavailable"
    pdfs = {}
    for f in ("NNPDF23_lo_as_0130_qed_0000.dat", "NNPDF31_lo_as_0118_0000.dat", "cteq6l1.tbl"):
        fp = os.path.join(pdf_dir, f)
        if os.path.exists(fp):
            pdfs[f] = sha256(fp)
    samples = {}
    for f in sorted(os.listdir(ev_dir)) if os.path.isdir(ev_dir) else []:
        if f.endswith(".npz"):
            d = np.load(os.path.join(ev_dir, f), allow_pickle=True)
            samples[f[:-4]] = dict(sha256=sha256(os.path.join(ev_dir, f)), n_events=int(len(d["ev"])), sigma_pb=float(d["sigma_pb"]))
    out = dict(python=platform.python_version(), packages=vers, pdf_files_sha256=pdfs,
               pythia_processes=P10.PROCESSES, selection_source=inspect.getsource(P10.fcc_monojet),
               fcc_edges=list(P10.FCC_EDGES), bkg_uncertainties=A.BKG_UNC, samples=samples)
    if extra:
        out.update(extra)
    return out


def write_json(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
