"""Round 9: the experimental design.

For each machine (FCC-hh 85 TeV, 100 TeV; HE-LHC 27 TeV for comparison) and each model-U threshold M_ext in the window:
* signal: the Round 7 lab-frame generator; selection "monojet + invisible remnant": leading object pT > p_cut, MET > p_cut, MET/S_T > 0.8;
  efficiency from Monte Carlo;
* background: Z(-> nu nu) + jet with MET > p_cut (order-of-magnitude estimate of Round 7, x 3 for W(-> l nu) with a lost lepton);
* discovery luminosity: Z_Asimov = sqrt(2((s+b) ln(1+s/b) - s)) = 5, plus s >= 5 events;
* class luminosity: enough selected events to separate model U from a GR hole at 5 sigma (Round 7: ~3 events) and to measure the onset;
* remnant-mass luminosity: if a fraction f_ch of remnants is charged, 10 heavy stable charged tracks with |eta| < 2.5 (mass from time of
  flight and momentum, M = p/(beta gamma)).
``decision_tree`` maps outcomes onto what they imply for the interior class and for the persistence of the daughter (closure lemma).
"""

from __future__ import annotations

from math import pi

import numpy as np
from scipy.optimize import brentq

from . import collider6l as C
from . import pipeline6o as P
from . import predictions7 as Q
from . import unified6l as U


def z_asimov(s, b):
    if b <= 0:
        return np.sqrt(2 * s) if s > 0 else 0.0
    return float(np.sqrt(max(2 * ((s + b) * np.log(1 + s / b) - s), 0.0)))


def selection_efficiency(raw: dict, p_cut=1000.0, frac_cut=0.8) -> float:
    MET = np.array(raw["MET"]); ST = np.array(raw["ST"]); lead = np.array(raw["lead"]); mf = np.array(raw["met_frac"])
    return float(np.mean((MET > p_cut) & (lead > p_cut) & (mf > frac_cut)))


def background_pb(lumi: C.Lumi, sqrt_s, p_cut_GeV) -> float:
    return 3.0 * Q.znunu_met_tail_pb(lumi, sqrt_s, p_cut_GeV / 1e3)


def design_point(lumi: C.Lumi, sqrt_s=85.0, n=6, M_D_TeV=8.0, c_l=0.6, kind="hayward", n_events=3000, seed=17, f_ch=0.3) -> dict:
    g = P.Generator(lumi, sqrt_s, n, M_D_TeV, c_l, kind, seed=seed)
    # the generator keeps a sliver below tau_max for sampling; a threshold beyond the kinematic limit is not reachable
    if g.sigma_pb <= 0 or (g.M_th / g.y) ** 2 / (sqrt_s * 1e3) ** 2 >= lumi.tau_max:
        return dict(sqrt_s=sqrt_s, c_l=c_l, M_ext_TeV=g.M_ext / 1e3, sigma_pb=0.0, reachable=False)
    lab = Q.LabGenerator(g, f_charged=f_ch, seed=seed)
    raw = lab.run(n_events)
    best = None
    for p_cut in (500.0, 1000.0, 2000.0, 3000.0, 5000.0, 8000.0, 12000.0, 16000.0):
        eff = selection_efficiency(raw, p_cut)
        sig = g.sigma_pb * eff
        bkg = background_pb(lumi, sqrt_s, p_cut)
        if sig <= 0:
            continue
        # luminosity (ab^-1) for Z = 5 and s >= 5
        f = lambda L: z_asimov(sig * L * 1e6, bkg * L * 1e6) - 5.0
        Lmin = 5.0 / (sig * 1e6)
        L5 = brentq(f, 1e-12, 1e6) if f(1e6) > 0 else np.inf
        L5 = max(L5, Lmin)
        if best is None or L5 < best["L_discovery_ab"]:
            best = dict(p_cut_GeV=p_cut, eff=eff, sigma_sel_fb=sig * 1e3, bkg_fb=bkg * 1e3, L_discovery_ab=float(L5))
    hscp_frac = float(np.mean(np.array(raw["hscp"]) > 0))
    L_class = 10.0 / (g.sigma_pb * max(best["eff"], 1e-9) * 1e6)                 # ~10 selected events: class + onset
    L_mass = 10.0 / (g.sigma_pb * max(hscp_frac, 1e-9) * 1e6) if hscp_frac > 0 else np.inf
    return dict(sqrt_s=sqrt_s, n=n, M_D_TeV=M_D_TeV, c_l=c_l, kind=kind, M_ext_TeV=g.M_ext / 1e3, sigma_pb=g.sigma_pb, reachable=True,
                **best, L_class_ab=float(L_class), hscp_fraction=hscp_frac, L_remnant_mass_ab=float(L_mass))


def luminosity_requirements(lumi: C.Lumi, machines=((85.0, 20.0), (100.0, 30.0), (27.0, 15.0)), n=6, M_D_TeV=8.0, kind="hayward",
                            c_ls=None, n_events=2000) -> list:
    c_ls = np.linspace(0.58, 0.76, 10) if c_ls is None else c_ls
    rows = []
    for rs, L_ab in machines:
        for cl in c_ls:
            d = design_point(lumi, rs, n, M_D_TeV, float(cl), kind, n_events)
            d["L_available_ab"] = L_ab
            d["discoverable"] = bool(d.get("reachable") and d["L_discovery_ab"] <= L_ab)
            rows.append(d)
    return rows


def reach(rows: list, sqrt_s, key="L_discovery_ab", L_available=None) -> float:
    """Largest M_ext (TeV) whose required luminosity ``key`` is within the available one: log-linear interpolation between grid points."""
    pts = sorted((r["M_ext_TeV"], r.get(key, np.inf)) for r in rows if r["sqrt_s"] == sqrt_s and r.get("reachable"))
    L_av = L_available if L_available is not None else next(r["L_available_ab"] for r in rows if r["sqrt_s"] == sqrt_s)
    best = 0.0
    for (m0, L0), (m1, L1) in zip(pts, pts[1:] + [(None, np.inf)]):
        if L0 > L_av:
            break
        best = m0
        if m1 is not None and L1 > L_av:
            if np.isfinite(L1):
                best = m0 + (m1 - m0) * (np.log(L_av) - np.log(L0)) / (np.log(L1) - np.log(L0))
            break
    return float(best)


def reach_table(rows: list) -> list:
    out = []
    for rs in sorted({r["sqrt_s"] for r in rows}):
        out.append(dict(sqrt_s=rs, L_available_ab=next(r["L_available_ab"] for r in rows if r["sqrt_s"] == rs),
                        discovery_TeV=reach(rows, rs, "L_discovery_ab"), class_TeV=reach(rows, rs, "L_class_ab"),
                        remnant_mass_TeV=reach(rows, rs, "L_remnant_mass_ab")))
    return out


def detector_requirements(M_ext_TeV=11.5, beta=0.39, path_m=10.0) -> dict:
    """Specifications that follow from the predicted signatures."""
    dt_ns = path_m / 0.2998 * (1 / beta - 1)
    return dict(
        met_trigger="single MET > 1 TeV trigger (unprescaled); jet pT > 1 TeV",
        coverage="hermetic calorimetry to |eta| ~ 6 (FCC-hh reference) for MET",
        jet_energy_scale="1-2 % at 1-10 TeV",
        met_resolution="~5 % at multi-TeV MET",
        tof_delay_ns_at_10m=float(dt_ns),
        timing="~30-50 ps timing layers: the slow remnant (beta ~ 0.4) arrives ~%.0f ns late at 10 m, separable from the bunch crossing" % dt_ns,
        dEdx="highly ionising track (dE/dx ~ 1/beta^2 ~ %.0f x MIP)" % (1 / beta**2),
        remnant_mass="M = p / (beta gamma) for charged remnants: a direct measurement of M_ext (~%.1f TeV)" % M_ext_TeV,
    )


def decision_tree() -> list:
    return [
        dict(step=1, observation="monojet / monophoton excess consistent with KK gravitons", implies="P1: low-scale gravity, n and M_D measured"),
        dict(step=2, observation="no KK excess at FCC-hh", implies="P1 excluded up to the reach: no collider black holes of any interior class"),
        dict(step=3, observation="onset of monojet + invisible-remnant events (MET/S_T -> 1) at M_th >> M_D",
             implies="regular-core (model-U) class; M_th = M_ext gives c_l; the remnant implies the trapped region never closes, so in this class the daughter persists (closure lemma)"),
        dict(step=4, observation="onset of high-multiplicity balanced events at a few M_D, complete evaporation",
             implies="GR-like exterior: either a singular interior, an FMM-type daughter (unobservable), or a closing trapped region (no daughter); undecided by experiment"),
        dict(step=5, observation="heavy stable charged tracks with M = p/(beta gamma) equal to the onset mass",
             implies="direct confirmation that remnants carry M_ext: the zero-temperature endpoint of model U"),
        dict(step=6, observation="KK excess but no black-hole onset up to the reach",
             implies="model U with c_l above the window (M_ext beyond reach), or no horizon formation at these energies"),
    ]
