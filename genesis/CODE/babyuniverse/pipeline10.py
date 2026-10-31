"""Round 10A: event-simulation pipeline (signal, Standard-Model backgrounds, detector response, pile-up).

Generation
* Backgrounds with Pythia 8.3 (pythia8mc): Z(-> nu nu)+jet, W(-> l nu)+jet (l = e, mu, tau), ttbar and QCD 2 -> 2, in pT-hat slices,
  with parton shower, hadronisation and decays (MPI off in the hard sample; the underlying activity is supplied by the pile-up overlay).
* Signal: model-U holes from the Round 6o/7 generator (geometric cross-section x NNPDF2.3, trapped mass, evaporation to the remnant).
  The production system (black hole + gravitational radiation, invisible, mass sqrt(s_hat)) is written to a Les Houches event file with
  the sampled x1, x2 and flavours mapped to a colour-singlet initial state (q qbar or g g, see ``_singlet_flavours``) and showered by Pythia:
  initial-state radiation, its hadronisation and the recoil of the system come from the parton shower, not from a leading-log formula.
  Hawking quanta (and, optionally, colour-discharge partons and spin-down quanta: P-R9-b, P-R6n-b) are added at object level in the hole
  frame and boosted with the showered system.
Detector (parametric, after the official Delphes FCC-hh card, M. Selvaggi): calorimetry to |eta| = 6, HCAL sigma/E = 50 %/sqrt(E) + 3 %
(|eta| < 1.7), 60 %/sqrt(E) + 3 % (to 4), 100 %/sqrt(E) + 10 % (to 6); ECAL 10 %/sqrt(E) + 0.7 %; tracking to |eta| = 6 (85-95 %);
anti-kT R = 0.4 (FastJet core in Pythia's SlowJet).  Jet response: Gaussian core plus a low-side tail (probability p_tail, a nuisance)
for QCD mismeasurement.  An ATLAS-like configuration (|eta| < 4.9, HCAL 50 %/sqrt(E) + 3 %) serves the 13 TeV validation.
Pile-up: a library of minimum-bias events (Pythia SoftQCD:nonDiffractive) overlaid with Poisson(mu) per crossing; charged pile-up removed
with efficiency eps_CHS (vertexing / 4D tracking), neutral pile-up scaled by a residual factor r_PU (PUPPI-like); the remainder enters
MET and adds noise to jets in proportion to their area.
Selections: the exact ATLAS monojet selection of arXiv:2102.10874 (as implemented in CheckMATE) and the FCC-hh model-U selection.
"""

from __future__ import annotations

import os
import tempfile
from math import pi

import numpy as np

from . import collider6l as C
from . import pipeline6o as P
from . import predictions7 as Q

INVISIBLE = {12, 14, 16, 1000022}
BH_ID = 1000022                        # an invisible, colour-neutral placeholder for the black-hole system in the LHE record


def _pythia(cmds, seed):
    import pythia8mc as p8
    p = p8.Pythia("", False)
    for c in list(cmds) + ["Print:quiet = on", "Next:numberCount = 0", "Random:setSeed = on", "Random:seed = %d" % (seed % 900000000)]:
        p.readString(c)
    if not p.init():
        raise RuntimeError("Pythia init failed: %s" % cmds)
    return p


def led_cmds(n=6, MD_PDG_TeV=21.2):
    """ADD real-graviton emission (Giudice-Rattazzi-Wells; Pythia ExtraDimensionsLED, M_D in the PDG/GRW convention)."""
    return ["ExtraDimensionsLED:monojet = on", "ExtraDimensionsLED:n = %d" % n, "ExtraDimensionsLED:MD = %g" % (MD_PDG_TeV * 1e3),
            "5000039:m0 = 5000", "5000039:mMin = 1", "5000039:mMax = %g" % (2 * MD_PDG_TeV * 1e3), "5000039:mWidth = 100000"]


def dm_cmds(M_med_TeV=20.0, m_chi_TeV=8.0, g_q=0.25, g_chi=1.0):
    """Vector-mediator dark matter, Z'(-> chi chi) + jet (Pythia DM module)."""
    return ["DM:ffbar2Zp2XXj = on", "DM:qg2Zp2XXj = on", "55:m0 = %g" % (M_med_TeV * 1e3), "52:m0 = %g" % (m_chi_TeV * 1e3),
            "Zp:gZp = 1.0", "Zp:vu = %g" % g_q, "Zp:vd = %g" % g_q, "Zp:au = 0", "Zp:ad = 0", "Zp:vX = %g" % g_chi, "Zp:aX = 0"]


PROCESSES = {
    "znunu": ["WeakBosonAndParton:qqbar2gmZg = on", "WeakBosonAndParton:qg2gmZq = on", "WeakZ0:gmZmode = 2", "23:onMode = off",
              "23:onIfAny = 12 14 16"],
    "wlnu": ["WeakBosonAndParton:qqbar2Wg = on", "WeakBosonAndParton:qg2Wq = on", "24:onMode = off", "24:onIfAny = 11 13 15"],
    "ttbar": ["Top:gg2ttbar = on", "Top:qqbar2ttbar = on"],
    "qcd": ["HardQCD:all = on"],
}


# --------------------------------------------------------------------------------------------------------------- detector

class Detector:
    def __init__(self, kind="fcchh", mu_pu=0.0, pu_lib=None, eps_chs=0.95, r_pu=0.3, p_tail=1e-3, jes=1.0, seed=1):
        self.kind = kind
        self.eta_cal = 6.0 if kind == "fcchh" else 4.9
        self.eta_trk = 6.0 if kind == "fcchh" else 2.5
        self.mu_pu, self.pu_lib, self.eps_chs, self.r_pu, self.p_tail, self.jes = mu_pu, pu_lib, eps_chs, r_pu, p_tail, jes
        self.rng = np.random.default_rng(seed)

    def jet_sigma(self, E, eta):
        a = np.abs(eta)
        if self.kind == "fcchh":
            s, c = np.where(a <= 1.7, 0.50, np.where(a <= 4.0, 0.60, 1.00)), np.where(a <= 4.0, 0.03, 0.10)
        else:
            s, c = np.where(a <= 3.2, 0.50, 1.00), np.where(a <= 3.2, 0.03, 0.10)
        return np.sqrt((c * E) ** 2 + s**2 * np.maximum(E, 0.0))

    def smear_jets(self, pt, eta):
        E = pt * np.cosh(eta)
        r = 1.0 + self.rng.normal(0.0, 1.0, len(pt)) * self.jet_sigma(E, eta) / np.maximum(E, 1e-9)
        tail = self.rng.random(len(pt)) < self.p_tail
        r = np.where(tail, r * self.rng.uniform(0.2, 0.8, len(pt)), r)
        return np.maximum(pt * r * self.jes, 0.0)

    def pileup(self):
        """Residual pile-up transverse-momentum vector and per-unit-area ET^2 sum (for jet noise)."""
        if not self.mu_pu or self.pu_lib is None:
            return np.zeros(2), 0.0
        n = self.rng.poisson(self.mu_pu)
        idx = self.rng.integers(0, len(self.pu_lib["neu"]), n)
        neu = self.pu_lib["neu"][idx].sum(0); chg = self.pu_lib["chg"][idx].sum(0)
        vec = self.r_pu * neu + (1.0 - self.eps_chs) * chg
        et2 = (self.r_pu**2 * self.pu_lib["neu_et2"][idx].sum() + (1.0 - self.eps_chs) ** 2 * self.pu_lib["chg_et2"][idx].sum())
        return vec, et2


def minbias_library(sqrt_s_TeV, n=3000, seed=11, eta_cal=6.0, eta_trk=6.0):
    """Per minimum-bias event: summed (px, py) of neutral and of tracked charged visible particles inside the calorimeter, and ET^2 sums
    per unit eta-phi area (jet noise)."""
    p = _pythia(["Beams:eCM = %g" % (sqrt_s_TeV * 1e3), "SoftQCD:nonDiffractive = on"], seed)
    neu, chg, ne2, ce2 = [], [], [], []
    area = 2 * eta_cal * 2 * pi
    for _ in range(n):
        if not p.next():
            continue
        a = np.array([(q.px(), q.py(), q.eta(), float(q.isCharged())) for q in p.event if q.isFinal() and q.isVisible()])
        if len(a) == 0:
            neu.append([0, 0]); chg.append([0, 0]); ne2.append(0.0); ce2.append(0.0); continue
        ins = np.abs(a[:, 2]) < eta_cal
        c = (a[:, 3] > 0) & (np.abs(a[:, 2]) < eta_trk) & ins
        nn = ins & ~c
        neu.append(a[nn, :2].sum(0)); chg.append(a[c, :2].sum(0))
        pt2 = a[:, 0] ** 2 + a[:, 1] ** 2
        ne2.append(pt2[nn].sum() / area); ce2.append(pt2[c].sum() / area)
    return dict(sqrt_s=sqrt_s_TeV, neu=np.array(neu), chg=np.array(chg), neu_et2=np.array(ne2), chg_et2=np.array(ce2))


# --------------------------------------------------------------------------------------------------------------- reconstruction

def _dphi(a, b):
    d = np.abs(a - b) % (2 * pi)
    return np.minimum(d, 2 * pi - d)


def reconstruct(p, sj, det: Detector, extra=None):
    """Objects from a Pythia event (and optional extra particle-level objects: list of (kind, px, py, pz, E))."""
    fin = [(q.px(), q.py(), q.pz(), q.e(), q.id(), q.isCharged(), q.isVisible(), q.eta(), q.phi(), q.pT()) for q in p.event if q.isFinal()]
    a = np.array(fin) if fin else np.zeros((0, 10))
    vis = (a[:, 6] > 0) & (np.abs(a[:, 7]) < det.eta_cal) if len(a) else np.zeros(0, bool)
    # isolated leptons and photons
    lep_e, lep_m, pho = [], [], []
    if len(a):
        cand = np.where(vis & ((np.abs(a[:, 4]) == 11) | (np.abs(a[:, 4]) == 13) | (a[:, 4] == 22)) & (a[:, 9] > 7.0))[0]
        for i in cand:
            dr = np.hypot(a[vis, 7] - a[i, 7], _dphi(a[vis, 8], a[i, 8]))
            iso = (a[vis, 9][(dr < 0.3)].sum() - a[i, 9]) < 0.15 * a[i, 9]
            if not iso:
                continue
            pid, pt, eta = abs(int(a[i, 4])), a[i, 9], a[i, 7]
            if pid == 11 and abs(eta) < (4.0 if det.kind == "fcchh" else 2.47) and det.rng.random() < 0.95:
                lep_e.append((pt, eta, a[i, 8]))
            elif pid == 13 and abs(eta) < (4.0 if det.kind == "fcchh" else 2.5) and det.rng.random() < 0.95:
                lep_m.append((pt, eta, a[i, 8]))
            elif pid == 22 and pt > 10.0 and abs(eta) < (4.0 if det.kind == "fcchh" else 2.37) and det.rng.random() < 0.9:
                pho.append((pt, eta, a[i, 8]))
    # jets (anti-kT 0.4 over visible final particles)
    sj.analyze(p.event)
    J = np.array([(sj.pT(k), sj.y(k), sj.phi(k)) for k in range(sj.sizeJet())]) if sj.sizeJet() else np.zeros((0, 3))
    # truth hadronic taus (for the tau veto): visible tau products are inside jets
    taus = [(q.eta(), q.phi(), q.pT()) for q in p.event if abs(q.id()) == 15 and q.isLastCopy() if hasattr(q, "isLastCopy")] if False else []
    true_vis = a[vis, :2].sum(0) if len(a) else np.zeros(2)
    if extra:
        for (k, px, py, pz, E) in extra:
            pt = np.hypot(px, py)
            if pt < 1e-9:
                continue
            eta = np.arcsinh(pz / pt)
            if k == "invisible" or abs(eta) > det.eta_cal:
                continue
            true_vis = true_vis + np.array([px, py])
            if k == "jet":
                J = np.vstack([J, [pt, eta, np.arctan2(py, px)]]) if len(J) else np.array([[pt, eta, np.arctan2(py, px)]])
            elif k == "lepton" and abs(eta) < 2.5:
                lep_e.append((pt, eta, np.arctan2(py, px)))
            elif k == "photon" and abs(eta) < 2.37:
                pho.append((pt, eta, np.arctan2(py, px)))
    # jet smearing + pile-up noise; MET = -(visible truth) - (jet mismeasurement) - pile-up residual
    pu_vec, pu_et2 = det.pileup()
    if len(J):
        pts = det.smear_jets(J[:, 0], J[:, 1])
        if pu_et2 > 0:
            pts = np.maximum(pts + det.rng.normal(0.0, np.sqrt(0.5 * pu_et2), len(pts)), 0.0)
        dpx = ((pts - J[:, 0]) * np.cos(J[:, 2])).sum(); dpy = ((pts - J[:, 0]) * np.sin(J[:, 2])).sum()
        J = np.column_stack([pts, J[:, 1], J[:, 2]])
        J = J[np.argsort(-J[:, 0])]
    else:
        dpx = dpy = 0.0
    met_vec = -(true_vis + np.array([dpx, dpy]) + pu_vec)
    return dict(jets=J, e=lep_e, mu=lep_m, gamma=pho, met=float(np.hypot(*met_vec)), met_phi=float(np.arctan2(met_vec[1], met_vec[0])))


def summarise(obj, det_kind="fcchh"):
    J = obj["jets"]
    if det_kind == "fcchh":
        sel = J[(J[:, 0] > 100.0) & (np.abs(J[:, 1]) < 4.0)] if len(J) else J
    else:
        sel = J[(J[:, 0] > 30.0) & (np.abs(J[:, 1]) < 2.8)] if len(J) else J
    lead = sel[0] if len(sel) else (0.0, 9.9, 0.0)
    dphimin = float(_dphi(sel[:4, 2], obj["met_phi"]).min()) if len(sel) else pi
    ST = float(sel[:, 0].sum() + sum(x[0] for x in obj["e"] + obj["mu"] + obj["gamma"]))
    return (obj["met"], float(lead[0]), float(lead[1]), len(sel), len(obj["e"]), len(obj["mu"]), len(obj["gamma"]), dphimin, ST)


FIELDS = ["met", "lead_pt", "lead_eta", "n_jets", "n_e", "n_mu", "n_gamma", "dphi_min", "ST", "weight_pb"]


# --------------------------------------------------------------------------------------------------------------- backgrounds

def background(proc, sqrt_s_TeV, pthat_lo, pthat_hi, n, det: Detector, seed=1, n_smear=1):
    """Events of one process in one pT-hat slice -> structured array (FIELDS); weight = sigma / N_accepted in pb.
    n_smear > 1 re-smears each event (independent detector draws) to populate response tails (QCD), weights divided accordingly."""
    import pythia8mc as p8
    cmds = ["Beams:eCM = %g" % (sqrt_s_TeV * 1e3), "PartonLevel:MPI = off", "PhaseSpace:pTHatMin = %g" % pthat_lo] + (
        PROCESSES[proc] if isinstance(proc, str) else list(proc))
    if pthat_hi:
        cmds.append("PhaseSpace:pTHatMax = %g" % pthat_hi)
    p = _pythia(cmds, seed)
    sj = p8.SlowJet(-1, 0.4, 20.0, det.eta_cal, 2, 1)
    rows = []
    for _ in range(n):
        if not p.next():
            continue
        for _k in range(n_smear):
            rows.append(summarise(reconstruct(p, sj, det), det.kind))
    sig_pb = p.infoPython().sigmaGen() * 1e9
    out = np.zeros(len(rows), dtype=[(f, "f8") for f in FIELDS])
    for i, f in enumerate(FIELDS[:-1]):
        out[f] = [r[i] for r in rows]
    out["weight_pb"] = sig_pb / max(len(rows), 1)
    return out, sig_pb


# --------------------------------------------------------------------------------------------------------------- signal

def _singlet_flavours(pdf, x1, x2, Q, rng):
    """Sample (f1, f2) from the PDFs and map to a colour-singlet initial state for the shower: any gluon -> g g; otherwise q qbar with the
    first quark's flavour (the shower depends on the colour charge of the incoming partons; documented approximation)."""
    from scipy.interpolate import RectBivariateSpline
    k = pdf._block(Q)
    b = pdf.blocks[k]
    cache = pdf.__dict__.setdefault("_flav_cache", {})
    fl = [f for f in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21) if (f if f != 21 else (21 if 21 in b["flavours"] else 0)) in b["flavours"]]

    def xf(f, x):
        pid = f if f != 21 else (21 if 21 in b["flavours"] else 0)
        key = (k, pid)
        if key not in cache:
            j = b["flavours"].index(pid)
            cache[key] = RectBivariateSpline(np.log(b["x"]), np.log(b["Q"] ** 2), b["xf"][:, :, j], kx=3, ky=3)
        return max(float(cache[key].ev(np.log(x), np.log(Q**2))), 0.0)

    w1 = np.array([xf(f, x1) for f in fl])
    w2 = np.array([xf(f, x2) for f in fl])
    f1 = fl[rng.choice(len(fl), p=w1 / w1.sum())]; f2 = fl[rng.choice(len(fl), p=w2 / w2.sum())]
    if f1 == 21 or f2 == 21:
        return 21, 21
    return f1, -f1


def _lhe(events, sqrt_s_GeV, path):
    with open(path, "w") as fh:
        fh.write('<LesHouchesEvents version="1.0">\n<init>\n2212 2212 %.6e %.6e 0 0 247000 247000 -4 1\n1.0 0.0 1.0 1\n</init>\n'
                 % (sqrt_s_GeV / 2, sqrt_s_GeV / 2))
        for (f1, f2, x1, x2, M) in events:
            E1, E2 = x1 * sqrt_s_GeV / 2, x2 * sqrt_s_GeV / 2
            c = (501, 502, 502, 501) if f1 == 21 else ((501, 0, 0, 501) if f1 > 0 else (0, 501, 501, 0))
            fh.write("<event>\n3 1 1.0 %.6e 0.0078 0.1\n" % M)
            fh.write("%d -1 0 0 %d %d 0 0 %.8e %.8e 0 0 9\n" % (f1, c[0], c[1], E1, E1))
            fh.write("%d -1 0 0 %d %d 0 0 %.8e %.8e 0 0 9\n" % (f2, c[2], c[3], -E2, E2))
            fh.write("%d 1 1 2 0 0 0 0 %.8e %.8e %.8e 0 9\n" % (BH_ID, E1 - E2, E1 + E2, M))
            fh.write("</event>\n")
        fh.write("</LesHouchesEvents>\n")


def _boost(p4, beta):
    return Q.LabGenerator._boost(p4, beta)


VIS_KIND = {"quark": "jet", "gluon": "jet", "top": "jet", "W": "jet", "Z": "jet", "higgs": "jet", "charged_lepton": "lepton",
            "photon": "photon", "neutrino": "invisible", "final": "jet"}


def signal(lumi: C.Lumi, sqrt_s_TeV=100.0, n=6, M_D_TeV=8.0, c_l=0.6, kind="hayward", n_events=2000, det: Detector = None, seed=7,
           discharge=False, spin=False, alpha_s=0.09, f_charged=0.0, M_th_TeV=None, trapped="conservative"):
    """Model-U (or GR, kind='gr') signal through Pythia's shower and the detector.  discharge: add one colour-discharge parton of energy
    U(0.3, 1) x alpha_s C_F / r_+ in the hole frame (P-R9-b); spin: add the Myers-Perry rotational energy as visible quanta (P-R6n-b)."""
    import pythia8mc as p8
    from . import openprob10 as O
    det = det or Detector("fcchh", seed=seed)
    rng = np.random.default_rng(seed)
    g = P.Generator(lumi, sqrt_s_TeV, n, M_D_TeV, c_l, kind, trapped=trapped, seed=seed, M_th_TeV=M_th_TeV)
    g.keep_quanta = True
    lab = Q.LabGenerator(g, isr=False, seed=seed)
    pdf = lab.pdf
    s = sqrt_s_TeV * 1e3
    empty = np.zeros(0, dtype=[(f, "f8") for f in FIELDS] + [("beta_rem", "f8"), ("hscp", "f8")])
    if (g.M_th / g.y) ** 2 / s**2 >= lumi.tau_max:                       # kinematically closed: no black holes
        return empty, 0.0, dict(M_ext_TeV=g.M_ext / 1e3, M_th_TeV=g.M_th / 1e3)
    evs, info = [], []
    n_try = 0
    while len(evs) < n_events and n_try < 200 * n_events:
        n_try += 1
        e = g.event()
        if not e["horizon"]:
            continue
        tau = e["tau"]; y = lab._rapidity(tau, e["M"])
        x1, x2 = np.sqrt(tau) * np.exp(y), np.sqrt(tau) * np.exp(-y)
        if x1 >= 1 or x2 >= 1:
            continue
        f1, f2 = _singlet_flavours(pdf, x1, x2, e["sqrt_shat"], rng) if pdf is not None else (21, 21)
        evs.append((f1, f2, x1, x2, e["sqrt_shat"])); info.append(e)
    tmp = tempfile.NamedTemporaryFile(suffix=".lhe", delete=False); tmp.close()
    _lhe(evs, s, tmp.name)
    p = _pythia(["Beams:frameType = 4", "Beams:LHEF = %s" % tmp.name, "PartonLevel:MPI = off", "%d:mayDecay = off" % BH_ID,
                 "%d:mMin = 0" % BH_ID, "%d:mMax = 200000" % BH_ID, "Check:event = off"], seed)
    sj = p8.SlowJet(-1, 0.4, 20.0, det.eta_cal, 2, 1)
    rows, betas, hscp = [], [], []
    rp_l = None
    for e in info:
        if not p.next():
            continue
        sysP = None
        for q in p.event:
            if q.id() == BH_ID and q.isFinal():
                sysP = np.array([q.e(), q.px(), q.py(), q.pz()])
                break
        if sysP is None:
            continue
        beta = sysP[1:] / sysP[0]
        extra = []
        tot = np.zeros(4)
        for (name, E, vis, px, py, pz) in e.get("quanta", []):
            q4 = _boost(np.array([E, px, py, pz]), beta)
            tot += q4
            extra.append((VIS_KIND.get(name, "jet") if vis else "invisible", q4[1], q4[2], q4[3], q4[0]))
        if discharge and kind != "gr" and e["remnant"] > 0:
            if rp_l is None:
                rp_l = O.remnant_discharge(n, M_D_TeV, c_l, kind)["E_discharge_GeV"]
            Ed = rng.uniform(0.3, 1.0) * rp_l
            ct = rng.uniform(-1, 1); ph = rng.uniform(0, 2 * pi); st = np.sqrt(1 - ct**2)
            q4 = _boost(np.array([Ed, Ed * st * np.cos(ph), Ed * st * np.sin(ph), Ed * ct]), beta)
            tot += q4; extra.append(("jet", q4[1], q4[2], q4[3], q4[0]))
        if spin and kind != "gr":
            Er = float(O.rotational_fraction(e["a_star"], n)) * e["M"]
            nq = max(1, int(round(Er / 1500.0)))
            for _ in range(nq):
                Eq = Er / nq
                ct = rng.uniform(-0.5, 0.5); ph = rng.uniform(0, 2 * pi); st = np.sqrt(1 - ct**2)   # superradiance favours the equator
                q4 = _boost(np.array([Eq, Eq * st * np.cos(ph), Eq * st * np.sin(ph), Eq * ct]), beta)
                tot += q4; extra.append(("jet", q4[1], q4[2], q4[3], q4[0]))
        rem = sysP - tot
        b_rem = float(np.linalg.norm(rem[1:]) / max(rem[0], 1e-9)) if e["remnant"] > 0 else 0.0
        betas.append(b_rem)
        hscp.append(bool(e["remnant"] > 0 and rng.random() < f_charged))
        rows.append(summarise(reconstruct(p, sj, det, extra), det.kind))
    os.unlink(tmp.name)
    out = np.zeros(len(rows), dtype=[(f, "f8") for f in FIELDS] + [("beta_rem", "f8"), ("hscp", "f8")])
    for i, f in enumerate(FIELDS[:-1]):
        out[f] = [r[i] for r in rows]
    out["weight_pb"] = g.sigma_pb / max(len(rows), 1)
    out["beta_rem"] = betas; out["hscp"] = hscp
    return out, g.sigma_pb, dict(M_ext_TeV=g.M_ext / 1e3, M_th_TeV=g.M_th / 1e3)


# --------------------------------------------------------------------------------------------------------------- selections

ATLAS_EM_EDGES = [200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1e9]


def atlas_monojet(ev):
    """arXiv:2102.10874 signal region: MET > 200, 1-4 jets (pT > 30, |eta| < 2.8), lepton/photon veto, dphi(jet, MET) > 0.4 (0.6 below
    250 GeV), leading jet pT > 150 with |eta| < 2.4.  (The tau veto is approximated by the lepton veto plus the jet multiplicity.)"""
    dcut = np.where(ev["met"] < 250, 0.6, 0.4)
    return ((ev["met"] > 200) & (ev["n_jets"] >= 1) & (ev["n_jets"] <= 4) & (ev["n_e"] + ev["n_mu"] == 0) & (ev["n_gamma"] == 0)
            & (ev["dphi_min"] > dcut) & (ev["lead_pt"] > 150) & (np.abs(ev["lead_eta"]) < 2.4))


FCC_EDGES = [1000, 2000, 3000, 4000, 6000, 8000, 11000, 15000, 1e9]


def fcc_monojet(ev, lead_min=500.0):
    """FCC-hh model-U selection: MET > 1 TeV, leading jet pT > lead_min (|eta| < 4), lepton and photon veto, dphi(jets, MET) > 0.4."""
    return ((ev["met"] > 1000) & (ev["n_jets"] >= 1) & (ev["n_e"] + ev["n_mu"] == 0) & (ev["n_gamma"] == 0) & (ev["dphi_min"] > 0.4)
            & (ev["lead_pt"] > lead_min))


def binned(ev, mask, edges, L_fb):
    """Expected events per MET bin (L in fb^-1) and the MC statistical uncertainty."""
    w = ev["weight_pb"] * 1e3 * L_fb
    idx = np.digitize(ev["met"][mask], edges) - 1
    n = np.bincount(idx[(idx >= 0) & (idx < len(edges) - 1)], weights=w[mask][(idx >= 0) & (idx < len(edges) - 1)], minlength=len(edges) - 1)
    e2 = np.bincount(idx[(idx >= 0) & (idx < len(edges) - 1)], weights=(w[mask] ** 2)[(idx >= 0) & (idx < len(edges) - 1)], minlength=len(edges) - 1)
    return n, np.sqrt(e2)
