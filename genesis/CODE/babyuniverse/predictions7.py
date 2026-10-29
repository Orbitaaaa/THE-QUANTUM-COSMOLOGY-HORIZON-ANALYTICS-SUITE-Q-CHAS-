"""Round 7: experimentally distinguishable predictions of model U.

Lab-frame events: the Round 6o generator plus (i) the hole's rapidity from the parton momentum fractions (x1 drawn from F(x1) F(tau/x1)/x1
with the NNPDF grid; symmetric fallback without it), (ii) hard initial-state radiation, P(pT,ISR > p) = 1 - exp(-A ln^2(sqrt(s_hat)/p)) with
A = C_A alpha_s/(2 pi) ~ 0.043 (leading-log estimate), against which the hole recoils, and (iii) a generic detector: |eta| < 4,
pT > 50 GeV; jets, leptons and photons visible; W, Z, H, t visible at 80 % of their energy; neutrinos and the remnant invisible (neutral
remnant; a charged fraction f_ch gives heavy stable charged tracks with beta = p/E).

Observables: S_T (scalar sum of visible pT), missing transverse momentum MET, number of visible objects, leading-object pT, MET/S_T, the
remnant's velocity.  Discrimination: expected log-likelihood ratio per event (Kullback-Leibler divergence of binned distributions) between
model U and a GR hole at the same threshold and cross-section, and the number of events for a 5 sigma separation (n KL = 12.5).

Backgrounds (order of magnitude, stated as such): QCD 2 -> 2 with sigma_hat = K pi alpha_s^2 / s_hat (K = 10, |eta| < 2.5) for the S_T tail;
Z(-> nu nu) + jet as 2 x 10^-3 of the QCD tail at the same pT for the MET tail.
"""

from __future__ import annotations

import os
from math import pi

import numpy as np

from . import collider6l as C
from . import pipeline6o as P

VISIBLE_FRAC = {"quark": 1.0, "top": 0.8, "gluon": 1.0, "W": 0.8, "Z": 0.8, "photon": 1.0, "higgs": 0.8, "charged_lepton": 1.0,
                "neutrino": 0.0, "final": 1.0}


class LabGenerator:
    def __init__(self, gen: P.Generator, pdf_path=None, isr=True, f_charged=0.0, eta_max=4.0, pt_min=50.0, seed=5):
        self.g = gen
        self.g.keep_quanta = True
        self.rng = np.random.default_rng(seed)
        path = pdf_path or C.PDF_SEARCH[0]
        self.pdf = C.PDF(path) if os.path.exists(path) else None
        self.isr, self.f_ch, self.eta_max, self.pt_min = isr, f_charged, eta_max, pt_min

    def _rapidity(self, tau, Q):
        lim = 0.5 * np.log(1 / tau)
        if self.pdf is None:
            return self.rng.uniform(-1, 1) * lim * self.rng.random()
        lx = np.linspace(np.log(tau), 0, 200)
        x1 = np.exp(lx)
        w = self.pdf.xF(x1, Q) / x1 * self.pdf.xF(np.minimum(tau / x1, 1.0), Q) / (tau / x1) * x1 / x1
        w = np.maximum(w, 0)
        c = np.cumsum(w)
        c = c / c[-1]
        x = float(np.exp(np.interp(self.rng.random(), c, lx)))
        return 0.5 * np.log(x * x / tau)

    def _isr(self, shat):
        if not self.isr:
            return 0.0
        A = 0.043
        u = self.rng.random()
        # invert P(> p) = 1 - exp(-A ln^2(shat/p)) for p in [50 GeV, 0.3 shat]
        L = np.sqrt(max(-np.log(1 - u), 0) / A)
        p = shat * np.exp(-L)
        return float(p) if 50.0 < p < 0.3 * shat else 0.0

    @staticmethod
    def _boost(p4, beta):
        b2 = beta @ beta
        if b2 < 1e-20:
            return p4
        g = 1 / np.sqrt(1 - b2)
        bp = beta @ p4[1:]
        E = g * (p4[0] + bp)
        pv = p4[1:] + ((g - 1) * bp / b2 + g * p4[0]) * beta
        return np.concatenate([[E], pv])

    def event(self):
        e = self.g.event()
        shat = e["sqrt_shat"]
        M = e["M"]
        y = self._rapidity(e["tau"], M)
        ptI = self._isr(shat)
        phiI = 2 * pi * self.rng.random()
        objs = []
        if ptI > 0:
            etaI = self.rng.uniform(-3, 3)
            objs.append(("isr_jet", ptI, etaI, phiI))
        # hole four-momentum in the lab: transverse recoil against ISR, longitudinal rapidity y
        pt_bh = np.array([-ptI * np.cos(phiI), -ptI * np.sin(phiI)])
        mT = np.sqrt(M**2 + pt_bh @ pt_bh)
        P_bh = np.array([mT * np.cosh(y), pt_bh[0], pt_bh[1], mT * np.sinh(y)])
        beta = P_bh[1:] / P_bh[0]
        vis_sum = np.zeros(4)
        for (name, E, vis, px, py, pz) in e.get("quanta", []):
            p4 = self._boost(np.array([E, px, py, pz]), beta)
            frac = VISIBLE_FRAC.get(name, 1.0)
            if frac <= 0:
                continue
            pt = np.hypot(p4[1], p4[2]) * frac
            pabs = np.linalg.norm(p4[1:])
            eta = np.arctanh(np.clip(p4[3] / max(pabs, 1e-12), -0.999999, 0.999999))
            if abs(eta) < self.eta_max and pt > self.pt_min:
                objs.append((name, pt, eta, np.arctan2(p4[2], p4[1])))
            vis_sum += p4 * frac
        # remnant
        P_rem = P_bh - sum((self._boost(np.array([E, px, py, pz]), beta) for (_, E, _, px, py, pz) in e.get("quanta", [])),
                           np.zeros(4)) if e.get("quanta") else P_bh
        rem_mass = e.get("remnant", 0.0)
        beta_rem = float(np.linalg.norm(P_rem[1:]) / max(P_rem[0], 1e-9)) if rem_mass > 0 else 0.0
        charged = bool(rem_mass > 0 and self.rng.random() < self.f_ch)
        px = sum(o[1] * np.cos(o[3]) for o in objs)
        py = sum(o[1] * np.sin(o[3]) for o in objs)
        ST = sum(o[1] for o in objs)
        MET = float(np.hypot(px, py))
        lead = max((o[1] for o in objs), default=0.0)
        return dict(M=M, N_obj=len(objs), ST=ST, MET=MET, lead=lead, met_frac=MET / ST if ST > 0 else 0.0, isr=ptI > 0, N_quanta=e.get("N", 0),
                    remnant=rem_mass, beta_remnant=beta_rem, hscp=charged and abs(np.arctanh(np.clip(P_rem[3] / max(np.linalg.norm(P_rem[1:]), 1e-9), -0.9999, 0.9999))) < 2.5)

    def run(self, n=3000) -> dict:
        ev = [self.event() for _ in range(n)]
        a = lambda k: np.array([x[k] for x in ev], dtype=float)
        return dict(kind=self.g.kind, sigma_pb=self.g.sigma_pb, M_ext_TeV=self.g.M_ext / 1e3, M_th_TeV=self.g.M_th / 1e3, n=n,
                    M=a("M").tolist(), N_obj=a("N_obj").tolist(), ST=a("ST").tolist(), MET=a("MET").tolist(), lead=a("lead").tolist(),
                    met_frac=a("met_frac").tolist(), isr=a("isr").tolist(), N_quanta=a("N_quanta").tolist(), remnant=a("remnant").tolist(),
                    beta_remnant=a("beta_remnant").tolist(), hscp=a("hscp").tolist())


def summary(r: dict) -> dict:
    ST, MET, N = np.array(r["ST"]), np.array(r["MET"]), np.array(r["N_obj"])
    return dict(kind=r["kind"], sigma_pb=r["sigma_pb"], frac_no_visible=float(np.mean(N == 0)), frac_isr_only=float(np.mean((N == 1) & (np.array(r["isr"]) > 0) & (np.array(r["N_quanta"]) == 0))),
                mean_N_obj=float(N.mean()), median_ST_TeV=float(np.median(ST) / 1e3), mean_MET_TeV=float(MET.mean() / 1e3),
                frac_MET_gt_1TeV=float(np.mean(MET > 1e3)), frac_ST_gt_5TeV=float(np.mean(ST > 5e3)),
                mean_met_frac=float(np.mean(np.array(r["met_frac"])[ST > 0])) if np.any(ST > 0) else 0.0,
                mean_beta_remnant=float(np.mean(np.array(r["beta_remnant"])[np.array(r["remnant"]) > 0])) if np.any(np.array(r["remnant"]) > 0) else 0.0)


def kl_separation(rA: dict, rB: dict) -> dict:
    """Binned (N_obj, MET/S_T) distributions; expected log-likelihood ratio per event and events for 5 sigma (n KL = 12.5)."""
    def hist(r):
        N = np.minimum(np.array(r["N_obj"]), 6)
        f = np.array(r["met_frac"])
        fb = np.digitize(f, [0.05, 0.2, 0.5, 0.8])
        H = np.zeros((7, 5))
        for a, b in zip(N.astype(int), fb):
            H[a, b] += 1
        H = H + 0.5
        return H / H.sum()
    A, B = hist(rA), hist(rB)
    kl_ab = float(np.sum(A * np.log(A / B)))
    kl_ba = float(np.sum(B * np.log(B / A)))
    return dict(kl_A_vs_B=kl_ab, kl_B_vs_A=kl_ba, events_5sigma=float(12.5 / max(min(kl_ab, kl_ba), 1e-12)))


# --------------------------------------------------------------------------- backgrounds (order of magnitude)


def qcd_tail_pb(lumi: C.Lumi, sqrt_s_TeV, ST_min_TeV, K=10.0, alpha_s=0.08) -> float:
    """sigma(QCD 2 -> 2, m > ST_min) with sigma_hat = K pi alpha_s^2 / s_hat (central jets); S_T ~ m for central 2 -> 2."""
    s = (sqrt_s_TeV * 1e3) ** 2
    tau_min = (ST_min_TeV * 1e3) ** 2 / s
    if tau_min >= lumi.tau_max:
        return 0.0
    taus = np.geomspace(tau_min, lumi.tau_max, 400)
    shat = taus * s
    val = lumi(taus, np.sqrt(shat)) * K * pi * alpha_s**2 / shat * taus
    return float(np.trapezoid(val, np.log(taus)) * C.GEV2_TO_PB)


def znunu_met_tail_pb(lumi: C.Lumi, sqrt_s_TeV, MET_min_TeV) -> float:
    return 2e-3 * qcd_tail_pb(lumi, sqrt_s_TeV, 2 * MET_min_TeV)


def predictions_table(lumi: C.Lumi, n_events=3000, seed=21) -> dict:
    """Benchmarks at 100 TeV (n = 6, M_D = 8 TeV): model U (Hayward c_l = 0.6, 0.7; Dymnikova c_l = 0.7) and GR at the same threshold."""
    out = {}
    for name, kind, cl, mth in (("U-hayward-0.6", "hayward", 0.6, None), ("U-hayward-0.7", "hayward", 0.7, None),
                                ("U-dymnikova-0.7", "dymnikova", 0.7, None), ("GR-11.5", "gr", 0.6, 11.5), ("GR-33.8", "gr", 0.7, 33.8)):
        g = P.Generator(lumi, 100.0, 6, 8.0, cl, kind, seed=seed, M_th_TeV=mth)
        lab = LabGenerator(g, seed=seed)
        r = lab.run(n_events)
        out[name] = dict(summary=summary(r), raw=r)
    out["separation"] = {"U-hayward-0.6 vs GR-11.5": kl_separation(out["U-hayward-0.6"]["raw"], out["GR-11.5"]["raw"]),
                         "U-dymnikova-0.7 vs GR-11.5": kl_separation(out["U-dymnikova-0.7"]["raw"], out["GR-11.5"]["raw"]),
                         "U-hayward-0.6 vs U-dymnikova-0.7": kl_separation(out["U-hayward-0.6"]["raw"], out["U-dymnikova-0.7"]["raw"]),
                         "U-hayward-0.7 vs GR-33.8": kl_separation(out["U-hayward-0.7"]["raw"], out["GR-33.8"]["raw"])}
    return out
