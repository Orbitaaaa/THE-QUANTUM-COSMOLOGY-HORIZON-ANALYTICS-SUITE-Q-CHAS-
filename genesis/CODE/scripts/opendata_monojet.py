"""Round 10D: the model-U monojet selection on CMS Open Data (NanoAOD), for anyone with the files.

The CMS 2016 collision data are public at opendata.cern.ch (e.g. the /MET/Run2016*-UL2016_MiniAODv2_NanoAODv9-v*/NANOAOD primary datasets).
This script applies a monojet + invisible-remnant selection close to ATLAS arXiv:2102.10874 (MET > 200 GeV, 1-4 jets with pT > 30 GeV and
|eta| < 2.8 (2.4 for CMS tracking), leading jet pT > 150 GeV, lepton and photon veto, dphi(jet, MET) > 0.4) and writes the MET spectrum in
the 13 ATLAS-style bins, ready to compare with the model-U 13 TeV templates of round10.py (DATA/round10.json, lhc_reinterpretation).

Status: the Open Data portal was not reachable from the environment in which this was written, so the script has been exercised only
on synthetic NanoAOD-like arrays (``--selftest``), not on real files.

Usage:  python opendata_monojet.py file1.root [file2.root ...] [--out spectrum.json]
        python opendata_monojet.py --selftest
"""
import json
import sys

import numpy as np

EDGES = [200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1e9]


def select(arr):
    """arr: dict-like of jagged (awkward) or per-event arrays with NanoAOD branch names."""
    import awkward as ak
    jpt, jeta, jphi = arr["Jet_pt"], arr["Jet_eta"], arr["Jet_phi"]
    good = (jpt > 30) & (abs(jeta) < 2.4) & (arr["Jet_jetId"] >= 2)
    jpt, jeta, jphi = jpt[good], jeta[good], jphi[good]
    njet = ak.num(jpt)
    met, mphi = arr["MET_pt"], arr["MET_phi"]
    lead = ak.fill_none(ak.firsts(jpt), 0.0)
    lead_eta = ak.fill_none(ak.firsts(jeta), 9.9)
    dphi = abs((jphi[:, :4] - mphi + np.pi) % (2 * np.pi) - np.pi)
    dmin = ak.fill_none(ak.min(dphi, axis=1), np.pi)
    nlep = ak.sum((arr["Electron_pt"] > 10) & (abs(arr["Electron_eta"]) < 2.5) & (arr["Electron_cutBased"] >= 1), axis=1) + \
           ak.sum((arr["Muon_pt"] > 10) & (abs(arr["Muon_eta"]) < 2.4) & arr["Muon_looseId"], axis=1)
    npho = ak.sum((arr["Photon_pt"] > 15) & (abs(arr["Photon_eta"]) < 2.5) & (arr["Photon_cutBased"] >= 1), axis=1)
    m = (met > 200) & (njet >= 1) & (njet <= 4) & (lead > 150) & (abs(lead_eta) < 2.4) & (nlep == 0) & (npho == 0) & (dmin > 0.4)
    return np.asarray(met[m])


BRANCHES = ["Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId", "MET_pt", "MET_phi", "Electron_pt", "Electron_eta", "Electron_cutBased",
            "Muon_pt", "Muon_eta", "Muon_looseId", "Photon_pt", "Photon_eta", "Photon_cutBased"]


def run(files):
    import uproot
    h = np.zeros(len(EDGES) - 1)
    for f in files:
        for arr in uproot.iterate(f + ":Events", BRANCHES, step_size="200 MB"):
            h += np.histogram(select(arr), EDGES)[0]
    return h


def selftest(n=20000, seed=1):
    import awkward as ak
    rng = np.random.default_rng(seed)
    nj = rng.integers(0, 6, n)
    def jag(f):
        return ak.unflatten(f(nj.sum()), nj)
    arr = {"Jet_pt": jag(lambda k: rng.exponential(200, k) + 30), "Jet_eta": jag(lambda k: rng.normal(0, 1.5, k)),
           "Jet_phi": jag(lambda k: rng.uniform(-np.pi, np.pi, k)), "Jet_jetId": jag(lambda k: np.full(k, 6)),
           "MET_pt": rng.exponential(150, n), "MET_phi": rng.uniform(-np.pi, np.pi, n)}
    z = ak.unflatten(np.zeros(0), np.zeros(n, int))
    for b in ("Electron_pt", "Electron_eta", "Electron_cutBased", "Muon_pt", "Muon_eta", "Photon_pt", "Photon_eta", "Photon_cutBased"):
        arr[b] = z
    arr["Muon_looseId"] = ak.unflatten(np.zeros(0, bool), np.zeros(n, int))
    met = select(arr)
    return dict(n_events=n, n_selected=int(len(met)), spectrum=np.histogram(met, EDGES)[0].tolist())


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        print(json.dumps(selftest(), indent=1))
    else:
        files = [a for a in sys.argv[1:] if a.endswith(".root")]
        out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else "opendata_spectrum.json"
        h = run(files)
        json.dump(dict(edges=EDGES, counts=h.tolist(), files=files), open(out, "w"), indent=1)
        print("wrote", out)
