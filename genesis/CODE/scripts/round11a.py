"""Round 11 (G11-10, G11-11): systematic-variation samples at 100 TeV, cached in DATA/round10_events (not committed).
Background variations (Z(nunu)+jet, W(l nu)+jet): renormalisation and factorisation scales x0.5 and x2; PDF sets CTEQ6L1 (pSet 8) and
NNPDF3.1 LO (pSet 17) instead of the default NNPDF2.3 QCD+QED LO (pSet 13).  Detector variations (backgrounds and the model-U signal):
jet-response tail probability x10; pile-up suppression degraded (eps_CHS 0.95 -> 0.80, neutral residual 0.3 -> 1.0).
Usage: python round11a.py [--workers 3]"""
import os
import sys
import time
import warnings
import multiprocessing as mp

import numpy as np

from _common import DATA
from babyuniverse import pipeline10 as P10
from babyuniverse import collider6l as C
import round10a as R

warnings.simplefilter("ignore")
EV = R.EV
GEN_VARIANTS = {"scaleUp": ["SigmaProcess:renormMultFac = 2.0", "SigmaProcess:factorMultFac = 2.0"],
                "scaleDn": ["SigmaProcess:renormMultFac = 0.5", "SigmaProcess:factorMultFac = 0.5"],
                "pdfCTEQ6L1": ["PDF:pSet = 8"], "pdfNNPDF31": ["PDF:pSet = 17"]}
DET_VARIANTS = {"tail10": dict(p_tail=1e-2), "puWorse": dict(eps_chs=0.80, r_pu=1.0)}


def det_variant(seed, **over):
    lib = R.pu_lib(100.0)
    kw = dict(mu_pu=1000.0, pu_lib=lib, eps_chs=0.95, r_pu=0.3, seed=seed)
    kw.update(over)
    return P10.Detector("fcchh", **kw)


def job(spec):
    name, kind, args = spec
    out = os.path.join(EV, name + ".npz")
    if os.path.exists(out):
        return name, "cached"
    t = time.time()
    if kind == "bkg":
        cmds, lo, hi, n, seed, dv = args
        ev, sig = P10.background(cmds, 100.0, lo, hi, n, det_variant(seed, **dv), seed=seed)
        np.savez_compressed(out, ev=ev, sigma_pb=sig)
    else:
        cl, nev, seed, dv = args
        L = C.load_lumi(os.path.join(DATA, "round6l_lumi.json"))
        ev, sig, info = P10.signal(L, 100.0, 6, 8.0, cl, "hayward", nev, det_variant(seed, **dv), seed=seed)
        np.savez_compressed(out, ev=ev, sigma_pb=sig, M_ext_TeV=info["M_ext_TeV"], M_th_TeV=info["M_th_TeV"])
    return name, "%.0f s" % (time.time() - t)


def specs():
    S = []
    for var, extra in GEN_VARIANTS.items():
        for proc in ("znunu", "wlnu"):
            for k, (lo, hi) in enumerate(R.VJ_FCC):
                S.append(("b100%s_%s_%d" % (var, proc, k), "bkg", (P10.PROCESSES[proc] + extra, lo, hi, 3000, 1200 + k, {})))
    for var, dv in DET_VARIANTS.items():
        for proc in ("znunu", "wlnu"):
            for k, (lo, hi) in enumerate(R.VJ_FCC):
                S.append(("b100%s_%s_%d" % (var, proc, k), "bkg", (P10.PROCESSES[proc], lo, hi, 3000, 1300 + k, dv)))
        for cl in (0.58, 0.62, 0.66, 0.68, 0.69, 0.70, 0.71, 0.72):
            S.append(("s100%s_U_%.2f" % (var, cl), "sig", (cl, 1500, 1400 + int(cl * 100), dv)))
    return S


if __name__ == "__main__":
    nw = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 3
    R.pu_lib(100.0)
    S = specs()
    print(len(S), "jobs", flush=True)
    with mp.Pool(nw) as pool:
        for name, status in pool.imap_unordered(job, S):
            print(name, status, flush=True)
