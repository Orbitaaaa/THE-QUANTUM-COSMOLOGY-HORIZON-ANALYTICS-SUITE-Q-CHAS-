"""Round 10A: produce the event samples (cached in DATA/round10_events, not committed): 13 TeV ATLAS-validation backgrounds with mu = 34
pile-up, 85 and 100 TeV FCC-hh backgrounds with mu = 1000, model-U signal grids and variants, and alternative hypotheses.
Usage: python round10a.py [--only 13|100|85|signal|alt|lhc] [--workers 4]"""
import os
import sys
import time
import pickle
import warnings
import multiprocessing as mp

import numpy as np

from _common import DATA
from babyuniverse import pipeline10 as P10
from babyuniverse import collider6l as C

warnings.simplefilter("ignore")
EV = os.path.join(DATA, "round10_events")
os.makedirs(EV, exist_ok=True)
VJ_FCC = [(800, 2000), (2000, 4000), (4000, 8000), (8000, 0)]


def pu_lib(rs):
    f = os.path.join(EV, "minbias_%g.pkl" % rs)
    if not os.path.exists(f):
        lib = P10.minbias_library(rs, 3000, 11, eta_cal=6.0 if rs > 20 else 4.9, eta_trk=6.0 if rs > 20 else 2.5)
        pickle.dump(lib, open(f, "wb"))
    return pickle.load(open(f, "rb"))


def det_for(rs, seed, **kw):
    if rs < 20:
        return P10.Detector("atlas", mu_pu=34.0, pu_lib=pu_lib(rs), eps_chs=0.9, r_pu=0.5, seed=seed, **kw)
    return P10.Detector("fcchh", mu_pu=1000.0, pu_lib=pu_lib(rs), eps_chs=0.95, r_pu=0.3, seed=seed, **kw)


def job(spec):
    name, kind, args = spec
    out = os.path.join(EV, name + ".npz")
    if os.path.exists(out):
        return name, "cached"
    t = time.time()
    if kind == "bkg":
        proc, rs, lo, hi, n, nsm, seed = args
        ev, sig = P10.background(proc, rs, lo, hi, n, det_for(rs, seed), seed=seed, n_smear=nsm)
        np.savez_compressed(out, ev=ev, sigma_pb=sig)
    elif kind == "sig":
        rs, n, MD, cl, kd, nev, seed, extra = args
        L = C.load_lumi(os.path.join(DATA, "round6l_lumi.json"))
        ev, sig, info = P10.signal(L, rs, n, MD, cl, kd, nev, det_for(rs, seed), seed=seed, **extra)
        np.savez_compressed(out, ev=ev, sigma_pb=sig, M_ext_TeV=info["M_ext_TeV"], M_th_TeV=info["M_th_TeV"])
    return name, "%.0f s" % (time.time() - t)


def specs(which):
    S = []
    if which in ("all", "13"):
        for proc, sl, n, nsm in (("znunu", [(150, 300), (300, 600), (600, 1000), (1000, 0)], 5000, 1),
                                 ("wlnu", [(150, 300), (300, 600), (600, 1000), (1000, 0)], 5000, 1),
                                 ("ttbar", [(100, 300), (300, 600), (600, 0)], 3000, 1),
                                 ("qcd", [(250, 500), (500, 1000), (1000, 0)], 2000, 10)):
            for k, (lo, hi) in enumerate(sl):
                S.append(("b13_%s_%d" % (proc, k), "bkg", (proc, 13.0, lo, hi, n, nsm, 100 + k)))
    for rs in (100.0, 85.0):
        if which in ("all", "%g" % rs):
            for proc, sl, n, nsm in (("znunu", VJ_FCC, 4000, 1), ("wlnu", VJ_FCC, 4000, 1), ("ttbar", [(800, 3000), (3000, 0)], 3000, 1),
                                     ("qcd", [(1500, 4000), (4000, 8000), (8000, 0)], 2000, 10)):
                for k, (lo, hi) in enumerate(sl):
                    S.append(("b%g_%s_%d" % (rs, proc, k), "bkg", (proc, rs, lo, hi, n, nsm, 200 + k)))
    if which in ("all", "signal"):
        for rs in (100.0, 85.0):
            for cl in (0.58, 0.62, 0.64, 0.66, 0.68, 0.69, 0.70, 0.71, 0.72, 0.74):
                S.append(("s%g_U_%.2f" % (rs, cl), "sig", (rs, 6, 8.0, cl, "hayward", 1500, 300 + int(cl * 100), {})))
            for cl in (0.62, 0.70):
                S.append(("s%g_Udis_%.2f" % (rs, cl), "sig", (rs, 6, 8.0, cl, "hayward", 1500, 400 + int(cl * 100), dict(discharge=True))))
                S.append(("s%g_Uspin_%.2f" % (rs, cl), "sig", (rs, 6, 8.0, cl, "hayward", 1500, 500 + int(cl * 100), dict(spin=True))))
            for cl in (0.73, 0.82):
                S.append(("s%g_Udym_%.2f" % (rs, cl), "sig", (rs, 6, 8.0, cl, "dymnikova", 1500, 600 + int(cl * 100), {})))
    if which in ("all", "alt"):
        L = None
        for cl, mth in ((0.62, 14.36), (0.70, 33.76)):          # GR holes at the same threshold as model U (M_ext of the c_l point)
            S.append(("a100_GR_%.2f" % cl, "sig", (100.0, 6, 8.0, cl, "gr", 1500, 700 + int(cl * 100), dict(M_th_TeV=mth))))
        for k, (lo, hi) in enumerate([(1000, 3000), (3000, 6000), (6000, 0)]):
            S.append(("a100_LED_%d" % k, "bkg", (P10.led_cmds(6, 21.2), 100.0, lo, hi, 3000, 1, 800 + k)))
        for Mm, mx in ((20.0, 8.0), (40.0, 15.0), (12.0, 5.0)):
            for k, (lo, hi) in enumerate([(1000, 4000), (4000, 0)]):
                S.append(("a100_DM_%g_%g_%d" % (Mm, mx, k), "bkg", (P10.dm_cmds(Mm, mx), 100.0, lo, hi, 2500, 1, 900 + k)))
    if which in ("all", "lhc"):
        for n_, MD in ((6, 2.0), (6, 4.0), (4, 4.0)):
            for cl in np.linspace(0.30, 0.80, 11):
                for trap in ("conservative", "optimistic"):
                    S.append(("l13_n%d_MD%g_%.2f_%s" % (n_, MD, cl, trap[:3]), "sig",
                              (13.0, n_, MD, float(cl), "hayward", 800, 1000 + int(cl * 100), dict(trapped=trap))))
    return S


if __name__ == "__main__":
    which = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else "all"
    nw = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 4
    for rs in ((13.0,) if which == "13" else (100.0, 85.0, 13.0) if which in ("all", "lhc") else (float(which),) if which in ("100", "85") else (100.0, 85.0)):
        pu_lib(rs)
    S = specs(which)
    print(len(S), "jobs", flush=True)
    with mp.Pool(nw) as pool:
        for name, status in pool.imap_unordered(job, S):
            print(name, status, flush=True)
