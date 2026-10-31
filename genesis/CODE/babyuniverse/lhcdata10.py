"""Published LHC numbers used in Round 10 (item 4: reinterpretation with the exact model-U signature).

ATLAS, "Search for new phenomena in events with an energetic jet and missing transverse momentum in pp collisions at sqrt(s) = 13 TeV",
Phys. Rev. D 103, 112006 (2021), arXiv:2102.10874, 139 fb^-1: observed events, post-fit Standard-Model prediction and its uncertainty in
the 13 exclusive (EM0-EM12) and 13 inclusive (IM0-IM12) MET regions, and the model-independent 95 % CL limits on signal events (S95),
as transcribed in the CheckMATE 2 analysis database (data/analysis_info/ATLAS_13TeV/atlas_2102_10874_ref.dat,
github.com/CheckMATE2/checkmate2).  (HEPData itself was not reachable from this environment.)
"""

ATLAS_2102_10874 = {
    # SR: (MET_lo, MET_hi, obs, bkg, bkg_err, S95_obs, S95_exp)
    "EM00": (200, 250, 1791624, 1783000, 26000, 57820, 52296), "EM01": (250, 300, 752328, 753000, 9000, 16946, 17450),
    "EM02": (300, 350, 313912, 314000, 3500, 6776, 7100), "EM03": (350, 400, 141036, 140100, 1600, 3827, 3128),
    "EM04": (400, 500, 102888, 101600, 1200, 3428, 2384), "EM05": (500, 600, 29458, 29200, 400, 1031, 862),
    "EM06": (600, 700, 10203, 10000, 180, 555, 404), "EM07": (700, 800, 3986, 3870, 80, 290, 197),
    "EM08": (800, 900, 1663, 1640, 40, 127, 111), "EM09": (900, 1000, 738, 754, 20, 58, 66),
    "EM10": (1000, 1100, 413, 359, 10, 92, 44), "EM11": (1100, 1200, 187, 182, 6, 34, 31),
    "EM12": (1200, 1e9, 207, 218, 9, 28, 35),
    "IM00": (200, 1e9, 3148643, 3120000, 40000, 100380, 75549), "IM04": (400, 1e9, 149743, 146400, 2300, 7072, 4443),
    "IM08": (800, 1e9, 3208, 3180, 130, 289, 264), "IM10": (1000, 1e9, 807, 720, 60, 191, 131),
    "IM12": (1200, 1e9, 207, 223, 19, 35, 46),
}
LUMI_FB = 139.0
EM = ["EM%02d" % i for i in range(13)]


def em_arrays():
    import numpy as np
    rows = [ATLAS_2102_10874[k] for k in EM]
    return dict(lo=np.array([r[0] for r in rows]), hi=np.array([r[1] for r in rows]), obs=np.array([r[2] for r in rows], float),
                bkg=np.array([r[3] for r in rows], float), err=np.array([r[4] for r in rows], float))
