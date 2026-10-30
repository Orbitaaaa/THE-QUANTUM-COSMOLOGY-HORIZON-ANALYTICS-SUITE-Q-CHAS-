"""Round 8: model U against current and future data.

* ``convention``  The fundamental-scale conventions.  This project uses the Dimopoulos-Landsberg form of the horizon radius, equivalent to
  G_D = 1/M_D^{n+2} (checked symbolically in Round 6p).  Experiments quote M_D in the PDG / Giudice-Rattazzi-Wells convention,
  G_D = (2 pi)^n / (8 pi M_D^{n+2}).  Hence M_D(DL) = M_D(PDG) [8 pi / (2 pi)^n]^{1/(n+2)} = 0.89, 0.66, 0.50, 0.43, 0.38 M_D(PDG) for
  n = 2..6.  [Correction to Rounds 6k-6l, which compared the two directly (F64).]
* ``monojet_bounds``  CMS (arXiv:1712.02345, 35.9 fb^-1) observed M_D(PDG) > 9.9, 7.5, 6.3, 5.7, 5.3 TeV for n = 2..6, i.e. M_D(DL) > 8.8,
  4.9, 3.2, 2.4, 2.0 TeV.
* ``lhc_model_u``  Model-U holes at 13 TeV: kinematically absent for M_ext > y_n x 13 TeV (conservative trapping, 7.6 TeV at n = 6); with full
  trapping their events would sit in the monojet high-recoil region, bounded by the ATLAS model-independent limit sigma x A x eps < 0.3 fb in
  the tightest inclusive region (arXiv:2102.10874, 139 fb^-1).  The multi-object black-hole searches (high S_T, >= 3 objects) are blind to
  near-threshold model-U holes.
* ``remnant_safety``  [Proposition-N 8-1] A remnant captured in a neutron star accretes at ~n sigma v ~ 1e13 GeV/s; any excess above
  extremality is re-radiated (T > 0), and the equilibrium excess is ~1e-9 of the remnant mass at a temperature ~0.1 GeV: no runaway growth.
* ``future_reach``  HL-LHC, FCC-hh 85 TeV and a 100 TeV collider: the range of M_ext with >= 10 (discovery-level, background-free topology)
  or >= 3 (evidence) model-U events.
"""

from __future__ import annotations

import os
from math import pi

import numpy as np
from scipy.optimize import brentq

from . import collider6l as C
from . import unified6l as U

CMS_MONOJET_PDG_TeV = {2: 9.9, 3: 7.5, 4: 6.3, 5: 5.7, 6: 5.3}
ATLAS_MONOJET_LIMIT_FB = 0.3
HBAR_GEV_S = 6.582119569e-25
HBARC_GEV_M = 1.97327e-16


def dl_over_pdg(n):
    return (8 * pi / (2 * pi) ** n) ** (1 / (n + 2))


def convention() -> list:
    return [dict(n=n, MD_DL_over_MD_PDG=float(dl_over_pdg(n))) for n in range(1, 8)]


def monojet_bounds() -> list:
    return [dict(n=n, MD_PDG_TeV=v, MD_DL_TeV=float(v * dl_over_pdg(n))) for n, v in CMS_MONOJET_PDG_TeV.items()]


def lhc_model_u(lumi: C.Lumi, n=6, M_D_TeV=8.0, kinds=("hayward", "dymnikova"), sqrt_s=13.0, L_fb=139.0) -> list:
    rows = []
    y = C.trapped_fraction_headon(n)
    for kind in kinds:
        for cl in np.linspace(0.3, 0.8, 26):
            Mext = C.thresholds(n, M_D_TeV, cl, kind)["model_U_TeV"]
            s_cons = C.sigma_pp(lumi, sqrt_s, n, M_D_TeV, Mext, "conservative")
            s_opt = C.sigma_pp(lumi, sqrt_s, n, M_D_TeV, Mext, "optimistic")
            rows.append(dict(kind=kind, c_l=float(cl), M_ext_TeV=float(Mext), sigma_cons_fb=s_cons * 1e3, sigma_opt_fb=s_opt * 1e3,
                             kinematic_limit_cons_TeV=float(y * sqrt_s), excluded_opt=bool(s_opt * 1e3 * 0.5 > ATLAS_MONOJET_LIMIT_FB)))
    return rows


def lhc_lower_edge(lumi: C.Lumi, n=6, M_D_TeV=8.0, kind="hayward", eff=0.5, trapped="optimistic") -> dict:
    """Smallest M_ext allowed by the monojet limit, with an assumed efficiency of the high-recoil region for model-U events
    (eff ~ 0.5: most events have an ISR jet and MET ~ pT(jet); conservatively bounded)."""
    def excess(M):
        return np.log(max(C.sigma_pp(lumi, 13.0, n, M_D_TeV, M, trapped) * 1e3 * eff, 1e-30)) - np.log(ATLAS_MONOJET_LIMIT_FB)
    if excess(13.0) > 0:
        return dict(M_ext_min_TeV=13.0)
    if excess(1.0) < 0:
        return dict(M_ext_min_TeV=1.0)
    m = brentq(excess, 1.0, 13.0)
    return dict(n=n, M_D_TeV=M_D_TeV, kind=kind, eff=eff, trapped=trapped, M_ext_min_TeV=float(m))


def multiobject_search_acceptance(raw: dict, ST_min_TeV=5.0, N_min=3) -> float:
    ST = np.array(raw["ST"]); N = np.array(raw["N_obj"])
    return float(np.mean((ST > ST_min_TeV * 1e3) & (N >= N_min)))


def remnant_safety(n=6, M_D_TeV=8.0, c_l=0.6, kind="hayward", n_nucleon=2.4e44, v_rel=0.3, g_star=10.0) -> dict:
    """Neutron-star test: geometric accretion (capture radius r_c = ((n+3)/2)^{1/(n+1)} sqrt((n+3)/(n+1)) r_+) of nucleons at the Fermi
    velocity against Hawking re-emission of the excess above extremality (g_star ~ 10 light dof at T << 100 GeV)."""
    D = n + 4
    MD = M_D_TeV * 1e3
    l = c_l / MD
    me = U.extremal_mu(D, 1.0, kind)
    rh = U.static_horizons(me * 1.000001, D, 1.0, kind)[1] * l          # GeV^-1
    rc_fac = ((n + 3) / 2) ** (1 / (n + 1)) * np.sqrt((n + 3) / (n + 1))
    sigma_m2 = pi * (rc_fac * rh * HBARC_GEV_M) ** 2
    acc_GeV_s = n_nucleon * sigma_m2 * v_rel * 2.998e8 * 0.94
    from .pipeline6o import a_n
    A = a_n(n, c_l) * MD                                                 # M = A mu_hat
    M_ext = A * me

    def lum(delta):
        tp = U.temperature(me * (1 + delta), D, 1.0, kind)
        if tp["r_plus"] is None or tp["T"] <= 0:
            return 1e-300
        T = tp["T"] / l
        rc = rc_fac * tp["r_plus"] * l
        return (pi**2 / 120) * g_star * pi * rc**2 * T**4 / HBAR_GEV_S       # GeV/s

    d_eq = brentq(lambda d: np.log(lum(d)) - np.log(acc_GeV_s), 1e-11, 1.0)
    T_eq = U.temperature(me * (1 + d_eq), D, 1.0, kind)["T"] / l
    return dict(n=n, M_D_TeV=M_D_TeV, c_l=c_l, kind=kind, remnant_TeV=float(M_ext / 1e3), r_horizon_m=float(rh * HBARC_GEV_M),
                accretion_GeV_per_s=float(acc_GeV_s), equilibrium_excess_fraction=float(d_eq), equilibrium_T_GeV=float(T_eq),
                verdict="accreted energy is re-radiated at a sub-GeV temperature just above extremality: no growth")


def future_reach(lumi: C.Lumi, n=6, kind="hayward") -> list:
    """Largest M_ext with >= 10 and >= 3 model-U black holes (conservative trapping) for each machine."""
    rows = []
    for name, (rs, L_ab) in C.MACHINES.items():
        for N_req in (3.0, 10.0):
            def excess(M):
                return np.log(max(C.events(C.sigma_pp(lumi, rs, n, 8.0, M, "conservative"), L_ab), 1e-30)) - np.log(N_req)
            lo, hi = 1.0, 0.6 * rs
            if excess(lo) < 0:
                m = 0.0
            elif excess(hi) > 0:
                m = hi
            else:
                m = brentq(excess, lo, hi)
            rows.append(dict(machine=name, sqrt_s_TeV=rs, L_ab=L_ab, N_required=N_req, M_ext_reach_TeV=float(m)))
    return rows


def allowed_region(lumi: C.Lumi, n=6, kind="hayward", MDs=None) -> list:
    """(M_D(DL), c_l) plane: monojet bound, LHC lower edge of M_ext, and the 85/100 TeV reach."""
    MDs = np.linspace(1.5, 15, 28) if MDs is None else MDs
    mono = {b["n"]: b["MD_DL_TeV"] for b in monojet_bounds()}.get(n, 0.0)
    rows = []
    kap = U.model_u_thresholds(n, 1.0, kind)["M_ext_over_MD"]
    for md in MDs:
        lo = lhc_lower_edge(lumi, n, md, kind)["M_ext_min_TeV"]
        c_lo = (lo / (kap * md)) ** (1 / (n + 1))
        c_hi100 = C.c_l_max_reachable(lumi, n, md, 100.0, 30.0, 10.0, kind)
        c_hi85 = C.c_l_max_reachable(lumi, n, md, 85.0, 20.0, 10.0, kind)
        rows.append(dict(M_D_DL_TeV=float(md), monojet_allowed=bool(md > mono), c_l_min_lhc=float(c_lo), c_l_max_85=float(c_hi85),
                         c_l_max_100=float(c_hi100)))
    return rows
