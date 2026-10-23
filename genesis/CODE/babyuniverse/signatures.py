"""Problem C (Stage 6 tasks 40-46): what an exterior observer can and cannot see.

* State level (Theorem V-2, Birkhoff): after the event the exterior is Schwarzschild of mass M on both branches.
* Event level (Proposition C1, this module): the two branches differ in the history of the region 2GM < r < R1.
  On the time-symmetric slice at R1 and on the daughter slice at R2, the geometry for r > R1 is the same Schwarzschild.
  On the daughter branch the wall is in region III and no light from it reaches region I.  Radial photons emitted by the
  wall at a constant proper rate Gamma and rest-frame frequency nu_e therefore arrive at infinity as follows:
    collapse branch:  rate Gamma dtau/du and frequency nu_e dtau/du with du/dtau = (A_out - Rdot)/f_out, an exponential tail
                      ~ exp(-u/4GM) in rate; the total number arriving after the turning point is Gamma * tau_collapse;
    daughter branch:  emission stops at the turning point; the last photons arrive with finite redshift, nu = nu_e sqrt(f(R1)).
  A detector that counts late photons distinguishes the branches; it cannot see the daughter itself.
* Dynamical assistance (Stage 4, problem B): the Euclidean duration of the barrier crossing sets the frequency a driving
  field needs, omega_c ~ 1 / Delta tau_E.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import quad

from .model import ThinWallModel
from . import tunnelling as T


def _rstar(model, R):
    rs = 2 * model.G * model.M
    return R + rs * np.log(R / rs - 1)


def light_curve(model: ThinWallModel, n: int = 400, R_start_frac: float = 0.6, eps: float = 1e-9) -> dict:
    """Retarded arrival time u (= t - r*, up to a constant) and the observed photon rate / frequency factor dtau/du for a
    wall emitting radially outward.  Expanding leg from R_start = 2GM + R_start_frac (R1 - 2GM) to R1 (common to both
    branches), then the collapsing leg from R1 towards 2GM (collapse branch only)."""
    if model.eta != 0:
        raise ValueError("Schwarzschild exterior required")
    R1 = model.turning_points()[0]
    rs = 2 * model.G * model.M
    Ra = rs + R_start_frac * (R1 - rs)

    def leg(R_end, sgn):
        """Integrate from the turning point R1 (tau = t = 0) down to R_end.  sgn = +1: expanding leg, traversed backwards
        in time (tau, t < 0); sgn = -1: collapsing leg (tau, t > 0)."""
        s = np.linspace(0, 1, n)
        R = R1 - (R1 - R_end) * s**2
        us, rates, taus = [], [], []
        tau = t = 0.0
        for i, Ri in enumerate(R):
            if i > 0:
                a, b = Ri, R[i - 1]
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    tau += quad(lambda r: 1 / np.sqrt(max(-float(model.V(r)), 1e-300)), a, b, limit=200)[0]
                    t += quad(lambda r: float(model.A_out(r) / model.f_out(r)) / np.sqrt(max(-float(model.V(r)), 1e-300)),
                              a, b, limit=200)[0]
            rd = sgn * np.sqrt(max(-float(model.V(Ri)), 0.0))
            du_dtau = (float(model.A_out(Ri)) - rd) / float(model.f_out(Ri))
            tt = -t if sgn > 0 else t
            us.append(tt - _rstar(model, Ri))
            rates.append(1 / du_dtau)
            taus.append(-tau if sgn > 0 else tau)
        return np.array(us), np.array(rates), np.array(taus)

    ue, re, te = leg(Ra, +1)
    uc, rc, tc = leg(rs * (1 + eps), -1)
    order = np.argsort(te)
    ue, re, te = ue[order], re[order], te[order]
    return dict(R1=float(R1), rs=float(rs), u_expanding=ue, rate_expanding=re, tau_expanding=te,
                u_collapse=uc, rate_collapse=rc, tau_collapse=tc,
                last_redshift_factor=float(np.sqrt(model.f_out(R1))))


def late_photon_fraction(model: ThinWallModel) -> dict:
    """Photons (per unit proper emission rate) arriving after the turning-point signal: collapse branch = proper time of the
    collapsing leg outside the horizon; daughter branch = 0.  Also the e-folding of the collapse tail in u (-> 4GM)."""
    R1 = model.turning_points()[0]
    rs = 2 * model.G * model.M
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tau_c = quad(lambda r: 1 / np.sqrt(max(-float(model.V(r)), 1e-300)), rs, R1, limit=400)[0]
    lc = light_curve(model)
    u, r = lc["u_collapse"], lc["rate_collapse"]
    k = slice(int(0.85 * len(u)), len(u) - 1)
    slope = np.polyfit(u[k], np.log(r[k]), 1)[0]
    return dict(N_late_collapse=float(tau_c), N_late_daughter=0.0, tail_efold=float(-1 / slope), four_GM=float(2 * rs),
                last_redshift_factor=lc["last_redshift_factor"])


def euclidean_crossing_time(model: ThinWallModel) -> float:
    """Proper Euclidean duration of the barrier crossing, Delta tau_E = int dR / sqrt(V) over (R1, R2)."""
    R1, R2 = T.barrier(model)

    def g(s):
        R = R1 + (R2 - R1) * (1 - np.cos(s)) / 2
        dR = (R2 - R1) * np.sin(s) / 2
        return dR / np.sqrt(max(float(model.V(R)), 1e-300))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(quad(g, 0, np.pi, limit=400)[0])


def assistance_frequency(E_vac_GeV: float, E_wall_GeV: float, M_planck: float) -> dict:
    """Frequency omega_c = 1/Delta tau_E a weak driving field must exceed to assist the laboratory transition (in GeV),
    with the thin-barrier estimate Delta tau_E ~ (pi/3) R_CdL, R_CdL = 3 sigma / rho."""
    from .observability import planck_params, M_PL_GEV
    p = planck_params(E_vac_GeV, E_wall_GeV)
    chi, k = p["chi"], p["k"]
    m = ThinWallModel(1.0, k, M_planck * chi)
    dtau = euclidean_crossing_time(m) / chi                      # Planck units
    R0 = 3 * p["sigma"] / p["rho"]
    return dict(delta_tau_E_lP=float(dtau), thin_estimate_lP=float(np.pi / 3 * R0), omega_c_GeV=float(M_PL_GEV / dtau),
                omega_c_over_Evac=float(M_PL_GEV / dtau / E_vac_GeV))
