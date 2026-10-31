"""Round 11 (G11-17 .. G11-23): GENESIS Timing Demonstrator, simulation support for the engineering design.

Sensor model (50 um LGAD, the ATLAS HGTD / CMS ETL class; 1.3 x 1.3 mm^2 pads):
* a minimum-ionising particle deposits charge in 1 um slices with Landau (Moyal-approximated) fluctuations, mean ~73 e-h pairs/um;
* electrons drift (v_e ~ 100 um/ns) to the gain layer and are multiplied by G; the multiplied holes drift back across the bulk
  (v_h ~ 70 um/ns) and dominate the induced current (Ramo: i = q v / d);
* front end: transimpedance R_T (UCSC single-channel board + second stage, ~5 kOhm total in published LGAD studies) with a second-order
  low-pass of bandwidth BW; band-limited white noise sigma_n; sampling f_s; constant-fraction discrimination with linear interpolation.
Outputs: per-plane time resolution (jitter (+) Landau term) versus gain, noise, bandwidth and sampling; two-plane time of flight;
the three-sensor method for unbiased single-plane resolutions; coincidence rates for cosmic rays and a Sr-90 source; beam-test
time-of-flight separation.  All numbers are model estimates to be validated by measurement (the design's acceptance tests).
"""

from __future__ import annotations

import numpy as np

Q_E_FC = 1.602e-4             # fC per electron


def moyal(rng, n, loc=0.0, scale=1.0):
    """Moyal approximation to the Landau distribution (heavy right tail)."""
    u = rng.uniform(1e-12, 1.0, n)
    # inverse of the Moyal CDF via the chi-square(1) representation: x = loc - scale ln(z), z ~ chi2_1
    z = rng.chisquare(1, n)
    return loc - scale * np.log(np.maximum(z, 1e-30))


LANDAU_SCALE = 3.0          # ionisation-fluctuation scale, calibrated so that the noiseless (Landau) term is ~29 ps for 50 um, the
                            # level reported for 50 um LGADs; the mean charge is renormalised so the scale does not change the mean signal
_KNORM = {}


def _knorm(ls):
    if ls not in _KNORM:
        m = moyal(np.random.default_rng(12345), 200000)
        _KNORM[ls] = float(np.mean(np.maximum(1.0 + ls * (m - 1.27) / 2.0, 0.0)))
    return _KNORM[ls]


def waveform(rng, gain=20.0, d_um=50.0, v_e=100.0, v_h=70.0, R_T=5000.0, bw_GHz=1.5, noise_mV=1.5, fs_GSps=10.0, t_window_ns=4.0,
             mean_eh_per_um=73.0, t0=1.0, dt_fine=0.002, landau_scale=None):
    """One MIP signal (mV versus ns) sampled at fs.  Returns sample times, samples and the true arrival time t0."""
    nz = int(d_um)
    z = np.arange(nz) + 0.5                                         # depth (um) from the gain layer
    ls = LANDAU_SCALE if landau_scale is None else landau_scale
    q = np.maximum(mean_eh_per_um * (1.0 + ls * (moyal(rng, nz) - 1.27) / 2.0), 0.0) / _knorm(ls)   # heavy delta-ray tail
    t = np.arange(0.0, t_window_ns, dt_fine)
    i = np.zeros_like(t)                                            # induced current in electrons/ns
    # primary electrons (small, unamplified) and amplified holes
    for zk, qk in zip(z, q):
        te = zk / v_e                                               # electron reaches the gain layer
        i += qk * v_e / d_um * ((t >= t0) & (t < t0 + te))
        th = d_um / v_h
        i += gain * qk * v_h / d_um * ((t >= t0 + te) & (t < t0 + te + th))
    I_uA = i * 1.602e-19 * 1e9 * 1e6                                # electrons/ns -> uA
    v = I_uA * 1e-6 * R_T * 1e3                                     # mV
    # second-order low-pass (two identical poles) at bw
    tau = 1.0 / (2 * np.pi * bw_GHz)
    a = dt_fine / (tau + dt_fine)
    y1 = np.zeros_like(v); y2 = np.zeros_like(v)
    for k in range(1, len(v)):
        y1[k] = y1[k - 1] + a * (v[k] - y1[k - 1])
        y2[k] = y2[k - 1] + a * (y1[k] - y2[k - 1])
    ts = np.arange(0.0, t_window_ns, 1.0 / fs_GSps)
    s = np.interp(ts, t, y2)
    # band-limited noise: white noise filtered by the same response, normalised to sigma_n
    wn = rng.normal(0, 1, len(t))
    n1 = np.zeros_like(wn); n2 = np.zeros_like(wn)
    for k in range(1, len(wn)):
        n1[k] = n1[k - 1] + a * (wn[k] - n1[k - 1])
        n2[k] = n2[k - 1] + a * (n1[k] - n2[k - 1])
    n2 *= noise_mV / max(np.std(n2[len(n2) // 4:]), 1e-12)
    s = s + np.interp(ts, t, n2)
    return ts, s, t0


def cfd_time(ts, s, frac=0.5):
    k = int(np.argmax(s)); amp = s[k]
    thr = frac * amp
    j = np.where(s[:k + 1] < thr)[0]
    if len(j) == 0:
        return np.nan, amp
    j = j[-1]
    if j + 1 >= len(s):
        return np.nan, amp
    return ts[j] + (thr - s[j]) * (ts[j + 1] - ts[j]) / (s[j + 1] - s[j]), amp


def plane_resolution(n=1500, seed=1, frac=0.5, **kw):
    rng = np.random.default_rng(seed)
    dts, amps = [], []
    for _ in range(n):
        ts, s, t0 = waveform(rng, **kw)
        tc, a = cfd_time(ts, s, frac)
        if np.isfinite(tc):
            dts.append(tc - t0); amps.append(a)
    dts = np.array(dts); amps = np.array(amps)
    # robust width (central 90 %), in ps
    lo, hi = np.quantile(dts, [0.05, 0.95])
    core = dts[(dts > lo) & (dts < hi)]
    return dict(sigma_ps=float(np.std(core) * 1e3), mean_amp_mV=float(np.mean(amps)), median_amp_mV=float(np.median(amps)), n=len(dts),
                params={k: v for k, v in kw.items()})


def jitter_estimate(noise_mV, amp_mV, rise_ns):
    """Analytic jitter sigma_t = sigma_n / (dV/dt) with dV/dt ~ 0.8 amp / rise (10-90 %)."""
    return 1e3 * noise_mV / (0.8 * amp_mV / rise_ns)


def three_plane(sig_a, sig_b, sig_c):
    """Three-sensor method: from the measured pairwise widths s_ab, s_bc, s_ac recover each plane's resolution."""
    s_ab2, s_bc2, s_ac2 = sig_a**2 + sig_b**2, sig_b**2 + sig_c**2, sig_a**2 + sig_c**2
    a = np.sqrt(max((s_ab2 + s_ac2 - s_bc2) / 2, 0)); b = np.sqrt(max((s_ab2 + s_bc2 - s_ac2) / 2, 0)); c = np.sqrt(max((s_ac2 + s_bc2 - s_ab2) / 2, 0))
    return dict(s_ab=float(np.sqrt(s_ab2)), s_bc=float(np.sqrt(s_bc2)), s_ac=float(np.sqrt(s_ac2)), recovered=(float(a), float(b), float(c)))


def cosmic_coincidence_rate(pad_mm=1.3, baseline_m=1.0, flux_per_cm2_min=1.0):
    """Vertical-muon coincidence rate through two equal pads a baseline apart (small-angle: R = I0 A Omega, I0 ~ flux/(pi) per sr
    for a cos^2 distribution integrated ~ 1 /cm^2/min).  Returns per day."""
    A = (pad_mm * 0.1) ** 2                                         # cm^2
    omega = A / (baseline_m * 100.0) ** 2                           # sr
    I0 = flux_per_cm2_min / (np.pi * 2 / 3)                         # /cm^2/min/sr for a cos^2 law (integral of cos^2 over hemisphere = 2 pi/3)
    return float(I0 * A * omega * 60 * 24)


def sr90_rate(activity_kBq=37.0, distance_mm=10.0, pad_mm=1.3, f_mip=0.3):
    """Rate of beta electrons through a pad at a distance (isotropic source; fraction f_mip energetic enough to traverse a stack)."""
    A = (pad_mm) ** 2
    frac = A / (4 * np.pi * distance_mm**2)
    return float(activity_kBq * 1e3 * frac * f_mip)                 # per second


def beam_tof(p_list=(0.5, 1.0, 2.0, 3.0), L=1.0, sigma_tof_ps=50.0):
    c = 0.299792458
    out = []
    for p in p_list:
        masses = {"pi": 0.13957, "K": 0.49368, "p": 0.93827}
        t = {k: L / c * np.hypot(p, m) / p for k, m in masses.items()}
        out.append(dict(p_GeV=p, dt_K_pi_ps=float((t["K"] - t["pi"]) * 1e3), dt_p_pi_ps=float((t["p"] - t["pi"]) * 1e3),
                        sep_K_pi_sigma=float((t["K"] - t["pi"]) * 1e3 / sigma_tof_ps), sep_p_pi_sigma=float((t["p"] - t["pi"]) * 1e3 / sigma_tof_ps)))
    return out


def design_scan(n=400):
    """Per-plane resolution versus gain, noise, bandwidth and sampling (calibrated model); the basis of the front-end requirements."""
    rows = []
    base = dict(gain=20.0, noise_mV=1.5, bw_GHz=1.5, fs_GSps=10.0)
    rows.append(dict(case="Landau floor (noiseless, 20 GS/s)", **plane_resolution(n, gain=20.0, noise_mV=0.01, fs_GSps=20.0)))
    for g in (8.0, 15.0, 20.0, 30.0):
        rows.append(dict(case="gain %g" % g, **plane_resolution(n, **{**base, "gain": g})))
    for nm in (0.8, 3.0, 5.0):
        rows.append(dict(case="noise %g mV" % nm, **plane_resolution(n, **{**base, "noise_mV": nm})))
    for bw in (0.5, 1.0, 2.5):
        rows.append(dict(case="bandwidth %g GHz" % bw, **plane_resolution(n, **{**base, "bw_GHz": bw})))
    for fs in (2.0, 5.0, 20.0):
        rows.append(dict(case="sampling %g GS/s" % fs, **plane_resolution(n, **{**base, "fs_GSps": fs})))
    return rows


# --------------------------------------------------------------------------------------------- DAQ / analysis chain (software skeleton)

def simulate_bench_run(n_events=3000, planes=3, seed=7, noise_mV=(1.5, 1.5, 1.5), gain=(20.0, 15.0, 25.0), offsets_ns=(0.0, 0.137, -0.081),
                       trigger_rate_hz=15.0, clock_jitter_ps=2.0):
    """A Sr-90 stack run: per event, each plane records a waveform; timestamps come from a common digitiser clock (DRS4-like), with
    fixed cable/board offsets and a small clock jitter.  Returns a list of events with per-plane (t_cfd, amplitude) and the trigger time."""
    rng = np.random.default_rng(seed)
    t_trig = np.cumsum(rng.exponential(1.0 / trigger_rate_hz, n_events))
    events = []
    for k in range(n_events):
        rec = []
        for p in range(planes):
            ts, s, t0 = waveform(rng, gain=gain[p], noise_mV=noise_mV[p], fs_GSps=5.0)
            tc, a = cfd_time(ts, s, 0.5)
            rec.append((tc - t0 + offsets_ns[p] + rng.normal(0, clock_jitter_ps * 1e-3), a))
        events.append(dict(t_trigger_s=float(t_trig[k]), planes=rec))
    return events


def analyze_run(events, amp_cut_mV=10.0, walk_correct=True):
    """Event selection (all planes above an amplitude cut), offset calibration, optional amplitude-walk correction, pairwise time
    differences and the three-sensor solution."""
    sel = [e for e in events if all(np.isfinite(t) and a > amp_cut_mV for t, a in e["planes"])]
    T = np.array([[t for t, a in e["planes"]] for e in sel]); A = np.array([[a for t, a in e["planes"]] for e in sel])
    T = T - np.median(T, axis=0)                                     # offset calibration
    if walk_correct:                                                 # residual dependence on amplitude (linear in 1/A), fitted per plane
        for p in range(T.shape[1]):
            ref = np.mean(np.delete(T, p, axis=1), axis=1)
            c = np.polyfit(1.0 / A[:, p], T[:, p] - ref, 1)
            T[:, p] -= np.polyval(c, 1.0 / A[:, p]) - np.polyval(c, np.median(1.0 / A[:, p]))
    def width(x):
        lo, hi = np.quantile(x, [0.05, 0.95]); core = x[(x > lo) & (x < hi)]; return float(np.std(core) * 1e3)
    s_ab, s_bc, s_ac = width(T[:, 0] - T[:, 1]), width(T[:, 1] - T[:, 2]), width(T[:, 0] - T[:, 2])
    a = np.sqrt(max((s_ab**2 + s_ac**2 - s_bc**2) / 2, 0)); b = np.sqrt(max((s_ab**2 + s_bc**2 - s_ac**2) / 2, 0)); c = np.sqrt(max((s_ac**2 + s_bc**2 - s_ab**2) / 2, 0))
    # statistical uncertainty of a width estimate ~ sigma / sqrt(2 N)
    N = len(sel)
    return dict(n_selected=N, pairwise_ps=(s_ab, s_bc, s_ac), plane_ps=(float(a), float(b), float(c)),
                stat_unc_ps=float(np.mean([s_ab, s_bc, s_ac]) / np.sqrt(2 * N)))


def tof_beta_resolution(sigma_plane_ps=35.0, L=1.0, betas=(0.3, 0.5, 0.7, 0.9, 0.99)):
    c = 0.299792458
    s_tof = np.sqrt(2) * sigma_plane_ps * 1e-3
    return [dict(beta=b, tof_ns=float(L / (c * b)), sigma_beta=float(b**2 * c * s_tof / L)) for b in betas]
