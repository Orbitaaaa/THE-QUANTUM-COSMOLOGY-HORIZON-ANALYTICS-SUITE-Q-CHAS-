"""Real gravitational-wave strain: LOSC HDF5 loading, whitening, ringdown fits, GWOSC catalogue.

Pipeline for the four events of the LIGO open-data tutorial (32 s of 4096 Hz strain per
detector, GWOSC release): Welch PSD (4 s Hann segments, 50 % overlap, mean-averaged as in
the LOSC tutorial) interpolated onto the FFT grid, frequency-domain whitening of the
Tukey-windowed strain, zero-phase Butterworth band-pass to the event's analysis band, and
the H1-L1 lag from the cross-correlation of the last loud cycles within the 10 ms light
travel time (sign flip allowed: the two detectors are nearly anti-aligned).

The remnant of each event follows from the template parameters (detector-frame masses,
aligned spins) through :mod:`qchas.horizon.kerr`. A single damped sinusoid is fitted to the
whitened H1 strain, band-passed 20-1000 Hz as in [QG] Sec. 19.2, and compared with the Kerr
(2,2,0) prediction. Its uncertainty combines the spread over start times 1-6 ms after the
peak with refits of the best fit injected into off-source noise; fits to pure off-source
noise give the false-alarm probability of the mode, so weak ringdowns are flagged rather
than quoted.

The [QG] Sec. 19 echo search and ringdown tables are quoted, never recomputed here.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
from scipy import signal
from scipy.optimize import least_squares

from ..horizon import kerr
from ..horizon.page import page_time_fraction
from .fetch import DEFAULT_RAW_DIR

TUTORIAL_EVENTS_FILE = "BBH_events_v3.json"
GWTC_FILE = "gwtc_allevents.json"
PSD_SEGMENT_S = 4.0
CROP_S = 0.2
MAX_LAG_S = 0.010
#: stretch (s before, s after tevent) used for the H1-L1 cross-correlation: the loud last cycles
LAG_WINDOW_S = (0.10, 0.02)
EXPORT_FS = 2048
#: ringdown band of [QG] Sec. 19.2 (extended to 2 f220 for light remnants)
RINGDOWN_BAND = (20.0, 1000.0)

#: [QG] Sec. 19.1 echo search in public LIGO-Virgo strain (R90 = None means "> 1": unconstrained)
_ECHO_FIELDS = ("event", "snr", "mf", "af", "f220", "dt_ms", "gamma", "echo_snr", "p", "R90")
PAPER_ECHO_TABLE = [dict(zip(_ECHO_FIELDS, row)) for row in [
    ("GW150914", 23.8, 68.7, 0.68, 248, 293, 0.50, 3.07, 0.83, 0.67),
    ("GW151012", 9.7, 47.4, 0.67, 355, 199, 0.50, 3.08, 0.73, None),
    ("GW151226", 12.7, 21.4, 0.79, 875, 100, 0.46, 3.10, 0.51, None),
    ("GW170104", 13.3, 58.3, 0.64, 282, 241, 0.51, 3.33, 0.51, None),
    ("GW170814", 16.9, 59.4, 0.72, 296, 262, 0.49, 3.17, 0.70, 0.95),
    ("GW190521", 24.3, 85.6, 0.68, 199, 366, 0.50, 2.83, 0.97, 0.63),
    ("GW200129", 25.4, 71.4, 0.73, 248, 318, 0.49, 3.64, 0.30, 0.74),
    ("GW241011", 34.2, 21.0, 0.90, 1020, 122, 0.35, 3.41, 0.30, None),
    ("GW250114", 77.8, 68.6, 0.66, 244, 289, 0.50, 3.83, 0.18, 0.28),
]]
PAPER_ECHO_SUMMARY = {
    "combined_R90": 0.28, "sum_snr2": 97.3, "p_combined": 0.78, "p_fisher": 0.81, "A90_gamma09": 0.08,
    "dt_gw150914_ms": {"ours": 293.1, "abedi": 292.5},
    "masses": "detector frame, from the paper's own IMRPhenomD fits",
    "source": "[QG] Sec. 19.1",
}
#: [QG] Sec. 19.2 ringdown table (detector-frame masses)
PAPER_RINGDOWN_TABLE = [
    {"event": "GW150914", "imr_mf": [67.4, 1.9], "imr_af": [0.64, 0.05], "kerr_f220": [246, 4], "kerr_tau_ms": [3.99, 0.18],
     "free_f": [249, 14], "free_tau_ms": [3.77, 1.50], "p_free": 0.96, "p_free_range": [0.23, 0.98],
     "kerr_221_mf": [59.1, 5.9], "kerr_221_af": [0.42, 0.18], "p_kerr": 0.40, "p_kerr_range": [0.30, 0.55]},
    {"event": "GW250114", "imr_mf": [68.5, 0.3], "imr_af": [0.66, 0.01], "kerr_f220": [245, 1], "kerr_tau_ms": [4.08, 0.03],
     "free_f": [243, 4], "free_tau_ms": [4.89, 0.44], "p_free": 0.09, "p_free_range": [0.00, 0.72],
     "kerr_221_mf": [72.2, 2.9], "kerr_221_af": [0.68, 0.06], "p_kerr": 0.01, "p_kerr_range": [0.00, 0.89]},
]
#: [QG] Sec. 18 prediction 2: |R|^2 <= e^{-S}, 10^{-1.75e80} for the 62 Msun GW150914 remnant
PAPER_ECHO_BOUND = {"m_sun": 62.0, "log10_R2": -1.75e80, "source": "[QG] Sec. 18 prediction 2"}


# --- input -----------------------------------------------------------------------------

def read_losc_hdf5(path: str | Path) -> tuple[np.ndarray, float, float]:
    """Strain, GPS start time and sampling rate of a LOSC/GWOSC HDF5 file."""
    import h5py

    with h5py.File(path, "r") as f:
        dset = f["strain/Strain"]
        strain = np.asarray(dset[()], dtype=float)
        gps_start = float(f["meta/GPSstart"][()])
        fs = 1.0 / float(dset.attrs["Xspacing"])
    return strain, gps_start, fs


def load_tutorial_events(raw_dir: str | Path | None = None) -> dict[str, dict]:
    """The LOSC tutorial event list (template masses and aligned spins are detector frame)."""
    path = Path(raw_dir or DEFAULT_RAW_DIR) / TUTORIAL_EVENTS_FILE
    with open(path) as f:
        return json.load(f)


# --- conditioning ------------------------------------------------------------------------

def psd_welch(strain: np.ndarray, fs: float, seg_s: float = PSD_SEGMENT_S) -> tuple[np.ndarray, np.ndarray]:
    """One-sided Welch PSD with Hann segments of ``seg_s`` seconds and 50 % overlap."""
    nper = int(seg_s * fs)
    return signal.welch(strain, fs=fs, window="hann", nperseg=nper, noverlap=nper // 2)


def whiten(strain: np.ndarray, fs: float, psd: tuple[np.ndarray, np.ndarray] | None = None) -> np.ndarray:
    """Divide the Tukey-windowed strain by the amplitude spectral density (unit-variance noise)."""
    freqs_psd, pxx = psd if psd is not None else psd_welch(strain, fs)
    n = strain.size
    dt = 1.0 / fs
    freqs = np.fft.rfftfreq(n, dt)
    asd = np.sqrt(np.interp(freqs, freqs_psd, pxx))
    hf = np.fft.rfft(strain * signal.windows.tukey(n, alpha=1.0 / 8.0))
    white = hf / asd * math.sqrt(2.0 * dt)
    return np.fft.irfft(white, n=n)


def bandpass(x: np.ndarray, fs: float, band: tuple[float, float], order: int = 4) -> np.ndarray:
    """Zero-phase Butterworth band-pass (second-order sections, forward-backward)."""
    sos = signal.butter(order, [band[0], band[1]], btype="bandpass", fs=fs, output="sos")
    return signal.sosfiltfilt(sos, x)


def detector_lag(h1: np.ndarray, l1: np.ndarray, fs: float, max_lag_s: float = MAX_LAG_S) -> dict:
    """Lag of L1 relative to H1 maximising |normalised cross-correlation| within +-max_lag_s.

    A positive lag means the signal reaches L1 later than H1; ``sign`` is -1 when L1 is
    inverted relative to H1.
    """
    kmax = int(round(max_lag_s * fs))
    a = (h1 - h1.mean()) / (h1.std() * h1.size)
    b = (l1 - l1.mean()) / l1.std()
    lags = np.arange(-kmax, kmax + 1)
    cc = np.array([np.dot(a[max(0, -k): a.size - max(0, k)], b[max(0, k): b.size - max(0, -k)]) for k in lags])
    i = int(np.argmax(np.abs(cc)))
    return {"lag_s": float(lags[i] / fs), "sign": int(np.sign(cc[i])), "corr": float(cc[i])}


# --- ringdown ----------------------------------------------------------------------------

def damped_sinusoid(t: np.ndarray, amp: float, tau: float, f: float, phi: float) -> np.ndarray:
    return amp * np.exp(-t / tau) * np.cos(2.0 * math.pi * f * t + phi)


def _fit_damped(t: np.ndarray, y: np.ndarray, f0: float, tau0: float) -> dict:
    """Least-squares fit of A e^{-t/tau} cos(2 pi f t + phi); A and phi enter linearly."""

    def basis(f: float, tau: float) -> np.ndarray:
        env = np.exp(-t / tau)
        return np.stack([env * np.cos(2 * math.pi * f * t), env * np.sin(2 * math.pi * f * t)], axis=1)

    def resid(p: np.ndarray) -> np.ndarray:
        b = basis(p[0], math.exp(p[1]))
        c, *_ = np.linalg.lstsq(b, y, rcond=None)
        return b @ c - y

    # global search on a (f, tau) grid with the linear amplitudes solved in closed form
    fg = np.linspace(max(20.0, 0.5 * f0), min(2000.0, 2.0 * f0), 151)
    tg = np.geomspace(max(2e-4, 0.25 * tau0), min(0.1, 8.0 * tau0), 41)
    env = np.exp(-t[None, :] / tg[:, None])[:, None, :]
    ph = 2 * math.pi * fg[:, None] * t[None, :]
    b1, b2 = env * np.cos(ph)[None], env * np.sin(ph)[None]
    g11, g12, g22 = (b1 * b1).sum(-1), (b1 * b2).sum(-1), (b2 * b2).sum(-1)
    r1, r2 = (b1 * y).sum(-1), (b2 * y).sum(-1)
    det = g11 * g22 - g12**2
    explained = (g22 * r1**2 - 2 * g12 * r1 * r2 + g11 * r2**2) / np.where(det > 0, det, np.inf)
    i, j = np.unravel_index(int(np.argmax(explained)), explained.shape)
    best = least_squares(resid, [fg[j], math.log(tg[i])], bounds=([20.0, math.log(2e-4)], [2000.0, math.log(0.1)]))
    f, tau = float(best.x[0]), float(math.exp(best.x[1]))
    b = basis(f, tau)
    c, *_ = np.linalg.lstsq(b, y, rcond=None)
    return {"f": f, "tau": tau, "amp": float(math.hypot(*c)), "phi": float(math.atan2(-c[1], c[0])),
            "rms_resid": float(np.sqrt(2.0 * best.cost / y.size))}


def ringdown_fit(white_bp: np.ndarray, fs: float, t_peak_idx: int, f_pred: float, tau_pred: float,
                 noise: np.ndarray | None = None, starts_ms=None, window_s: float | None = None,
                 n_noise: int = 40, seed: int = 0) -> dict:
    """Single-mode ringdown fit to whitened, band-passed strain after the peak.

    Fits are repeated for every start time 1-6 ms after the peak (the start-time spread is
    the systematic). With off-source ``noise`` of the same detector and filtering, the fit
    is (i) refitted after injecting the median-start fit into ``n_noise`` noise stretches
    (the statistical spread), and (ii) run on ``n_noise`` pure-noise stretches, whose mode
    SNRs give the false-alarm probability ``p_noise`` of the on-source mode SNR. Quoted f and
    tau are medians over start times with both spreads (16-84 % half widths) in quadrature.
    """
    starts_ms = np.arange(1.0, 6.0 + 1e-9, 1000.0 / fs) if starts_ms is None else np.asarray(starts_ms)
    window_s = window_s or max(0.02, 6.0 * tau_pred)
    nwin = int(round(window_s * fs))
    t = np.arange(nwin) / fs
    fits = []
    for s in starts_ms:
        i0 = t_peak_idx + int(round(s * 1e-3 * fs))
        fits.append(_fit_damped(t, white_bp[i0:i0 + nwin], f_pred, tau_pred))
    f_all = np.array([r["f"] for r in fits])
    tau_all = np.array([r["tau"] for r in fits])
    f_med, tau_med = float(np.median(f_all)), float(np.median(tau_all))
    k = int(np.argmin(np.abs(f_all - f_med) / f_med + np.abs(tau_all - tau_med) / tau_med))
    ref = fits[k]
    model = damped_sinusoid(t, ref["amp"], ref["tau"], ref["f"], ref["phi"])

    def half_width(x) -> float:
        q16, q84 = np.percentile(np.asarray(x), [16, 84])
        return float(0.5 * (q84 - q16))

    out_noise: dict = {"f_stat_Hz": None, "tau_stat_ms": None, "snr": None, "p_noise": None, "sigma": None}
    if noise is not None and n_noise > 0:
        sigma = float(np.std(noise))
        snr = float(np.sqrt(np.sum(model**2)) / sigma)
        rng = np.random.default_rng(seed)
        fn, tn, null_snr = [], [], []
        for _ in range(n_noise):
            j = int(rng.integers(0, noise.size - nwin))
            seg = noise[j:j + nwin]
            r = _fit_damped(t, model + seg, f_pred, tau_pred)
            fn.append(r["f"])
            tn.append(r["tau"])
            j2 = int(rng.integers(0, noise.size - nwin))
            r0 = _fit_damped(t, noise[j2:j2 + nwin], f_pred, tau_pred)
            null_snr.append(np.sqrt(np.sum(damped_sinusoid(t, r0["amp"], r0["tau"], r0["f"], r0["phi"]) ** 2)) / sigma)
        p_noise = (1 + int(np.sum(np.array(null_snr) >= snr))) / (n_noise + 1)
        out_noise = {"f_stat_Hz": half_width(fn), "tau_stat_ms": 1e3 * half_width(tn), "snr": snr,
                     "p_noise": p_noise, "sigma": sigma}
    sys_f, sys_tau = half_width(f_all), 1e3 * half_width(tau_all)
    err_f = math.hypot(sys_f, out_noise["f_stat_Hz"] or 0.0)
    err_tau = math.hypot(sys_tau, out_noise["tau_stat_ms"] or 0.0)
    return {
        "f_Hz": f_med, "f_err_Hz": err_f, "tau_ms": 1e3 * tau_med, "tau_err_ms": err_tau,
        "f_sys_Hz": sys_f, "tau_sys_ms": sys_tau, **out_noise,
        "significant": out_noise["p_noise"] is not None and out_noise["p_noise"] < 0.05,
        "t0_ms_after_peak": float(np.median(starts_ms)), "t0_range_ms": [float(starts_ms[0]), float(starts_ms[-1])],
        "window_ms": 1e3 * window_s,
        "per_start": {"t0_ms": [float(s) for s in starts_ms], "f_Hz": f_all.tolist(),
                      "tau_ms": (1e3 * tau_all).tolist()},
        "n_noise": n_noise if noise is not None else 0, "computed": True,
    }


# --- Kerr and horizon numbers of a remnant ------------------------------------------------

def remnant_summary(mf: float, af: float, mf_thermo: float | None = None) -> dict:
    """Ringdown, echo delay, entropy, echo bound, Hawking temperature and lifetimes of a remnant.

    ``mf`` (detector frame) sets the observed f220, tau220 and echo delay. The thermodynamic
    quantities use ``mf_thermo`` (source frame) when given, else ``mf``. Lifetimes assume
    dM/dt = -alpha/M^2 with the massless-scalar alpha = 7.44e-5 ([BH] Sec. 17) on a
    Schwarzschild hole of that mass: an order-of-magnitude scale, not a full Kerr spin-down.
    """
    f220, tau220 = kerr.qnm(mf, af)
    m_t = mf if mf_thermo is None else mf_thermo
    bh = kerr.KerrHorizon(m_t, af)
    t_evap = kerr.KerrHorizon(m_t, 0.0).evaporation_time_years()
    return {
        "f220_Hz": f220, "tau220_ms": 1e3 * tau220, "echo_delay_ms": 1e3 * kerr.echo_delay(mf, af),
        "S_BH_nats": bh.entropy, "S_BH_log10": bh.log10_entropy,
        "S_BH_bits_log10": math.log10(bh.entropy / math.log(2)),
        "echo_bound_log10": kerr.echo_reflectivity_bound(bh.entropy),
        "T_H_K": bh.hawking_temperature_K, "r_plus_km": bh.r_plus_km,
        "t_evap_log10_yr": math.log10(t_evap), "page_time_log10_yr": math.log10(page_time_fraction() * t_evap),
        "lifetime_model": "massless scalar, alpha = 7.44e-5 / M^2 ([BH] Sec. 17), Schwarzschild",
    }


# --- one event end to end -------------------------------------------------------------------

def _decimate(x: np.ndarray, factor: int) -> np.ndarray:
    return x[::factor] if factor > 1 else x


def analyse_event(name: str, params: dict, raw_dir: str | Path | None = None, n_noise: int = 40,
                  redshift: float | None = None) -> dict:
    """Whitened strain, detector lag, remnant and ringdown fit for one tutorial event."""
    raw = Path(raw_dir or DEFAULT_RAW_DIR)
    m1, m2, a1, a2 = (float(params[k]) for k in ("m1", "m2", "a1", "a2"))
    mf, af = kerr.remnant_aligned(m1, m2, a1, a2)
    mf_src = None if redshift is None else mf / (1.0 + redshift)
    out = {"name": name, "source": "LOSC tutorial (GWOSC)", "tevent": float(params["tevent"]),
           "m1": m1, "m2": m2, "a1": a1, "a2": a2, "mass_frame": "detector", "approximant": params.get("approx"),
           "mf": mf, "af": af, "mf_source": mf_src, "redshift": redshift,
           "thermo_frame": "detector" if mf_src is None else "source",
           **remnant_summary(mf, af, mf_src), "strain": None, "ringdown_fit": None, "warnings": []}
    paths = {det: raw / params[f"fn_{det}"] for det in ("H1", "L1")}
    missing = [str(p.name) for p in paths.values() if not p.exists()]
    if missing:
        out["warnings"].append(f"{name}: strain files not available ({', '.join(missing)})")
        return out
    band = tuple(float(b) for b in params["fband"])
    data = {}
    for det, p in paths.items():
        strain, gps0, fs = read_losc_hdf5(p)
        white = whiten(strain, fs)
        data[det] = {"white": white, "bp": bandpass(white, fs, band), "gps0": gps0, "fs": fs}
    fs = data["H1"]["fs"]
    gps0 = data["H1"]["gps0"]
    t = gps0 + np.arange(data["H1"]["bp"].size) / fs
    i_ev = int(round((out["tevent"] - gps0) * fs))
    ncrop = int(round(CROP_S * fs))
    sl = slice(i_ev - ncrop, i_ev + ncrop + 1)
    h1, l1 = data["H1"]["bp"], data["L1"]["bp"]
    sl_lag = slice(i_ev - int(round(LAG_WINDOW_S[0] * fs)), i_ev + int(round(LAG_WINDOW_S[1] * fs)))
    lag = detector_lag(h1[sl_lag], l1[sl_lag], fs)
    factor = max(1, int(round(fs / EXPORT_FS)))
    out["strain"] = {
        "t": np.round(_decimate(t[sl] - out["tevent"], factor), 6).tolist(),
        "H1": np.round(_decimate(h1[sl], factor), 4).tolist(),
        "L1": np.round(_decimate(l1[sl], factor), 4).tolist(),
        "fs": fs, "fs_export": fs / factor, "band": list(band),
        "units": "whitened (noise sigma = 1 before band-pass)",
        "lag_ms": 1e3 * lag["lag_s"], "lag_sign": lag["sign"], "lag_corr": lag["corr"],
        "lag_convention": "positive: L1 later than H1; sign -1: L1 inverted",
        "lag_window_s": list(LAG_WINDOW_S),
    }
    # peak of the chirp in the event band; ringdown in the 20-1000 Hz band of [QG] Sec. 19.2
    half = int(0.05 * fs)
    envelope = np.abs(signal.hilbert(h1[i_ev - 4 * half: i_ev + 4 * half]))  # padded against edge effects
    i_peak = i_ev - half + int(np.argmax(envelope[3 * half: 5 * half]))
    rd_band = (20.0, min(0.45 * fs, max(RINGDOWN_BAND[1], 2.0 * out["f220_Hz"])))
    h1_rd = bandpass(data["H1"]["white"], fs, rd_band)
    noise = np.concatenate([h1_rd[int(2 * fs): i_ev - int(fs)], h1_rd[i_ev + int(fs): -int(2 * fs)]])
    fit = ringdown_fit(h1_rd, fs, i_peak, out["f220_Hz"], out["tau220_ms"] * 1e-3, noise=noise, n_noise=n_noise)
    fit["t_peak_gps"] = float(t[i_peak])
    fit["band"] = list(rd_band)
    fit["kerr_f220_Hz"] = out["f220_Hz"]
    fit["kerr_tau220_ms"] = out["tau220_ms"]
    fit["f_pull"] = (fit["f_Hz"] - out["f220_Hz"]) / fit["f_err_Hz"] if fit["f_err_Hz"] > 0 else None
    fit["tau_pull"] = (fit["tau_ms"] - out["tau220_ms"]) / fit["tau_err_ms"] if fit["tau_err_ms"] > 0 else None
    out["ringdown_fit"] = fit
    return out


# --- GWOSC catalogue ------------------------------------------------------------------------

_GWTC_FIELDS = {
    "name": "commonName", "gps": "GPS", "m1_source": "mass_1_source", "m2_source": "mass_2_source",
    "mf_source": "final_mass_source", "chi_eff": "chi_eff", "redshift": "redshift",
    "distance_Mpc": "luminosity_distance", "snr": "network_matched_filter_snr", "catalog": "catalog.shortName",
}


def _num(x) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def parse_gwosc_allevents(data: dict | str | Path) -> list[dict]:
    """Rows of the GWOSC event API JSON (``/eventapi/json/allevents/``), one per event.

    Robust to null and missing fields. When several versions of an event are listed, the
    highest version is kept. Where the final mass and chi_eff are known, the remnant spin is
    estimated with the aligned-spin fit (both spins set to chi_eff, an approximation) and
    the Kerr numbers are added.
    """
    if not isinstance(data, dict):
        with open(data) as f:
            data = json.load(f)
    events = data.get("events", data)
    best: dict[str, tuple[int, dict]] = {}
    for key, ev in events.items():
        if not isinstance(ev, dict):
            continue
        row = {k: ev.get(v) for k, v in _GWTC_FIELDS.items()}
        if row["catalog"] is None and isinstance(ev.get("catalog"), dict):
            row["catalog"] = ev["catalog"].get("shortName")
        row["name"] = row["name"] or re.sub(r"-v\d+$", "", str(key))
        for k in ("gps", "m1_source", "m2_source", "mf_source", "chi_eff", "redshift", "distance_Mpc", "snr"):
            row[k] = _num(row[k])
        m = re.search(r"-v(\d+)$", str(key))
        version = int(ev.get("version") or (m.group(1) if m else 0) or 0)
        if row["name"] not in best or version > best[row["name"]][0]:
            best[row["name"]] = (version, row)
    rows = []
    for version, row in best.values():
        row["version"] = version
        row.update({"af_estimate": None, "f220_Hz": None, "echo_delay_ms": None, "S_BH_log10": None})
        if row["m1_source"] and row["m2_source"] and row["chi_eff"] is not None:
            mf_fit, af = kerr.remnant_aligned(row["m1_source"], row["m2_source"], row["chi_eff"], row["chi_eff"])
            mf_src = row["mf_source"] or mf_fit
            z = row["redshift"] or 0.0
            row["af_estimate"] = af
            row["f220_Hz"] = kerr.qnm(mf_src * (1 + z), af)[0]
            row["echo_delay_ms"] = 1e3 * kerr.echo_delay(mf_src * (1 + z), af)
            row["S_BH_log10"] = kerr.KerrHorizon(mf_src, af).log10_entropy
        rows.append(row)
    rows.sort(key=lambda r: (r["gps"] is None, r["gps"] or 0.0))
    return rows


# --- bundle ---------------------------------------------------------------------------------

def to_json(raw_dir: str | Path | None = None, n_noise: int = 40, events=None) -> dict:
    """gw.json of the data contract (docs/ARCHITECTURE.md)."""
    raw = Path(raw_dir or DEFAULT_RAW_DIR)
    warnings: list[str] = []
    try:
        tutorial = load_tutorial_events(raw)
    except FileNotFoundError:
        tutorial = {}
        warnings.append(f"{TUTORIAL_EVENTS_FILE} not available; no LIGO strain analysed")
    gwtc = None
    if (raw / GWTC_FILE).exists():
        try:
            gwtc = parse_gwosc_allevents(raw / GWTC_FILE)
        except (ValueError, KeyError, TypeError) as exc:
            warnings.append(f"GWOSC catalogue unreadable: {exc}")
    else:
        warnings.append("GWOSC event API not reachable from build host; gwtc is null")
    z_by_name = {r["name"]: r["redshift"] for r in gwtc or [] if r.get("redshift") is not None}
    out_events = []
    for name, params in tutorial.items():
        if events is not None and name not in events:
            continue
        z = z_by_name.get(name) or z_by_name.get(name.replace("LVT", "GW"))
        ev = analyse_event(name, params, raw, n_noise=n_noise, redshift=z)
        warnings.extend(ev.pop("warnings"))
        out_events.append(ev)
    check = kerr.echo_reflectivity_bound(kerr.KerrHorizon(PAPER_ECHO_BOUND["m_sun"], 0.0).entropy)
    paper_rows = []
    for row in PAPER_ECHO_TABLE:
        f220, _ = kerr.qnm(row["mf"], row["af"])
        paper_rows.append({**row, "f220_ours": f220, "dt_ms_ours": 1e3 * kerr.echo_delay(row["mf"], row["af"])})
    return {
        "events": out_events,
        "paper_events": paper_rows,
        "paper_echo_summary": PAPER_ECHO_SUMMARY,
        "paper_ringdown": {"rows": PAPER_RINGDOWN_TABLE, "source": "[QG] Sec. 19.2"},
        "echo_bound_check": {**PAPER_ECHO_BOUND, "log10_R2_ours": check, "computed": True},
        "gwtc": gwtc,
        "combined_R90": PAPER_ECHO_SUMMARY["combined_R90"],
        "source": "[QG] Sec. 19",
        "warnings": warnings,
    }
