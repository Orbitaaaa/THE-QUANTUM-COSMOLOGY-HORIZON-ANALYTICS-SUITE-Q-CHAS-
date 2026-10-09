"""Real LIGO strain pipeline and GWOSC catalogue parsing, qchas.data.gw ([QG] Sec. 19)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from qchas.data import gw
from qchas.data.fetch import DEFAULT_RAW_DIR
from qchas.horizon import kerr

RAW = Path(DEFAULT_RAW_DIR)
H1_150914 = RAW / "H-H1_LOSC_4_V2-1126259446-32.hdf5"
needs_strain = pytest.mark.skipif(not H1_150914.exists() or not (RAW / gw.TUTORIAL_EVENTS_FILE).exists(),
                                  reason="LIGO tutorial strain not downloaded (run `qchas fetch`)")


# --- synthetic checks (no data needed) -----------------------------------------------------

def test_detector_lag_recovers_shift_and_sign():
    fs = 4096.0
    rng = np.random.default_rng(0)
    t = np.arange(int(0.3 * fs)) / fs
    chirp = np.sin(2 * np.pi * (60 * t + 300 * t**2)) * np.exp(-((t - 0.2) / 0.05) ** 2)
    shift = 29  # samples: L1 later by 7.08 ms
    h1 = chirp + 0.05 * rng.standard_normal(t.size)
    l1 = -np.roll(chirp, shift) + 0.05 * rng.standard_normal(t.size)
    r = gw.detector_lag(h1, l1, fs)
    assert r["lag_s"] == pytest.approx(shift / fs)
    assert r["sign"] == -1 and r["corr"] < -0.8


def test_damped_sinusoid_fit_recovers_parameters():
    fs = 4096.0
    t = np.arange(int(0.03 * fs)) / fs
    rng = np.random.default_rng(1)
    y = gw.damped_sinusoid(t, 3.0, 0.004, 250.0, 0.7) + 0.1 * rng.standard_normal(t.size)
    r = gw._fit_damped(t, y, 230.0, 0.005)
    assert r["f"] == pytest.approx(250.0, abs=2.0)
    assert r["tau"] == pytest.approx(0.004, rel=0.05)
    assert r["amp"] == pytest.approx(3.0, rel=0.05)


def test_ringdown_fit_noise_bootstrap_and_significance():
    fs = 4096.0
    rng = np.random.default_rng(2)
    noise = rng.standard_normal(int(20 * fs))
    data = rng.standard_normal(int(1 * fs))
    i_peak = 1000
    t = np.arange(data.size - i_peak) / fs
    data[i_peak:] += gw.damped_sinusoid(t, 6.0, 0.004, 250.0, 0.3)
    r = gw.ringdown_fit(data, fs, i_peak, 245.0, 0.004, noise=noise, n_noise=20, seed=3)
    assert abs(r["f_Hz"] - 250.0) < 3 * r["f_err_Hz"]
    assert abs(r["tau_ms"] - 4.0) < 3 * r["tau_err_ms"]
    assert r["significant"] and r["p_noise"] < 0.05
    pure = gw.ringdown_fit(rng.standard_normal(int(1 * fs)), fs, i_peak, 245.0, 0.004, noise=noise,
                           n_noise=20, seed=4)
    assert pure["p_noise"] > 0.05


def test_whiten_white_noise_has_unit_variance():
    fs = 4096.0
    rng = np.random.default_rng(5)
    x = 1e-21 * rng.standard_normal(int(32 * fs))
    w = gw.whiten(x, fs)
    core = w[int(4 * fs): -int(4 * fs)]
    assert np.std(core) == pytest.approx(1.0, rel=0.05)


def _gwosc_sample() -> dict:
    return {"events": {
        "GW150914-v3": {"commonName": "GW150914", "GPS": 1126259462.4, "version": 3, "mass_1_source": 35.6,
                        "mass_2_source": 30.6, "final_mass_source": 63.1, "chi_eff": -0.01, "redshift": 0.09,
                        "luminosity_distance": 440.0, "network_matched_filter_snr": 24.4,
                        "catalog.shortName": "GWTC-1-confident"},
        "GW150914-v1": {"commonName": "GW150914", "GPS": 1126259462.4, "version": 1, "mass_1_source": 99.0,
                        "mass_2_source": None, "catalog.shortName": "O1_O2-Preliminary"},
        "GW190425-v2": {"commonName": "GW190425", "GPS": "1240215503.0", "mass_1_source": 2.0,
                        "mass_2_source": 1.4, "final_mass_source": None, "chi_eff": None, "redshift": None,
                        "luminosity_distance": None, "network_matched_filter_snr": "NaN",
                        "catalog.shortName": "GWTC-2"},
        "S200105ae-v1": {"commonName": None, "GPS": None},
    }}


def test_parse_gwosc_allevents_robust():
    rows = gw.parse_gwosc_allevents(_gwosc_sample())
    names = [r["name"] for r in rows]
    assert names == ["GW150914", "GW190425", "S200105ae"]
    g = rows[0]
    assert g["version"] == 3 and g["m1_source"] == 35.6 and g["catalog"] == "GWTC-1-confident"
    assert g["f220_Hz"] == pytest.approx(kerr.qnm(63.1 * 1.09, g["af_estimate"])[0])
    assert 0.6 < g["af_estimate"] < 0.75
    n = rows[1]
    assert n["gps"] == pytest.approx(1240215503.0) and n["snr"] is None and n["f220_Hz"] is None
    assert rows[2]["gps"] is None
    json.dumps(rows)


def test_to_json_without_data(tmp_path):
    d = gw.to_json(tmp_path)
    assert d["events"] == [] and d["gwtc"] is None
    assert d["combined_R90"] == 0.28
    assert any("GWOSC" in w for w in d["warnings"])
    assert len(d["paper_events"]) == 9
    for row in d["paper_events"]:
        assert row["f220_ours"] == pytest.approx(row["f220"], rel=0.01)
        assert row["dt_ms_ours"] == pytest.approx(row["dt_ms"], rel=0.015)
    assert d["echo_bound_check"]["log10_R2_ours"] / 1e80 == pytest.approx(-1.75, abs=0.005)
    labels = {r["event"]: r["R90_label"] for r in d["paper_events"]}
    assert labels["GW150914"] == "0.67" and labels["GW151012"] == "> 1" and labels["GW250114"] == "0.28"


def test_malformed_catalogue_is_reported(tmp_path):
    (tmp_path / gw.GWTC_FILE).write_text("[1, 2, 3]")
    d = gw.to_json(tmp_path)
    assert d["gwtc"] is None and any("unreadable" in w for w in d["warnings"])


def test_missing_strain_is_reported(tmp_path):
    params = {"tevent": 1.0, "m1": 36.0, "m2": 29.0, "a1": 0.0, "a2": 0.0, "fband": [43, 300],
              "fn_H1": "nope_H1.hdf5", "fn_L1": "nope_L1.hdf5"}
    ev = gw.analyse_event("TEST", params, tmp_path)
    assert ev["strain"] is None and ev["ringdown_fit"] is None
    assert ev["warnings"] and ev["f220_Hz"] > 0


def test_unreadable_strain_is_reported(tmp_path):
    for det in ("H1", "L1"):
        (tmp_path / f"bad_{det}.hdf5").write_bytes(b"not an HDF5 file")
    params = {"tevent": 1.0, "m1": 36.0, "m2": 29.0, "a1": 0.0, "a2": 0.0, "fband": [43, 300],
              "fn_H1": "bad_H1.hdf5", "fn_L1": "bad_L1.hdf5"}
    ev = gw.analyse_event("TEST", params, tmp_path)
    assert ev["strain"] is None and ev["ringdown_fit"] is None
    assert any("unreadable" in w for w in ev["warnings"])
    (tmp_path / gw.TUTORIAL_EVENTS_FILE).write_text(json.dumps({"TEST": params, "BROKEN": {"m1": 30.0}}))
    d = gw.to_json(tmp_path)
    assert [e["name"] for e in d["events"]] == ["TEST"]
    assert any("BROKEN" in w for w in d["warnings"]) and any("unreadable" in w for w in d["warnings"])


def test_catalogue_neutron_star_events_get_no_bbh_remnant():
    rows = gw.parse_gwosc_allevents({"events": {
        "GW200115_042309-v2": {"commonName": "GW200115_042309", "GPS": 1263097407.7, "mass_1_source": 5.9,
                               "mass_2_source": 1.44, "chi_eff": -0.19, "redshift": 0.06},
        "GW170104-v2": {"commonName": "GW170104", "GPS": 1167559936.6, "mass_1_source": 30.8,
                        "mass_2_source": 20.0, "chi_eff": -0.04, "redshift": 0.20}}})
    nsbh, bbh = rows[1], rows[0]
    assert nsbh["possible_neutron_star"] and nsbh["f220_Hz"] is None and nsbh["af_estimate"] is None
    assert not bbh["possible_neutron_star"] and bbh["f220_Hz"] > 0


def test_remnant_summary_units():
    s = gw.remnant_summary(68.0, 0.68, mf_thermo=62.0)
    assert s["S_BH_log10"] == pytest.approx(math.log10(kerr.KerrHorizon(62.0, 0.68).entropy))
    assert s["echo_bound_log10"] == pytest.approx(-kerr.KerrHorizon(62.0, 0.68).entropy / math.log(10))
    assert s["page_time_log10_yr"] < s["t_evap_log10_yr"]
    assert 70 < s["t_evap_log10_yr"] < 74


# --- the real strain -----------------------------------------------------------------------

@needs_strain
def test_read_losc_hdf5():
    strain, gps0, fs = gw.read_losc_hdf5(H1_150914)
    assert fs == 4096.0 and gps0 == 1126259446.0 and strain.size == 32 * 4096
    assert 1e-24 < np.std(strain) < 1e-16


@needs_strain
def test_tutorial_events_and_gw150914_remnant():
    ev = gw.load_tutorial_events()
    assert set(ev) == {"GW150914", "LVT151012", "GW151226", "GW170104"}
    p = ev["GW150914"]
    mf, af = kerr.remnant_aligned(p["m1"], p["m2"], p["a1"], p["a2"])
    assert 62.0 < mf < 68.5
    assert 0.64 < af < 0.70
    s = gw.remnant_summary(mf, af)
    assert s["f220_Hz"] == pytest.approx(250, abs=8)
    assert s["tau220_ms"] == pytest.approx(4.0, abs=0.3)
    assert s["echo_delay_ms"] == pytest.approx(290, abs=10)


@needs_strain
def test_gw150914_whitened_strain_lag_and_ringdown():
    ev = gw.load_tutorial_events()
    r = gw.analyse_event("GW150914", ev["GW150914"], n_noise=20)
    st = r["strain"]
    assert len(st["t"]) == len(st["H1"]) == len(st["L1"])
    assert st["fs_export"] == 2048.0
    assert st["t"][0] == pytest.approx(-0.2, abs=1e-3) and st["t"][-1] == pytest.approx(0.2, abs=1e-3)
    # L1 sees the signal ~7 ms earlier and inverted (detectors nearly anti-aligned)
    assert -8.0 < st["lag_ms"] < -6.0 and st["lag_sign"] == -1 and st["lag_corr"] < -0.7
    # the whitened, band-passed chirp stands well above the noise away from the merger
    h1, t = np.array(st["H1"]), np.array(st["t"])
    assert np.max(np.abs(h1)) > 6.0 * np.std(h1[np.abs(t) > 0.15])
    fit = r["ringdown_fit"]
    assert fit["significant"]
    assert abs(fit["f_Hz"] - r["f220_Hz"]) < 2.5 * fit["f_err_Hz"]
    assert abs(fit["tau_ms"] - r["tau220_ms"]) < 2.5 * fit["tau_err_ms"]
    assert 230 < fit["f_Hz"] < 265
    json.dumps(r)


@pytest.mark.slow
@needs_strain
def test_gw_to_json_full():
    d = gw.to_json(n_noise=40)
    assert [e["name"] for e in d["events"]] == ["GW150914", "LVT151012", "GW151226", "GW170104"]
    text = json.dumps(d)
    assert "NaN" not in text and "Infinity" not in text
    for e in d["events"]:
        assert e["strain"] is not None and len(e["strain"]["H1"]) > 800
        assert e["ringdown_fit"]["p_noise"] is not None
    assert d["events"][0]["ringdown_fit"]["significant"]
