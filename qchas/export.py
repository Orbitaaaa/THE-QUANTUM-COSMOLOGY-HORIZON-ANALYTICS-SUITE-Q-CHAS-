"""Build every dashboard bundle of the data contract (docs/ARCHITECTURE.md) from real data.

Each bundle is built in isolation: a failure (for example an unreachable archive) is
recorded in ``manifest.json`` (``warnings`` / ``errors``) and the other bundles are still
written. Nothing here invents data; a missing input produces ``null`` and a warning.
"""

from __future__ import annotations

import datetime as _dt
import json
import math
import time
import traceback
from pathlib import Path
from typing import Callable

import numpy as np

from . import __version__
from .data.fetch import DEFAULT_RAW_DIR, load_provenance


# --- JSON helpers ----------------------------------------------------------------------

def _clean(obj):
    """Make ``obj`` strict-JSON serialisable: numpy -> python, NaN/inf -> None."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _clean(obj.tolist())
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        f = float(obj)
        return f if math.isfinite(f) else None
    if isinstance(obj, Path):
        return str(obj)
    return obj


def _write_json(path: Path, obj) -> int:
    text = json.dumps(_clean(obj), separators=(",", ":"), allow_nan=False)
    path.write_text(text)
    return len(text)


class _Builder:
    def __init__(self, out_dir: Path):
        self.out = out_dir
        self.bundles: list[str] = []
        self.warnings: list[str] = []
        self.errors: list[dict] = []
        self.timings: dict[str, float] = {}
        self.results: dict[str, object] = {}

    def run(self, name: str, fn: Callable[[], object], filename: str | None = None):
        t0 = time.perf_counter()
        print(f"[build] {name} ...", flush=True)
        try:
            obj = fn()
        except Exception as exc:  # recorded, not fatal: other bundles still build
            self.errors.append({"bundle": name, "error": f"{type(exc).__name__}: {exc}",
                                "trace": traceback.format_exc(limit=6)})
            print(f"[build] {name} FAILED: {exc}", flush=True)
            return None
        finally:
            self.timings[name] = round(time.perf_counter() - t0, 2)
        self.results[name] = obj
        if filename and obj is not None:
            size = _write_json(self.out / filename, obj)
            self.bundles.append(filename)
            print(f"[build] {name} -> {filename} ({size / 1e3:.0f} kB, {self.timings[name]} s)", flush=True)
        return obj


# --- bundle builders ---------------------------------------------------------------------

def _topology() -> dict:
    from .topology import circles, manifolds, vacuum

    out = manifolds.to_json()
    out["circles"] = circles.to_json()
    out["vacuum"] = vacuum.to_json()
    return out


def _lattice(fast: bool) -> dict:
    from .topology import lattice, shapeflow

    return {**lattice.to_json(), **shapeflow.to_json(n=24 if fast else 240)}


def _sky(b: _Builder, raw_dir: Path):
    from .cmb import maps

    sky = maps.prepare_sky(raw_dir=raw_dir)
    meta, arrays = maps.sky_map_export(sky)
    for fname, arr in arrays.items():
        (b.out / fname).write_bytes(np.ascontiguousarray(arr).tobytes())
        b.bundles.append(fname)
    b.warnings += list(sky.meta.get("warnings", []))
    return sky, meta


def _cmb(sky, raw_dir: Path, fast: bool, likelihood_part: dict | None) -> dict:
    from .cmb import covariance, resonance, stats

    wmap = stats.load_wmap_tt(raw_dir)
    planck = stats.load_planck_tt(raw_dir)
    out = stats.cmb_fields(sky, wmap, planck, n_realizations=20_000 if fast else 100_000)
    out["torus_D"] = covariance.torus_D_bundle(
        models=(("FCC", 0.985), ("cubic", 0.985), ("BCC", 0.985), ("FCC", 0.70), ("FCC", 1.05)))
    out["kl"] = covariance.kl_table()
    out["sec8"] = covariance.sec8_table()
    out["ceiling"] = resonance.PAPER_CEILING
    out["likelihood"] = None if likelihood_part is None else likelihood_part.get("likelihood")
    if likelihood_part is not None and "likelihood_meta" in likelihood_part:
        out["likelihood_meta"] = likelihood_part["likelihood_meta"]
    return out


def _likelihood(sky, fast: bool, null_sims: int) -> dict | None:
    try:
        from .cmb import likelihood
    except ImportError:
        return None
    return likelihood.export_fields(sky, fast=fast, null_sims=null_sims)


def _resonance(fast: bool, likelihood_part: dict | None) -> dict:
    from .cmb import resonance

    out = resonance.resonance_bundle(benchmark_orientations=5 if fast else 20)
    out["search"] = None if likelihood_part is None else likelihood_part.get("search")
    out["source"] = "[TOP] Secs. 11.6-11.8"
    return out


def _spherex(n_frames: int, out_dir: Path) -> dict:
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    from .data import spherex

    frames = spherex.sample_frames(n_frames)
    if not frames:
        raise RuntimeError("no SPHEREx frames could be read")
    ra = np.array([[c[0] for c in f.corners] + [f.ra] for f in frames])
    dec = np.array([[c[1] for c in f.corners] + [f.dec] for f in frames])
    gal = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs").galactic
    gl, gb = gal.l.deg, gal.b.deg
    rows = []
    for i, f in enumerate(frames):
        d = f.as_dict()
        d.pop("corners")
        d.update({"l": gl[i, 4], "b": gb[i, 4], "corners_gal": [[gl[i, j], gb[i, j]] for j in range(4)],
                  "corners_eq": [list(c) for c in f.corners]})
        rows.append(d)
    thumbs = []
    thumb_dir = out_dir / "spherex"
    thumb_dir.mkdir(exist_ok=True)
    picked = {}
    for f, row in zip(frames, rows):
        picked.setdefault(f.detector, (f, row))
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import image as mpimg

    for det, (f, row) in sorted(picked.items()):
        try:
            th = spherex.thumbnail(f.key)
        except Exception as exc:  # pragma: no cover - network dependent
            print(f"[build] SPHEREx thumbnail {f.key} failed: {exc}")
            continue
        img = np.asarray(th["pixels"], dtype=np.uint8).reshape(th["size"], th["size"])
        name = f"spherex/{f.obsid}_D{det}.png"
        mpimg.imsave(out_dir / name, img, cmap="inferno", vmin=0, vmax=255, origin="lower")
        thumbs.append({"obsid": f.obsid, "detector": det, "file": name, "l": row["l"], "b": row["b"],
                       "band_um": list(f.band_um), "unit": th["unit"], "p1": th["p1"], "p995": th["p995"],
                       "key": f.key, "date": f.date})
    return {"release": "qr3", "bucket": "nasa-irsa-spherex", "n_frames": len(rows),
            "bands": {str(k): list(v) for k, v in spherex.SPHEREX_BANDS.items()},
            "frames": rows, "thumbnails": thumbs,
            "note": ("Frames are a uniform random sample of the released Quick-Release Level-2 images, "
                     "located from their own FITS WCS headers (TAN-SIP)."),
            "source": "NASA/IPAC IRSA, s3://nasa-irsa-spherex"}


def _horizon(fast: bool) -> dict:
    from . import claims
    from .horizon import atmosphere, capacity, evaporation, kerr, page, qes

    qes_rows = []
    for m, c, rb, depth, lag, tscr, ds in qes.PAPER_TABLE:
        r = qes.solve_qes(m, c, rb).as_dict()
        r["paper"] = {"depth": depth, "lag": lag, "t_scr": tscr, "s_gen_minus_s_bh": ds}
        qes_rows.append(r)
    atm = atmosphere.atmosphere_profile(ell_max=600 if fast else 3000)
    atm["paper_share"] = [atmosphere.PAPER_SHARE.get(g) for g in atm["gaps"]]
    atm["continuum"] = {"area_coeff": 0.2954, "f1": 0.1483, "f2": 0.0305, "f3": 0.0133, "source": "[BH] Sec. 15"}
    m_ex = 10.0
    beta = 8 * math.pi * m_ex
    windows = [capacity.capacity_window(m_ex, x / beta) for x in (0.01, 0.1, 1.0, 3.0, 10.0, 30.0, 100.0)]
    evap = evaporation.to_json(n_qubits=12 if fast else 20)
    return {
        "page": page.curves(),
        "capacity": {
            "ladder": capacity.ALGEBRA_LADDER,
            "de_sitter": capacity.de_sitter_entropy(),
            "min_area_lp2": capacity.minimal_area(),
            "window_examples": windows,
            "window_mass_planck": m_ex,
        },
        "qes": {"rows": qes_rows,
                "laws": {"depth_coeff": 0.0399, "depth_err": 0.0007, "lag_offset_rs": 4.24, "source": "[BH] Sec. 9"}},
        "atmosphere": atm,
        "evaporation": evap["evaporation"],
        "decoding": evap["decoding"],
        "kerr_spin_down": kerr.to_json(),
        "claims": evap.get("claims", []) + claims.horizon_claims(),
    }


def _bounce(fast: bool) -> dict:
    from .horizon import bounce

    return bounce.to_json(fast=fast)


def _gw(b: _Builder, raw_dir: Path) -> dict:
    from .data import gw

    out = gw.to_json(raw_dir)
    b.warnings += list(out.pop("warnings", []) or [])
    return out


# --- orchestration ---------------------------------------------------------------------

def build_all(out_dir: Path | None = None, raw_dir: Path | None = None, *, fast: bool = False,
              spherex_frames: int = 240, null_sims: int = 100, only: list[str] | None = None) -> dict:
    """Run every analysis and write the bundles; returns the manifest."""
    out_dir = Path(out_dir or Path(__file__).resolve().parents[1] / "web" / "public" / "data")
    raw_dir = Path(raw_dir or DEFAULT_RAW_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    b = _Builder(out_dir)
    want = (lambda name: True) if not only else (lambda name: name in only)

    if want("topology"):
        b.run("topology", _topology, "topology.json")
    if want("lattice"):
        b.run("lattice", lambda: _lattice(fast), "lattice.json")

    sky = None
    if want("sky") or want("cmb") or want("resonance"):
        for stale in ("sky_map_full.bin", "sky_map_lowl.bin", "sky_mask.bin"):
            (out_dir / stale).unlink(missing_ok=True)
        res = b.run("sky", lambda: _sky(b, raw_dir))
        if res is not None:
            sky, sky_meta = res
            _write_json(out_dir / "sky_map.json", sky_meta)
            b.bundles.append("sky_map.json")
    lik = None
    if sky is not None and (want("cmb") or want("resonance")):
        lik = b.run("likelihood", lambda: _likelihood(sky, fast, null_sims))
        if lik is None:
            b.warnings.append("Torus likelihood and resonance search not computed (qchas.cmb.likelihood unavailable or failed)")
    if want("cmb"):
        b.run("cmb", lambda: _cmb(sky, raw_dir, fast, lik), "cmb.json")
    if want("resonance"):
        b.run("resonance", lambda: _resonance(fast, lik), "resonance.json")
    if want("spherex") and spherex_frames > 0:
        if b.run("spherex", lambda: _spherex(spherex_frames, out_dir), "spherex.json") is None:
            b.warnings.append("SPHEREx frames could not be read from the IRSA S3 bucket on this host")
    if want("horizon"):
        b.run("horizon", lambda: _horizon(fast), "horizon.json")
    if want("bounce"):
        b.run("bounce", lambda: _bounce(fast), "bounce.json")
    if want("gw"):
        b.run("gw", lambda: _gw(b, raw_dir), "gw.json")
    if want("validation"):
        from .validation import collect

        b.run("validation", lambda: collect(out_dir), "validation.json")

    prov = load_provenance(raw_dir)
    manifest = {
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "qchas_version": __version__,
        "cmb_map_source": None if sky is None else sky.source_name,
        "fast_build": fast,
        "datasets": prov,
        "bundles": sorted(set(b.bundles + _existing(out_dir))),
        "warnings": sorted(set(b.warnings)),
        "errors": b.errors,
        "timings_s": b.timings,
    }
    _write_json(out_dir / "manifest.json", manifest)
    return manifest


def _existing(out_dir: Path) -> list[str]:
    """Bundles already on disk from earlier partial builds (``--only``).

    Binary maps count only when the current ``sky_map.json`` still references them, so a stale
    mask from an earlier build is never advertised.
    """
    names = ["topology.json", "lattice.json", "sky_map.json", "cmb.json", "resonance.json", "spherex.json",
             "horizon.json", "bounce.json", "gw.json", "validation.json"]
    found = [n for n in names if (out_dir / n).exists()]
    sky = out_dir / "sky_map.json"
    if sky.exists():
        meta = json.loads(sky.read_text())
        refs = [m.get("file") for m in meta.get("maps", [])] + [meta.get("mask_file")]
        found += [f for f in refs if f and (out_dir / f).exists()]
    return found
