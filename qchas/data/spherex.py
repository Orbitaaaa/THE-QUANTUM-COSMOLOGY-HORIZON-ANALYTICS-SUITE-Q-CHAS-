"""SPHEREx Level-2 spectral images from the IRSA public S3 bucket.

The bucket ``nasa-irsa-spherex`` holds the Quick Release (QR) Level-2 calibrated
spectral images (one 2040x2040 MJy/sr frame per detector per exposure, ~65 MB each).
For the all-sky atlas we only need where each frame points, so we read the FITS
headers with HTTP range requests (~0.2 MB per frame) and build the TAN-SIP WCS
locally. A handful of full frames can also be pulled for thumbnails.
"""

from __future__ import annotations

import json
import re
import urllib.request
import warnings
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .registry import SPHEREX_S3

#: nominal SPHEREx detector bands (micron) and spectral resolving power, Crill et al. 2020 (SPIE 11443)
SPHEREX_BANDS = {
    1: (0.75, 1.09, 39), 2: (1.10, 1.62, 41), 3: (1.63, 2.41, 41),
    4: (2.42, 3.82, 35), 5: (3.83, 4.41, 112), 6: (4.42, 5.00, 128),
}
_HEADER_BYTES = 2880 * 80  # primary + IMAGE header (with SIP terms) fit well inside this


@dataclass
class Frame:
    key: str
    obsid: str
    detector: int
    date: str
    ra: float
    dec: float
    corners: list[tuple[float, float]]
    band_um: tuple[float, float]

    def as_dict(self) -> dict:
        return {
            "key": self.key, "obsid": self.obsid, "detector": self.detector, "date": self.date,
            "ra": self.ra, "dec": self.dec, "corners": self.corners, "band_um": list(self.band_um),
        }


def _get(url: str, byte_range: tuple[int, int] | None = None, timeout: float = 60.0) -> bytes:
    headers = {"User-Agent": "qchas/0.1"}
    if byte_range is not None:
        headers["Range"] = f"bytes={byte_range[0]}-{byte_range[1]}"
    req = urllib.request.Request(url, headers=headers)
    last = None
    for _ in range(3):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:  # transient proxy / S3 errors
            last = exc
    raise RuntimeError(f"{url}: {last}")


def list_keys(prefix: str = "qr3/level2/", max_keys: int = 4000, suffix: str = ".fits") -> list[str]:
    """List object keys under ``prefix`` (paginated ListObjectsV2)."""
    keys: list[str] = []
    token = None
    while len(keys) < max_keys:
        url = f"{SPHEREX_S3}/?list-type=2&max-keys=1000&prefix={prefix}"
        if token:
            url += "&continuation-token=" + urllib.request.quote(token, safe="")
        xml = _get(url).decode()
        keys += [k for k in re.findall(r"<Key>([^<]+)</Key>", xml) if k.endswith(suffix)]
        m = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", xml)
        if not m:
            break
        token = m.group(1)
    return keys[:max_keys]


def list_weeks(release: str = "qr3") -> list[str]:
    xml = _get(f"{SPHEREX_S3}/?list-type=2&delimiter=/&prefix={release}/level2/").decode()
    return re.findall(rf"<Prefix>({release}/level2/[^<]+/)</Prefix>", xml)


def _hdu_spans(raw: bytes) -> list[tuple[dict[str, str], int, int]]:
    """Walk the FITS headers present in a byte prefix of a file.

    Returns ``(cards, header_start, header_end)`` per HDU; the data of an HDU starts at
    ``header_end``. Stops when the prefix runs out.
    """
    spans, off = [], 0
    while off + 2880 <= len(raw):
        start, cards, ended = off, [], False
        while off + 2880 <= len(raw):
            blk = raw[off:off + 2880]
            off += 2880
            chunk = [blk[i:i + 80].decode("ascii", "replace") for i in range(0, 2880, 80)]
            cards += chunk
            if any(c.rstrip() == "END" for c in chunk):
                ended = True
                break
        if not ended:
            break
        hdr = {c[:8].strip(): c[10:].split(" /")[0].strip().strip("'").strip() for c in cards if c[8:10] == "= "}
        spans.append((hdr, start, off))
        naxis = int(hdr.get("NAXIS", 0))
        size = abs(int(hdr.get("BITPIX", 8))) // 8 if naxis else 0
        for ax in range(1, naxis + 1):
            size *= int(hdr.get(f"NAXIS{ax}", 0))
        off += -(-size // 2880) * 2880
    return spans


def _image_header(raw: bytes):
    """astropy Header, card dict and data offset of the IMAGE extension (HDU 1)."""
    from astropy.io import fits

    spans = _hdu_spans(raw)
    if len(spans) < 2:
        raise ValueError("IMAGE header not found in prefix")
    cards, start, end = spans[1]
    return fits.Header.fromstring(raw[start:end].decode("ascii", "replace")), cards, end


def frame_from_key(key: str) -> Frame:
    from astropy.wcs import WCS

    raw = _get(f"{SPHEREX_S3}/{key}", (0, _HEADER_BYTES - 1))
    header, cards, _ = _image_header(raw)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w = WCS(header, naxis=2, relax=True)
    nx, ny = int(cards["NAXIS1"]), int(cards["NAXIS2"])
    px = np.array([[0.5, 0.5], [nx + 0.5, 0.5], [nx + 0.5, ny + 0.5], [0.5, ny + 0.5], [(nx + 1) / 2, (ny + 1) / 2]])
    world = w.all_pix2world(px, 1)
    det = int(cards.get("DETECTOR", re.search(r"_(\d)D\d_", key).group(1) if re.search(r"_(\d)D\d_", key) else 0))
    lo, hi, _ = SPHEREX_BANDS.get(det, (np.nan, np.nan, 0))
    return Frame(
        key=key, obsid=cards.get("OBSID", ""), detector=det, date=cards.get("DATE-OBS", ""),
        ra=float(world[4, 0]), dec=float(world[4, 1]),
        corners=[(float(a), float(b)) for a, b in world[:4]], band_um=(lo, hi),
    )


def sample_frames(n_frames: int = 240, release: str = "qr3", workers: int = 12, seed: int = 0) -> list[Frame]:
    """Pick ``n_frames`` Level-2 frames spread over all released weeks and read their WCS."""
    weeks = list_weeks(release)
    per_week = max(1, n_frames // max(1, len(weeks)))
    rng = np.random.default_rng(seed)
    chosen: list[str] = []
    for wk in weeks:
        keys = list_keys(wk, max_keys=3000)
        if not keys:
            continue
        idx = rng.choice(len(keys), size=min(per_week, len(keys)), replace=False)
        chosen += [keys[i] for i in sorted(idx)]
    frames: list[Frame] = []
    with ThreadPoolExecutor(workers) as pool:
        for fr in pool.map(_safe_frame, chosen):
            if fr is not None:
                frames.append(fr)
    return frames


def _safe_frame(key: str) -> Frame | None:
    try:
        return frame_from_key(key)
    except Exception as exc:  # pragma: no cover - network dependent
        warnings.warn(f"SPHEREx header read failed for {key}: {exc}")
        return None


def thumbnail(key: str, size: int = 204) -> dict:
    """Download one full frame's IMAGE plane and return a robustly-scaled thumbnail (uint8)."""
    raw_head = _get(f"{SPHEREX_S3}/{key}", (0, _HEADER_BYTES - 1))
    _, cards, data_start = _image_header(raw_head)
    nx, ny = int(cards["NAXIS1"]), int(cards["NAXIS2"])
    nbytes = nx * ny * 4
    blob = _get(f"{SPHEREX_S3}/{key}", (data_start, data_start + nbytes - 1), timeout=300)
    img = np.frombuffer(blob, dtype=">f4").reshape(ny, nx).astype(np.float64)
    f = nx // size
    img = np.nanmedian(img[: f * size, : f * size].reshape(size, f, size, f).swapaxes(1, 2).reshape(size, size, -1), axis=2)
    finite = np.isfinite(img)
    lo, hi = np.percentile(img[finite], [1, 99.5]) if finite.any() else (0.0, 1.0)
    scaled = np.clip((np.arcsinh((img - lo) / max(hi - lo, 1e-12) * 10) / np.arcsinh(10)), 0, 1)
    scaled[~finite] = 0
    return {"key": key, "size": size, "unit": cards.get("BUNIT", "MJy/sr"), "p1": float(lo), "p995": float(hi),
            "pixels": (scaled * 255).astype(np.uint8).ravel().tolist()}


def save_frames(frames: list[Frame], path: Path) -> None:
    path.write_text(json.dumps([f.as_dict() for f in frames]))


def load_frames(path: Path) -> list[dict]:
    return json.loads(Path(path).read_text())


__all__ = ["SPHEREX_BANDS", "Frame", "list_keys", "list_weeks", "frame_from_key", "sample_frames", "thumbnail"]

