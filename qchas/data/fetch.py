"""Download public datasets with mirror fallback and a provenance record.

Uses only the standard library so it runs anywhere; ``urllib`` honours ``HTTPS_PROXY``
and the system CA bundle. Each successful download is recorded in
``<raw_dir>/provenance.json`` (source URL, size, SHA-256, UTC time), and the dashboard
shows that record next to every result derived from the file.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import shutil
import sys
import time
import urllib.error
import urllib.request
from dataclasses import asdict
from pathlib import Path

from .registry import DATASETS, Dataset

DEFAULT_RAW_DIR = Path(os.environ.get("QCHAS_RAW_DIR", Path(__file__).resolve().parents[2] / "data" / "raw"))
USER_AGENT = "qchas/0.1 (+https://github.com/Orbitaaaa/THE-QUANTUM-COSMOLOGY-HORIZON-ANALYTICS-SUITE-Q-CHAS-)"


class FetchError(RuntimeError):
    pass


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _provenance_path(raw_dir: Path) -> Path:
    return raw_dir / "provenance.json"


def load_provenance(raw_dir: Path | None = None) -> dict:
    raw_dir = Path(raw_dir or DEFAULT_RAW_DIR)
    p = _provenance_path(raw_dir)
    if p.exists():
        return json.loads(p.read_text())
    return {}


def _record(raw_dir: Path, ds: Dataset, url: str, path: Path) -> None:
    prov = load_provenance(raw_dir)
    prov[ds.key] = {
        "filename": ds.filename,
        "url": url,
        "bytes": path.stat().st_size,
        "sha256": _sha256(path),
        "retrieved_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "mission": ds.mission,
        "description": ds.description,
        "citation": ds.citation,
    }
    _provenance_path(raw_dir).write_text(json.dumps(prov, indent=2, sort_keys=True))


def _download(url: str, dest: Path, timeout: float, retries: int, quiet: bool) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "wb") as out:
                total = int(resp.headers.get("Content-Length") or 0)
                done = 0
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
                    done += len(chunk)
                    if not quiet and total > 50_000_000:
                        sys.stderr.write(f"\r  {done / 1e6:8.1f} / {total / 1e6:.1f} MB")
                if not quiet and total > 50_000_000:
                    sys.stderr.write("\n")
            shutil.move(tmp, dest)
            return
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            last = exc
            if tmp.exists():
                tmp.unlink()
            time.sleep(2 ** attempt)
    raise FetchError(f"{url}: {last}")


def fetch(key: str, raw_dir: Path | None = None, *, force: bool = False, timeout: float = 120.0,
          retries: int = 3, quiet: bool = False) -> Path:
    """Return the local path of dataset ``key``, downloading it from the first working mirror."""
    ds = DATASETS[key]
    raw_dir = Path(raw_dir or DEFAULT_RAW_DIR)
    raw_dir.mkdir(parents=True, exist_ok=True)
    dest = raw_dir / ds.filename
    if dest.exists() and dest.stat().st_size >= ds.min_bytes and not force:
        prov = load_provenance(raw_dir)
        if key not in prov:
            _record(raw_dir, ds, "(pre-existing local file)", dest)
        return dest
    errors = []
    for url in ds.mirrors:
        try:
            if not quiet:
                print(f"[fetch] {key}: {url}", file=sys.stderr)
            _download(url, dest, timeout, retries, quiet)
            if dest.stat().st_size < ds.min_bytes:
                head = dest.read_bytes()[:200]
                dest.unlink()
                raise FetchError(f"{url}: only {len(head)}+ bytes (archive error page?)")
            _record(raw_dir, ds, url, dest)
            return dest
        except FetchError as exc:
            errors.append(str(exc))
    raise FetchError(f"could not fetch {key} from any mirror:\n  " + "\n  ".join(errors))


def available(key: str, raw_dir: Path | None = None) -> Path | None:
    """Local path of ``key`` if it is already on disk, else ``None`` (never downloads)."""
    ds = DATASETS[key]
    p = Path(raw_dir or DEFAULT_RAW_DIR) / ds.filename
    return p if p.exists() and p.stat().st_size >= ds.min_bytes else None


def fetch_many(keys, raw_dir: Path | None = None, **kw) -> dict[str, Path | str]:
    """Fetch several datasets; failures are returned as error strings instead of raising."""
    out: dict[str, Path | str] = {}
    for k in keys:
        try:
            out[k] = fetch(k, raw_dir, **kw)
        except FetchError as exc:
            out[k] = str(exc)
    return out


def describe() -> list[dict]:
    return [asdict(d) for d in DATASETS.values()]
