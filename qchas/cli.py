"""Command line entry point: ``qchas fetch | build | info``."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .data.fetch import DEFAULT_RAW_DIR, describe, fetch_many, load_provenance
from .data.registry import DATASETS

ALWAYS = [k for k, d in DATASETS.items() if "large" not in d.tags]
PLANCK_LARGE = [k for k, d in DATASETS.items() if "large" in d.tags]


def _cmd_fetch(args: argparse.Namespace) -> int:
    keys = list(args.keys) if args.keys else ALWAYS + (PLANCK_LARGE if args.planck else [])
    unknown = [k for k in keys if k not in DATASETS]
    if unknown:
        print(f"unknown dataset keys: {unknown}; see 'qchas info'", file=sys.stderr)
        return 2
    results = fetch_many(keys, Path(args.raw_dir), timeout=args.timeout)
    failed = {k: v for k, v in results.items() if isinstance(v, str)}
    for k, v in results.items():
        status = "ok  " if not isinstance(v, str) else "FAIL"
        print(f"[{status}] {k:28s} {v if not isinstance(v, str) else v.splitlines()[0]}")
    if failed:
        print(f"\n{len(failed)} dataset(s) unreachable from this host; analyses fall back to the "
              "datasets that are present and the build records a warning.", file=sys.stderr)
    return 1 if args.strict and failed else 0


def _cmd_build(args: argparse.Namespace) -> int:
    from .export import build_all

    manifest = build_all(out_dir=Path(args.out), raw_dir=Path(args.raw_dir), fast=args.fast,
                         spherex_frames=args.spherex_frames, null_sims=args.null_sims)
    print(json.dumps({"bundles": manifest["bundles"], "warnings": manifest["warnings"]}, indent=2))
    return 0


def _cmd_info(args: argparse.Namespace) -> int:
    prov = load_provenance(Path(args.raw_dir))
    for d in describe():
        have = "local" if d["key"] in prov else "     "
        print(f"{have}  {d['key']:28s} {d['mission']:16s} {d['description']}")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="qchas", description="Quantum Cosmology & Horizon Analytics Suite")
    p.add_argument("--raw-dir", default=str(DEFAULT_RAW_DIR), help="where raw downloads live")
    sub = p.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="download the public datasets (mirrors tried in order)")
    f.add_argument("keys", nargs="*", help="dataset keys (default: everything except the 1.6 GB Planck map)")
    f.add_argument("--planck", action="store_true", help="also fetch the Planck PR3 Commander map (1.6 GB)")
    f.add_argument("--strict", action="store_true", help="exit non-zero if any dataset is unreachable")
    f.add_argument("--timeout", type=float, default=120.0)
    f.set_defaults(func=_cmd_fetch)

    b = sub.add_parser("build", help="run the analyses and write the dashboard bundles")
    b.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "web" / "public" / "data"))
    b.add_argument("--fast", action="store_true", help="smaller grids and fewer simulations (CI smoke build)")
    b.add_argument("--spherex-frames", type=int, default=240, help="SPHEREx frames to locate (0 to skip)")
    b.add_argument("--null-sims", type=int, default=100, help="isotropic simulations for likelihood calibration")
    b.set_defaults(func=_cmd_build)

    i = sub.add_parser("info", help="list the dataset registry and what is on disk")
    i.set_defaults(func=_cmd_info)

    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
