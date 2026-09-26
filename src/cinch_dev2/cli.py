"""Command-line interface for the Cinch dev2 research preview."""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

from . import __version__
from .workflow import load_workflow, verify_frozen_files


def _pipeline_path() -> Path:
    return Path(__file__).resolve().parent / "scripts" / "cinch_pipeline.py"


def _list_stages(as_json: bool) -> int:
    stages = load_workflow()["stages"]
    if as_json:
        print(json.dumps(stages, indent=2))
    else:
        for stage in stages:
            print(f"{stage['id']:30s} {stage['track']:10s} {stage['status']:12s} {stage['name']}")
    return 0


def _verify() -> int:
    rows = verify_frozen_files()
    for row in rows:
        print(f"{row['status']:8s} {row['path']}")
    return int(any(row["status"] != "OK" for row in rows))


def _doctor() -> int:
    modules = ["numpy", "pandas", "scipy", "sklearn", "matplotlib", "seaborn", "Bio", "numba", "umap", "yaml"]
    failed = False
    for module in modules:
        available = importlib.util.find_spec(module) is not None
        print(f"{'OK' if available else 'MISSING':8s} python:{module}")
        failed |= not available
    print("OPTIONAL external:configure and external:uberBlast are required only for the legacy WGS mapping stage")
    return int(failed)


def _run(extra: list[str]) -> int:
    if extra[:1] == ["--"]:
        extra = extra[1:]
    command = [sys.executable, str(_pipeline_path()), *extra]
    return subprocess.run(command, check=False).returncode


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cinch-dev2", description=__doc__)
    parser.add_argument("--version", action="version", version=f"Cinch dev2 {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)
    list_parser = sub.add_parser("list", help="list registered analysis stages")
    list_parser.add_argument("--json", action="store_true")
    sub.add_parser("verify", help="verify frozen source snapshots")
    sub.add_parser("doctor", help="check runtime dependencies")
    run = sub.add_parser("run", help="run the frozen staged workflow")
    run.add_argument("pipeline_args", nargs=argparse.REMAINDER, help="arguments passed to the pipeline runner")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "list":
        return _list_stages(args.json)
    if args.command == "verify":
        return _verify()
    if args.command == "doctor":
        return _doctor()
    return _run(args.pipeline_args)


if __name__ == "__main__":
    raise SystemExit(main())
