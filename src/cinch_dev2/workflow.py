"""Workflow metadata and frozen-file verification."""

from __future__ import annotations

import hashlib
import json
from importlib.resources import files
from pathlib import Path


def workflow_path() -> Path:
    return Path(str(files("cinch_dev2").joinpath("workflow.json")))


def load_workflow() -> dict:
    return json.loads(workflow_path().read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_frozen_files(package_root: Path | None = None) -> list[dict]:
    root = package_root or Path(__file__).resolve().parent
    results = []
    for item in load_workflow().get("frozen_files", []):
        path = root / item["path"]
        actual = sha256(path) if path.is_file() else None
        results.append(
            {
                "path": item["path"],
                "expected": item["sha256"],
                "actual": actual,
                "status": "OK" if actual == item["sha256"] else ("MISSING" if actual is None else "CHANGED"),
            }
        )
    return results

