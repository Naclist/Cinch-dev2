#!/usr/bin/env python3
"""Refresh package integrity hashes after an intentional source freeze."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "src" / "cinch_dev2"
WORKFLOW = PACKAGE / "workflow.json"
MANIFEST = ROOT / "SOURCE_SHA256.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def selected_files() -> list[Path]:
    paths = [PACKAGE / "core" / "statistics.py", PACKAGE / "scripts" / "cinch_pipeline.py"]
    paths.extend(sorted((PACKAGE / "scripts" / "snapshots").rglob("*.py")))
    return sorted(set(paths))


def main() -> None:
    data = json.loads(WORKFLOW.read_text(encoding="utf-8"))
    paths = selected_files()
    data["frozen_files"] = [
        {"path": str(path.relative_to(PACKAGE)), "sha256": sha256(path)} for path in paths
    ]
    WORKFLOW.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    repository_files = sorted(
        path for path in ROOT.rglob("*")
        if path.is_file()
        and ".git" not in path.parts
        and ".venv" not in path.parts
        and "__pycache__" not in path.parts
        and path != MANIFEST
    )
    rows = ["sha256\tpath"]
    rows.extend(f"{sha256(path)}\t{path.relative_to(ROOT)}" for path in repository_files)
    MANIFEST.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"Frozen {len(paths)} source files; indexed {len(repository_files)} repository files")


if __name__ == "__main__":
    main()

