"""Byte identities for verifying build inputs and installed artifacts."""

import hashlib
from pathlib import Path


def record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def check(item):
    if record(item["path"]) != item:
        raise ValueError("Build input/source identity changed: " + item["path"])
