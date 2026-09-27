"""Check the corrected QfO proteomes and their explicitly pinned metadata."""

import json
from pathlib import Path

from benchmark_tools.run_qfo_order_replay import record

STAGING_BYTES = 27130
STAGING_SHA = "07a890eb816944f946d46039559a6046c2a3664f033eaa0f30f35bb25b9c9ab8"


def verify_inventory(directory, fastas, metadata):
    directory = Path(directory)
    expected = sorted(fastas, key=lambda r: r["path"])
    records = [*expected, metadata]
    paths = [Path(r["path"]) for r in records]
    if (len(set(paths)) != len(paths) or any(p.parent != directory for p in paths)
            or Path(metadata["path"]).name != "staging_manifest.json"):
        raise ValueError("Invalid planned input inventory")
    if set(directory.iterdir()) != set(paths) or any(not p.is_file() for p in paths):
        raise ValueError("Missing or unexpected input directory entry")
    if any(record(r["path"]) != r for r in records):
        raise ValueError("Changed input or metadata content")
    staged = json.loads(Path(metadata["path"]).read_text())
    if sorted(staged["input_fastas"], key=lambda r: r["path"]) != expected:
        raise ValueError("Staging metadata disagrees with planned FASTAs")
    return records


def verify_qfo_inputs(plan):
    directory = Path(plan["input_directory"])
    metadata = dict(path=str(directory / "staging_manifest.json"),
                    bytes=STAGING_BYTES, sha256=STAGING_SHA)
    return verify_inventory(directory, plan["fastas"], metadata)
