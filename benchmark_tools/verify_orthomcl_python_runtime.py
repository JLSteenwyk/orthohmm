"""Enforce the pinned dedicated BPO Python runtime before and after execution."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.freeze_orthomcl_python_runtime import inspect
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.snapshot_runtime_trees import verify

RUNTIME_SHA = "3519dbf43a376c935cc76ab00559af045c3e94ce47bf421319b970e660c8cddc"


def compare(expected, observed):
    for key in ("python_version", "prefix", "base_prefix", "executable", "packages"):
        if observed[key] != expected["inspection"][key]:
            raise ValueError("Changed dedicated Python identity: " + key)
    files = {r["path"]: r for r in expected["runtime"]["records"] if r["kind"] == "file"}
    for mapped in observed["mapped_files"]:
        old = files.get(mapped["path"])
        if old is None or any(old[key] != mapped[key] for key in ("bytes", "sha256")):
            raise ValueError("Unbound or changed mapped Python runtime file: " + mapped["path"])


def verify_runtime(root):
    path = root / "benchmark_tools/results/orthomcl_python_runtime_20260918.json"
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != RUNTIME_SHA:
        raise ValueError("Changed dedicated Python runtime manifest")
    expected = json.loads(raw)
    if expected["status"] != "dedicated_bpo_python_runtime_inventoried":
        raise ValueError("Unexpected Python runtime snapshot status")
    observed = inspect()
    compare(expected, observed)
    checked = verify(expected["runtime"])
    return {"status": "dedicated_bpo_python_runtime_verified", "manifest": record(path),
            "runtime_records": checked["records"], "python_version": observed["python_version"],
            "packages": observed["packages"], "mapped_files": observed["mapped_files"],
            "limitations": ["Explicit runtime identity, not an OS-wide hermetic execution claim.",
                            "Helper source identity is bound separately by the frozen executor."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(verify_runtime(args.root.resolve()), sort_keys=True))
