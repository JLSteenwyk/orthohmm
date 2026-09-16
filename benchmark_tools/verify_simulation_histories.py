"""Verify frozen paired simulation histories after successful generation."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record
from benchmark_tools.run_simulation_generation import child_path, verify_file


def history_inventory(native):
    records = {}
    excluded = {"T/SpeciesTreeParameters.tsv", "G/GenomeParameters.tsv"}
    for stage in ("T", "G"):
        files = [p for p in sorted((native / stage).rglob("*"))
                 if p.is_file() and p.relative_to(native).as_posix() not in excluded]
        if not files:
            raise ValueError(f"No biological history files in {native / stage}")
        for path in files:
            record = file_record(path, native)
            records[record["path"]] = record
    return records


def compare_histories(first, second):
    a, b = history_inventory(first), history_inventory(second)
    differences = [name for name in sorted(a.keys() | b.keys()) if a.get(name) != b.get(name)]
    return {"matched": not differences, "first_files": len(a), "second_files": len(b),
            "differences": differences, "first_inventory": a, "second_inventory": b}


def verify_generation(panel, run, manifest_hash):
    status_path = panel / "execution" / run["label"] / "status.json"
    status = json.loads(status_path.read_text())
    if status["label"] != run["label"] or status["status"] != "complete" or not status["output_inventory_recorded"]:
        raise ValueError(f"Generation not verified complete: {run['label']}")
    if status["provenance"]["manifest"]["sha256"] != manifest_hash:
        raise ValueError("Generation used a different manifest")
    expected = [(c["stage"], c["argv"]) for c in run["commands"]]
    actual = [(s["stage"], s["command"]) for s in status["stages"]]
    if actual != expected or any(s["status"] != "complete" or s["exit_code"] != 0 for s in status["stages"]):
        raise ValueError("Generation stage evidence mismatch")
    recorded = {}
    for record in status["outputs"]:
        if record["path"] in recorded:
            raise ValueError("Duplicate output inventory path")
        verify_file(child_path(panel, record["path"]), record)
        recorded[record["path"]] = record
    native = child_path(panel, run["native_output"])
    for relative, record in history_inventory(native).items():
        key = (native / relative).relative_to(panel.resolve()).as_posix()
        if key not in recorded or any(record[field] != recorded[key][field] for field in ("sha256", "bytes")):
            raise ValueError("History absent from generation output inventory")
    return native


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    raw = args.manifest.read_bytes()
    if hashlib.sha256(raw).hexdigest() != args.manifest_sha256:
        raise ValueError("Manifest checksum mismatch")
    manifest = json.loads(raw)
    runs = manifest["simulation_runs"]
    checks = manifest["history_equivalence_checks"]
    if len(runs) != 40 or len({r["label"] for r in runs}) != 40 or len(checks) != 20:
        raise ValueError("Incomplete frozen panel")
    native = {r["label"]: verify_generation(args.panel.resolve(), r, args.manifest_sha256) for r in runs}
    results = [{**check, **compare_histories(native[check["first"]], native[check["second"]])} for check in checks]
    report = {"schema_version": 1, "manifest_sha256": args.manifest_sha256,
              "verifier": file_record(Path(__file__).resolve(), Path(__file__).resolve().parent),
              "all_histories_matched": all(r["matched"] for r in results), "checks": results,
              "method_accuracy_evaluated": False}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        handle.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not report["all_histories_matched"]:
        raise SystemExit("Paired histories differ; inspect preserved mismatch report")


if __name__ == "__main__":
    main()
