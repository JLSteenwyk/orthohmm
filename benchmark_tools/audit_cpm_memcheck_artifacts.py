"""Read-only post-failure artifact checks; never admits the failed diagnostic."""

import argparse
import json
from pathlib import Path

from benchmark_tools.audit_failed_recovery_refinement import coverage, record

STATUS_SHA = "33d543f8d47b1df8175030aa209a1ea19406545f37401ba1b7c3f467ebea3d3a"


def collect(value):
    if isinstance(value, dict):
        if set(value) == {"path", "bytes", "sha256"}:
            yield value
        else:
            for child in value.values():
                yield from collect(child)
    elif isinstance(value, list):
        for child in value:
            yield from collect(child)


def inspect(records):
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting retained identities")
        unique[item["path"]] = item
    mismatches = []
    for path, expected in unique.items():
        try:
            actual = record(path)
        except OSError as error:
            mismatches.append(dict(expected=expected, error=str(error)))
            continue
        if actual != expected:
            mismatches.append(dict(expected=expected, actual=actual))
    return dict(unique_files=len(unique), matching_files=len(unique) - len(mismatches), mismatches=mismatches)


def audit(root):
    directory = root / "benchmarks/results/qfo_cpm_refinement_memcheck_diagnostic_v1"
    native = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
    status_record = record(directory / "status.json")
    if status_record["sha256"] != STATUS_SHA:
        raise ValueError("Changed failed diagnostic status")
    failed = json.loads(Path(status_record["path"]).read_text())
    if failed["status"] != "memcheck_diagnostic_failed" or failed["returncode"] != 97:
        raise ValueError("Unexpected diagnostic outcome")
    retained = [status_record, *failed["checked_records"], *collect(failed["runtime_before"]),
                failed["child_log"], failed["memcheck_xml"]]
    before = inspect(retained)
    if before["mismatches"]:
        raise ValueError("Retained diagnostic inputs changed: " + json.dumps(before["mismatches"]))
    child_path, reference_path = directory / "refinement_repeat.json", native / "refinement.json"
    child_identity, reference_identity = record(child_path), record(reference_path)
    child, reference = json.loads(child_path.read_text()), json.loads(reference_path.read_text())
    if {k: v for k, v in child.items() if k != "output"} != {
            k: v for k, v in reference.items() if k != "output"}:
        raise ValueError("Refinement metadata differs")
    output = record(directory / "refinement_repeat.txt")
    if child["output"] != output or any(output[k] != reference["output"][k] for k in ("bytes", "sha256")):
        raise ValueError("Refinement output differs")
    names_path = native / "payload/gene_names.txt"
    names_identity = record(names_path)
    checked = [*retained, child_identity, reference_identity, output, names_identity,
               record(Path(__file__).resolve())]
    result = coverage(Path(output["path"]), names_path.read_text().splitlines(), reference["groups"])
    after = inspect(checked)
    if after["mismatches"]:
        raise ValueError("Artifacts changed during readback")
    return dict(status="failed_memcheck_artifacts_verified_not_admitted", original_returncode=97,
                original_record_entries=len(failed["checked_records"]), before=before, after=after,
                coverage=result, metadata_equal_except_output_path=True, checked_records=checked,
                seed_admitted=False, accuracy_evaluated=False, publication_ready=False,
                limitations=["Byte checks and partition coverage do not establish safe native execution.",
                    "Runtime file identities checked; no runtime probe or inference was rerun.",
                    "Original nonzero diagnostic and failed scientific admission remain failed."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = audit(args.root.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
