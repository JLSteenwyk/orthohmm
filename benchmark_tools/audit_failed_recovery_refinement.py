"""Read-only diagnosis of admission 22155; never authorizes recovered results."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess


def record(path):
    path = Path(path)
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            size += len(block)
            digest.update(block)
    return dict(path=str(path), bytes=size, sha256=digest.hexdigest())


def coverage(path, names, expected_groups):
    universe = set(names)
    if len(universe) != len(names) or not universe:
        raise ValueError("Invalid gene universe")
    seen = set()
    groups = 0
    with Path(path).open() as stream:
        for line in stream:
            genes = line.split()
            if not genes:
                continue
            members = set(genes)
            if len(members) != len(genes) or members & seen or not members <= universe:
                raise ValueError("Invalid partition membership")
            seen.update(members)
            groups += 1
    if seen != universe or groups != expected_groups:
        raise ValueError("Incomplete partition or wrong group count")
    return dict(genes=len(seen), groups=groups)


def audit(root):
    directory = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_admission_v1"
    native = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
    accounting = subprocess.check_output(["sacct", "-j", "22155", "--parsable2",
        "--format=JobID,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)
    rows = [row for row in csv.DictReader(io.StringIO(accounting), delimiter="|") if row["JobID"] == "22155"]
    if len(rows) != 1 or tuple(rows[0][key] for key in ("State", "ExitCode", "AllocCPUS", "ReqMem", "NodeList")) != (
            "FAILED", "1:0", "2", "64G", "bizon"):
        raise ValueError("Unexpected admission scheduler state")
    files = [directory / "status.json", directory / "refinement.log",
        native / "status.json", native / "refinement.json", native / "refinement_repeat.json",
        native / "payload/gene_names.txt", native / "orthogroups_profiles_refined.txt",
        native / "refinement_repeat.txt", directory / "refinement/refinement_repeat.txt",
        root / "benchmarks/work/qfo_cpm_checkpoint_admit_22155.log", Path(__file__).resolve()]
    records = [record(path) for path in files]
    failed = json.loads(files[0].read_text())
    if (failed["status"] != "recovery_admission_failed" or failed["error_type"] != "CalledProcessError"
            or "SIGSEGV" not in failed["error"]
            or any(failed[key] is not False for key in ("seed_admitted", "accuracy_evaluated", "publication_ready"))):
        raise ValueError("Unexpected failure report")
    if (directory / "refinement/refinement_repeat.json").exists():
        raise ValueError("Unexpected completed child report")
    log = files[1].read_text()
    if not all(token in log for token in ("Segmentation fault", "Garbage-collecting", "read_partition", "refinement_worker")):
        raise ValueError("Unexpected failure stack")
    children = [json.loads(path.read_text()) for path in files[3:5]]
    if any(child["genes"] != 984137 or child["groups"] != 390845 for child in children):
        raise ValueError("Unexpected native partition metadata")
    if any(tuple(item[key] for key in ("bytes", "sha256")) !=
           tuple(records[6][key] for key in ("bytes", "sha256")) for item in records[7:9]):
        raise ValueError("Failed-child output differs from native refinements")
    names = files[5].read_text().splitlines()
    if len(names) != 984137:
        raise ValueError("Wrong saved gene universe")
    observed = coverage(files[8], names, 390845)
    if [record(path) for path in files] != records:
        raise ValueError("Diagnostic inputs changed")
    return dict(status="failed_refinement_output_readback_verified_not_admitted", scheduler=rows[0],
        accounting=accounting, coverage=observed, byte_identical_to_two_native_refinements=True,
        child_metadata_missing=True, checked_records=records, source=records[-1],
        seed_admitted=False, retry_authorized=False, accuracy_evaluated=False, publication_ready=False,
        limitations=["The child died during post-write readback; matching bytes do not convert its failure to success.",
                     "Garbage collection is the observed stack location, not an established root cause.",
                     "No optimizer or refinement was rerun; no scientific native libraries are imported.",
                     "Independent recovery admission remains failed and candidate job 22156 must remain blocked."])


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
