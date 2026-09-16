"""Label-blind file and completion checks for the frozen YGOB inference run."""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record
from benchmark_tools.build_publication_runtime import verify_runtime


def require_completed_job(accounting, job_id):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    matching = [r for r in rows if r.get("JobIDRaw") == str(job_id)]
    if len(matching) != 1 or matching[0].get("State") != "COMPLETED" or matching[0].get("ExitCode") != "0:0":
        raise ValueError("Requested job is not uniquely confirmed COMPLETED with exit 0:0")
    return matching[0]


def read_metadata(path):
    result = {}
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) != 2 or not row[0] or row[0] in result:
                raise ValueError("Malformed or duplicate runner metadata")
            result[row[0]] = row[1]
    return result


def verify_records(records, base):
    """Verify all recorded bytes/hashes; reject duplicate and escaping paths."""
    if not records:
        raise ValueError("Empty artifact manifest")
    base = base.resolve()
    observed = set()
    for record in records:
        path = (base / record["path"]).resolve()
        if not path.is_relative_to(base) or path in observed:
            raise ValueError("Duplicate or escaping artifact path")
        observed.add(path)
        actual = file_record(path, base)
        if any(actual[k] != record[k] for k in ("bytes", "sha256")):
            raise ValueError(f"Artifact changed: {path}")
    return observed


def validate_runtime_probe(probe, frozen_root, binary_hash, source_hash):
    root = frozen_root.resolve()
    if probe.get("status") != "passed" or probe.get("exit_code") != 0 or probe.get("root") != str(root):
        raise ValueError("Missing successful exact-checkout runtime probe")
    library = probe.get("pair_align", {})
    if library.get("sha256") != binary_hash or library.get("path") != str(root / "orthohmm/search/csrc/pair_align.so"):
        raise ValueError("Profile probe used a different native library")
    if probe.get("profile_source_sha256") != source_hash or probe.get("profile_source") != str(root / "orthohmm/search/msa_profile.py"):
        raise ValueError("Profile probe used a different source")
    if probe.get("profile_length", 0) < 1:
        raise ValueError("Profile probe produced no profile")


def verify_native_evidence(out, frozen_root, expected_manifest):
    snapshot = out / "native_runtime.json"
    if snapshot.read_bytes() != expected_manifest.read_bytes():
        raise ValueError("YGOB native runtime differs from the frozen build manifest")
    build = verify_runtime(snapshot, frozen_root)
    binary = next(r for r in build["binaries"] if Path(r["path"]).name == "pair_align.so")
    source = file_record(frozen_root / "orthohmm/search/msa_profile.py", frozen_root)
    probes = []
    for method in ("high_sensitivity", "satellite_v2"):
        for stage in ("before", "after"):
            path = out / f"{method}_{stage}.runtime.json"
            validate_runtime_probe(json.loads(path.read_text()), frozen_root, binary["sha256"], source["sha256"])
            probes.append(file_record(path, out))
    return {"manifest": file_record(snapshot, out), "probes": probes}


def verify_run(root, frozen_root, job_id, accounting):
    scheduler = require_completed_job(accounting, job_id)
    out = root / "benchmarks/results/ygob_validation_v1"
    prepared = root / "benchmarks/work/ygob_validation_v1"
    metadata = read_metadata(out / "run_metadata.tsv")
    if metadata.get("job_id") != str(job_id) or metadata.get("runner_exit_code") != "0" or not metadata.get("finished"):
        raise ValueError("Runner metadata does not confirm successful requested job")
    native_runtime = verify_native_evidence(out, frozen_root,
                        root / "benchmark_tools/results/publication_native_runtime_20260916.json")
    source_commit = subprocess.check_output(["git", "-C", str(frozen_root), "rev-parse", "HEAD"], text=True).strip()
    if not source_commit.startswith("7f3a9e4") or (out / "inference_source_commit.txt").read_text().strip() != source_commit:
        raise ValueError("Inference source is not the frozen revision")
    subprocess.run(["git", "-C", str(frozen_root), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools/benchmark_production.py"], check=True, capture_output=True)
    for name in ("launcher_inputs.sha256", "orthofinder_entrypoint.sha256", "tool_entrypoints.sha256"):
        subprocess.run(["sha256sum", "--check", "--status", str(out / name)], check=True)
    manifest = json.loads((prepared / "manifest.json").read_text())
    frozen = json.loads((root / "benchmark_tools/results/ygob_validation_inputs_20260916.json").read_text())
    for key in ("inputs", "reference", "proteins", "species", "reference_genes", "reference_groups"):
        if manifest[key] != frozen[key]:
            raise ValueError(f"Prepared manifest differs from frozen specification: {key}")
    if (manifest["proteins"], len(manifest["species"]), manifest["reference_genes"], manifest["reference_groups"]) != (83404, 16, 83391, 10250):
        raise ValueError("Frozen input dimensions changed")
    input_paths = verify_records(manifest["inputs"], prepared)
    verify_records([manifest["reference"]], prepared)
    expected_inputs = {p.name: file_record(p, p.parent) for p in input_paths}
    verified_methods = {}
    for method in ("high_sensitivity", "satellite_v2"):
        path = out / f"{method}.json"
        metrics = json.loads(path.read_text())
        harness = metrics["harness"]
        if metrics.get("status") != "complete" or harness.get("exit_code") != 0 or harness.get("git_commit") != source_commit:
            raise ValueError(f"Incomplete or wrong-source OrthoHMM run: {method}")
        records = harness["input_manifest"]
        verify_records(records, prepared / "input")
        if {r["path"]: r for r in records} != expected_inputs:
            raise ValueError("OrthoHMM input manifest mismatch")
        sources = verify_records(harness["source_manifest"], frozen_root)
        expected_sources = set((frozen_root / "orthohmm").rglob("*.py")) | {frozen_root / "benchmark_tools/benchmark_production.py"}
        if sources != {p.resolve() for p in expected_sources}:
            raise ValueError("Incomplete OrthoHMM source manifest")
        outputs = verify_records(harness["output_manifest"], out / method)
        native = out / method / ("orthohmm_orthogroups.txt" if method == "high_sensitivity"
                                else "orthohmm_phylogeny/orthohmm_root_hogs.tsv")
        if native.resolve() not in outputs:
            raise ValueError("Required native output missing from verified manifest")
        verified_methods[method] = {"metrics": file_record(path, out), "native_output": str(native)}
    copies = list((out / "orthofinder/input").glob("*.fasta"))
    if {p.name: file_record(p, p.parent) for p in copies} != expected_inputs:
        raise ValueError("OrthoFinder input copies differ")
    time_log = (out / "orthofinder.time.log").read_text()
    if [line.strip() for line in time_log.splitlines() if "Exit status:" in line] != ["Exit status: 0"]:
        raise ValueError("OrthoFinder did not finish successfully")
    full = list((out / "orthofinder").glob("**/Orthogroups/Orthogroups.txt"))
    if len(full) != 1:
        raise ValueError("Expected exactly one finalized OrthoFinder orthogroup file")
    return {"schema_version": 1, "inference_files_verified": True, "accuracy_computed": False,
            "all_scoring_gates_verified": False, "scheduler": scheduler, "runner_metadata": metadata,
            "native_runtime": native_runtime,
            "source_commit": source_commit, "orthohmm": verified_methods,
            "orthofinder_final_groups": file_record(full[0], out),
            "remaining_checks": ["Exact inference commands and installed tool versions",
                                 "Overlap/reference-resource audits and native ID conversion",
                                 "Independent reference reconstruction and scoring"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--frozen-root", type=Path, required=True)
    parser.add_argument("--job-id", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Refusing to overwrite verification evidence")
    accounting = subprocess.check_output(["sacct", "-j", str(args.job_id), "--parsable2",
                                          "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    report = verify_run(args.root.resolve(), args.frozen_root.resolve(), args.job_id, accounting)
    report["accounting_raw"] = accounting
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
