"""Admit corrected candidate preparation before reconciliation or conversion."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_candidate_arm import audit as audit_arm
from benchmark_tools.checked_replay_payload_worker import corrected_evidence
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_orthobench_factorial import plan_cells
from benchmark_tools.prepare_qfo_corrected_factorial import validate_replay_admission, ADMISSION_EXECUTOR
from benchmark_tools.prepare_qfo_factorial import verify_inputs
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_qfo_replay_launcher import verify
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR = "0a028a743fca7626b86376475f4c7fd438093717"
ARMS = ("p0_c0", "p0_c1", "p1_c0", "p1_c1")
HELPERS = ("prepare_orthobench_factorial.py", "prepare_qfo_factorial.py", "replay_phylogeny.py",
           "audit_historical_profile_ablation.py", "audit_qfo_replay_inputs.py", "audit_accuracy_checkpoint.py",
           "verify_qfo_replay_launcher.py", "prepare_ob_candidate_neighborhood.py", "checked_replay_payload_worker.py",
           "audit_candidate_arm.py")
PARAMETERS = {"max_component_genes": 500, "max_satellite_genes": 12, "max_satellite_to_anchor_ratio": .75,
              "min_margin": 1.5, "iteration_margin_increment": .5, "max_satellites_per_anchor": 4,
              "max_species_overlap_fraction": 1., "min_avg_score": 0., "min_max_score": 0.,
              "min_coverage": .5, "min_norm": .03, "max_iterations": 2}


def validate_manifest(manifest, scheduler, root, executor, fasta):
    output = root / "benchmarks/results/qfo_corrected_factorial_v1"
    if (manifest["status"] != "corrected_qfo_four_candidate_arms_prepared_unscored"
            or manifest["accuracy_computed"] is not False or set(manifest["candidate_arms"]) != set(ARMS)):
        raise ValueError("Incomplete corrected candidate preparation")
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "2"
            or manifest["job_id"] != scheduler["JobIDRaw"]):
        raise ValueError("Preparation job identity/allocation differs")
    if (manifest["cells"] != plan_cells(output, fasta, executor / "benchmark_tools/replay_phylogeny.py", 32)
            or manifest["core_root"] != str(root / "benchmarks/work/publication_method_native_v2")
            or manifest["launcher_root"] != str(root / "benchmarks/work/publication_qfo_replay_native_v1")):
        raise ValueError("Candidate cell plan or frozen roots differ")
    expected_env = {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if manifest["environment_overrides"] != expected_env:
        raise ValueError("Preparation environment differs")
    for label, arm in manifest["candidate_arms"].items():
        working = output / "candidates" / label / "orthohmm_working_res"
        expanded = label.endswith("c1")
        if (arm["candidate_expansion"] is not expanded
                or arm["candidate_partition"]["path"] != str(working / "orthohmm_edges_clustered.txt")
                or expanded != ("membership_constraints" in arm)):
            raise ValueError("Candidate arm path or membership policy differs")
        if expanded and arm["expansion"]["parameters"] != PARAMETERS:
            raise ValueError("Candidate expansion settings differ from frozen satellite_v2")
    return output


def admit(root, manifest_path, manifest_sha, job, destination):
    if destination.exists():
        raise FileExistsError(destination)
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, job)
    expected_path = root / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
    if manifest_path != expected_path:
        raise ValueError("Unexpected corrected candidate manifest path")
    manifest = read_frozen(manifest_path, manifest_sha)
    executor = root / "benchmarks/work/publication_qfo_corrected_candidates_v1"
    admission_executor = root / "benchmarks/work/publication_qfo_corrected_replay_admission_v1"
    for directory, revision in ((executor, EXECUTOR), (admission_executor, ADMISSION_EXECUTOR)):
        if subprocess.check_output(["git", "-C", str(directory), "rev-parse", "HEAD"], text=True).strip() != revision:
            raise ValueError("Frozen executor revision changed")
        subprocess.run(["git", "-C", str(directory), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    source = record(executor / "benchmark_tools/prepare_qfo_corrected_factorial.py")
    helpers = [record(executor / "benchmark_tools" / name) for name in HELPERS]
    if manifest["source"] != source or manifest["helpers"] != helpers:
        raise ValueError("Candidate producer/helper inventory differs")
    check(manifest["admission"])
    admission = json.loads(Path(manifest["admission"]["path"]).read_text())
    plan, plan_record, native_admission_record, names_record = corrected_evidence(
        Path(manifest["plan"]["path"]), manifest["plan"]["sha256"])
    seeds = dict(validate_replay_admission(admission, plan_record,
        record(admission_executor / "benchmark_tools/admit_qfo_corrected_replay.py")))
    if manifest["plan"] != plan_record or manifest["native_partition_comparison"] != admission["native_partition_comparison"]:
        raise ValueError("Prepared plan or native/replay comparison differs")
    parents = {Path(item["path"]).parent for item in plan["input_fastas"]}
    if len(parents) != 1:
        raise ValueError("Ambiguous corrected FASTA directory")
    fasta = next(iter(parents))
    verify_inputs(manifest["input_fastas"], plan["input_fastas"])
    verify_inputs([record(p) for p in sorted(fasta.glob("*.fasta"))], manifest["input_fastas"])
    output = validate_manifest(manifest, scheduler, root, executor, fasta)
    prior = manifest["admission_scheduler"]
    observed_accounting = subprocess.check_output(["sacct", "-j", prior["JobIDRaw"], "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    if require_completed_job(observed_accounting, prior["JobIDRaw"]) != prior:
        raise ValueError("Recorded replay-admission accounting differs")
    native = json.loads(Path(native_admission_record["path"]).read_text())["content"]
    if manifest["species_code_mapping"] != native["species_ownership"]:
        raise ValueError("Candidate species ownership differs from native admission")
    if any(manifest["numeric_checkpoint"][k] != native["numeric_checkpoint"][k] for k in ("status", "summary", "manifest")):
        raise ValueError("Candidate numeric checkpoint differs")
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(Path(manifest["core_root"]), Path(manifest["launcher_root"]), runtime_path)
    if not runtime == plan["runtime"] == manifest["runtime_before"] == manifest["runtime_after"]:
        raise ValueError("Candidate runtime observations differ")
    names = Path(names_record["path"]).read_text().splitlines()
    universe = set(names)
    if len(names) != 984137 or len(universe) != len(names):
        raise ValueError("Wrong corrected gene universe")
    summaries = {}
    records = [record(manifest_path), source, *helpers, manifest["admission"], plan_record,
               native_admission_record, names_record, *admission["checked_records"], *plan["checked_records"],
               record(__file__), record(Path(__file__).with_name("audit_candidate_arm.py"))]
    for label in ARMS:
        arm = manifest["candidate_arms"][label]
        summary = audit_arm(arm, seeds[label.startswith("p1")], output / "candidates" / label, universe, label.endswith("c1"))
        if summary != arm["content_audit"]:
            raise ValueError("Independent candidate content check differs from producer")
        records.extend(summary["checked_records"])
        summaries[label] = summary
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting candidate provenance records")
        unique[item["path"]] = item
    for item in unique.values():
        check(item)
    if verify(Path(manifest["core_root"]), Path(manifest["launcher_root"]), runtime_path) != runtime:
        raise ValueError("Runtime changed during candidate admission")
    result = {"status": "corrected_qfo_candidates_admitted", "accuracy_evaluated": False, "publication_ready": False,
              "source": record(__file__), "prepared_manifest": record(manifest_path), "scheduler": scheduler,
              "accounting": accounting, "candidate_arms": summaries, "cells": manifest["cells"],
              "checked_records": list(unique.values()),
              "limitations": ["Recorded candidate consistency and provenance; search-support values are not independently recomputed.",
                  "Incremental shared-host preparation is not dedicated end-to-end timing.",
                  "Native reconciliation, conversion and benchmark assessment require separate admission."]}
    with destination.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "manifest", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--job", required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.manifest.resolve(), args.manifest_sha256, args.job, args.output.resolve())
