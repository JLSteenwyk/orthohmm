"""Prepare the corrected-release four candidate arms from independently admitted replay."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
ADMISSION_EXECUTOR = "09dec90118e9280990295f8aa9c1aa1a9171714e"


def validate_replay_admission(admission, plan_record, source_record):
    if (admission.get("status") != "corrected_checked_replay_admitted"
            or admission.get("accuracy_evaluated") is not False or admission.get("publication_ready") is not False
            or admission["source"] != source_record or admission["plan"] != plan_record):
        raise ValueError("Require independently admitted corrected replay")
    scheduler = admission["scheduler"]
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "32"):
        raise ValueError("Replay scheduler has not completed under frozen resources")
    stages = admission["coverage"]
    if [s["label"] for s in stages] != ["multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"]:
        raise ValueError("Require all four admitted replay partitions")
    checked = admission["checked_records"]
    if any(s["output"] not in checked for s in stages) or admission["source_report"] not in checked:
        raise ValueError("Replay outputs missing from checked provenance")
    if (len(admission["clustering"]) != 4
            or any(type(s["genes"]) is not int or s["genes"] != 984137 for s in admission["clustering"])):
        raise ValueError("Wrong corrected clustering universe")
    # Equality is an observation, never a prerequisite for retaining a replay.
    if type(admission["native_partition_comparison"]["partition_equal"]) is not bool:
        raise ValueError("Missing native/replay comparison")
    return [(False, stages[1]["output"]), (True, stages[3]["output"])]


def prepare(root, admission_path, admission_sha, admission_job, output):
    if output.exists():
        raise FileExistsError(output)
    environment = {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if any(os.environ.get(k) != v for k, v in environment.items()):
        raise ValueError("Require fixed hash seed and single-threaded numerical libraries")
    if not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "2":
        raise ValueError("Require scheduled 2-CPU preparation")
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.verify_ygob_validation import require_completed_job
    from benchmark_tools.checked_replay_payload_worker import corrected_evidence
    accounting = subprocess.check_output(["sacct", "-j", str(admission_job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, admission_job)
    admission = read_frozen(admission_path, admission_sha)
    admission_executor = root / "benchmarks/work/publication_qfo_corrected_replay_admission_v1"
    if subprocess.check_output(["git", "-C", str(admission_executor), "rev-parse", "HEAD"], text=True).strip() != ADMISSION_EXECUTOR:
        raise ValueError("Replay admission executor changed")
    subprocess.run(["git", "-C", str(admission_executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    plan, plan_record, _, _ = corrected_evidence(Path(admission["plan"]["path"]), admission["plan"]["sha256"])
    seeds = validate_replay_admission(admission, plan_record,
        record(admission_executor / "benchmark_tools/admit_qfo_corrected_replay.py"))
    for item in admission["checked_records"]:
        check(item)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    core = root / "benchmarks/work/publication_method_native_v2"
    sys.path.insert(0, str(launcher))
    from orthohmm.accuracy import load_accuracy_checkpoint
    from orthohmm.orthohmm import _expand_phylogeny_candidates
    sys.path.remove(str(launcher))
    if Path(sys.modules["orthohmm.orthohmm"].__file__).resolve() != launcher / "orthohmm/orthohmm.py":
        raise ValueError("Candidate engine imported outside frozen launcher")
    from Bio import SeqIO
    from benchmark_tools.prepare_orthobench_factorial import prepare_partition, plan_cells
    from benchmark_tools.prepare_qfo_factorial import verify_inputs
    from benchmark_tools.replay_phylogeny import load_membership_constraints
    from benchmark_tools.verify_qfo_replay_launcher import verify
    from benchmark_tools.audit_qfo_replay_inputs import verify_species_partition
    from benchmark_tools.audit_accuracy_checkpoint import audit
    parents = {Path(r["path"]).parent for r in plan["input_fastas"]}
    if len(parents) != 1:
        raise ValueError("Ambiguous corrected FASTA directory")
    fasta = next(iter(parents))
    fastas = [record(p) for p in sorted(fasta.glob("*.fasta"))]
    verify_inputs(fastas, plan["input_fastas"])
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    if runtime != plan["runtime"]:
        raise ValueError("Runtime differs from admitted corrected plan")
    checkpoint = Path(plan["checkpoint_manifest"]["path"]).parent
    numeric = audit(checkpoint, plan["checkpoint_manifest"]["sha256"])
    names, species, queries, targets, scores = load_accuracy_checkpoint(checkpoint, verify=False)
    if len(names) != 984137 or len(set(names)) != len(names):
        raise ValueError("Wrong corrected checkpoint universe")
    ownership = verify_species_partition(names, species,
        ((r["path"], (s.id for s in SeqIO.parse(r["path"], "fasta"))) for r in fastas))
    executor = Path(__file__).resolve().parent.parent
    helpers = [record(executor / "benchmark_tools" / name) for name in
        ("prepare_orthobench_factorial.py", "prepare_qfo_factorial.py", "replay_phylogeny.py",
         "audit_historical_profile_ablation.py", "audit_qfo_replay_inputs.py", "audit_accuracy_checkpoint.py",
         "verify_qfo_replay_launcher.py", "prepare_ob_candidate_neighborhood.py", "checked_replay_payload_worker.py")]
    cells = plan_cells(output, fasta, executor / "benchmark_tools/replay_phylogeny.py", 32)
    output.mkdir(parents=True, exist_ok=False)
    manifest = output / "manifest.json"
    report = {"status": "preparing", "accuracy_computed": False, "source": record(__file__),
        "job_id": os.environ["SLURM_JOB_ID"], "admission": record(admission_path), "admission_scheduler": scheduler,
        "plan": plan_record, "runtime_before": runtime, "numeric_checkpoint": numeric, "input_fastas": fastas,
        "species_code_mapping": ownership, "helpers": helpers, "core_root": str(core), "launcher_root": str(launcher),
        "candidate_arms": {}, "cells": cells, "environment_overrides": environment,
        "native_partition_comparison": admission["native_partition_comparison"],
        "limitations": ["Corrected-release factorial; original-release predictions and scores are not reused.",
            "Profile-off retains HMM initial search; profile-on includes downstream sequence refinement.",
            "Shared-host incremental preparation, not matched end-to-end timing.",
            "Candidate output admission, reconciliation and scoring are still required."]}
    try:
        for profile, seed in seeds:
            check(seed)
            for expand in (False, True):
                label = f"p{int(profile)}_c{int(expand)}"
                started = time.perf_counter()
                arm = prepare_partition(Path(seed["path"]), output / "candidates" / label, names, species,
                    (queries, targets, scores), expand, _expand_phylogeny_candidates, load_membership_constraints)
                arm["incremental_preparation_seconds"] = time.perf_counter() - started
                arm["output_files"] = [record(p) for p in sorted((output / "candidates" / label).rglob("*")) if p.is_file()]
                report["candidate_arms"][label] = arm
                manifest.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
            check(seed)
        verify_inputs([record(p) for p in sorted(fasta.glob("*.fasta"))], fastas)
        if audit(checkpoint, numeric["manifest"]["sha256"]) != numeric:
            raise ValueError("Checkpoint changed during preparation")
        report["runtime_after"] = verify(core, launcher, runtime_path)
        if report["runtime_after"] != runtime:
            raise ValueError("Frozen runtime changed")
        for item in helpers + [report["source"], report["admission"], plan_record]:
            check(item)
        for arm in report["candidate_arms"].values():
            for item in arm["output_files"]:
                check(item)
        report["status"] = "corrected_qfo_four_candidate_arms_prepared_unscored"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        manifest.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "admission", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--admission-job", required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.admission.resolve(), args.admission_sha256, args.admission_job, args.output.resolve())
