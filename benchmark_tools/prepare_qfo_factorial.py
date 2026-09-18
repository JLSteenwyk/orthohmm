"""Prepare frozen QfO candidate arms from the admitted recovered replay."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def select_seeds(stages):
    labels = [row["label"] for row in stages]
    if labels != ["multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"]:
        raise ValueError("Require all four admitted stages in order")
    return [(False, stages[1]["output"]), (True, stages[3]["output"])]


def verify_inputs(records, audited):
    if len(records) != 78 or len({r["path"] for r in records}) != 78:
        raise ValueError("Require 78 unique QfO FASTAs")
    if sorted(records, key=lambda r: r["path"]) != sorted(audited, key=lambda r: r["path"]):
        raise ValueError("FASTA inventory differs from admitted replay")


def prepare(root, output):
    if output.exists():
        raise FileExistsError(output)
    expected_env = {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if any(os.environ.get(k) != v for k, v in expected_env.items()):
        raise ValueError("Require fixed hash seed and single-threaded numerical libraries")
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    core = root / "benchmarks/work/publication_method_native_v2"
    # Load the frozen package before helper imports that transitively use it.
    sys.path.insert(0, str(launcher))
    from orthohmm.accuracy import load_accuracy_checkpoint
    from orthohmm.orthohmm import _expand_phylogeny_candidates
    sys.path.remove(str(launcher))
    if Path(sys.modules["orthohmm.orthohmm"].__file__).resolve() != launcher / "orthohmm/orthohmm.py":
        raise ValueError("Candidate engine imported outside frozen launcher")
    from Bio import SeqIO
    from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
    from benchmark_tools.prepare_qfo_recovered_pairs import ADMISSION_SHA, stage_outputs
    from benchmark_tools.prepare_orthobench_factorial import prepare_partition, plan_cells
    from benchmark_tools.replay_phylogeny import load_membership_constraints
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.verify_qfo_replay_launcher import verify
    from benchmark_tools.verify_ygob_validation import require_completed_job
    from benchmark_tools.audit_qfo_replay_inputs import verify_species_partition
    from benchmark_tools.audit_accuracy_checkpoint import audit

    admission_path = root / "benchmark_tools/results/qfo_checked_full_replay_recovered_verified_20260917.json"
    admission = read_frozen(admission_path, ADMISSION_SHA)
    stages = stage_outputs(admission)
    seeds = select_seeds(stages)
    accounting = subprocess.check_output(["sacct", "-j", "21480", "--parsable2",
                                         "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21480)
    for item in admission["provenance_checked"]:
        check(item)
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    input_link = admission["recovery"]["input_recheck"]
    check(input_link)
    inputs = json.loads(Path(input_link["path"]).read_text())
    fasta = root / "qfo_benchmark/input"
    fastas = [record(p) for p in sorted(fasta.glob("*.fasta"))]
    verify_inputs(fastas, inputs["input_fastas"])
    checkpoint = root / "qfo_benchmark/results/orthohmm_high_sensitivity_isolated/output/orthohmm_working_res/high_sensitivity_checkpoint"
    numeric = audit(checkpoint, inputs["numeric_checkpoint"]["manifest"]["sha256"])
    names, species, queries, targets, scores = load_accuracy_checkpoint(checkpoint, verify=False)
    if len(names) != 976504:
        raise ValueError("Unexpected QfO gene universe")
    ownership = verify_species_partition(names, species,
        ((r["path"], (s.id for s in SeqIO.parse(r["path"], "fasta"))) for r in fastas))
    executor = Path(__file__).resolve().parent.parent
    helper_names = ("prepare_orthobench_factorial.py", "replay_phylogeny.py", "audit_historical_profile_ablation.py",
                    "prepare_qfo_recovered_pairs.py", "audit_qfo_replay_inputs.py", "audit_accuracy_checkpoint.py",
                    "verify_qfo_replay_launcher.py", "prepare_ob_candidate_neighborhood.py")
    helpers = [record(executor / "benchmark_tools" / name) for name in helper_names]
    cells = plan_cells(output, fasta, executor / "benchmark_tools/replay_phylogeny.py", 32)
    output.mkdir(parents=True)
    manifest = output / "manifest.json"
    report = {"status": "preparing", "accuracy_computed": False, "source": record(__file__),
              "job_id": os.environ.get("SLURM_JOB_ID"), "admission": record(admission_path),
              "admission_scheduler": scheduler, "runtime_before": runtime, "numeric_checkpoint": numeric,
              "input_fastas": fastas, "species_code_mapping": ownership, "helpers": helpers,
              "core_root": str(core), "launcher_root": str(launcher), "candidate_arms": {}, "cells": cells,
              "environment_overrides": {k: os.environ.get(k) for k in ("PYTHONHASHSEED", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
              "limitations": ["New factorial on admitted recovered checkpoints; not equivalent to historical comparison rows.",
                "Profile-off retains initial HMM search. Profile branch includes its downstream sequence refinement.",
                "Preparation cost is incremental on a shared workstation, not controlled end-to-end timing.",
                "Reconciliation, output admission and all accuracy scoring remain separate pending steps."]}
    try:
        for profile, seed in seeds:
            check(seed)
            for expand in (False, True):
                label = f"p{int(profile)}_c{int(expand)}"
                started = time.perf_counter()
                arm = prepare_partition(Path(seed["path"]), output / "candidates" / label, names, species,
                    (queries, targets, scores), expand, _expand_phylogeny_candidates, load_membership_constraints)
                arm["incremental_preparation_seconds"] = time.perf_counter() - started
                report["candidate_arms"][label] = arm
                manifest.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
            check(seed)
        verify_inputs([record(p) for p in sorted(fasta.glob("*.fasta"))], fastas)
        after = audit(checkpoint, numeric["manifest"]["sha256"])
        if after != numeric:
            raise ValueError("Numeric checkpoint changed during preparation")
        report["runtime_after"] = verify(core, launcher, runtime_path)
        if report["runtime_after"] != runtime:
            raise ValueError("Frozen runtime changed")
        for item in helpers + [report["source"], report["admission"]]:
            check(item)
        report["status"] = "qfo_four_candidate_arms_prepared_unscored"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        manifest.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output.resolve())
