"""Prepare four frozen corrected-QfO candidate variants after a byte-equal control."""

import argparse
import importlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import ARMS, check, controlled_expansion, record
from benchmark_tools.run_simulation_methods import read_frozen

PLAN_SHA = "4f25aff958e327a3d133f08326260c600e51e87e1d56c81f114db9993930479b"
PROTOCOL_SHA = "7b66c2f29b0098f26a68cb9f230bd6260e04089db9d40c98f4ce0b9af1924044"
ENVIRONMENT = {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}


def require_environment(environment):
    if (not environment.get("SLURM_JOB_ID") or environment.get("SLURM_CPUS_PER_TASK") != "2"
            or environment.get("SLURM_JOB_NODELIST") != "bizon"
            or any(environment.get(k) != v for k, v in ENVIRONMENT.items())):
        raise ValueError("Require scheduled 2-CPU bizon job and frozen numerical environment")


def prepare_arms(engine, baseline, names, species, hits, output, report, auditor):
    universe = set(names)
    for label, _ in ARMS:
        directory = output / label
        working = directory / "orthohmm_working_res"
        working.mkdir(parents=True, exist_ok=False)
        partition = working / "orthohmm_edges_clustered.txt"
        check(baseline["seed_partition"])
        shutil.copyfile(baseline["seed_partition"]["path"], partition)
        row = {"label": label, "status": "preparing"}
        report["arms"].append(row)
        (output / "progress.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        started = time.monotonic()
        row.update(controlled_expansion(engine, baseline["expansion"]["parameters"], label,
                                        (str(directory), names, species, hits)))
        row["incremental_seconds"] = time.monotonic() - started
        arm = {"seed_partition": baseline["seed_partition"], "candidate_expansion": True,
               "candidate_partition": record(partition),
               "membership_constraints": record(working / "phylogeny_candidate_merges.json"),
               "expansion": row["engine_fixed_profile_report"],
               "output_files": [record(p) for p in sorted(directory.rglob("*")) if p.is_file()]}
        arm["content_audit"] = auditor(arm, baseline["seed_partition"], directory, universe, True)
        row["candidate_arm"] = arm
        if label == "control":
            for key in ("candidate_partition", "membership_constraints"):
                if any(arm[key][k] != baseline[key][k] for k in ("sha256", "bytes")):
                    raise ValueError("Unchanged corrected control differs: " + key)
            row["baseline_byte_equivalent"] = True
        row["status"] = "candidate_prepared_unscored"
        (output / "progress.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def prepare(root, output):
    if output.exists():
        raise FileExistsError(output)
    require_environment(os.environ)
    results = root / "benchmark_tools/results"
    plan_path = results / "qfo_parameter_neighborhood_plan_20260919.json"
    protocol_path = results / "QFO_PARAMETER_NEIGHBORHOOD_PROTOCOL_20260919.md"
    plan = read_frozen(plan_path, PLAN_SHA)
    if record(protocol_path)["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed neighborhood protocol")
    inputs = [record(plan_path), record(protocol_path)]
    evidence = []
    for item in plan["inputs"]:
        path = root / item["path"]
        if record(path)["sha256"] != item["sha256"]:
            raise ValueError("Changed baseline evidence: " + item["path"])
        inputs.append(record(path))
        if path.suffix == ".json":
            evidence.append(json.loads(path.read_text()))
    replay, manifest, native, _score = evidence
    if (replay["status"] != "corrected_checked_replay_admitted"
            or replay["native_partition_comparison"]["partition_equal"] is not True
            or native["cell"] != "p1_c1_r1"):
        raise ValueError("Wrong corrected baseline")
    from benchmark_tools.verify_qfo_replay_launcher import verify
    core, launcher = Path(manifest["core_root"]), Path(manifest["launcher_root"])
    runtime_path = results / "publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    if not runtime == manifest["runtime_before"] == manifest["runtime_after"]:
        raise ValueError("Changed scientific runtime")
    sys.path.insert(0, str(launcher))
    try:
        engine = importlib.import_module("orthohmm.orthohmm")
        accuracy = importlib.import_module("orthohmm.accuracy")
    finally:
        sys.path.remove(str(launcher))
    for module in (engine, accuracy):
        if Path(module.__file__).resolve().parent != launcher / "orthohmm":
            raise ValueError("Scientific module imported outside frozen launcher")
    from benchmark_tools.run_qfo_corrected_factorial_cell import verify_admission
    from benchmark_tools.audit_accuracy_checkpoint import audit
    from benchmark_tools.audit_candidate_arm import audit as audit_arm
    admission_record = native["candidate_admission"]
    admission_path = Path(admission_record["path"])
    check(admission_record)
    _, verified, _, _, _, _ = verify_admission(root, admission_path, admission_record["sha256"], "21758", 3)
    if verified != manifest:
        raise ValueError("Candidate admission refers to a different manifest")
    check(manifest["plan"])
    replay_plan = json.loads(Path(manifest["plan"]["path"]).read_text())
    checkpoint_record = replay_plan["checkpoint_manifest"]
    checkpoint = Path(checkpoint_record["path"]).parent
    numeric = audit(checkpoint, checkpoint_record["sha256"])
    if (numeric["summary"] != manifest["numeric_checkpoint"]["summary"]
            or numeric["summary"]["species"] != 78):
        raise ValueError("Numeric checkpoint differs from admitted candidate baseline")
    names, species, queries, targets, scores = accuracy.load_accuracy_checkpoint(checkpoint, verify=False)
    if len(names) != 984137 or len(set(names)) != len(names):
        raise ValueError("Wrong corrected gene universe")
    baseline = manifest["candidate_arms"]["p1_c1"]
    inputs.extend([admission_record, manifest["plan"], record(runtime_path), baseline["seed_partition"]])
    # Record every imported local analysis module, not only this entry point.
    helpers = [record(module.__file__) for name, module in sorted(sys.modules.items())
               if name.startswith("benchmark_tools.") and getattr(module, "__file__", None)]
    inputs.extend(helpers + [record(__file__)])
    for item in inputs:
        check(item)
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "preparing_unscored", "accuracy_evaluated": False, "publication_ready": False,
              "source": record(__file__), "inputs": inputs, "job_id": os.environ["SLURM_JOB_ID"],
              "environment": dict(ENVIRONMENT), "numeric_checkpoint": numeric,
              "runtime_before": runtime, "arms": [],
              "limitations": ["Control plus four candidate variants only; two CPM variants remain separate.",
                  "Content checks do not admit native phylogeny, pair conversion or official scores.",
                  "Incremental shared-host preparation is not controlled end-to-end timing."]}
    try:
        prepare_arms(engine, baseline, names, species, (queries, targets, scores), output, report, audit_arm)
        if audit(checkpoint, checkpoint_record["sha256"]) != numeric:
            raise ValueError("Numeric checkpoint changed")
        report["runtime_after"] = verify(core, launcher, runtime_path)
        if report["runtime_after"] != runtime:
            raise ValueError("Scientific runtime changed during preparation")
        verify_admission(root, admission_path, admission_record["sha256"], "21758", 3)
        for item in inputs:
            check(item)
        for row in report["arms"]:
            for item in row["candidate_arm"]["output_files"]:
                check(item)
        report["status"] = "corrected_qfo_five_candidate_neighborhood_arms_prepared_unscored"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output.resolve())
