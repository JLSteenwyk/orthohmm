"""Run admitted corrected-QfO candidate variants through frozen inferred phylogeny."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import ARMS, check, record
from benchmark_tools.run_ob_candidate_neighborhood import variant_cell
from benchmark_tools.run_qfo_corrected_factorial_cell import verify_admission
from benchmark_tools.run_qfo_factorial_cell import ENVIRONMENT_SHA, native_command
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment, execute
from benchmark_tools.validate_simulation_outputs import verify_process
from benchmark_tools.verify_ygob_validation import require_completed_job

ADMISSION_SHA = "77759ee0f733242d6d4c368f44c2c1fd84ef3b1ebbb679494acce7fab68f0141"
BASELINE_SHA = "c8a289ac8128711da6c4e93654ce6953e48854b8d21b2adfea074d3feb8c70aa"
ADMISSION_COMMIT = "6a0258c15595775752bff5528d7ba678a66725df"
VARIANTS = tuple(label for label, _ in ARMS if label != "control")


def select_arm(admission, index):
    if (type(index) is not int or not 0 <= index < len(VARIANTS)
            or admission["status"] != "corrected_qfo_candidate_neighborhood_admitted_unscored"
            or admission["accuracy_evaluated"] is not False or admission["publication_ready"] is not False
            or [r["label"] for r in admission["arms"]] != ["control", *VARIANTS]):
        raise ValueError("Require exact independently admitted candidate panel")
    row = admission["arms"][index + 1]
    if row["status"] != "candidate_prepared_unscored":
        raise ValueError("Incomplete candidate arm")
    arm = row["candidate_arm"]
    return {"label": row["label"], "partition": arm["candidate_partition"],
            "constraints": arm["membership_constraints"]}


def verify_sources(root, index):
    results = root / "benchmark_tools/results"
    path = results / "qfo_parameter_candidate_admission_21929.json"
    admission = read_frozen(path, ADMISSION_SHA)
    arm = select_arm(admission, index)
    accounting = subprocess.check_output(["sacct", "-j", "21929", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, "21929")
    if scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "2":
        raise ValueError("Wrong candidate admission allocation")
    executor = root / "benchmarks/work/publication_qfo_parameter_candidate_admission_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMISSION_COMMIT:
        raise ValueError("Candidate admission executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    if admission["source"] != record(executor / "benchmark_tools/admit_qfo_candidate_neighborhood.py"):
        raise ValueError("Candidate admission source differs")
    for item in [*admission["checked_records"], admission["preparation"], arm["partition"], arm["constraints"]]:
        check(item)
    manifest, original, launcher, prepared_executor, environment, records = verify_baseline(root)
    return arm, manifest, original, launcher, prepared_executor, environment, [record(path), *records], scheduler


def verify_baseline(root):
    """Verify the reusable full-pipeline baseline without selecting a parameter arm."""
    results = root / "benchmark_tools/results"
    baseline_path = results / "qfo_corrected_factorial_native_admission_21764.json"
    baseline = read_frozen(baseline_path, BASELINE_SHA)
    if baseline["status"] != "corrected_qfo_native_pair_output_verified" or baseline["cell"] != "p1_c1_r1":
        raise ValueError("Require complete corrected full-pipeline baseline")
    for item in baseline["checked_records"]:
        check(item)
    old = baseline["candidate_admission"]
    _, manifest, original, _, prepared_executor, _ = verify_admission(
        root, Path(old["path"]), old["sha256"], "21758", 3)
    launcher = Path(manifest["launcher_root"])
    argv, equivalence = native_command(original, launcher, prepared_executor)
    status_record = baseline["native_group_integrity"]["integrity"]["execution_status"]
    check(status_record)
    status = json.loads(Path(status_record["path"]).read_text())
    config = {"argv": argv, "output": argv[argv.index("--output-directory") + 1],
              "metrics": argv[argv.index("--json") + 1]}
    verify_process(config, status["methods"][original["label"]])
    environment_path = results / "publication_variable_native_methods_20260916.json"
    environment = read_frozen(environment_path, ENVIRONMENT_SHA)
    verify_environment(environment)
    records = [record(baseline_path), record(environment_path), status_record]
    return manifest, original, launcher, prepared_executor, environment, records


def run(root, index):
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "32"
            or os.environ.get("SLURM_JOB_NODELIST") != "bizon"
            or os.environ.get("SLURM_ARRAY_TASK_ID") != str(index)):
        raise ValueError("Require scheduled 32-CPU bizon array task")
    if type(index) is not int or not 0 <= index < len(VARIANTS):
        raise ValueError("Unknown candidate variant")
    output = root / "benchmarks/results/qfo_parameter_phylogeny_v1" / VARIANTS[index]
    if output.exists():
        raise FileExistsError(output)
    arm, manifest, original, launcher, prepared_executor, environment, records, scheduler = verify_sources(root, index)
    cell = variant_cell(original, arm, output)
    argv, equivalence = native_command(cell, launcher, prepared_executor)
    env, resolved = execution_environment(environment)
    env.update(manifest["environment_overrides"], PYTHONPATH=str(launcher))
    helpers = [record(module.__file__) for name, module in sorted(sys.modules.items())
               if name.startswith("benchmark_tools.") and getattr(module, "__file__", None)]
    provenance = {"source": record(__file__), "helpers": helpers, "inputs": records,
                  "candidate_admission_scheduler": scheduler, "arm": arm, "cell": cell,
                  "executed_argv": argv, "launcher_source_equivalence": equivalence,
                  "resolved_tools": resolved, "cwd": str(launcher),
                  "job_id": os.environ["SLURM_JOB_ID"], "array_task_id": str(index),
                  "scope": "Corrected candidate variant; independently inferred phylogeny; unscored incremental execution"}
    output.mkdir(parents=True, exist_ok=False)
    (output / "preflight.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    method = {"argv": argv, "output": str(output / "output"), "metrics": str(output / "metrics.json")}
    inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in manifest["input_fastas"]]}
    postflight = {"status": "running", "accuracy_evaluated": False, "native_outputs_validated": False}
    cwd = Path.cwd()
    try:
        os.chdir(launcher)
        execution = execute({"label": cell["label"], "methods": {cell["label"]: method}},
                            [cell["label"]], env, output / "execution", inputs, provenance)
        os.chdir(cwd)
        verify_sources(root, index)
        for item in [*helpers, *records, provenance["source"]]:
            check(item)
        for pair in equivalence:
            check(pair["prepared"])
            check(pair["executed"])
        if execution["failed_methods"]:
            raise RuntimeError("Candidate phylogeny failed; preserve outputs without retry")
        postflight.update(status="complete_pending_native_validation", cell=cell)
    except BaseException as error:
        postflight.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        os.chdir(cwd)
        (output / "postflight.json").write_text(json.dumps(postflight, indent=2, sort_keys=True) + "\n")
    return postflight


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--index", required=True, type=int, choices=range(4))
    args = parser.parse_args()
    run(args.root.resolve(), args.index)
