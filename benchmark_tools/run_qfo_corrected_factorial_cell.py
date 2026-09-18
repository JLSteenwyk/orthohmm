"""Run a corrected QfO reconciliation cell from independently admitted candidates."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_candidates import validate_manifest
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_qfo_factorial_cell import native_command, ENVIRONMENT_SHA
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment, execute
from benchmark_tools.verify_qfo_replay_launcher import verify
from benchmark_tools.verify_ygob_validation import require_completed_job

ADMISSION_EXECUTOR = "9082ccee291176b8883884d80e80ff4817053b86"


def select_cell(admission, manifest, root, index):
    if (admission["status"] != "corrected_qfo_candidates_admitted"
            or admission["accuracy_evaluated"] is not False or admission["publication_ready"] is not False):
        raise ValueError("Require admitted corrected candidate evidence")
    executor = root / "benchmarks/work/publication_qfo_corrected_candidates_v1"
    parents = {Path(r["path"]).parent for r in manifest["input_fastas"]}
    if len(parents) != 1:
        raise ValueError("Ambiguous corrected FASTA directory")
    output = validate_manifest(manifest, admission["scheduler"], root, executor, next(iter(parents)))
    if admission["cells"] != manifest["cells"] or set(admission["candidate_arms"]) != set(manifest["candidate_arms"]):
        raise ValueError("Admitted cell inventory differs")
    for label, arm in manifest["candidate_arms"].items():
        if admission["candidate_arms"][label] != arm["content_audit"]:
            raise ValueError("Admitted candidate content differs")
    if type(index) is not int or not 0 <= index < 4:
        raise ValueError("Unknown corrected reconciliation cell")
    return [c for c in manifest["cells"] if c["reconciliation"]][index], output, executor


def verify_admission(root, path, digest, job, index):
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, job)
    if scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "2":
        raise ValueError("Unexpected candidate-admission resources")
    admission = read_frozen(path, digest)
    executor = root / "benchmarks/work/publication_qfo_corrected_candidate_admission_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMISSION_EXECUTOR:
        raise ValueError("Candidate admission executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    if admission["source"] != record(executor / "benchmark_tools/admit_qfo_corrected_candidates.py"):
        raise ValueError("Wrong independent candidate admission source")
    expected = root / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
    prepared = admission["prepared_manifest"]
    if prepared["path"] != str(expected) or prepared not in admission["checked_records"]:
        raise ValueError("Prepared manifest is not admitted corrected output")
    for item in admission["checked_records"]:
        check(item)
    manifest = read_frozen(expected, prepared["sha256"])
    cell, output, prepared_executor = select_cell(admission, manifest, root, index)
    runtime = verify(Path(manifest["core_root"]), Path(manifest["launcher_root"]),
                     root / "benchmark_tools/results/publication_native_runtime_20260916.json")
    if not runtime == manifest["runtime_before"] == manifest["runtime_after"]:
        raise ValueError("Runtime differs from admitted candidate preparation")
    if record(path)["sha256"] != digest:
        raise ValueError("Candidate admission changed during verification")
    return admission, manifest, cell, output, prepared_executor, scheduler


def run(root, admission_path, admission_sha, admission_job, environment_path, index, check_only=False):
    if not check_only:
        if not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "32":
            raise ValueError("Require scheduled 32-CPU reconciliation")
        if os.environ.get("SLURM_ARRAY_TASK_ID", str(index)) != str(index):
            raise ValueError("Array task differs from requested cell")
    admission, manifest, cell, output, prepared_executor, scheduler = verify_admission(
        root, admission_path, admission_sha, admission_job, index)
    launcher = Path(manifest["launcher_root"])
    argv, sources = native_command(cell, launcher, prepared_executor)
    environment = read_frozen(environment_path, ENVIRONMENT_SHA)
    verify_environment(environment)
    env, resolved = execution_environment(environment)
    env.update(manifest["environment_overrides"])
    env["PYTHONPATH"] = str(launcher)
    if check_only:
        return {"status": "corrected_cell_preflight_verified", "cell": cell["label"], "accuracy_evaluated": False}
    method = {"argv": argv, "output": argv[argv.index("--output-directory") + 1],
              "metrics": argv[argv.index("--json") + 1]}
    evidence = output / "execution" / cell["label"]
    provenance = {"candidate_admission": record(admission_path), "admission_scheduler": scheduler,
        "prepared_manifest": admission["prepared_manifest"], "environment": record(environment_path),
        "resolved_tools": resolved, "planned_argv": cell["argv"], "executed_argv": argv,
        "launcher_source_equivalence": sources, "executor": record(__file__),
        "execution_helper": record(Path(__file__).with_name("run_simulation_methods.py")),
        "command_helper": record(Path(__file__).with_name("run_qfo_factorial_cell.py")),
        "slurm_job_id": os.environ["SLURM_JOB_ID"], "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "cell_index": index, "cwd": str(launcher), "scope": "corrected-release incremental reconciliation; no scoring"}
    inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in manifest["input_fastas"]]}
    cwd = Path.cwd()
    try:
        os.chdir(launcher)
        result = execute({"label": cell["label"], "methods": {cell["label"]: method}},
                         [cell["label"]], env, evidence, inputs, provenance)
    finally:
        os.chdir(cwd)
    verify_admission(root, admission_path, admission_sha, admission_job, index)
    verify_environment(environment)
    for key in ("candidate_admission", "prepared_manifest", "environment", "executor", "execution_helper", "command_helper"):
        check(provenance[key])
    for pair in sources:
        check(pair["prepared"])
        check(pair["executed"])
    postflight = {"status": "corrected_inputs_runtime_sources_reverified", "failed_methods": result["failed_methods"],
                  "native_outputs_validated": False, "accuracy_evaluated": False}
    with (evidence / "postflight.json").open("x") as stream:
        json.dump(postflight, stream, indent=2, sort_keys=True)
        stream.write("\n")
    if result["failed_methods"]:
        raise RuntimeError("Corrected reconciliation failed; artifacts retained without retry or scoring")
    return postflight


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "admission", "environment-manifest"):
        parser.add_argument("--" + name, required=True, type=Path)
    for name in ("admission-sha256", "admission-job"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run(args.root.resolve(), args.admission.resolve(), args.admission_sha256,
                         args.admission_job, args.environment_manifest.resolve(), args.index, args.check_only), sort_keys=True))
