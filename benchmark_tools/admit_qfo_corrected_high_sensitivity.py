"""Admit corrected native HMM evidence for downstream replay, not accuracy."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_corrected_primary import verify
from benchmark_tools.verify_ygob_validation import require_completed_job
from benchmark_tools.validate_high_sensitivity_outputs import validate

PLAN_SHA = "dbebd3a6915fddeb2b89e0e591a2ac5c798ee6bccec9925fef485260f41baa5a"
EXECUTOR = "9d8b608ab9b4a152d67ceacc04fe49b4d7595788"
METHOD = "orthohmm_high_sensitivity"


def validate_execution(plan, execution, scheduler, plan_record, runner_record):
    config = plan["methods"][METHOD]
    if execution.get("status") != "process_succeeded_pending_native_admission" or execution.get("exit_code") != 0:
        raise ValueError("Native execution has not succeeded")
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "32"
            or execution["job_id"] != scheduler["JobIDRaw"] or execution["node"] != "bizon"):
        raise ValueError("Wrong scheduler identity/allocation")
    if (execution["method"] != METHOD or type(execution["index"]) is not int or execution["index"] != 0
            or execution["array_task_id"] != "0" or not execution["array_job_id"]):
        raise ValueError("Wrong primary method/task")
    if execution["source"] != runner_record or execution["manifest"] != plan_record:
        raise ValueError("Source/plan binding differs")
    if execution["native_argv"] != config["native_argv"] or execution["cwd"] != config["cwd"]:
        raise ValueError("Native command/cwd differs")
    if not execution["finished_epoch"] >= execution["started_epoch"] > 0:
        raise ValueError("Invalid execution timestamps")
    if execution["accuracy_admitted"] is not False or execution["native_outputs_validated"] is not False:
        raise ValueError("Unexpected native admission state")
    directory = Path(plan["output_root"]) / "execution" / METHOD
    for key, name in (("log", "native.log"), ("timing", "time.txt")):
        if execution[key]["path"] != str(directory / name) or execution[key]["bytes"] <= 0:
            raise ValueError("Invalid execution log/timing")


def validate_metrics_execution(metrics, config, execution):
    argv = config["native_argv"]
    if argv[1:3] != ["-m", "orthohmm"]:
        raise ValueError("Unexpected native module invocation")
    expected = [argv[0], str(Path(config["cwd"]) / "orthohmm/__main__.py"), *argv[3:]]
    if metrics["command"] != expected or metrics["cwd"] != config["cwd"]:
        raise ValueError("Metrics command/cwd differs from execution")
    if not (execution["started_epoch"] <= metrics["started_at_epoch_s"]
            <= metrics["finished_at_epoch_s"] <= execution["finished_epoch"]):
        raise ValueError("Metrics timestamps outside native execution")


def admit(root, job, destination):
    if destination.exists():
        raise FileExistsError(destination)
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, job)
    plan_path = root / "benchmark_tools/results/qfo_corrected_primary_commands_20260918.json"
    plan = read_frozen(plan_path, PLAN_SHA)
    config = plan["methods"][METHOD]
    executor = root / "benchmarks/work/publication_qfo_corrected_primary_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Inference executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    status_path = Path(plan["output_root"]) / "execution" / METHOD / "status.json"
    status_record = record(status_path)
    execution = json.loads(status_path.read_text())
    validate_execution(plan, execution, scheduler, record(plan_path),
                       record(executor / "benchmark_tools/run_qfo_corrected_primary.py"))
    _, inputs = verify(plan)
    output = Path(config["output"])
    observed = [record(p) for p in sorted(output.rglob("*")) if p.is_file()]
    if not observed or observed != execution["outputs"]:
        raise ValueError("Native output inventory/content differs")
    # Native metrics are a sibling of the output directory, not in its inventory.
    metrics_path = Path(config["metrics"])
    metrics_record = record(metrics_path)
    metrics = json.loads(metrics_path.read_text())
    validate_metrics_execution(metrics, config, execution)
    checkpoint = output / "orthohmm_working_res/high_sensitivity_checkpoint"
    manifest = record(checkpoint / "manifest.json")
    checked = [status_record, record(plan_path), execution["source"], execution["log"], execution["timing"],
               metrics_record, *plan["inputs"], *inputs, *observed,
               *[record(Path(__file__).with_name(name)) for name in (
                   "run_qfo_corrected_primary.py", "validate_high_sensitivity_outputs.py",
                   "audit_accuracy_checkpoint.py", "audit_qfo_replay_inputs.py", "score_ygob_groups.py")]]
    for item in checked:
        check(item)
    content = validate(output, metrics_path, [Path(r["path"]) for r in inputs], manifest["sha256"])
    if content["genes"] != 984137 or len(content["species_ownership"]) != 78:
        raise ValueError("Corrected QfO input universe differs")
    for item in checked:
        check(item)
    if {str(p) for p in output.rglob("*") if p.is_file()} != {r["path"] for r in observed}:
        raise ValueError("Output membership changed during admission")
    report = {"status": "corrected_high_sensitivity_native_evidence_admitted", "source": record(__file__),
              "scheduler": scheduler, "accounting": accounting, "checked_records": checked,
              "content": content, "checkpoint_manifest": manifest, "accuracy_evaluated": False,
              "limitations": ["Native evidence admission for replay/conversion, not biological accuracy or search completeness.",
                              "Metrics hash acquired after completion; native runner did not inventory the sibling metrics file.",
                              "Runtime availability and runner checks are not a complete historical process execution trace.",
                              "Shared-host resources are not dedicated matched timing; conversion/scoring require separate admission."]}
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.job, args.output.resolve())
