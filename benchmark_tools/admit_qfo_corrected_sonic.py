"""Admit corrected SonicParanoid native pair tables before conversion."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_corrected_sonic import PLAN_SHA, verify
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job
from benchmark_tools.sonicparanoid_to_pairwise import find_pair_directory
from benchmark_tools.validate_sonicparanoid_tables import validate_tables

EXECUTOR = "d3f6dc0a7931579af779470b7f982e34a9d09b2b"


def validate_execution(plan, execution, scheduler, plan_record, runner_record):
    if execution.get("status") != "process_succeeded_pending_native_admission" or execution.get("exit_code") != 0:
        raise ValueError("Native execution has not succeeded")
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "32"
            or execution["job_id"] != scheduler["JobIDRaw"] or execution["node"] != "bizon"):
        raise ValueError("Wrong scheduler identity/allocation")
    if execution["source"] != runner_record or execution["plan"] != plan_record:
        raise ValueError("Source/plan binding differs")
    if execution["native_argv"] != plan["native_argv"] or execution["cwd"] != plan["cwd"]:
        raise ValueError("Command/cwd differs")
    if not execution["finished_epoch"] >= execution["started_epoch"] > 0:
        raise ValueError("Invalid execution timestamps")
    if not execution["runtime_before"] or execution["runtime_before"] != execution["runtime_after"]:
        raise ValueError("Runtime inventories differ")
    expected = [{**item, "path": str(Path(plan["copy_inputs_to"]) / Path(item["path"]).name)}
                for item in plan["input_fastas"]]
    if execution["copied_inputs"] != expected:
        raise ValueError("Copied input inventory differs")
    root = Path(plan["output_root"])
    for key, name in (("log", "native.log"), ("timing", "time.txt")):
        if execution[key]["path"] != str(root / name) or execution[key]["bytes"] <= 0:
            raise ValueError("Invalid execution log/timing")


def admit(root, job, destination):
    if destination.exists():
        raise FileExistsError(destination)
    plan_path = root / "benchmark_tools/results/qfo_corrected_sonic_commands_20260918.json"
    plan = read_frozen(plan_path, PLAN_SHA)
    run_root = Path(plan["output_root"])
    status_path = run_root / "execution.json"
    status_record = record(status_path)
    execution = json.loads(status_path.read_text())
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, job)
    executor = root / "benchmarks/work/publication_qfo_corrected_sonic_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Inference executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    validate_execution(plan, execution, scheduler, record(plan_path), record(executor / "benchmark_tools/run_qfo_corrected_sonic.py"))
    _, environment, runtime_inventory = verify(plan_path)
    if execution["environment"] != environment or execution["runtime_after"] != runtime_inventory:
        raise ValueError("Effective environment/runtime inventory differs")
    observed = [record(p) for p in sorted((run_root / "output").rglob("*")) if p.is_file()]
    if not observed or observed != execution["outputs"]:
        raise ValueError("Native output set/content differs")
    if {str(p) for p in (run_root / "input").iterdir()} != {r["path"] for r in execution["copied_inputs"]}:
        raise ValueError("Copied input directory membership differs")
    checked = [status_record, record(plan_path), execution["source"], execution["log"], execution["timing"],
               plan["source"], *plan["checked_records"], *execution["copied_inputs"], *observed,
               *[item["manifest"] for item in runtime_inventory],
               *[record(Path(__file__).with_name(name)) for name in (
                   "run_qfo_corrected_sonic.py", "validate_sonicparanoid_tables.py",
                   "sonicparanoid_to_pairwise.py", "fastoma_to_pairwise.py")]]
    for item in checked:
        check(item)
    directory = find_pair_directory(run_root / "output")
    inventory = validate_tables(directory, [Path(r["path"]) for r in plan["input_fastas"]])
    for item in checked:
        check(item)
    if {str(p) for p in (run_root / "output").rglob("*") if p.is_file()} != {r["path"] for r in observed}:
        raise ValueError("Output membership changed during validation")
    report = {"status": "corrected_sonic_native_tables_admitted_for_conversion", "source": record(__file__),
              "scheduler": scheduler, "accounting": accounting, "checked_records": checked,
              "runtime_inventory": runtime_inventory, "pair_directory": str(directory),
              "table_inventory": inventory, "accuracy_evaluated": False,
              "limitations": ["Representation/provenance validation, not biological accuracy or independent proof of algorithm completeness.",
                              "Native repeated relations use the frozen converter's deduplication behavior.",
                              "Conversion, reference filtering and assessment require separate validation.",
                              "Shared-host inference is not matched dedicated timing."]}
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
