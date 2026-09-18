"""Admit corrected Proteinortho native graph provenance before pair conversion."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_corrected_proteinortho import PLAN_SHA
from benchmark_tools.validate_proteinortho_graph import validate_graph

RUNTIME_SHA = "e0fc669cbbd805f9ec768ccd9088977609be3b4a88d9ff95bfb74f5b73f33ae5"


def validate_execution(plan, runtime, execution, plan_record, runtime_record):
    if execution.get("status") != "process_succeeded_pending_native_admission" or execution.get("exit_code") != 0:
        raise ValueError("Native execution has not succeeded")
    if runtime.get("status") != "proteinortho_runtime_frozen_unrun" or runtime["plan"] != plan_record:
        raise ValueError("Runtime plan binding differs")
    if execution["plan"] != plan_record or execution["runtime"] != runtime_record:
        raise ValueError("Execution manifest binding differs")
    if execution["source"]["sha256"] != runtime["source"]["sha256"]:
        raise ValueError("Execution source differs from frozen runner")
    if execution["native_argv"] != plan["native_argv"] or execution["cwd"] != plan["cwd"]:
        raise ValueError("Execution command/cwd differs")
    if execution["node"] != plan["resources"]["node"] or not execution.get("job_id"):
        raise ValueError("Missing/different scheduler provenance")
    if not execution["finished_epoch"] >= execution["started_epoch"] > 0:
        raise ValueError("Invalid execution timestamps")
    directory = Path(plan["cwd"])
    expected = [{**item, "path": str(directory / Path(item["path"]).name)} for item in plan["input_fastas"]]
    if execution["copied_inputs"] != expected:
        raise ValueError("Copied input inventory differs")
    outputs = execution["outputs"]
    output_map = {item["path"]: item for item in outputs}
    if len(output_map) != len(outputs) or any(not Path(p).is_relative_to(directory) for p in output_map):
        raise ValueError("Invalid output inventory")
    for item in expected:
        if output_map.get(item["path"]) != item:
            raise ValueError("Input absent/changed in final inventory")
    root = Path(plan["output_root"])
    for name in ("log", "timing"):
        expected_path = root / ("native.log" if name == "log" else "time.txt")
        if execution[name]["path"] != str(expected_path) or execution[name]["bytes"] <= 0:
            raise ValueError("Invalid execution log/timing record")
    native = {}
    for name, relative in plan["native_outputs"].items():
        item = output_map.get(str(root / relative))
        if item is None or item["bytes"] <= 0:
            raise ValueError("Missing/empty native output: " + name)
        native[name] = item
    return native


def admit(root, destination):
    if destination.exists():
        raise FileExistsError(destination)
    results = root / "benchmark_tools/results"
    plan_path = results / "qfo_corrected_proteinortho_commands_20260918.json"
    runtime_path = results / "qfo_corrected_proteinortho_runtime_20260918.json"
    plan = read_frozen(plan_path, PLAN_SHA)
    runtime = read_frozen(runtime_path, RUNTIME_SHA)
    execution_path = Path(plan["output_root"]) / "execution.json"
    execution_record = record(execution_path)
    execution = json.loads(execution_path.read_text())
    native = validate_execution(plan, runtime, execution, record(plan_path), record(runtime_path))
    actual_paths = {str(p) for p in Path(plan["cwd"]).rglob("*") if p.is_file()}
    if actual_paths != {r["path"] for r in execution["outputs"]}:
        raise ValueError("Native output tree differs from final execution inventory")
    records = [record(plan_path), record(runtime_path), execution_record, plan["source"],
               execution["source"], *plan["checked_records"], *execution["outputs"],
               execution["log"], execution["timing"],
               record(Path(__file__).with_name("validate_proteinortho_graph.py")),
               record(Path(__file__).with_name("fastoma_to_pairwise.py"))]
    for item in records:
        check(item)
    inventory = validate_graph(Path(native["pairs"]["path"]), [Path(r["path"]) for r in plan["input_fastas"]])
    for item in records:
        check(item)
    report = {"status": "corrected_proteinortho_native_graph_admitted_for_conversion",
              "source": record(__file__), "checked_records": records, "native_outputs": native,
              "graph_inventory": inventory, "accuracy_evaluated": False,
              "limitations": ["Graph representation/provenance admission, not independent proof of native algorithm completeness.",
                              "Group table is bound by hash, not independently validated for group-based scoring.",
                              "Pair conversion, reference filtering and assessment require separate validation.",
                              "Shared-host inference is not dedicated timing; current runtime availability is not a historical execution trace."]}
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve())
