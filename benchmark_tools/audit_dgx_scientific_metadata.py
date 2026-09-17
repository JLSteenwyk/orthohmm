"""Check completed DGX command metadata without admitting scientific timings."""

import argparse
from copy import deepcopy
import csv
import io
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.gnu_time_companion import command as time_command
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.snapshot_orthohmm_input_order import record

SPEC_SHA = "fe96893489935124b675e3fb06317870cc7c155d9765e08067ad5b38b83df362"
PLAN_SHA = "7096348236f8372ef7f6ad12a3e829eae3dfee812b17c3b0bdd582480538ff1b"
RECIPE_SHA = "1de0ec6320defadef9bd22fda16ca7cc2f6f0ef78ffdf14b82d095fc8793fcc2"
WRAPPER_SHA = "36243fcdfe6f49e73dbe1dd1b627bf4a330f19f3cd52dcf94a9618ddf0415cf8"
PROJECT = Path("/home/jlsteenwyk/projects/orthohmm-publication")


def runtime_rows(spec):
    rows = [{**item, "records": count, "status": "runtime_tree_identity_matches",
             "scientific_execution_authorized": False}
            for item, count in zip(spec["runtime_manifests"], (26673, 10066))]
    rows.append({"path": str(PROJECT / "runtime_inventory_v1/native_launcher_recipe_v2.json"),
                 "sha256": RECIPE_SHA, "records": 14, "status": "runtime_tree_identity_matches",
                 "scientific_execution_authorized": False})
    return rows


def verify_contract(spec, run, prepared, verified, measured, scheduler):
    if (scheduler[0] != f"21656_{run['index']}" or scheduler[2:5] != ["COMPLETED", "0:0", "20"]
            or scheduler[5] not in ("96G", "96Gn") or scheduler[6] != "spark-7ff0"
            or int(scheduler[1]) != measured["job_id"]):
        raise ValueError("Scheduler identity, state or allocation differs")
    order = [r for r in spec["orders"] if r["input_directory"] == run["dataset"]["input_directory"]]
    if len(order) != 1:
        raise ValueError("Missing or ambiguous input order")
    order = order[0]
    expected = deepcopy(run)
    expected["gnu_time"] = {"executable": "/usr/bin/time", "output": str(Path(run["measurement_directory"]).parent / "native.time.tsv")}
    argv = time_command(run["native_argv"], expected["gnu_time"]["output"])
    if (prepared["run"] != expected or prepared["status"] != "fresh_native_inputs_prepared"
            or prepared["inference_started"] is not False or prepared["measured_argv"] != argv
            or measured["command"] != argv or measured["cwd"] != run["cwd"]):
        raise ValueError("Recorded preparation or measured command differs")
    if prepared["original_inputs"] != order["inputs_in_native_order"]:
        raise ValueError("Prepared original input identities differ")
    expected_copies = []
    expected_order = order["native_order"]
    if run["native_method"] == "orthofinder_full":
        expected_copies = sorted([{**r, "path": str(Path(run["configuration"]["copy_inputs_to"]) / Path(r["path"]).name)}
                                  for r in run["dataset"]["inputs"]], key=lambda r: r["path"])
        expected_order = sorted(expected_order)
    if prepared["copied_inputs"] != expected_copies or prepared["expected_native_basename_order"] != expected_order:
        raise ValueError("Prepared input copies or native order differ")
    if (verified["status"] != "command_exited_zero" or verified["measurement"] != measured
            or verified["source_sha256"] != WRAPPER_SHA
            or measured["status"] != "command_exited_zero" or measured["exit_code"] != 0
            or measured["timed_out"] is not False or "measurement_error" in measured):
        raise ValueError("Incomplete or inconsistent verified measurement")
    for side in ("before", "after"):
        observed = verified[side]
        if (observed["native_order"] != order["native_order"]
                or observed["original_inputs"] != order["inputs_in_native_order"]
                or observed["runtime"] != runtime_rows(spec)):
            raise ValueError("Runtime or original input verification differs")
        if side == "after" and run["native_method"] == "orthofinder_full" and observed.get("copied_inputs") != expected_copies:
            raise ValueError("Post-run copied input verification differs")
    expected_controls = {"requested_cpus": 20, "requested_memory_bytes": 96 * 1024**3,
                         "timeout_s": 85800, "interval_s": 1, "host_interval_s": 30}
    if any(measured[key] != value for key, value in expected_controls.items()):
        raise ValueError("Measurement controls differ")
    if measured["clock_domain"]["hostname"] != "spark-7ff0":
        raise ValueError("Wrong measurement host")
    for value in (prepared["preparation_wall_s"], verified["before_check_wall_s"], verified["after_check_wall_s"]):
        if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValueError("Invalid separate preparation or verification duration")


def audit(root, results, indices):
    if not indices or len(set(indices)) != len(indices) or any(i not in range(27) for i in indices):
        raise ValueError("Require unique explicit scientific indices")
    spec_path = results / "dgx_scientific_execution_20260917.json"
    plan_path = results / "dgx_scaling_commands_20260917.json"
    spec, plan = read_pinned(spec_path, SPEC_SHA), read_pinned(plan_path, PLAN_SHA)
    if spec["purpose"] != "scientific_scaling" or spec["runs"] != plan["runs"] or len(spec["runs"]) != 27:
        raise ValueError("Scientific specification differs from frozen panel")
    command = ["sacct", "-j", "21656", "--parsable2", "--noheader", "-X",
               "--format=JobID,JobIDRaw,State,ExitCode,AllocCPUS,ReqMem,NodeList"]
    accounting = subprocess.check_output(command, text=True)
    jobs = {}
    for row in csv.reader(io.StringIO(accounting), delimiter="|"):
        if row and row[0] in {f"21656_{i}" for i in indices}:
            if row[0] in jobs or len(row) != 7:
                raise ValueError("Duplicate or malformed accounting row")
            jobs[row[0]] = row
    records = []
    for index in indices:
        paths = [root / f"run_{index:02d}" / name for name in ("preparation.json", "verification.json", "results.json")]
        evidence = [record(path) for path in paths]
        prepared, verified, measured = [json.loads(p.read_text()) for p in paths]
        run = spec["runs"][index]
        verify_contract(spec, run, prepared, verified, measured, jobs[f"21656_{index}"])
        if evidence != [record(path) for path in paths]:
            raise ValueError("Metadata changed during audit")
        records.append({"index": index, "method": run["native_method"], "proteomes": run["proteomes"],
                        "repeat": run["repeat"], "scheduler": jobs[f"21656_{index}"], "evidence": evidence,
                        "status": "recorded_contract_matches", "host_evidence_status": measured["host_workload"]["status"]})
    return {"status": "selected_completed_metadata_contracts_verified_not_admitted", "runs": records,
            "specification": record(spec_path), "plan": record(plan_path), "source": record(__file__),
            "accounting_command": command, "accounting_snapshot": accounting,
            "scientific_timings_admitted": 0, "controlled_workload_verified": False,
            "limitations": ["Checks retained metadata assertions, not current remote bytes or complete runtime trees.",
                            "Does not independently replay resource samples or validate native outputs.",
                            "Host evidence status is retained, not upgraded to quiet-host certification.",
                            "Explicit subset only; pending, running and unselected tasks are not admitted or treated as failures."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "results", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--indices", nargs="+", type=int, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.root, args.results, args.indices)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
