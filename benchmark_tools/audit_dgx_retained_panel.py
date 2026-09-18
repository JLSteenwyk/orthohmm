"""Review transferred DGX metadata and host evidence without admitting timings."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_dgx_scientific_metadata import verify_contract, SPEC_SHA, PLAN_SHA
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.review_dgx_host_intervals import review
from benchmark_tools.gnu_time_companion import parse as parse_time


def scheduler_rows(snapshot, final):
    rows = snapshot["runs"]
    if len(rows) != 27 or {r["JobID"] for r in rows} != {f"21656_{i}" for i in range(27)}:
        raise ValueError("Incomplete or duplicate retained scheduler panel")
    result = {}
    for row in rows:
        if row["JobID"] == "21656_26":
            observed = final["preceding_successful_controller_observation"]
            if (observed["exit_code"] != 0 or observed["ArrayJobId"] != "21656"
                    or observed["ArrayTaskId"] != "26" or observed["JobState"] != "COMPLETED"
                    or observed["ExitCode"] != "0:0" or observed["JobId"] != row["JobIDRaw"]
                    or observed["NodeList"] != row["NodeList"] or observed["NumCPUs"] != 20
                    or observed["AllocTRES"] != "cpu=20,mem=96G,node=1,billing=20"
                    or row["AllocCPUS"] != "20" or row["ReqMem"] not in ("96G", "96Gn")):
                raise ValueError("Final controller evidence does not bind retained task")
            state, exit_code = observed["JobState"], observed["ExitCode"]
            origin = "transcribed_terminal_controller_fields_not_fresh_sacct"
        else:
            state, exit_code = row["State"], row["ExitCode"]
            origin = "retained_sacct_snapshot"
        if state != "COMPLETED" or exit_code != "0:0":
            raise ValueError("Nonterminal or failed retained task")
        result[row["JobID"]] = {"origin": origin, "contract_row": [row["JobID"], row["JobIDRaw"],
            state, exit_code, row["AllocCPUS"], row["ReqMem"], row["NodeList"]]}
    return result


def audit(root, results):
    spec_path, plan_path = (results / name for name in
        ("dgx_scientific_execution_20260917.json", "dgx_scaling_commands_20260917.json"))
    spec, plan = read_pinned(spec_path, SPEC_SHA), read_pinned(plan_path, PLAN_SHA)
    if spec["purpose"] != "scientific_scaling" or spec["runs"] != plan["runs"] or len(spec["runs"]) != 27:
        raise ValueError("Frozen scientific panel differs")
    snapshot_path = results / "dgx_scheduler_progress_20260918.json"
    final_path = results / "dgx_final_task_controller_20260918.json"
    provenance = [record(p) for p in (spec_path, plan_path, snapshot_path, final_path)]
    jobs = scheduler_rows(json.loads(snapshot_path.read_text()), json.loads(final_path.read_text()))
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "audit_dgx_scientific_metadata.py", "review_dgx_host_intervals.py", "command_host_monitor.py",
        "observe_host_competition.py", "gnu_time_companion.py")]
    inventory = [record(p) for p in sorted(root.rglob("*")) if p.is_file()]
    rows = []
    for run in spec["runs"]:
        index = run["index"]
        directory = root / f"run_{index:02d}"
        scheduler = jobs[f"21656_{index}"]
        row = {"index": index, "method": run["native_method"], "proteomes": run["proteomes"],
               "repeat": run["repeat"], "scheduler_evidence": scheduler}
        try:
            prepared = json.loads((directory / "preparation.json").read_text())
            verified = json.loads((directory / "verification.json").read_text())
            measurement_path = directory / "measurement/results.json"
            measured = json.loads(measurement_path.read_text())
            verify_contract(spec, run, prepared, verified, measured, scheduler["contract_row"])
            raw_records = []
            for name in ("samples.jsonl", "host_samples.jsonl", "command.log"):
                actual = record(directory / "measurement" / name)
                expected = measured[name]
                if (actual["sha256"] != expected["sha256"]
                        or expected["path"] != str(Path(run["measurement_directory"]) / name)):
                    raise ValueError("Retained measurement payload hash/path differs")
                raw_records.append(actual)
            timing = parse_time((directory / "native.time.tsv").read_text())
            if timing["exit_status"] != measured["exit_code"]:
                raise ValueError("GNU-time exit status differs")
            row.update(status="metadata_and_payload_hashes_verified_not_admitted", raw_records=raw_records,
                       gnu_time=timing, host_review=review(directory / "measurement"))
        except Exception as error:
            row.update(status="review_failed", error_type=type(error).__name__, error=str(error))
        rows.append(row)
    for item in [*provenance, *helpers, *inventory]:
        check(item)
    if {str(p) for p in root.rglob("*") if p.is_file()} != {r["path"] for r in inventory}:
        raise ValueError("Transferred file inventory changed")
    return {"status": "retained_panel_review_complete_not_timing_admission", "source": record(__file__),
        "provenance": provenance, "helpers": helpers, "transferred_inventory": inventory, "runs": rows,
        "review_failures": sum(r["status"] == "review_failed" for r in rows),
        "scientific_timings_admitted": 0, "controlled_workload_verified": False,
        "limitations": ["Final-task scheduler fields are transcribed from a controller response; fresh sacct unavailable.",
            "Recorded commands/runtime assertions checked, not independently proven historical execution.",
            "Host replay retains inconclusive evidence; process names do not authenticate kernel threads.",
            "Resource-sample accounting replay and native output validation remain separate requirements.",
            "GNU-time is retained with its native scope, not substituted for aggregate memory accounting."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "results", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.root.resolve(), args.results.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"runs": len(result["runs"]), "review_failures": result["review_failures"],
                      "scientific_timings_admitted": 0}))
