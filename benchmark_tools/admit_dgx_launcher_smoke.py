"""Independently admit the three native launcher smokes, not scientific timings."""

import argparse
from copy import deepcopy
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_slurm_measurement import audit
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.monitor_slurm_resources import source_record
from benchmark_tools.validate_scaling_outputs import validate

SPEC_SHA = "8af8a049c69b7e479e8443e1e0e14b74c78d67b345c03d7f557768488d057fa2"
RECIPE_SHA = "ee62c74d1c9710967b00b417040c9d6a4a8b5b495c5f8afd4184542b3d76f595"
PROJECT = Path("/home/jlsteenwyk/projects/orthohmm-publication")


def verify_run(run, order, prepared, verified, measurement):
    expected = deepcopy(run)
    expected["gnu_time"] = {"executable": "/usr/bin/time",
                            "output": str(Path(run["measurement_directory"]).parent / "native.time.tsv")}
    if prepared["run"] != expected or prepared["status"] != "fresh_native_inputs_prepared":
        raise ValueError("Prepared scientific command differs from frozen smoke specification")
    if verified["status"] != "command_exited_zero" or verified["measurement"] != measurement:
        raise ValueError("Unsuccessful or inconsistent verified measurement")
    for side in ("before", "after"):
        observed = verified[side]
        if (observed["native_order"] != order["native_order"]
                or observed["original_inputs"] != order["inputs_in_native_order"]):
            raise ValueError("Input identity/enumeration verification differs")
        rows = observed["runtime"]
        expected_runtimes = [
            ("trees.json", "2f38fc57683e51a7b6768b4709db16293590983640c41cfe12a56ed6036750ae", 26673),
            ("system_trees.json", "4083d15c0ffc40013756c4588763aa8af24ca39165622665ee5074662e4b46ce", 10066),
            ("native_launcher_recipe_v1.json", RECIPE_SHA, 12)]
        expected_rows = [{"path": str(PROJECT / "runtime_inventory_v1" / name), "sha256": sha,
                          "records": count, "scientific_execution_authorized": False,
                          "status": "runtime_tree_identity_matches"} for name, sha, count in expected_runtimes]
        if rows != expected_rows:
            raise ValueError("Missing or changed runtime/recipe verification")
    if verified["source_sha256"] != "36243fcdfe6f49e73dbe1dd1b627bf4a330f19f3cd52dcf94a9618ddf0415cf8":
        raise ValueError("Unexpected verified wrapper source")
    if (measurement["requested_cpus"] != 20 or measurement["requested_memory_bytes"] != 96 * 1024 ** 3
            or measurement["timeout_s"] != 900 or measurement["interval_s"] != 1 or measurement["host_interval_s"] != 30):
        raise ValueError("Changed smoke measurement configuration")
    return expected


def report(root, results, original_inputs):
    spec_path = results / "dgx_launcher_smoke_spec_20260917.json"
    spec = read_pinned(spec_path, SPEC_SHA)
    accounting = subprocess.check_output(["sacct", "-j", "21653", "--parsable2", "--noheader",
        "--format=JobID,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)
    jobs = {}
    for row in csv.reader(io.StringIO(accounting), delimiter="|"):
        if row and row[0] in {"21653_0", "21653_1", "21653_2"}:
            if row[0] in jobs:
                raise ValueError("Duplicate scheduler record")
            jobs[row[0]] = row
    if len(jobs) != 3 or any(row[1:3] != ["COMPLETED", "0:0"] or row[4] != "20"
            or row[5] not in {"96G", "96Gn"} or row[6] != "spark-7ff0" for row in jobs.values()):
        raise ValueError("All three matched-resource smokes must complete successfully")
    maps = {PROJECT / "launcher_smoke_v1": root, PROJECT / "smoke_missing20_20261101_input": original_inputs}
    rows = []
    for run in spec["runs"]:
        directory = root / f"run_{run['index']:02d}"
        paths = [directory / name for name in ("preparation.json", "verification.json", "measurement/results.json")]
        evidence = [source_record(path) for path in paths]
        prepared, verified, measured = [json.loads(path.read_text()) for path in paths]
        expected = verify_run(run, spec["orders"][0], prepared, verified, measured)
        replay = audit(directory / "measurement", evidence[2]["sha256"])
        native = validate(expected, measured, evidence_roots=maps)
        if native["input_genes"] != 645:
            raise ValueError("Wrong smoke input universe")
        if evidence != [source_record(path) for path in paths]:
            raise ValueError("Smoke evidence changed during admission")
        rows.append({"index": run["index"], "method": run["native_method"],
                     "scheduler": jobs[f"21653_{run['index']}"], "evidence": evidence,
                     "collector_replay": replay, "native_outputs": native,
                     "native_wall_s_diagnostic_only": measured["command_wall_s"]})
    return {"status": "all_native_launcher_smokes_admitted", "runs": rows, "specification": source_record(spec_path),
            "source": source_record(__file__), "accounting": accounting,
            "scientific_scaling_runs_completed": 0, "scientific_execution_authorized": False,
            "limitations": ["Three small engineering smokes, not scientific scaling measurements or accuracy estimates.",
                            "Evidence is relocated for file access only; original remote native commands are checked unchanged.",
                            "Short host observation spans do not establish comparative timing quality.",
                            "Scientific launch requires a separate pinned authorization for the unchanged27-run plan."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "results", "original-inputs", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = report(args.root.resolve(), args.results.resolve(), args.original_inputs.resolve())
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
