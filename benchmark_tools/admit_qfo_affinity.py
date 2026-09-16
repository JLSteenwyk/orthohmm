"""Independently admit the completed saved-graph QfO CPU-affinity panel."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.repeat_qfo_saved_graph import affinity_arms, check_worker
from benchmark_tools.run_qfo_publication_replay import compare_partition
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job

REPORT_SHA = "e9311a5862b1db605b7d741003b2033de6f41f2c6db7e726316c181d00105de6"
IDENTITY_KEYS = ("modules", "native_libraries", "python", "versions", "environment", "platform", "host")


def validate_panel(report):
    if (report["status"] != "affinity_panel_complete" or report["accuracy_evaluated"] is not False
            or report["job_id"] != "21311" or report["affinity_panel"] is not True or len(report["repeats"]) != 4):
        raise ValueError("Incomplete or wrong affinity experiment")
    expected = affinity_arms(report["parent_cpu_affinity"])
    if report["planned_arms"] != [{"name": name, "cpu_affinity": cpus} for name, cpus in expected]:
        raise ValueError("Changed prespecified arm inventory")
    first = report["repeats"][0]["worker"]
    for index, (row, (label, cpus)) in enumerate(zip(report["repeats"], expected)):
        worker = row["worker"]
        if (row["index"] != index or row["arm"] != label or row["execution"]["exit_code"] != 0
                or worker["cpu_affinity"] != cpus or worker["requested_cpu_affinity"] != cpus
                or worker["inherited_cpu_affinity"] != report["parent_cpu_affinity"]
                or worker["inputs"] != report["graph_inputs"]
                or any(worker[key] != first[key] for key in IDENTITY_KEYS)):
            raise ValueError("Worker arm, graph or software identity differs")
        if index:
            if row["software_identity_excluding_affinity_equal"] is not True:
                raise ValueError("Incorrect software identity flag")
            if row["software_identity_equal"] != (worker["cpu_affinity"] == first["cpu_affinity"]):
                raise ValueError("Incorrect full identity flag")


def admit(root, output):
    if output.exists():
        raise FileExistsError(output)
    path = root / "benchmarks/results/qfo_saved_graph_affinity_v1/results.json"
    report = read_frozen(path, REPORT_SHA)
    validate_panel(report)
    accounting = subprocess.check_output(["sacct", "-j", "21311", "--parsable2", "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21311)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    overrides = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                 "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    records = [report["source"], *report["graph_inputs"]]
    comparisons = []
    universe = set(Path(report["graph_inputs"][0]["path"]).read_text().splitlines())
    if len(universe) != 976504:
        raise ValueError("Changed gene universe")
    for row in report["repeats"]:
        worker = row["worker"]
        directory = Path(worker["metadata"]["output_directory"])
        payload = directory / "payload"
        if json.loads((payload / "worker_before.json").read_text()) != worker:
            raise ValueError("Stored worker snapshot differs from native file")
        if json.loads((directory / "execution.json").read_text()) != row["execution"]:
            raise ValueError("Stored execution record differs")
        check_worker(worker, launcher, payload, overrides)
        records += [*worker["modules"].values(), *worker["native_libraries"], worker["python"], worker["observer"], row["partition"]]
        for key in ("versus_capture", "versus_diagnostic", "versus_first_repeat", "versus_previous_same_affinity"):
            if key not in row:
                continue
            comparison = row[key]
            fresh = compare_partition(Path(comparison["expected"]["path"]), Path(row["partition"]["path"]), universe)
            if fresh != comparison:
                raise ValueError("Complete native partition comparison differs: " + row["arm"] + ":" + key)
            records += [comparison["expected"], comparison["observed"]]
            comparisons.append({"arm": row["arm"], "comparison": key, "result": fresh})
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting identity for one file")
        unique[item["path"]] = item
    for item in unique.values():
        verify_file(Path(item["path"]), item)
    one = compare_partition(Path(report["repeats"][0]["partition"]["path"]), Path(report["repeats"][2]["partition"]["path"]), universe)
    many = compare_partition(Path(report["repeats"][1]["partition"]["path"]), Path(report["repeats"][3]["partition"]["path"]), universe)
    result = {"status": "affinity_panel_native_evidence_verified", "accuracy_evaluated": False, "publication_ready": False,
              "scheduler": scheduler, "source_report": file_provenance(path), "source": file_provenance(Path(__file__)),
              "verified_unique_files": len(unique), "verified_comparisons": comparisons,
              "one_cpu_repeat": one, "thirty_two_cpu_repeat": many, "native_report": report,
              "limitations": ["Same-affinity disagreement disproves general repeatability for this recorded configuration, not the validity of every earlier observation.",
                  "This bounded experiment does not isolate the cause of native partition variability.",
                  "No partition was selected by accuracy or closeness to a historical output; all outputs remain preserved.",
                  "Matching source/library/input hashes do not record every aspect of native runtime state.",
                  "Shared-node durations are diagnostic, not controlled efficiency measurements."]}
    output.mkdir(parents=True)
    (output / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
