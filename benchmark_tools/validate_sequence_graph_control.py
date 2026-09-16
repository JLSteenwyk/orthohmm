"""Admit both completed sequence-search graph controls before accuracy scoring."""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.audit_historical_profile_ablation import read_partition, verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.run_sequence_graph_control import check_profile_off, graph_command
from benchmark_tools.verify_qfo_replay_launcher import verify

EXECUTOR = "588d92fcf26195cc3b9a17ec3f0755e78f4438e0"
VARIANTS = ("all_hits", "top100")


def completed_tasks(accounting):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    selected = [r for r in rows if r["JobID"] in {"21293_0", "21293_1"}]
    if len(selected) != 2 or len({r["JobID"] for r in selected}) != 2:
        raise ValueError("Missing or duplicate graph-control scheduler tasks")
    if any(r["State"] != "COMPLETED" or r["ExitCode"] != "0:0" for r in selected):
        raise ValueError("Both graph controls must complete successfully before admission")
    return {VARIANTS[int(r["JobID"].split("_")[1])]: r for r in selected}


def check_records(preflight, result, metrics, label, scheduler, command, launcher):
    if preflight["variant"] != label or preflight["job_id"] != scheduler["JobIDRaw"]:
        raise ValueError("Graph-control scheduler or variant identity differs")
    if preflight["accuracy_evaluated"] is not False or result["accuracy_evaluated"] is not False:
        raise ValueError("Unexpected pre-admission accuracy evaluation")
    if result["status"] != "graph_complete_pending_scoring" or result["exit_code"] != 0:
        raise ValueError("Graph inference or its postflight failed")
    if preflight["command"] != command or metrics["command"] != command or metrics["cwd"] != str(launcher):
        raise ValueError("Graph replay command or working directory differs")
    expected = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if preflight["environment_overrides"] != expected:
        raise ValueError("Graph environment differs")
    check_profile_off(metrics)


def validate(root):
    root = root.resolve()
    accounting = subprocess.check_output(["sacct", "-j", "21293", "--parsable2",
                                         "--format=JobID,JobIDRaw,State,ExitCode,Elapsed"], text=True)
    tasks = completed_tasks(accounting)
    executor = root / "benchmarks/work/publication_ob_sequence_graph_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Graph executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    frozen = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(frozen, launcher, runtime_path)
    conversion_path = root / "benchmarks/results/ob_sequence_numeric_v1/manifest.json"
    conversion_record = file_provenance(conversion_path)
    conversion = json.loads(conversion_path.read_text())
    if (conversion["status"] != "numeric_checkpoints_verified" or conversion["accuracy_evaluated"] is not False
            or set(conversion["variants"]) != set(VARIANTS) or conversion["scheduler"]["JobIDRaw"] != "21291"):
        raise ValueError("Wrong conversion status or inventory")
    for record in (conversion["source"], conversion["execution"]):
        verify_file(Path(record["path"]), record)
    reports, tracked = {}, [conversion_record]
    for label in VARIANTS:
        directory = root / "benchmarks/results/ob_sequence_graph_v1" / label
        records = {name: file_provenance(directory / name) for name in ("preflight.json", "verification.json", "replay.json")}
        tracked.extend(records.values())
        preflight, result, metrics = [json.loads((directory / name).read_text()) for name in records]
        variant = conversion["variants"][label]
        if variant["cap"] != (None if label == "all_hits" else 100):
            raise ValueError("Hit cap differs")
        checkpoint = Path(variant["checkpoint"])
        command = graph_command(launcher, checkpoint, variant["manifest"]["sha256"], directory)
        check_records(preflight, result, metrics, label, tasks[label], command, launcher)
        if preflight["conversion"] != conversion_record or preflight["runtime"] != runtime:
            raise ValueError("Preflight conversion or runtime changed")
        if preflight["source"] != file_provenance(executor / "benchmark_tools/run_sequence_graph_control.py"):
            raise ValueError("Graph executor source differs")
        if metrics["source"] != file_provenance(launcher / "benchmark_tools/replay_high_sensitivity.py"):
            raise ValueError("Replay source differs")
        numeric = audit(checkpoint, variant["manifest"]["sha256"])
        if metrics["input"]["manifest"] != numeric["manifest"] or metrics["input"]["summary"] != numeric["summary"]:
            raise ValueError("Replay checkpoint admission differs")
        universe = set((checkpoint / "gene_names.txt").read_text().splitlines())
        if len(universe) != 251378 or result["genes"] != len(universe):
            raise ValueError("Incomplete input universe")
        prediction = directory / "replay/orthogroups_multipass_refined.txt"
        groups = read_partition(prediction, universe)
        if result["prediction"] != file_provenance(prediction) or result["groups"] != len(groups):
            raise ValueError("Final prediction changed")
        if result["metrics"] != records["replay.json"]:
            raise ValueError("Native metrics changed")
        for stage in metrics["stages"]:
            path = Path(stage["output"]["path"])
            expected = directory / "replay" / ("orthogroups_" + stage["label"] + ".txt")
            if path != expected or stage["clusters"] != len(read_partition(path, universe)):
                raise ValueError("Stage output path or cluster count differs")
            verify_file(path, stage["output"])
            tracked.append(stage["output"])
        reports[label] = {"scheduler": tasks[label], "native_records": records, "prediction": result["prediction"],
                          "genes": len(universe), "groups": len(groups), "numeric_admission": numeric}
    for record in tracked:
        verify_file(Path(record["path"]), record)
    if verify(frozen, launcher, runtime_path) != runtime:
        raise ValueError("Runtime changed during validation")
    return {"status": "native_graph_controls_verified", "accuracy_evaluated": False,
            "variants": reports, "conversion": conversion_record, "runtime": runtime,
            "executor_commit": EXECUTOR, "verifier": file_provenance(Path(__file__))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = validate(args.root)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
