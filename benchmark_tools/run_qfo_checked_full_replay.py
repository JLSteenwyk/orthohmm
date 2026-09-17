"""Run one complete frozen cached QfO replay with every clustering payload checked."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


def replay_worker(root, output):
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    executor = Path(__file__).resolve().parent.parent
    sys.path[0] = str(launcher)
    from benchmark_tools import replay_high_sensitivity as replay
    from orthohmm import externals
    if Path(replay.__file__).resolve() != launcher / "benchmark_tools/replay_high_sensitivity.py":
        raise ValueError("Wrong frozen replay import")
    sys.path.insert(0, str(executor / "benchmark_tools"))
    from checked_replay_interceptor import CheckedReplaySubprocess
    from validate_checked_replay_payload import validate
    from prepare_ob_candidate_neighborhood import record
    command = json.loads((output / "replay_command.json").read_text())
    if command[:2] != [sys.executable, str(launcher / "benchmark_tools/replay_high_sensitivity.py")]:
        raise ValueError("Unexpected replay command")
    original = externals.subprocess
    proxy = CheckedReplaySubprocess(original, root, output / "clustering", executor / "benchmark_tools/checked_replay_payload_worker.py",
        lambda payload, manifest: validate(payload, manifest, root, executor))
    report = {"status": "running", "accuracy_evaluated": False, "source": record(__file__),
              "replay_source": record(replay.__file__), "command": command, "calls": proxy.calls,
              "helpers": [record(executor / "benchmark_tools" / name) for name in
                          ("checked_replay_interceptor.py", "validate_checked_replay_payload.py", "checked_replay_payload_worker.py")]}
    externals.subprocess = proxy
    try:
        sys.argv = command[1:]
        replay.main(command[2:])
        if len(proxy.calls) != 4 or any(row["status"] != "checked" for row in proxy.calls):
            raise ValueError("Replay did not complete exactly four checked clustering calls")
        report["status"] = "checked_full_replay_returned"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        externals.subprocess = original
        (output / "checked_worker.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def run(root, output):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.verify_qfo_replay_launcher import verify
    from benchmark_tools.run_qfo_publication_replay import command_for, compare_partition
    from benchmark_tools.checked_replay_payload_worker import ADMISSION_SHA
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    report = {"status": "preflight", "accuracy_evaluated": False, "job_id": os.environ.get("SLURM_JOB_ID"),
              "source": record(__file__), "limitations": ["One checked cached replay; not an end-to-end timing or general determinism experiment.",
                  "No accuracy-based selection, default promotion or silent replacement of historical outputs."]}
    executor = Path(__file__).resolve().parent.parent
    report["executor_commit"] = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    report["executor_sources"] = [record(executor / "benchmark_tools" / name) for name in
        ("run_qfo_checked_full_replay.py", "checked_replay_interceptor.py", "validate_checked_replay_payload.py",
         "checked_replay_payload_worker.py", "checked_python_pair_worker.py", "repeat_qfo_saved_graph.py", "probe_leiden_boundary.py")]
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    try:
        admission_path = root / "benchmark_tools/results/qfo_checked_repeats_verified_20260917.json"
        admission = read_frozen(admission_path, ADMISSION_SHA)
        if admission["status"] != "checked_repeats_verified" or not admission["all_three_partitions_equal"]:
            raise ValueError("Initial checked repeats not admitted")
        for item in admission["provenance_checked"]:
            check(item)
        report["admission"] = record(admission_path)
        report["runtime"] = verify(core, launcher, runtime_path)
        def audit(label):
            path = output / (label + ".json")
            subprocess.run([sys.executable, str(executor / "benchmark_tools/audit_qfo_replay_inputs.py"),
                            "--root", str(root), "--output", str(path)], check=True)
            return json.loads(path.read_text())
        inputs = audit("inputs_before")
        report["inputs_before"] = record(output / "inputs_before.json")
        command = command_for(root, launcher, output)
        (output / "replay_command.json").write_text(json.dumps(command) + "\n")
        env = {**os.environ, "PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
               "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        worker_command = [sys.executable, str(Path(__file__).resolve()), "--root", str(root), "--output", str(output), "--replay-worker"]
        report.update(status="running", worker_command=worker_command, replay_command=command)
        with (output / "replay.log").open("x") as log:
            result = subprocess.run(worker_command, cwd=launcher, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = result.returncode
        if result.returncode:
            raise RuntimeError("Full checked replay failed; retained without retry")
        checked = json.loads((output / "checked_worker.json").read_text())
        replay = json.loads((output / "replay.json").read_text())
        if (checked["status"] != "checked_full_replay_returned" or len(checked["calls"]) != 4
                or replay["counts"]["genes"] != 976504 or replay["counts"].get("profiles_built", 0) <= 0
                or [s["label"] for s in replay["stages"]] != ["multipass", "multipass_refined", "profiles", "profiles_refined"]):
            raise ValueError("Incomplete full replay or missing profile construction")
        names = Path(admission["native_report"]["graph_inputs"][0]["path"]).read_text().splitlines()
        universe = set(names)
        coverage = []
        for stage in replay["stages"]:
            path = Path(stage["output"]["path"])
            check(stage["output"])
            coverage.append({"stage": stage["label"], **compare_partition(path, path, universe)})
        first = checked["calls"][0]
        for actual, expected in zip(json.loads(Path(first["manifest"]["path"]).read_text())["inputs"][:4], admission["native_report"]["graph_inputs"]):
            if any(actual[k] != expected[k] for k in ("bytes", "sha256")):
                raise ValueError("Regenerated initial payload changed")
        report["initial_versus_checked_repeat"] = compare_partition(
            Path(admission["native_report"]["repeats"][0]["partition"]["path"]), Path(first["partition"]["path"]), universe)
        if audit("inputs_after") != inputs or verify(core, launcher, runtime_path) != report["runtime"]:
            raise ValueError("Inputs or frozen native runtime changed")
        for item in admission["provenance_checked"]:
            check(item)
        report.update(status="full_checked_replay_complete_unscored", worker=record(output / "checked_worker.json"),
            replay=record(output / "replay.json"), inputs_after=record(output / "inputs_after.json"), coverage=coverage)
        for item in report["executor_sources"]:
            check(item)
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--replay-worker", action="store_true")
    args = parser.parse_args()
    (replay_worker if args.replay_worker else run)(args.root.resolve(), args.output.resolve())
