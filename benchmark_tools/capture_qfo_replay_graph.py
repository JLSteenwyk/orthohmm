"""Capture frozen replay's initial graph path in a fresh, instrumented process."""

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


class CaptureComplete(Exception):
    """Expected stop before second clustering, refinement or profiles."""


def fingerprint(array):
    import numpy as np
    value = np.ascontiguousarray(array)
    return {"dtype": str(value.dtype), "shape": list(value.shape),
            "sha256": hashlib.sha256(memoryview(value).cast("B")).hexdigest()}


class Recorder:
    def __init__(self, module, output):
        self.module, self.output = module, output
        self.cluster = module.execute_leiden
        self.singletons = module.build_singleton_assignment_edges
        self.calls = 0
        self.report = {"status": "capturing", "accuracy_evaluated": False}

    def cluster_call(self, resolution, output_directory, edges=None, include_isolates=False, seed=0):
        self.calls += 1
        if resolution != .1 or seed != 4 or include_isolates is not True:
            raise ValueError("Unexpected initial clustering parameters")
        record = {key: fingerprint(getattr(edges, key)) for key in ("sources", "targets", "weights")}
        if self.calls == 1:
            import numpy as np
            self.report["rbnh_arrays"] = record
            for key in record:
                np.save(self.output / ("rbnh_" + key + ".npy"), getattr(edges, key), allow_pickle=False)
            with (self.output / "gene_names.txt").open("x") as handle:
                handle.write("\n".join(edges.gene_names) + "\n")
            self.cluster(resolution, output_directory, edges=edges, include_isolates=include_isolates, seed=seed)
            shutil.copyfile(Path(output_directory) / "orthohmm_working_res/orthohmm_edges_clustered.txt",
                            self.output / "initial_partition.txt")
        elif self.calls == 2:
            if "singleton_arrays" not in self.report:
                raise ValueError("Second clustering reached without singleton assignment")
            self.report["multipass_arrays"] = record
            self.report["status"] = "captured_before_second_clustering"
            raise CaptureComplete()
        else:
            raise ValueError("Unexpected extra clustering call")

    def singleton_call(self, *args, **kwargs):
        if self.calls != 1 or "singleton_arrays" in self.report:
            raise ValueError("Unexpected singleton-assignment order")
        result = self.singletons(*args, **kwargs)
        self.report["singleton_edges"] = len(result)
        self.report["singleton_arrays"] = {key: fingerprint(getattr(result, key)) for key in ("sources", "targets", "weights")}
        return result


def worker(root, output):
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    # Start imports exactly in the launcher namespace, not the wrapper checkout.
    sys.path[0] = str(launcher)
    from benchmark_tools import replay_high_sensitivity as replay
    from benchmark_tools.orthobench_stage_diagnostics import file_provenance
    if Path(replay.__file__).resolve() != launcher / "benchmark_tools/replay_high_sensitivity.py":
        raise ValueError("Wrong replay import")
    command = json.loads((output / "replay_command.json").read_text())
    if command[:2] != [sys.executable, str(launcher / "benchmark_tools/replay_high_sensitivity.py")]:
        raise ValueError("Wrong frozen replay command")
    recorder = Recorder(replay, output)
    recorder.report.update(command=command, replay_source=file_provenance(Path(replay.__file__)),
                           cluster_source=file_provenance(Path(sys.modules[recorder.cluster.__module__].__file__)),
                           singleton_source=file_provenance(Path(sys.modules[recorder.singletons.__module__].__file__)))
    replay.execute_leiden = recorder.cluster_call
    replay.build_singleton_assignment_edges = recorder.singleton_call
    try:
        replay.main(command[2:])
    except CaptureComplete:
        pass
    except BaseException as error:
        recorder.report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "capture.json").write_text(json.dumps(recorder.report, indent=2, sort_keys=True) + "\n")
    if recorder.calls != 2 or recorder.report["status"] != "captured_before_second_clustering":
        raise ValueError("Replay did not reach the planned capture stop")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    if args.worker:
        worker(root, output)
        return
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from benchmark_tools.verify_qfo_replay_launcher import verify
    from benchmark_tools.orthobench_stage_diagnostics import file_provenance
    from benchmark_tools.audit_historical_profile_ablation import verify_file
    from benchmark_tools.run_qfo_publication_replay import compare_partition, command_for
    if output.exists():
        raise FileExistsError(output)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    frozen = root / "benchmarks/work/publication_method_native_v2"
    runtime = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    before = verify(frozen, launcher, runtime)
    reference_path = root / "benchmarks/results/qfo_initial_graph_diagnostic_v1/diagnostic.json"
    reference_record = file_provenance(reference_path)
    reference = json.loads(reference_path.read_text())
    if reference["status"] != "diagnostic_complete" or not reference["repeat_partition_comparison"]["partition_equal"]:
        raise ValueError("Initial graph diagnostic did not complete")
    for record in [*reference["graph_arrays"], *[r["partition"] for r in reference["repeats"]]]:
        verify_file(Path(record["path"]), record)
    output.mkdir(parents=True)
    audit_script = Path(__file__).resolve().parent / "audit_qfo_replay_inputs.py"
    def audit(name):
        path = output / name
        subprocess.run([sys.executable, str(audit_script), "--root", str(root), "--output", str(path)], check=True)
        return json.loads(path.read_text())
    inputs = audit("inputs_before.json")
    packages = sorted((d.metadata["Name"], d.version) for d in importlib.metadata.distributions())
    (output / "replay_command.json").write_text(json.dumps(command_for(root, launcher, output)) + "\n")
    env = os.environ.copy()
    env.update(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    command = [sys.executable, str(Path(__file__).resolve()), "--root", str(root), "--output", str(output), "--worker"]
    with (output / "worker.log").open("x") as handle:
        subprocess.run(command, cwd=launcher, env=env, stdout=handle, stderr=subprocess.STDOUT, check=True)
    capture = json.loads((output / "capture.json").read_text())
    universe = set((output / "gene_names.txt").read_text().splitlines())
    comparison = compare_partition(Path(reference["repeats"][0]["partition"]["path"]), output / "initial_partition.txt", universe)
    if verify(frozen, launcher, runtime) != before:
        raise ValueError("Runtime changed during capture")
    if audit("inputs_after.json") != inputs:
        raise ValueError("Inputs changed during capture")
    if sorted((d.metadata["Name"], d.version) for d in importlib.metadata.distributions()) != packages:
        raise ValueError("Installed package versions changed during capture")
    verify_file(reference_path, reference_record)
    report = {"status": "capture_complete", "accuracy_evaluated": False, "job_id": os.environ.get("SLURM_JOB_ID"),
              "source": file_provenance(Path(__file__)), "runtime": before, "reference": reference_record,
              "capture": file_provenance(output / "capture.json"), "command": command,
              "packages": packages, "environment_overrides": {key: env[key] for key in
                  ("PYTHONPATH", "PYTHONHASHSEED", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
              "inputs_before": file_provenance(output / "inputs_before.json"),
              "inputs_after": file_provenance(output / "inputs_after.json"),
              "rbnh_arrays_equal": all(capture["rbnh_arrays"][key] == value["frozen"] for key, value in reference["rbnh_comparison"]["arrays"].items()),
              "initial_partition_comparison": comparison, "singleton_edges": capture["singleton_edges"],
              "singleton_arrays_equal": capture["singleton_arrays"] == reference["repeats"][0]["singleton_arrays"],
              "limitations": "Instrumentation adds observation and an early stop; no profiles, second clustering, refinement or accuracy run."}
    (output / "comparison.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
