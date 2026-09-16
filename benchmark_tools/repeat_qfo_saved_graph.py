"""Three instrumented fresh-worker repeats of the preserved QfO RBNH graph."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

COMPARISON_SHA = "ae225a72f4fcfc24b36d2b022dfdf1a386e5cdb2f82ab48ae40db49a6edd29c8"
CAPTURE_SHA = "3b087dfe4ddf73273c2f66ce10e94eeadad16f5b4c38ccf8836dff267fd172fd"
ENV_KEYS = ("PYTHONPATH", "PYTHONHASHSEED", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS", "LD_LIBRARY_PATH", "LD_PRELOAD")


def record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def mapped_libraries(text):
    paths = set()
    for line in text.splitlines():
        fields = line.split(maxsplit=5)
        if len(fields) != 6 or not fields[5].startswith("/"):
            continue
        path = fields[5]
        if ".so" in Path(path).name:
            if path.endswith(" (deleted)"):
                raise ValueError("A loaded native library has been deleted")
            paths.add(Path(path))
    return sorted(paths)


def worker(launcher, payload):
    # Import in the frozen worker namespace; do not load the development core.
    sys.path[0] = str(launcher)
    from orthohmm import leiden_worker
    import igraph
    import leidenalg
    import importlib.metadata
    import numpy
    if Path(leiden_worker.__file__).resolve() != launcher / "orthohmm/leiden_worker.py":
        raise ValueError("Wrong frozen worker import")
    metadata = json.loads((payload / "metadata.json").read_text())
    if metadata["cpm_resolution"] != .1 or metadata["seed"] != 4 or metadata["include_isolates"] is not True:
        raise ValueError("Changed clustering parameters")
    modules = {name: record(module.__file__) for name, module in list(sys.modules.items())
               if name.split(".")[0] in {"numpy", "igraph", "leidenalg", "orthohmm"}
               and getattr(module, "__file__", None)}
    libraries = [record(path) for path in mapped_libraries(Path("/proc/self/maps").read_text())]
    snapshot = {"status": "before_native_clustering", "accuracy_evaluated": False,
                "modules": modules, "native_libraries": libraries, "python": record(sys.executable),
                "versions": {name: importlib.metadata.version(name) for name in ("numpy", "igraph", "leidenalg")},
                "environment": {key: os.environ.get(key) for key in ENV_KEYS},
                "cpu_affinity": sorted(os.sched_getaffinity(0)), "platform": platform.platform(),
                "host": platform.node(), "cwd": str(Path.cwd()), "metadata": metadata,
                "inputs": [record(payload / name) for name in ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")],
                "observer": record(__file__)}
    (payload / "worker_before.json").write_text(json.dumps(snapshot, indent=2, sort_keys=True) + "\n")
    # The original worker flushes and calls os._exit(0); the parent verifies its result.
    leiden_worker.main([str(payload)])
    raise RuntimeError("Frozen worker unexpectedly returned instead of exiting")


def check_worker(snapshot, launcher, payload, overrides):
    expected = {"cpm_resolution": .1, "seed": 4, "include_isolates": True,
                "output_directory": str(payload.parent)}
    if (snapshot["status"] != "before_native_clustering" or snapshot["accuracy_evaluated"] is not False
            or snapshot["metadata"] != expected or snapshot["cwd"] != str(launcher)):
        raise ValueError("Unexpected worker parameters, status or working directory")
    if any(snapshot["environment"][key] != value for key, value in overrides.items()):
        raise ValueError("Worker environment overrides differ")
    modules = snapshot["modules"]
    for name in ("orthohmm.leiden_worker", "orthohmm.externals", "orthohmm.helpers"):
        if Path(modules[name]["path"]) != launcher / (name.replace(".", "/") + ".py"):
            raise ValueError("Worker imported a different scientific source")
    if not snapshot["native_libraries"] or not snapshot["cpu_affinity"]:
        raise ValueError("Missing native library or affinity evidence")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--worker-payload", type=Path)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    if args.worker_payload is not None:
        worker(launcher, args.worker_payload.resolve())
        return
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import numpy as np
    from benchmark_tools.audit_historical_profile_ablation import verify_file
    from benchmark_tools.capture_qfo_replay_graph import fingerprint
    from benchmark_tools.run_qfo_publication_replay import compare_partition
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.verify_qfo_replay_launcher import verify
    from benchmark_tools.verify_ygob_validation import require_completed_job
    if output.exists():
        raise FileExistsError(output)
    accounting = subprocess.check_output(["sacct", "-j", "21305", "--parsable2", "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21305)
    results = root / "benchmark_tools/results"
    comparison = read_frozen(results / "qfo_replay_initial_capture_comparison_20260916.json", COMPARISON_SHA)
    capture = read_frozen(results / "qfo_replay_initial_capture_20260916.json", CAPTURE_SHA)
    if comparison["status"] != "capture_complete" or comparison["rbnh_arrays_equal"] is not True:
        raise ValueError("Saved graph capture was not admitted")
    reference_records = [comparison["capture"], comparison["reference"], comparison["inputs_before"], comparison["inputs_after"],
                         comparison["initial_partition_comparison"]["expected"], comparison["initial_partition_comparison"]["observed"]]
    for item in reference_records:
        verify_file(Path(item["path"]), item)
    saved = root / "benchmarks/results/qfo_replay_initial_capture_v1"
    names = saved / "gene_names.txt"
    if record(names)["sha256"] != "246d02da9635576f04e94a51dbec4093d33f3481bf7f761e03afb3a53d51a0d3":
        raise ValueError("Saved gene order changed")
    sources = {"gene_names.txt": names}
    for key in ("sources", "targets", "weights"):
        path = saved / ("rbnh_" + key + ".npy")
        if fingerprint(np.load(path, mmap_mode="r", allow_pickle=False)) != capture["rbnh_arrays"][key]:
            raise ValueError("Saved graph array changed")
        sources[key + ".npy"] = path
    source_records = [record(path) for path in sources.values()]
    universe = set(names.read_text().splitlines())
    frozen = root / "benchmarks/work/publication_method_native_v2"
    runtime = results / "publication_native_runtime_20260916.json"
    before = verify(frozen, launcher, runtime)
    output.mkdir(parents=True)
    overrides = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                 "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    env = os.environ.copy()
    env.update(overrides)
    report = {"status": "running", "accuracy_evaluated": False, "job_id": os.environ.get("SLURM_JOB_ID"),
              "capture_scheduler": scheduler, "source": record(__file__), "runtime": before,
              "graph_inputs": source_records, "repeats": [], "limitations": [
                  "Instrumentation imports native modules before the frozen worker and records their loaded libraries.",
                  "Binary identity is recorded for these repeats, not retroactively proven for historical jobs.",
                  "Three repeats cannot establish general determinism; every partition is retained without accuracy selection.",
                  "No HMM searches, profile expansion, graph rebuilding or accuracy scoring are performed."]}
    try:
        for index in range(3):
            directory = output / f"repeat_{index}"
            payload = directory / "payload"
            payload.mkdir(parents=True)
            (directory / "orthohmm_working_res").mkdir()
            for name, path in sources.items():
                (payload / name).symlink_to(path)
            (payload / "metadata.json").write_text(json.dumps({"cpm_resolution": .1, "seed": 4,
                "include_isolates": True, "output_directory": str(directory)}) + "\n")
            command = [sys.executable, str(Path(__file__).resolve()), "--root", str(root), "--output", str(output),
                       "--worker-payload", str(payload)]
            started = time.monotonic()
            with (directory / "worker.log").open("x") as log:
                run = subprocess.run(command, cwd=launcher, env=env, stdout=log, stderr=subprocess.STDOUT)
            execution = {"command": command, "exit_code": run.returncode, "wall_s": time.monotonic() - started}
            (directory / "execution.json").write_text(json.dumps(execution, indent=2) + "\n")
            if run.returncode != 0:
                raise RuntimeError(f"Repeat {index} failed with exit {run.returncode}; preserved without retry")
            snapshot = json.loads((payload / "worker_before.json").read_text())
            check_worker(snapshot, launcher, payload, overrides)
            if snapshot["inputs"] != source_records:
                raise ValueError("Worker did not load the intended saved graph")
            for item in [*snapshot["modules"].values(), *snapshot["native_libraries"], snapshot["python"], snapshot["observer"], *source_records]:
                verify_file(Path(item["path"]), item)
            partition = directory / "orthohmm_working_res/orthohmm_edges_clustered.txt"
            row = {"index": index, "execution": execution, "worker": snapshot, "partition": record(partition),
                   "versus_capture": compare_partition(saved / "initial_partition.txt", partition, universe),
                   "versus_diagnostic": compare_partition(Path(comparison["initial_partition_comparison"]["expected"]["path"]), partition, universe)}
            if report["repeats"]:
                first = report["repeats"][0]
                row["versus_first_repeat"] = compare_partition(Path(first["partition"]["path"]), partition, universe)
                row["software_identity_equal"] = all(snapshot[key] == first["worker"][key]
                    for key in ("modules", "native_libraries", "python", "versions", "environment", "platform", "host", "cpu_affinity"))
            report["repeats"].append(row)
            (output / "progress.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        if verify(frozen, launcher, runtime) != before:
            raise ValueError("Frozen runtime changed during repeats")
        for item in [*source_records, *reference_records]:
            verify_file(Path(item["path"]), item)
        report["status"] = "three_repeats_complete"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
