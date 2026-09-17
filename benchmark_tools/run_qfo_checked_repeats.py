"""Three fresh Python-pair replays with full native graph integrity gates."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.checked_python_pair_worker import require_admission
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_leiden_boundary import saved_fingerprint
from benchmark_tools.repeat_qfo_saved_graph import check_worker
from benchmark_tools.run_qfo_publication_replay import compare_partition
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_qfo_replay_launcher import verify

ADMISSION_SHA = "bfb9f49e1eb32afb31ab898d0101cc785da1b54110c70d7ea06ad0259fb87631"


def constructor_digest(payload):
    import numpy as np
    sources = np.load(payload / "sources.npy", mmap_mode="r")
    targets = np.load(payload / "targets.npy", mmap_mode="r")
    if sources.dtype != np.int32 or targets.dtype != np.int32 or sources.ndim != 1 or sources.shape != targets.shape:
        raise ValueError("Unexpected saved constructor arrays")
    digest = hashlib.sha256()
    for start in range(0, len(sources), 100000):
        digest.update(np.column_stack((sources[start:start + 100000], targets[start:start + 100000])).tobytes(order="C"))
    return digest.hexdigest()


def check_gate(boundary, adapter, saved, input_digest):
    arguments = {"initial_membership": None, "weights": "weight", "n_iterations": 2,
                 "max_comm_size": 0, "seed": 4, "kwargs": {"resolution_parameter": .1},
                 "partition_type": "leidenalg.VertexPartition.CPMVertexPartition"}
    if boundary["accuracy_evaluated"] is not False or len(boundary["calls"]) != 1:
        raise ValueError("Require one unscored optimizer observation")
    call = boundary["calls"][0]
    if (call["status"] != "optimizer_returned" or call["arguments"] != arguments
            or any(call[key] != saved for key in ("before", "after", "saved"))):
        raise ValueError("Native integrity or optimizer arguments differ")
    if adapter["format"] != "python_pairs" or adapter["accuracy_evaluated"] is not False or len(adapter["calls"]) != 1:
        raise ValueError("Require exactly one Python-pair constructor")
    construction = adapter["calls"][0]
    if (construction["status"] != "constructor_returned" or construction["dtype"] != "int32"
            or construction["shape"] != [saved["edges"], 2] or construction["n"] != saved["vertices"]
            or construction["directed"] is not False
            or construction["ordered_input_bytes_sha256"] != input_digest):
        raise ValueError("Unexpected constructor conversion")


def run(root, output):
    if output.exists():
        raise FileExistsError(output)
    admission_path = root / "benchmark_tools/results/qfo_constructor_formats_verified_20260916.json"
    admission = read_frozen(admission_path, ADMISSION_SHA)
    require_admission(admission)
    for item in admission["provenance_checked"]:
        check(item)
    inputs = admission["native_report"]["graph_inputs"]
    names = Path(inputs[0]["path"]).read_text().splitlines()
    if len(names) != 976504 or len(set(names)) != len(names):
        raise ValueError("Changed QfO gene universe")
    universe = set(names)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    frozen = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(frozen, launcher, runtime_path)
    worker_path = Path(__file__).with_name("checked_python_pair_worker.py")
    overrides = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                 "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    env = {**os.environ, **overrides}
    output.mkdir(parents=True)
    report = {"status": "running", "accuracy_evaluated": False, "job_id": os.environ.get("SLURM_JOB_ID"),
              "source": record(__file__), "worker_source": record(worker_path), "admission": record(admission_path),
              "runtime": runtime, "graph_inputs": inputs, "repeats": [], "planned_repeats": 3,
              "limitations": ["Initial saved graph only; not full HMM/profile replay or historical equivalence.",
                  "Three repeats are bounded evidence, not general determinism; no partition is selected by accuracy.",
                  "Python-pair conversion and instrumentation change allocation/timing; shared-node costs are diagnostic."]}
    try:
        for index in range(3):
            directory = output / f"repeat_{index}"
            payload = directory / "payload"
            payload.mkdir(parents=True)
            (directory / "orthohmm_working_res").mkdir()
            for item in inputs:
                source = Path(item["path"])
                (payload / source.name.removeprefix("rbnh_")).symlink_to(source)
            (payload / "metadata.json").write_text(json.dumps({"cpm_resolution": .1, "seed": 4,
                "include_isolates": True, "output_directory": str(directory)}) + "\n")
            command = [sys.executable, str(worker_path.resolve()), "--root", str(root), "--payload", str(payload),
                       "--admission", str(admission_path), "--admission-sha256", ADMISSION_SHA]
            started = time.monotonic()
            with (directory / "worker.log").open("x") as log:
                process = subprocess.run(command, cwd=launcher, env=env, stdout=log, stderr=subprocess.STDOUT)
            execution = {"command": command, "exit_code": process.returncode, "wall_s": time.monotonic() - started}
            (directory / "execution.json").write_text(json.dumps(execution, indent=2) + "\n")
            if process.returncode:
                raise RuntimeError(f"Checked worker {index} failed; preserved without retry")
            snapshot = json.loads((payload / "worker_before.json").read_text())
            check_worker(snapshot, launcher, payload, overrides)
            if snapshot["inputs"] != inputs or len(snapshot["cpu_affinity"]) != 1:
                raise ValueError("Worker graph or CPU differs")
            boundary = json.loads((payload / "native_boundary.json").read_text())
            adapter = json.loads((payload / "constructor_adapter.json").read_text())
            check_gate(boundary, adapter, saved_fingerprint(payload), constructor_digest(payload))
            provenance = json.loads((payload / "checked_worker_provenance.json").read_text())
            expected_helpers = [record(worker_path.with_name(name)) for name in
                                ("repeat_qfo_saved_graph.py", "probe_leiden_boundary.py")]
            if (provenance["source"] != report["worker_source"] or provenance["admission"] != report["admission"]
                    or provenance["inputs"] != inputs or provenance["helpers"] != expected_helpers):
                raise ValueError("Checked worker provenance changed")
            evidence = [record(payload / name) for name in ("worker_before.json", "native_boundary.json",
                        "constructor_adapter.json", "checked_worker_provenance.json", "metadata.json")]
            for item in [*snapshot["modules"].values(), *snapshot["native_libraries"], snapshot["python"], snapshot["observer"], *inputs]:
                check(item)
            partition = directory / "orthohmm_working_res/orthohmm_edges_clustered.txt"
            # Even the first output must have complete, nonduplicated gene coverage.
            comparison = compare_partition(Path(report["repeats"][0]["partition"]["path"]) if index else partition,
                                           partition, universe)
            if index and any(snapshot[key] != report["repeats"][0]["worker"][key] for key in
                    ("modules", "native_libraries", "python", "versions", "environment", "platform", "host", "cpu_affinity")):
                raise ValueError("Worker identity changed across repeats")
            report["repeats"].append({"index": index, "execution": execution, "worker": snapshot,
                "partition": record(partition), "versus_first": comparison, "native_boundary": boundary,
                "constructor_adapter": adapter, "evidence": evidence})
            (output / "progress.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        if verify(frozen, launcher, runtime_path) != runtime:
            raise ValueError("Frozen runtime changed")
        for item in admission["provenance_checked"]:
            check(item)
        report["status"] = "three_checked_repeats_complete_unscored"
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
    args = parser.parse_args()
    run(args.root.resolve(), args.output.resolve())
