"""Independently audit the fixed three-repeat checked QfO optimizer diagnostic."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.repeat_qfo_saved_graph import check_worker
from benchmark_tools.run_qfo_publication_replay import compare_partition
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_qfo_replay_launcher import verify
from benchmark_tools.verify_ygob_validation import require_completed_job

ADMISSION_SHA = "bfb9f49e1eb32afb31ab898d0101cc785da1b54110c70d7ea06ad0259fb87631"


def check_panel(report):
    if (report["status"] != "three_checked_repeats_complete_unscored"
            or report["accuracy_evaluated"] is not False or report["job_id"] != "21328"
            or report["planned_repeats"] != 3
            or [row["index"] for row in report["repeats"]] != [0, 1, 2]
            or any(row["execution"]["exit_code"] != 0 for row in report["repeats"])):
        raise ValueError("Require exact complete unscored three-repeat panel 21328")


def reconstruct(payload):
    import numpy as np
    arrays = [np.load(payload / (name + ".npy"), mmap_mode="r", allow_pickle=False)
              for name in ("sources", "targets", "weights")]
    if (any(a.ndim != 1 or a.shape != arrays[0].shape for a in arrays)
            or [a.dtype for a in arrays] != [np.dtype("int32"), np.dtype("int32"), np.dtype("float64")]):
        raise ValueError("Invalid saved graph arrays")
    names = (payload / "gene_names.txt").read_text().splitlines()
    if not names or len(set(names)) != len(names) or any(not name for name in names):
        raise ValueError("Invalid saved gene names")
    oriented, endpoints, weights = hashlib.sha256(), hashlib.sha256(), hashlib.sha256()
    for start in range(0, len(arrays[0]), 65536):
        left, right, values = [a[start:start + 65536] for a in arrays]
        pairs = np.stack((left, right), axis=1)
        if np.any(pairs < 0) or np.any(pairs >= len(names)) or not np.isfinite(values).all():
            raise ValueError("Invalid endpoint or weight")
        oriented.update(pairs.tobytes(order="C"))
        canonical = np.stack((np.minimum(left, right), np.maximum(left, right)), axis=1).astype("<i8")
        endpoints.update(canonical.tobytes(order="C"))
        weights.update(values.astype("<f8").tobytes(order="C"))
    return ({"vertices": len(names), "edges": len(arrays[0]), "directed": False,
             "ordered_endpoints_sha256": endpoints.hexdigest(), "ordered_weights_sha256": weights.hexdigest()},
            oriented.hexdigest(), set(names))


def check_native(boundary, adapter, saved, oriented):
    if boundary["accuracy_evaluated"] is not False or len(boundary["calls"]) != 1:
        raise ValueError("Unexpected optimizer inventory")
    call = boundary["calls"][0]
    expected = {"initial_membership": None, "weights": "weight", "n_iterations": 2,
                "max_comm_size": 0, "seed": 4, "kwargs": {"resolution_parameter": .1},
                "partition_type": "leidenalg.VertexPartition.CPMVertexPartition"}
    if (call["status"] != "optimizer_returned" or call["arguments"] != expected
            or any(call[key] != saved for key in ("before", "saved", "after"))):
        raise ValueError("Native graph or optimizer settings differ")
    expected_constructor = {"status": "constructor_returned", "dtype": "int32",
                            "shape": [saved["edges"], 2], "n": saved["vertices"], "directed": False,
                            "ordered_input_bytes_sha256": oriented}
    if adapter != {"format": "python_pairs", "accuracy_evaluated": False, "calls": [expected_constructor]}:
        raise ValueError("Constructor record differs")


def admit(root, output, digest):
    if output.exists():
        raise FileExistsError(output)
    path = root / "benchmarks/results/qfo_checked_repeats_v1/results.json"
    report = read_frozen(path, digest)
    check_panel(report)
    accounting = subprocess.check_output(["sacct", "-j", "21328", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,MaxRSS"], text=True)
    scheduler = require_completed_job(accounting, 21328)
    prior_path = root / "benchmark_tools/results/qfo_constructor_formats_verified_20260916.json"
    prior = read_frozen(prior_path, ADMISSION_SHA)
    if report["admission"] != record(prior_path) or report["graph_inputs"] != prior["native_report"]["graph_inputs"]:
        raise ValueError("Prior admission or graph inputs changed")
    executor = root / "benchmarks/work/publication_qfo_checked_repeats_v1"
    revision = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    expected_revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", "f6ad87c^{commit}"], text=True).strip()
    if revision != expected_revision:
        raise ValueError("Changed frozen executor revision")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    for key, name in (("source", "run_qfo_checked_repeats.py"), ("worker_source", "checked_python_pair_worker.py")):
        if report[key] != record(executor / "benchmark_tools" / name):
            raise ValueError("Executor source record differs")
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    if verify(core, launcher, runtime_path) != report["runtime"]:
        raise ValueError("Changed native runtime")
    records = [record(path), report["source"], report["worker_source"], record(prior_path), *prior["provenance_checked"]]
    overrides = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                 "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    partitions, first = [], None
    for row in report["repeats"]:
        directory = path.parent / f"repeat_{row['index']}"
        payload = directory / "payload"
        filenames = ("worker_before.json", "native_boundary.json", "constructor_adapter.json", "checked_worker_provenance.json", "metadata.json")
        evidence = [record(payload / name) for name in filenames]
        if evidence != row["evidence"]:
            raise ValueError("Preserved worker evidence differs")
        worker, boundary, adapter, provenance, metadata = [json.loads((payload / name).read_text()) for name in filenames]
        if worker != row["worker"] or boundary != row["native_boundary"] or adapter != row["constructor_adapter"]:
            raise ValueError("Native files disagree with parent report")
        check_worker(worker, launcher, payload, overrides)
        if worker["inputs"] != report["graph_inputs"] or len(worker["cpu_affinity"]) != 1 or worker["metadata"] != metadata:
            raise ValueError("Worker inputs, affinity or metadata differ")
        helpers = [record(executor / "benchmark_tools" / name) for name in ("repeat_qfo_saved_graph.py", "probe_leiden_boundary.py")]
        if (provenance["source"] != report["worker_source"] or provenance["admission"] != report["admission"]
                or provenance["inputs"] != report["graph_inputs"] or provenance["helpers"] != helpers
                or worker["observer"] != helpers[0]):
            raise ValueError("Worker provenance differs")
        if first and any(worker[k] != first[k] for k in ("modules", "native_libraries", "python", "versions", "environment", "platform", "host", "cpu_affinity")):
            raise ValueError("Worker runtime changed across repeats")
        first = first or worker
        inputs = [record(payload / name) for name in ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")]
        if inputs != report["graph_inputs"]:
            raise ValueError("Saved payload differs")
        saved, oriented, universe = reconstruct(payload)
        if saved["vertices"] != 976504 or saved["edges"] != 24148515:
            raise ValueError("Unexpected QfO graph size")
        check_native(boundary, adapter, saved, oriented)
        execution = json.loads((directory / "execution.json").read_text())
        command = [worker["python"]["path"], report["worker_source"]["path"], "--root", str(root), "--payload", str(payload),
                   "--admission", str(prior_path), "--admission-sha256", ADMISSION_SHA]
        # sys.executable may name a symlink; its resolved executable must match the snapshot.
        actual = execution["command"]
        if execution != row["execution"] or not actual or [str(Path(actual[0]).resolve()), *actual[1:]] != command:
            raise ValueError("Worker command record differs")
        partition = directory / "orthohmm_working_res/orthohmm_edges_clustered.txt"
        if record(partition) != row["partition"]:
            raise ValueError("Partition bytes changed")
        comparison = compare_partition(partitions[0] if partitions else partition, partition, universe)
        if comparison != row["versus_first"]:
            raise ValueError("Recomputed partition comparison differs")
        partitions.append(partition)
        records.extend([*evidence, *helpers, *inputs, record(partition), record(directory / "execution.json"),
                        record(directory / "worker.log"), *worker["modules"].values(), *worker["native_libraries"], worker["python"]])
    comparisons = [{"left": i, "right": j, **compare_partition(partitions[i], partitions[j], universe)}
                   for i, j in itertools.combinations(range(3), 2)]
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting provenance records")
        unique[item["path"]] = item
        check(item)
    if verify(core, launcher, runtime_path) != report["runtime"]:
        raise ValueError("Native runtime changed during audit")
    result = {"status": "checked_repeats_verified", "publication_ready": False, "accuracy_evaluated": False,
              "scheduler": scheduler, "scheduler_accounting": accounting, "source": record(__file__),
              "source_report": record(path), "executor_commit": revision, "native_report": report,
              "provenance_checked": list(unique.values()), "pairwise_comparisons": comparisons,
              "all_three_partitions_equal": all(r["partition_equal"] for r in comparisons),
              "limitations": ["Bounded three-repeat initial-graph diagnostic; no general determinism or causal library diagnosis.",
                  "Preserved native observations, not a retrospective live-memory inspection.",
                  "No historical-equivalence claim, full HMM replay, accuracy selection, or controlled timing comparison."]}
    output.mkdir(parents=True)
    (output / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report-sha256", required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve(), args.report_sha256)
