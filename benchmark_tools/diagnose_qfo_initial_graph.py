"""Localize QfO replay drift without rerunning search or profile expansion."""

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import types

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.run_qfo_publication_replay import compare_partition
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.verify_qfo_replay_launcher import verify

OLD = "694a77fe56167754ca949751bca88aa7d11353dc"
CORE = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
CHECKPOINT_SHA = "b90c787f050a087adeeb81b1cace18c9cbb30e1724926a40a1df6d7e49bd549c"


def load_revision(root, revision, suffix):
    source = subprocess.check_output(["git", "-C", str(root), "show", revision + ":orthohmm/accuracy.py"])
    name = "orthohmm._diagnostic_accuracy_" + suffix
    module = types.ModuleType(name)
    module.__package__ = "orthohmm"
    sys.modules[name] = module
    exec(compile(source, revision + ":orthohmm/accuracy.py", "exec"), module.__dict__)
    return module, {"revision": revision, "path": "orthohmm/accuracy.py",
                    "sha256": hashlib.sha256(source).hexdigest()}


def array_evidence(array):
    value = np.ascontiguousarray(array)
    return {"shape": list(value.shape), "dtype": str(value.dtype),
            "sha256": hashlib.sha256(memoryview(value).cast("B")).hexdigest()}


def compare_edges(first, second):
    fields = ("sources", "targets", "weights")
    arrays = {key: {"historical": array_evidence(getattr(first, key)),
                    "frozen": array_evidence(getattr(second, key))} for key in fields}
    names_equal = list(first.gene_names) == list(second.gene_names)
    return {"gene_order_equal": names_equal, "arrays": arrays,
            "byte_equal": names_equal and all(v["historical"] == v["frozen"] for v in arrays.values())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    frozen = root / "benchmarks/work/publication_method_native_v2"
    runtime = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    before = verify(frozen, launcher, runtime)
    # Bind imported graph helpers and isolated workers to the verified launcher.
    sys.path.insert(0, str(launcher))
    os.chdir(launcher)
    os.environ["PYTHONPATH"] = str(launcher)
    from orthohmm import helpers, externals
    for module, relative in ((helpers, "orthohmm/helpers.py"), (externals, "orthohmm/externals.py")):
        expected = subprocess.check_output(["git", "-C", str(root), "show", CORE + ":" + relative])
        if Path(module.__file__).read_bytes() != expected:
            raise ValueError("Imported graph helper differs from frozen source")
        if subprocess.check_output(["git", "-C", str(root), "show", OLD + ":" + relative]) != expected:
            raise ValueError("Historical graph helper differs; comparison needs separate environments")
    historical, old_record = load_revision(root, OLD, "historical")
    current, current_record = load_revision(root, CORE, "frozen")
    checkpoint = root / "qfo_benchmark/results/orthohmm_high_sensitivity_isolated/output/orthohmm_working_res/high_sensitivity_checkpoint"
    manifest = file_provenance(checkpoint / "manifest.json")
    if manifest["sha256"] != CHECKPOINT_SHA:
        raise ValueError("Changed QfO checkpoint")
    names, species, q, t, scores = current.load_accuracy_checkpoint(checkpoint, verify=True)
    metrics_path = root / "qfo_benchmark/results/orthohmm_high_sensitivity_isolated/metrics.json"
    metrics_record = file_provenance(metrics_path)
    if metrics_record["sha256"] != "fb6b8d7e6824e1801c7cdb794e319d30ee9c16d9eada9213ef345b8de54d8893":
        raise ValueError("Historical metrics changed")
    old_counts = json.loads(metrics_path.read_text())["counts"]
    output.mkdir(parents=True)
    result = {"status": "running", "accuracy_evaluated": False,
              "job_id": os.environ.get("SLURM_JOB_ID"), "source": file_provenance(Path(__file__)),
              "checkpoint": manifest, "historical_metrics": metrics_record,
              "accuracy_sources": [old_record, current_record], "runtime": before,
              "packages": {name: importlib.metadata.version(name) for name in ("numpy", "igraph", "leidenalg")},
              "repeats": [], "limitations": [
                  "Historical intermediate partitions and graph hashes were not retained.",
                  "Historical package/native versions are not proven by matching Python source.",
                  "This diagnostic does not rerun search, profiles, or benchmark accuracy.",
                  "Equal edge counts alone do not establish graph equality."]}
    try:
        old_edges = historical.build_rbnh_edges(names, species, q, t, scores)
        edges = current.build_rbnh_edges(names, species, q, t, scores)
        result["rbnh_comparison"] = compare_edges(old_edges, edges)
        result["historical_recorded_rbnh_edges"] = old_counts["high_sensitivity_rbnh_edges"]
        result["recomputed_rbnh_edges"] = len(edges)
        del old_edges
        if not result["rbnh_comparison"]["byte_equal"]:
            result["status"] = "rbnh_source_difference_identified"
            return
        for key in ("sources", "targets", "weights"):
            np.save(output / ("rbnh_" + key + ".npy"), getattr(edges, key), allow_pickle=False)
        result["graph_arrays"] = [file_provenance(output / ("rbnh_" + key + ".npy")) for key in ("sources", "targets", "weights")]
        for index in range(2):
            directory = output / f"repeat_{index}"
            (directory / "orthohmm_working_res").mkdir(parents=True)
            externals.execute_leiden(.1, str(directory), edges=edges, include_isolates=True, seed=4)
            partition = directory / "orthohmm_working_res/orthohmm_edges_clustered.txt"
            clusters = current.read_index_clusters(str(partition), names)
            singleton = current.build_singleton_assignment_edges(names, clusters, q, t, scores)
            result["repeats"].append({"index": index, "partition": file_provenance(partition),
                "groups": len(clusters), "singleton_edges": len(singleton),
                "singleton_arrays": {key: array_evidence(getattr(singleton, key)) for key in ("sources", "targets", "weights")},
                "matches_historical_singleton_count": len(singleton) == old_counts["high_sensitivity_singleton_edges"]})
        result["repeat_partition_comparison"] = compare_partition(
            Path(result["repeats"][0]["partition"]["path"]), Path(result["repeats"][1]["partition"]["path"]), set(names))
        result["status"] = "diagnostic_complete"
        if verify(frozen, launcher, runtime) != before:
            raise ValueError("Runtime changed")
        verify_file(metrics_path, metrics_record)
        verify_file(checkpoint / "manifest.json", manifest)
    except Exception as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "diagnostic.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
