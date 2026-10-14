"""Isolate graph and singleton order effects using frozen accepted hits."""

import argparse
import gc
import hashlib
import importlib
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

PROTOCOL = "hillclimb_canonical_graph_probe_20261014.json"
REVISION = "402d32323f7166c6135f4c81b849b0abfc97e55c"


def selection_trace(n_genes, clusters, queries, targets, scores):
    cluster_by_gene = np.full(n_genes, -1, dtype=np.int32)
    singleton = np.ones(n_genes, dtype=bool)
    for number, members in enumerate(clusters):
        if len(members) > 1:
            cluster_by_gene[members] = number
            singleton[members] = False
    eligible = singleton[queries] & (cluster_by_gene[targets] >= 0) & np.isfinite(scores) & (scores > 0)
    rows = np.flatnonzero(eligible)
    best = np.full(n_genes, -np.inf, dtype=np.float64)
    np.maximum.at(best, queries[rows], scores[rows])
    winners = rows[scores[rows] == best[queries[rows]]]
    chosen = np.full(n_genes, -1, dtype=np.int32)
    q, first = np.unique(queries[winners], return_index=True)
    chosen[q] = cluster_by_gene[targets[winners[first]]]
    # Distinct winning clusters, not the number of tied rows or genes.
    packed = queries[winners].astype(np.int64) * len(clusters) + cluster_by_gene[targets[winners]]
    distinct = np.unique(packed)
    counts = np.bincount(distinct // len(clusters), minlength=n_genes).astype(np.int32) if len(clusters) else np.zeros(n_genes, dtype=np.int32)
    return best, chosen, counts


def trace_comparison(left, right):
    best_a, chosen_a, ties_a = left
    best_b, chosen_b, ties_b = right
    if not np.array_equal(best_a, best_b) or not np.array_equal(ties_a, ties_b):
        raise ValueError("Matched singleton score or equally-best cluster distributions differ")
    changed = np.flatnonzero(chosen_a != chosen_b)
    if np.any(ties_a[changed] < 2):
        raise ValueError("A changed singleton choice has no multi-cluster exact tie")
    return changed


def verify_trace_edges(n_genes, clusters, trace, edges):
    cluster_by_gene = np.full(n_genes, -1, dtype=np.int32)
    for number, members in enumerate(clusters):
        if len(members) > 1:
            cluster_by_gene[members] = number
    left = cluster_by_gene[edges.sources] < 0
    right = cluster_by_gene[edges.targets] < 0
    if not np.all(left ^ right):
        raise ValueError("Singleton output edge does not connect a singleton to an anchor")
    queries = np.where(left, edges.sources, edges.targets)
    anchors = np.where(left, edges.targets, edges.sources)
    if (not np.array_equal(np.unique(queries), np.flatnonzero(trace[1] >= 0))
            or np.any(trace[1][queries] != cluster_by_gene[anchors])):
        raise ValueError("Independent selected-cluster trace differs from actual singleton edges")


def bind_science(source):
    sys.path.insert(0, str(source))
    from orthohmm import accuracy, externals
    if (Path(inspect.getfile(accuracy)).resolve() != source / "orthohmm/accuracy.py"
            or Path(inspect.getfile(externals)).resolve() != source / "orthohmm/externals.py"):
        raise ValueError("Loaded graph implementation is not the frozen scientific source")
    return accuracy, externals


def save_edges(edges, folder):
    from benchmark_tools.prepare_ob_candidate_neighborhood import record
    folder.mkdir(parents=True, exist_ok=False)
    pins = {}
    for name, values in (("sources", edges.sources), ("targets", edges.targets), ("weights", edges.weights)):
        path = folder / (name + ".npy")
        np.save(path, values, allow_pickle=False)
        pins[name] = record(path)
    return {"count": len(edges.sources), "arrays": pins}


def same_arrays(left, right):
    return left["count"] == right["count"] and all(
        (left["arrays"][key]["bytes"], left["arrays"][key]["sha256"]) ==
        (right["arrays"][key]["bytes"], right["arrays"][key]["sha256"])
        for key in ("sources", "targets", "weights"))


def prepare(repo, source, output):
    science = repo / "benchmarks/work/canonical_downstream_r1_v1_20261014"
    bind_science(science)
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.run_full_search_hillclimb import runtime
    from benchmark_tools.run_membership_hillclimb import write_json
    root = repo / "benchmarks/results/canonical_downstream_r1_v1_20261014"
    parent_path = root / "execution_plan.json"
    parent = json.loads(parent_path.read_text())
    report_path = root / "results.json"
    report = json.loads(report_path.read_text())
    if (output.exists() or source == repo or parent["source_revision"] != REVISION
            or report["source_revision"] != REVISION or Path(parent["source_root"]).resolve() != science):
        raise ValueError("Require fresh output, frozen runner and original scientific plan")
    case = next(c for c in parent["cases"] if c["dataset"] == "QfO")
    ordered = next(r for r in report["arms"] if (r["dataset"], r["label"]) == ("QfO", "canonical"))["ordering"]
    if ordered["whole_block_multiset_and_internal_rows_preserved"] is not True:
        raise ValueError("Canonical checkpoint lacks the whole-block preservation gate")
    check(case["checkpoint_manifest"])
    check(ordered["manifest"])
    checkpoints = [case["checkpoint"], ordered["checkpoint"]]
    graphs = [root / "QfO" / label / "graph_metrics.json" for label in ("arrival", "canonical")]
    graph_reports = [json.loads(path.read_text()) for path in graphs]
    if graph_reports[0]["parameters"] != graph_reports[1]["parameters"]:
        raise ValueError("Original graph settings differ")
    if (graph_reports[0]["parameters"].get("cpm_resolution") != .1
            or graph_reports[0]["parameters"].get("leiden_seed") != 4):
        raise ValueError("Require original CPM-0.1, seed-4 graph settings")
    inputs = [record(parent_path), record(report_path), record(source / "benchmark_tools" / PROTOCOL),
              *parent["sources"], *parent["runtime"]["files"], *[record(p) for p in graphs]]
    inputs.extend(record(p) for checkpoint in checkpoints for p in sorted(Path(checkpoint).iterdir()) if p.is_file())
    stages = [next(s for s in r["stages"] if s["label"] == "multipass") for r in graph_reports]
    inputs.extend(stage["output"] for stage in stages)
    drivers = [record(source / "benchmark_tools/probe_canonical_initial_graph.py"),
               record(source / "benchmark_tools/probe_canonical_initial_graph.slurm")]
    for name in ("numpy", "igraph", "leidenalg"):
        package = Path(importlib.import_module(name).__file__).parent
        inputs.extend(record(p) for p in sorted(package.rglob("*.so")))
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    for pin in [*drivers, record(source / "benchmark_tools" / PROTOCOL)]:
        relative = Path(pin["path"]).relative_to(source)
        blob = subprocess.check_output(["git", "-C", str(source), "show", revision + ":" + str(relative)])
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Graph probe driver or protocol differs from committed revision")
    for pin in (*inputs, *drivers):
        check(pin)
    if runtime() != parent["runtime"]:
        raise ValueError("Original scientific runtime changed")
    output.mkdir(parents=True)
    plan = {"status": "prepared_graph_only_diagnostic", "repo": str(repo), "driver_source": str(source),
            "science": str(science), "output": str(output), "inputs": inputs, "drivers": drivers,
            "runtime": parent["runtime"], "checkpoints": checkpoints, "original_stages": stages,
            "expected_multipass_edges": [r["counts"]["multipass_edges"] for r in graph_reports],
            "original_rbnh_edges": [r["counts"]["rbnh_edges"] for r in graph_reports],
            "parameters": graph_reports[0]["parameters"], "scientific_revision": REVISION,
            "runner_revision": revision,
            "cpu": 1, "defaults_promoted": False, "goal_completed": False}
    write_json(output / "execution_plan.json", plan)
    return plan


def run(path):
    plan = json.loads(path.read_text())
    science, output = Path(plan["science"]), Path(plan["output"])
    accuracy, externals = bind_science(science)
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.run_full_search_hillclimb import runtime
    from benchmark_tools.run_membership_hillclimb import write_json
    from benchmark_tools.audit_historical_profile_ablation import read_partition
    from benchmark_tools.run_search_recovery import host_activity
    if Path(__file__).resolve() != Path(plan["driver_source"]) / "benchmark_tools/probe_canonical_initial_graph.py" or (output / "results.json").exists():
        raise ValueError("Require frozen driver and preserve prior attempts")
    checked = [record(path), *plan["inputs"], *plan["drivers"]]
    result = {"status": "checking_graph_probe_inputs", "execution_plan": record(path),
              "job_id": os.environ.get("SLURM_JOB_ID"), "host_before": host_activity(),
              "accuracy_evaluated": False, "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", result)
    try:
        for pin in checked:
            check(pin)
        if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
            raise ValueError("Runtime or CUDA visibility changed")
        os.environ["PYTHONPATH"] = str(science)
        os.chdir(science)
        data = [accuracy.load_accuracy_checkpoint(p, verify=True) for p in plan["checkpoints"]]
        if data[0][0] != data[1][0] or not np.array_equal(data[0][1], data[1][1]) or len(data[0][0]) != 984137:
            raise ValueError("Checkpoint gene/species tables differ")
        names = data[0][0]
        if any(len(values[-1]) != 90687327 or len(np.unique(values[1])) != 78 for values in data):
            raise ValueError("Wrong original corrected QfO checkpoint inventory")
        universe = set(names)
        result["status"] = "building_matched_rbnh_graphs"
        write_json(output / "results.json", result)
        graphs = []
        for label, (_, species, q, t, scores) in zip(("arrival", "canonical"), data):
            edges = accuracy.build_rbnh_edges(names, species, q, t, scores)
            graphs.append(save_edges(edges, output / label / "rbnh"))
            del edges
            gc.collect()
        if [g["count"] for g in graphs] != plan["original_rbnh_edges"]:
            raise ValueError("RBNH counts do not reproduce original inference")
        result.update(rbnh=graphs, rbnh_exact_arrays_equal=same_arrays(*graphs))
        if not result["rbnh_exact_arrays_equal"]:
            result["status"] = "rbnh_graph_difference_verified"
        else:
            from orthohmm.helpers import IndexedEdges
            arrays = [np.load(graphs[0]["arrays"][k]["path"], mmap_mode="r", allow_pickle=False) for k in ("sources", "targets", "weights")]
            rbnh = IndexedEdges(names, *arrays)
            initial_folder = output / "initial"
            externals.execute_leiden(.1, str(initial_folder), edges=rbnh, include_isolates=True, seed=4)
            initial_path = initial_folder / "orthohmm_working_res/orthohmm_edges_clustered.txt"
            clusters = accuracy.read_index_clusters(str(initial_path), names)
            read_partition(initial_path, universe)
            checked.append(record(initial_path))
            result["initial_clusters"] = len(clusters)
            result["status"] = "probing_fixed_partition_singleton_selection"
            write_json(output / "results.json", result)
            traces, singleton_records, combined_records, reproduced = [], [], [], []
            for label, (_, _, q, t, scores), original in zip(("arrival", "canonical"), data, plan["original_stages"]):
                trace = selection_trace(len(names), clusters, q, t, scores)
                traces.append(trace)
                edges = accuracy.build_singleton_assignment_edges(names, clusters, q, t, scores)
                verify_trace_edges(len(names), clusters, trace, edges)
                singleton_records.append(save_edges(edges, output / label / "singletons"))
                combined = accuracy.combine_edges(rbnh, edges)
                combined_records.append(save_edges(combined, output / label / "multipass"))
                folder = output / label / "clustered"
                externals.execute_leiden(.1, str(folder), edges=combined, include_isolates=True, seed=4)
                generated = folder / "orthohmm_working_res/orthohmm_edges_clustered.txt"
                a = {frozenset(g) for g in read_partition(generated, universe)}
                b = {frozenset(g) for g in read_partition(Path(original["output"]["path"]), universe)}
                reproduced.append(a == b)
                checked.append(record(generated))
                trace_folder = output / label / "trace"
                trace_folder.mkdir()
                for name, values in zip(("maximum_score", "selected_cluster", "tied_clusters"), trace):
                    p = trace_folder / (name + ".npy")
                    np.save(p, values, allow_pickle=False)
                    checked.append(record(p))
                del edges, combined, a, b
                gc.collect()
            changed = trace_comparison(*traces)
            result.update(singletons=singleton_records, multipass=combined_records,
                singleton_exact_arrays_equal=same_arrays(*singleton_records),
                generated_multipass_counts=[g["count"] for g in combined_records],
                original_multipass_partitions_reproduced=reproduced,
                changed_singleton_choices=len(changed),
                changed_choice_samples=[{"gene": names[int(i)], "maximum_score": float(traces[0][0][i]),
                    "equally_best_clusters": int(traces[0][2][i]), "arrival_cluster": int(traces[0][1][i]),
                    "canonical_cluster": int(traces[1][1][i])} for i in changed[:20]],
                every_changed_choice_has_multicluster_exact_tie=True)
            if result["generated_multipass_counts"] != plan["expected_multipass_edges"] or not all(reproduced):
                result["status"] = "fixed_partition_probe_does_not_reproduce_original_multipass"
            elif len(changed) and not result["singleton_exact_arrays_equal"]:
                result["status"] = "matched_graph_singleton_tie_diagnosis_reproduces_original_multipass"
            else:
                result["status"] = "matched_graph_probe_reproduced_without_singleton_choice_difference"
        checked.extend(pin for graph in graphs for pin in graph["arrays"].values())
        for key in ("singletons", "multipass"):
            checked.extend(pin for graph in result.get(key, []) for pin in graph["arrays"].values())
        for pin in checked:
            check(pin)
        if runtime() != plan["runtime"]:
            raise ValueError("Scientific runtime changed during graph diagnosis")
        result["checked"] = checked
    except Exception as error:
        result.update(status="failed_graph_diagnostic", error_type=type(error).__name__, error=str(error))
    result["host_after"] = host_activity()
    write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    for name in ("repo", "source", "output"):
        prep.add_argument("--" + name, type=Path, required=True)
    worker = sub.add_parser("run")
    worker.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve()) if args.command == "prepare"
              else run(args.execution_plan.resolve()))
    print(result["status"])
    raise SystemExit(1 if result["status"] == "failed_graph_diagnostic" else 0)
