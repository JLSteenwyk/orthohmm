"""Locate admitted functional pair changes within label-independent seed structure."""

import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.compare_qfo_scored_pairs import read_scores
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from qfo_benchmark.og_to_pairwise import _strip_to_uniprot


def partition_membership(path, names):
    groups = read_partition(path, set(names))
    lookup = {gene: index for index, gene in enumerate(names)}
    membership = np.full(len(names), -1, dtype=np.int32)
    sizes = np.array([len(group) for group in groups], dtype=np.int32)
    for index, group in enumerate(groups):
        indices = np.fromiter((lookup[gene] for gene in group), dtype=np.int64)
        membership[indices] = index
    if np.any(membership < 0):
        raise ValueError("Structural partition lost genes")
    return membership, sizes


def direct_support(wanted, sources, targets, weights, n_genes, chunk_size=1_000_000):
    arrays = [np.asarray(a) for a in (wanted, sources, targets, weights)]
    wanted, sources, targets, weights = arrays
    if (type(n_genes) is not int or not 0 < n_genes <= 3_037_000_499
            or type(chunk_size) is not int or chunk_size <= 0
            or any(a.ndim != 1 for a in arrays)
            or any(not np.issubdtype(a.dtype, np.integer) for a in arrays[:3])
            or not np.issubdtype(weights.dtype, np.floating)
            or len(sources) != len(targets) or len(sources) != len(weights)):
        raise ValueError("Invalid direct graph arrays or indexing limits")
    if len(wanted) and (np.any(wanted < 0) or np.any(wanted >= n_genes**2) or np.any(wanted[1:] <= wanted[:-1])):
        raise ValueError("Watched keys must be unique sorted valid integer pairs")
    support = np.zeros(len(wanted), dtype=np.float64)
    for offset in range(0, len(sources), chunk_size):
        a = sources[offset:offset + chunk_size].astype(np.int64)
        b = targets[offset:offset + chunk_size].astype(np.int64)
        values = weights[offset:offset + chunk_size]
        if (np.any(a < 0) or np.any(b < 0) or np.any(a >= n_genes) or np.any(b >= n_genes)
                or np.any(~np.isfinite(values)) or np.any(values <= 0)):
            raise ValueError("Invalid graph endpoint or nonpositive/nonfinite score")
        if not len(wanted):
            continue
        # Widen endpoints before multiplication; the full QfO key exceeds int32.
        keys = np.minimum(a, b) * n_genes + np.maximum(a, b)
        positions = np.searchsorted(wanted, keys)
        matches = (positions < len(wanted)) & (wanted[np.minimum(positions, len(wanted) - 1)] == keys)
        np.maximum.at(support, positions[matches], values[matches])
    return support


def feature_key(a, b, parent, seed, candidate, retained, score):
    return (bool(parent[a] == parent[b]), bool(seed[a] == seed[b]), bool(candidate[a] == candidate[b]),
            int(retained[a]) + int(retained[b]),
            "absent" if score == 0 else "at_least_1_5" if score >= 1.5 else "positive_below_1_5")


def run(repo, output):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Require new CPU-only descriptive diagnostic namespace")
    admission_path = repo / "benchmarks/results/component_recovery_admission_r1_v1_20261013/report.json"
    admission = json.loads(admission_path.read_text())
    if admission["status"] != "terminal_component_candidates_and_raw_scores_independently_verified":
        raise ValueError("Component result has not been independently admitted")
    scheduler = subprocess.check_output(["sacct", "-j", "24222", "-X",
        "--format=JobID,State,ExitCode", "--parsable2"], text=True)
    if scheduler.splitlines() != ["JobID|State|ExitCode", "24222|COMPLETED|0:0"]:
        raise ValueError("Independent component admission is not successfully terminal")
    plan_path = Path(admission["execution_plan"]["path"])
    check(admission["execution_plan"])
    plan = json.loads(plan_path.read_text())
    producer_path = Path(admission["original_failed_producer"]["path"])
    check(admission["original_failed_producer"])
    producer = json.loads(producer_path.read_text())
    table_path = repo / "benchmark_tools/results/hillclimb_comparison_readback_v11_20261013/summary.json"
    table = json.loads(table_path.read_text())
    control = next(c for m in table["methods"] if m["key"] == "copy_branch:control" for c in m["cells"] if c["dataset"] == "QfO")
    gene_path = Path(plan["data"]["checkpoint"]) / "gene_names.txt"
    original_pins = {Path(pin["path"]): pin for pin in plan["inputs"]}
    check(original_pins[gene_path])
    names = gene_path.read_text().splitlines()
    if len(names) != 984137 or names != sorted(set(names)):
        raise ValueError("Original checkpoint gene order or full universe changed")
    accessions = {_strip_to_uniprot(gene): index for index, gene in enumerate(names)}
    if len(accessions) != len(names):
        raise ValueError("Non-injective original gene accession normalization")
    row = producer["structural_controls"]["components500"]
    checked = [record(admission_path), admission["execution_plan"], admission["original_failed_producer"],
               record(table_path), original_pins[gene_path], record(Path(__file__)),
               record(Path(__file__).with_name("hillclimb_component_pair_structure_20261014.json")),
               plan["data"]["partition"], row["seed"], row["copy_split_trace"], row["candidate"]["candidate_partition"]]
    for pin in checked:
        check(pin)
    parent, _ = partition_membership(Path(plan["data"]["partition"]["path"]), names)
    seed, seed_sizes = partition_membership(Path(row["seed"]["path"]), names)
    candidate, _ = partition_membership(Path(row["candidate"]["candidate_partition"]["path"]), names)
    trace = json.loads(Path(row["copy_split_trace"]["path"]).read_text())
    eligible = np.zeros(len(names), dtype=bool)
    for event in trace["events"]:
        if event["large_copy_eligible"]:
            members = np.asarray(event["members"], dtype=np.int64)
            if np.any(members < 0) or np.any(members >= len(names)) or len(np.unique(members)) != len(members) or np.any(eligible[members]):
                raise ValueError("Invalid executed large-copy trace membership")
            eligible[members] = True
    retained = eligible & (seed_sizes[seed] > 1) & (seed_sizes[seed] <= 500)
    if int(eligible.sum()) != admission["component_counts"]["genes"] or int(retained.sum()) != admission["component_counts"]["bounded_component_genes"]:
        raise ValueError("Structural activation differs from admitted component counts")
    folders = [repo / "benchmarks/results/deferred_qfo_inferred_r1_v2_20261013/control/qfo_scoring",
               repo / "benchmarks/results/component_scoring_recovery_r1_v1_20261013/qfo_scoring"]
    values = {}
    functional_receipt_path = repo / "benchmark_tools/results/large_copy_components_functional_diagnostic_receipt_20261013.json"
    functional_receipt = json.loads(functional_receipt_path.read_text())
    check(functional_receipt["diagnostic"])
    functional = json.loads(Path(functional_receipt["diagnostic"]["path"]).read_text())
    functional_pins = {Path(pin["path"]): pin for pin in functional["evidence"]}
    checked.extend((record(functional_receipt_path), functional_receipt["diagnostic"]))
    for metric in ("GO", "EC"):
        values[metric] = []
        for cell, folder in zip((control, admission), folders):
            paths = list((folder / "results" / metric).glob("*raw.txt.gz"))
            if len(paths) != 1:
                raise ValueError("Ambiguous retained functional scores")
            check(functional_pins[paths[0]])
            checked.append(functional_pins[paths[0]])
            scores = read_scores(paths[0], metric)
            mean = sum(scores.values()) / len(scores) / 1e6
            if len(scores) != cell["details"][metric]["assessed_relations"] or abs(mean - cell["scores"]["QfO: " + metric]) > 5.051e-7:
                raise ValueError("Raw functional values differ from admitted endpoint")
            values[metric].append(scores)
    watched = set().union(*(set(scores) for pair in values.values() for scores in pair))
    pair_indices = {pair: tuple(sorted((accessions[pair[0]], accessions[pair[1]]))) for pair in watched}
    keys = np.array(sorted(a * len(names) + b for a, b in pair_indices.values()), dtype=np.int64)
    payload = Path(plan["data"]["graph"]["payload"])
    graph_names = payload / "gene_names.txt"
    if graph_names.read_text().splitlines() != names:
        raise ValueError("Graph and checkpoint gene orders differ")
    graph_paths = [payload / (name + ".npy") for name in ("sources", "targets", "weights")]
    graph_pins = [original_pins[p] for p in (graph_names, *graph_paths)]
    for pin in graph_pins:
        check(pin)
    checked.extend(graph_pins)
    arrays = [np.load(p, mmap_mode="r", allow_pickle=False) for p in graph_paths]
    if len(arrays[0]) != plan["data"]["graph"]["fingerprint"]["edges"]:
        raise ValueError("Original graph edge count changed")
    support = direct_support(keys, *arrays, len(names))
    scores_by_key = dict(zip(keys.tolist(), support.tolist()))
    rows = []
    for metric, (left, right) in values.items():
        buckets = defaultdict(lambda: [0, 0])
        for pair in left.keys() | right.keys():
            if pair in left and pair in right:
                if left[pair] != right[pair]:
                    raise ValueError("Shared functional annotation score changed")
                category, value = "shared", right[pair]
            elif pair in right:
                category, value = "added", right[pair]
            else:
                category, value = "removed", left[pair]
            a, b = pair_indices[pair]
            features = feature_key(a, b, parent, seed, candidate, retained, scores_by_key[a * len(names) + b])
            if category != "removed" and not features[2]:
                raise ValueError("Admitted pair crosses final native candidate families")
            bucket = buckets[(category, *features)]
            bucket[0] += 1
            bucket[1] += value
        for key, (count, total) in sorted(buckets.items()):
            category, same_parent, same_seed, same_candidate, endpoints, evidence = key
            rows.append({"metric": metric, "category": category, "within_parent_graph_family": same_parent,
                "within_component_seed": same_seed, "within_final_candidate": same_candidate,
                "retained_component_endpoints": endpoints, "direct_graph_support": evidence,
                "scored_relations": count, "score_sum_millionths": total, "serialized_mean": total / count / 1e6})
    for pin in checked:
        check(pin)
    fields = dict(line.split(":", 1) for line in Path("/proc/self/status").read_text().splitlines() if ":" in line)
    peak = int(fields["VmHWM"].split()[0])
    if peak > 16 * 1024**2:
        raise MemoryError("Prespecified post-exec memory budget exceeded")
    result = {"status": "admitted_functional_pair_changes_structurally_localized", "checked": checked,
        "job_id": os.environ.get("SLURM_JOB_ID"),
        "genes": len(names), "watched_pairs": len(watched), "eligible_large_copy_genes": int(eligible.sum()),
        "retained_component_genes": int(retained.sum()), "rows": rows, "linux_postexec_peak_rss_kib": peak,
        "scope": "Development-exposed descriptive structure; direct HMM graph support is homology evidence, not orthology truth or causal attribution.",
        "native_inference_rerun": False, "new_scoring": False, "defaults_promoted": False, "goal_completed": False}
    write_json(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.output.resolve())["status"])
