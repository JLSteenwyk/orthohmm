"""Trace admitted native graph and partition stages for the complete VGNC cohort."""

import argparse
from collections import Counter
import csv
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np

from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_leiden_boundary import saved_fingerprint
from benchmark_tools.score_ygob_groups import membership, read_predictions

PROTOCOL = "hillclimb_qfo_graph_localization_20261013.json"
PARTITIONS = ("initial", "multipass", "multipass_refined", "profile_base", "profile_expanded",
              "profiles_refined", "candidate", "root")
GRAPHS = ("initial", "multipass", "profile_base", "profile_expanded")
EXPECTED_GENES = 984137


def watched_graph(payload, names, wanted, expected, chunk_size=1000000):
    if chunk_size < 1:
        raise ValueError("Positive chunk size required")
    if (payload / "gene_names.txt").read_text().splitlines() != names:
        raise ValueError("Saved graph gene order differs from checkpoint")
    if saved_fingerprint(payload) != expected:
        raise ValueError("Saved graph differs from admitted native fingerprint")
    sources, targets, weights = [np.load(payload / (name + ".npy"), mmap_mode="r", allow_pickle=False)
                                 for name in ("sources", "targets", "weights")]
    if sources.dtype != np.dtype("int32") or targets.dtype != np.dtype("int32") or weights.dtype != np.dtype("float64"):
        raise ValueError("Unexpected saved graph dtypes")
    mask = np.zeros(len(names), dtype=bool)
    mask[list({i for pair in wanted for i in pair})] = True
    found, previous = {}, -1
    for start in range(0, len(sources), chunk_size):
        a, b, w = (values[start:start + chunk_size] for values in (sources, targets, weights))
        if (np.any(a < 0) or np.any(b >= len(names)) or np.any(a >= b)
                or not np.isfinite(w).all() or np.any(w <= 0)):
            raise ValueError("Invalid saved graph endpoints or weights")
        packed = a.astype(np.int64) * len(names) + b
        if len(packed) and (packed[0] <= previous or np.any(packed[1:] <= packed[:-1])):
            raise ValueError("Duplicate or unordered saved graph edges")
        if len(packed):
            previous = int(packed[-1])
        selected = mask[a] & mask[b]
        for i, j, weight in zip(a[selected], b[selected], w[selected]):
            pair = int(i), int(j)
            if pair in wanted:
                found[pair] = float(weight)
    return found


def partition_memberships(path, universe, wanted, expected_groups, root=False):
    if root:
        groups = read_predictions(path, "root_hogs")
        full = membership(groups)
        if set(full) != universe:
            raise ValueError("Root-HOG gene universe differs")
        selected = {gene: full[gene] for gene in wanted}
    else:
        groups = read_partition(path, universe)
        selected = {gene: i for i, group in enumerate(groups) for gene in group if gene in wanted}
    if len(groups) != expected_groups or set(selected) != wanted:
        raise ValueError("Partition group count or watched universe differs")
    return selected


def extend_rows(rows, index, graphs, partitions):
    result = []
    seen = set()
    for row in rows:
        pair = row["left"], row["right"]
        if pair[0] >= pair[1] or pair in seen or row["category"] not in ("TP", "FN"):
            raise ValueError("Invalid retained cohort row")
        seen.add(pair)
        if row["reference_is_cross_species"] not in ("True", "False") or row["same_root_hog"] not in ("True", "False"):
            raise ValueError("Malformed retained cohort boolean")
        left, right = row["native_gene_left"], row["native_gene_right"]
        watched = tuple(sorted((index[left], index[right])))
        extended = dict(row)
        for name in GRAPHS:
            weight = graphs[name].get(watched)
            if weight is not None and (not math.isfinite(weight) or weight <= 0):
                raise ValueError("Invalid watched graph weight")
            extended[name + "_edge"] = weight is not None
            extended[name + "_edge_weight"] = "NA" if weight is None else weight
        for name in PARTITIONS:
            extended["same_" + name] = partitions[name][left] == partitions[name][right]
        if extended["same_root"] != (row["same_root_hog"] == "True"):
            raise ValueError("Root membership differs from retained trace")
        cross = row["reference_is_cross_species"] == "True"
        same_candidate = row["category"] == "TP" or row["stage"] == "within_candidate_not_predicted"
        if cross and extended["same_candidate"] != same_candidate:
            raise ValueError("Candidate membership differs from retained trace")
        if not cross and (row["category"] != "FN" or row["stage"] != "reference_within_species"):
            raise ValueError("Within-species reference eligibility differs")
        result.append(extended)
    return result


def accounting(rows):
    stages, transitions = Counter(), Counter()
    for row in rows:
        cohort = row["category"] if row["reference_is_cross_species"] == "True" else "within_species_FN"
        for name in PARTITIONS:
            stages[f'{cohort}/{name}/same_{row["same_" + name]}'] += 1
        for left, right in zip(PARTITIONS, PARTITIONS[1:]):
            transitions[f'{cohort}/{left}:{row["same_" + left]}/{right}:{row["same_" + right]}'] += 1
        for name in GRAPHS:
            stages[f'{cohort}/{name}/edge_{row[name + "_edge"]}'] += 1
    return {"stages": dict(stages), "transitions": dict(transitions)}


def prepare(repo, trace_path, output):
    if output.exists():
        raise FileExistsError(output)
    trace = json.loads(trace_path.read_text())
    if trace["status"] != "native_reference_pair_stage_trace_verified" or trace["reference_pairs"] != 23934:
        raise ValueError("Need complete verified stage trace")
    check(trace["checked"][0])
    previous = json.loads(Path(trace["checked"][0]["path"]).read_text())
    admission_path = repo / "benchmark_tools/results/qfo_corrected_replay_admission_21757.json"
    admission = json.loads(admission_path.read_text())
    factorial = json.loads(Path(previous["inputs"]["factorial"]["path"]).read_text())
    if admission["status"] != "corrected_checked_replay_admitted" or admission["plan"] != factorial["plan"]:
        raise ValueError("Saved graph ancestry differs from candidate factorial")
    checked = [record(trace_path), record(admission_path), trace["table"], *trace["checked"],
               record(repo / "benchmark_tools" / PROTOCOL), admission["plan"]]
    graphs, partitions = {}, {}
    for stage in admission["clustering"]:
        label = stage["stage"]
        if label not in GRAPHS or label in graphs:
            raise ValueError("Unexpected native graph inventory")
        parent = Path(stage["partition"]["path"]).parent
        manifest_record = next(item for item in admission["checked_records"]
                               if item["path"] == str(parent / "payload_manifest.json"))
        check(manifest_record)
        manifest = json.loads(Path(manifest_record["path"]).read_text())
        if manifest["stage"] != label:
            raise ValueError("Payload stage label differs")
        checked.extend([manifest_record, *manifest["inputs"], stage["partition"]])
        graphs[label] = {"payload": str(parent / "payload"), "fingerprint": stage["saved_graph"]}
        partitions[label] = {"file": stage["partition"], "groups": stage["groups"]}
    for label, old in (("multipass_refined", "multipass_refined"), ("profiles_refined", "strict_profiles_refined")):
        stage = next(item for item in admission["coverage"] if item["label"] == old)
        partitions[label] = {"file": stage["output"], "groups": stage["groups"]}
        checked.append(stage["output"])
    partitions["candidate"] = {"file": previous["inputs"]["candidates"], "groups": 351739}
    partitions["root"] = {"file": previous["inputs"]["roots"], "groups": 366068}
    if set(graphs) != set(GRAPHS) or set(partitions) != set(PARTITIONS):
        raise ValueError("Incomplete stage inventory")
    sources = [record(module.__file__) for module in (sys.modules[__name__],
        sys.modules["benchmark_tools.audit_historical_profile_ablation"],
        sys.modules["benchmark_tools.prepare_ob_candidate_neighborhood"],
        sys.modules["benchmark_tools.probe_leiden_boundary"], sys.modules["benchmark_tools.score_ygob_groups"])]
    checked = list({item["path"]: item for item in checked}.values())
    for item in [*checked, *sources]:
        check(item)
    plan = {"status": "prepared_unrun", "checked": checked, "sources": sources,
        "source_revision": subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
        "cohort": trace["table"], "checkpoint_names": previous["inputs"]["checkpoint_gene_names.txt"],
        "graphs": graphs, "partitions": partitions, "expected_reference_pairs": 23934,
        "expected_raw_counts": trace["raw_categories"]["counts"]}
    with output.open("x") as stream:
        json.dump(plan, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return plan


def run(plan_path, output):
    if output.exists():
        raise FileExistsError(output)
    plan = json.loads(plan_path.read_text())
    checked = [record(plan_path), *plan["checked"], *plan["sources"]]
    for item in checked:
        check(item)
    started, cpu = time.monotonic(), time.process_time()
    names = Path(plan["checkpoint_names"]["path"]).read_text().splitlines()
    if len(names) != EXPECTED_GENES or len(set(names)) != len(names):
        raise ValueError("Unexpected checkpoint gene universe")
    index = {gene: i for i, gene in enumerate(names)}
    with Path(plan["cohort"]["path"]).open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    counts = Counter(row["category"] for row in rows)
    if len(rows) != plan["expected_reference_pairs"] or any(counts[label] != plan["expected_raw_counts"][label] for label in ("TP", "FN")):
        raise ValueError("Reference cohort or native labels differ")
    wanted_genes = {row[key] for row in rows for key in ("native_gene_left", "native_gene_right")}
    wanted_pairs = {tuple(sorted((index[row["native_gene_left"]], index[row["native_gene_right"]]))) for row in rows}
    graphs = {name: watched_graph(Path(case["payload"]), names, wanted_pairs, case["fingerprint"])
              for name, case in plan["graphs"].items()}
    partitions = {name: partition_memberships(Path(case["file"]["path"]), set(names), wanted_genes,
                                             case["groups"], name == "root")
                  for name, case in plan["partitions"].items()}
    extended = extend_rows(rows, index, graphs, partitions)
    totals = accounting(extended)
    for item in checked:
        check(item)
    output.mkdir(parents=True, exist_ok=False)
    table = output / "reference_graph_stages.tsv"
    with table.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(extended[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(extended)
    result = {"status": "admitted_native_graph_partition_stages_traced", "checked": checked,
        "source_revision": plan["source_revision"], "reference_pairs": len(extended), "counts": dict(counts),
        "accounting": totals, "table": record(table), "job_id": os.environ.get("SLURM_JOB_ID"),
        "resources": {"wall_s": time.monotonic() - started, "cpu_s": time.process_time() - cpu,
                      "process_peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
        "defaults_changed": False, "goal_completed": False,
        "limitations": ["Development-exposed saved graph/partition trace; no fresh search or inference was run.",
            "Graph edges and co-membership are not biological orthology truth or evidence that all split pairs should merge.",
            "Stage changes localize separation, not all threshold, optimization, or refinement decisions that caused it.",
            "Within-species asserted rows and native pair labels remain unchanged.",
            "Shared-host diagnostic timing is not inference timing."]}
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="operation", required=True)
    pre = commands.add_parser("prepare")
    pre.add_argument("--repo", required=True, type=Path)
    pre.add_argument("--trace-report", required=True, type=Path)
    pre.add_argument("--output", required=True, type=Path)
    execution = commands.add_parser("run")
    execution.add_argument("--plan", required=True, type=Path)
    execution.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.operation == "prepare":
        prepare(args.repo.resolve(), args.trace_report.resolve(), args.output.resolve())
    else:
        result = run(args.plan.resolve(), args.output.resolve())
        print(json.dumps(result["accounting"], sort_keys=True))


if __name__ == "__main__":
    main()
