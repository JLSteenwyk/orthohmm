"""Join complete paired VGNC outcomes to their admitted baseline graph context."""

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import resource
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.trace_deferred_reference_pairs import paired_transitions
from benchmark_tools.trace_qfo_graph_stages import GRAPHS, PARTITIONS

PROTOCOL = "hillclimb_deferred_graph_context_20261013.json"
EXPECTED_PAIRS = 23934


def exec_peak_rss_kib(status):
    fields = dict(line.split(":", 1) for line in status.splitlines() if ":" in line)
    if "VmHWM" not in fields:
        raise ValueError("Require Linux post-exec RSS high-water accounting")
    value = fields["VmHWM"].split()
    if len(value) != 2 or value[1] != "kB" or not value[0].isdigit():
        raise ValueError("Malformed Linux RSS high-water field")
    return int(value[0])


def diagnose(control, variant, baseline):
    paired = paired_transitions(control, variant)
    if len(baseline) != len(control):
        raise ValueError("Baseline graph context must cover the complete paired cohort")
    groups = {}
    for left, right, graph in zip(control, variant, baseline):
        if any(graph.get(key) != value for key, value in left.items()):
            raise ValueError("Baseline graph row differs from the exact control trace")
        for row in (left, right):
            allowed = ({"predicted"} if row["category"] == "TP" else
                       {"reference_within_species"} if row["reference_is_cross_species"] == "False" else
                       {"different_candidate_families", "within_candidate_not_predicted"})
            if row["stage"] not in allowed or row["initial_hit_directions"] not in {"0", "1", "2"}:
                raise ValueError("Malformed native outcome stage or hit count")
        for stage in GRAPHS:
            present, weight = graph[stage + "_edge"], graph[stage + "_edge_weight"]
            if present not in {"True", "False"}:
                raise ValueError("Malformed baseline graph boolean")
            if present == "True":
                if not math.isfinite(float(weight)) or float(weight) <= 0:
                    raise ValueError("Present baseline edge needs a finite positive weight")
            elif weight != "NA":
                raise ValueError("Absent baseline edge must have explicit NA weight")
        if any(graph["same_" + stage] not in {"True", "False"} for stage in PARTITIONS):
            raise ValueError("Malformed baseline partition boolean")
        if graph["same_root"] != left["same_root_hog"]:
            raise ValueError("Baseline RootHOG context differs from control")
        if left["reference_is_cross_species"] == "True" and (
                (graph["same_candidate"] == "True") != (left["stage"] != "different_candidate_families")):
            raise ValueError("Baseline candidate context differs from control")
        keys = ("transition:" + left["category"] + "->" + right["category"],
                "variant:" + right["category"] + "/" + right["stage"])
        for key in keys:
            group = groups.setdefault(key, {"pairs": 0, "initial_hit_counts": Counter(),
                "baseline_edges": Counter(), "baseline_partitions": Counter(), "baseline_partition_boundaries": Counter()})
            group["pairs"] += 1
            group["initial_hit_counts"][left["initial_hit_directions"]] += 1
            for stage in GRAPHS:
                group["baseline_edges"][stage + "/present_" + graph[stage + "_edge"]] += 1
            for stage in PARTITIONS:
                group["baseline_partitions"][stage + "/same_" + graph["same_" + stage]] += 1
            for a, b in zip(PARTITIONS, PARTITIONS[1:]):
                group["baseline_partition_boundaries"][a + ":" + graph["same_" + a] + "/" +
                                                       b + ":" + graph["same_" + b]] += 1
    return {"paired_transitions": paired, "strata": groups}


def run(source, graph_receipt_path, paired_receipt_path, output):
    if output.exists():
        raise FileExistsError("Require an untouched diagnostic output")
    if Path(__file__).resolve() != source / "benchmark_tools/diagnose_deferred_graph_context.py":
        raise ValueError("Require explicitly bound source root")
    started, cpu = time.monotonic(), time.process_time()
    rusage_at_start = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    checked = {}

    def consume(item):
        if item["path"] in checked and checked[item["path"]] != item:
            raise ValueError("Conflicting consumed input identities")
        check(item)
        checked[item["path"]] = item
        return Path(item["path"])

    graph_receipt = json.loads(consume(record(graph_receipt_path)).read_text())
    paired_receipt = json.loads(consume(record(paired_receipt_path)).read_text())
    if (graph_receipt["status"] != "complete_native_graph_stage_table_readback_verified" or
            paired_receipt["status"] != "independently_recounted_paired_tables_and_raw_vgnc_verified"):
        raise ValueError("Require independently admitted graph and paired diagnostic receipts")
    graph = json.loads(consume(graph_receipt["checked"][0]).read_text())
    paired = json.loads(consume(paired_receipt["checked"][0]).read_text())
    if graph["status"] != "admitted_native_graph_partition_stages_traced" or paired["status"] != "complete_paired_reference_localization_verified":
        raise ValueError("Require completed graph and paired producers")
    graph_plan = json.loads(consume(graph["checked"][0]).read_text())
    paired_plan = json.loads(consume(paired["execution_plan"]).read_text())
    if (tuple(arm["label"] for arm in paired["arms"]) != ("control", "deferred500") or
            tuple(arm["label"] for arm in paired_receipt["arms"]) != ("control", "deferred500") or
            tuple(case["label"] for case in paired_plan["cases"]) != ("control", "deferred500") or
            any(arm["status"] != "complete_reference_trace_verified" for arm in paired["arms"])):
        raise ValueError("Require both completed ordered native traces")
    baseline_path = consume(graph["table"])
    pair_paths = [consume(arm["table"]) for arm in paired["arms"]]
    consume(graph_plan["cohort"])
    consume(paired_plan["original_control_table"])
    if not (graph_plan["cohort"]["sha256"] == paired_plan["original_control_table"]["sha256"] == paired["arms"][0]["table"]["sha256"]):
        raise ValueError("Different original control cohorts")
    tables = []
    for path in [*pair_paths, baseline_path]:
        with path.open() as stream:
            tables.append(list(csv.DictReader(stream, delimiter="\t")))
    if any(len(rows) != EXPECTED_PAIRS for rows in tables) or graph["reference_pairs"] != EXPECTED_PAIRS:
        raise ValueError("Require the full asserted reference cohort")
    result = diagnose(*tables)
    if result["paired_transitions"] != paired["paired_transitions"]:
        raise ValueError("Paired accounting differs from admitted producer")
    for rows, arm, receipt_arm in zip(tables, paired["arms"], paired_receipt["arms"]):
        counts = Counter(row["category"] for row in rows)
        if (arm["label"] != receipt_arm["label"] or any(counts[k] != arm["raw_categories"]["counts"][k] or
                counts[k] != receipt_arm["raw_categories"][k] for k in ("TP", "FN"))):
            raise ValueError("Physical table accounting differs from native raw labels")
    if sum(row["reference_is_cross_species"] == "False" for row in tables[0]) != 2:
        raise ValueError("Both within-species reference assertions must remain explicit")
    protocol = json.loads(consume(record(source / "benchmark_tools" / PROTOCOL)).read_text())
    for module in (__name__, "benchmark_tools.prepare_ob_candidate_neighborhood", "benchmark_tools.run_membership_hillclimb",
                   "benchmark_tools.trace_deferred_reference_pairs", "benchmark_tools.trace_qfo_graph_stages"):
        consume(record(Path(sys.modules[module].__file__)))
    for item in checked.values():
        check(item)
    peak = exec_peak_rss_kib(Path("/proc/self/status").read_text())
    if peak > protocol["resource_limits"]["memory_gib"] * 1024 ** 2:
        raise ValueError("Post-exec diagnostic memory exceeds prespecified budget")
    if time.monotonic() - started > protocol["resource_limits"]["wall_limit_minutes"] * 60:
        raise ValueError("Diagnostic wall time exceeds prespecified budget")
    result.update(status="complete_paired_baseline_graph_context_verified", checked=list(checked.values()),
        source_revision=subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        resources={"wall_s": time.monotonic() - started, "cpu_s": time.process_time() - cpu,
                   "linux_post_exec_peak_rss_kib": peak, "rusage_maxrss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   "rusage_maxrss_at_diagnostic_start_kib": rusage_at_start,
                   "rss_scope": "Linux /proc/self/status VmHWM covers this exec, including imports. getrusage maxrss can retain pre-exec launch history; it is not the diagnostic's attributable peak."},
        limitations=protocol["limitations"], defaults_promoted=False, goal_completed=False)
    write_json(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source-root", "graph-receipt", "paired-receipt", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    result = run(args.source_root.resolve(), args.graph_receipt.resolve(), args.paired_receipt.resolve(), args.output.resolve())
    print(json.dumps({"status": result["status"], "paired_transitions": result["paired_transitions"]}))
