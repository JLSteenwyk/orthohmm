"""Read back all five candidate partitions/traces without recalculating scores."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.link_factorial_scaling_resources import compare, partition
from benchmark_tools.diagnose_candidate_trace_variation import trace_comparison
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_ob_candidate_order_scores import LABELS
from benchmark_tools.replay_native_candidate_trace import replay


def collect(plan_ref, report_ref):
    inputs = {}

    def bind(ref):
        check(ref)
        inputs[ref["path"]] = ref
        return Path(ref["path"])

    plan = json.loads(bind(plan_ref).read_text())
    report = json.loads(bind(report_ref).read_text())
    if (plan.get("schema") != "native_candidate_fixed_seed_factorial_plan_v1"
            or plan.get("labels") != list(LABELS)
            or report.get("schema") != "native_candidate_fixed_seed_factorial_result_v1"
            or report.get("plan") != plan_ref or report.get("status") != "controls_reproduced"
            or report.get("accuracy_scored") is not False or report.get("frozen_method_modified") is not False
            or [r["label"] for r in report["rows"]] != list(LABELS)):
        raise ValueError("Require bound successful five-arm diagnostic, not benchmark accuracy")
    output = Path(plan["output"])
    if (output / "failure.json").exists():
        raise ValueError("Retained failure; no successful diagnostic readback")
    started = json.loads(bind(record(output / "started.json")).read_text())
    if started["plan"] != plan_ref or started["runtime"] != plan["runtime"]:
        raise ValueError("Startup bindings differ")
    for ref in plan["checked_records"]:
        bind(ref)
    seed = partition(bind(plan["seed"]), "space_separated_groups")
    groups, rows, traces = {}, [], {}
    for row in report["rows"]:
        label = row["label"]
        saved_row = json.loads(bind(record(output / (label + "_result.json"))).read_text())
        if saved_row != row:
            raise ValueError("Arm report differs from aggregate")
        predicted = partition(bind(row["partition"]), "space_separated_groups")
        trace_ref = record(output / (label + "_merges.json"))
        trace = json.loads(bind(trace_ref).read_text())
        for event in trace:
            if event["margin"] == "positive_infinity":
                event["margin"] = float("inf")
        if (predicted[0] != seed[0] or row["groups"] != len(predicted[1])
                or row["merges"] != len(trace) or len(seed[1]) - len(trace) != len(predicted[1])
                or replay(seed[1], trace) != predicted):
            raise ValueError("Accepted trace does not explain whole arm partition/counts")
        groups[label] = predicted
        traces[label] = trace
        rows.append({"label": label, "partition": row["partition"], "trace": trace_ref,
                     "genes": len(predicted[0]), "groups": len(predicted[1]),
                     "accepted_merges": len(trace), "whole_partition_reconstructed": True})
    comparisons = {a + "__vs__" + b: compare(groups[a], groups[b], set())
                   for i, a in enumerate(LABELS) for b in LABELS[i + 1:]}
    controls = {"historical": groups[LABELS[0]] == partition(bind(plan["original_prediction"]), "space_separated_groups"),
                "fresh_full": groups[LABELS[-1]] == partition(bind(plan["native_prediction"]), "space_separated_groups")}
    if controls != report["controls"] or not all(controls.values()) or comparisons != report["comparisons"]:
        raise ValueError("Complete comparisons or controls disagree")
    trace_comparisons = {label: trace_comparison(traces[LABELS[0]], traces[label],
                        plan["parameters"].get("max_satellites_per_anchor", 4)) for label in LABELS[1:]}
    for name in ("readback_native_candidate_factorial.py", "replay_native_candidate_trace.py",
                 "diagnose_candidate_trace_variation.py", "link_factorial_scaling_resources.py",
                 "prepare_ob_candidate_neighborhood.py", "score_ygob_groups.py"):
        bind(record(Path(__file__).with_name(name)))
    for ref in inputs.values():
        check(ref)
    return {"schema": "native_candidate_factorial_readback_v1", "plan": plan_ref, "report": report_ref,
            "rows": rows, "controls": controls, "comparisons": comparisons,
            "accepted_trace_comparisons_against_historical": trace_comparisons,
            "checked_records": sorted(inputs.values(), key=lambda r: r["path"]),
            "candidate_scores_recalculated": False, "accuracy_scored": False,
            "limitations": [
                "Complete partition/accepted-union consistency is checked; rejected candidates and score computation are not independently replayed.",
                "This readback pins retained merge traces now; their digests were not saved in the runner's per-arm rows.",
                "Fixed-seed interventions support boundary-specific effects, not universal determinism, full upstream historical causality or biological accuracy gains.",
                "No native inference restart, timing correction, job release or full runtime closure claim."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "report", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("plan-sha256", "report-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    refs = [record(args.plan), record(args.report)]
    if [ref["sha256"] for ref in refs] != [args.plan_sha256, args.report_sha256]:
        raise ValueError("Requested plan/report checksum differs")
    if args.output.exists():
        raise FileExistsError(args.output)
    value = collect(*refs)
    with args.output.open("x") as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(record(args.output), sort_keys=True))
