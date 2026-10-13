"""Reuse the tested VGNC stage tracer for independently admitted deferral cells."""

import argparse
from collections import Counter
from copy import deepcopy
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import trace_qfo_reference_stages as tracing
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import PRODUCTS, argument, write_json

PROTOCOL = "hillclimb_deferred_pair_localization_20261013.json"
LABELS = ("control", "deferred500")
IMMUTABLE_FIELDS = ("left", "right", "reference_label", "native_gene_left", "native_gene_right",
    "reference_is_cross_species", "initial_hit_directions", "forward_normalized_score", "reverse_normalized_score")


def paired_transitions(left, right):
    if len(left) != len(right) or not left:
        raise ValueError("Require complete paired reference cohorts")
    counts, stages, seen = Counter(), Counter(), set()
    for a, b in zip(left, right):
        pair = a["left"], a["right"]
        if pair[0] >= pair[1] or pair in seen or any(a[key] != b[key] for key in IMMUTABLE_FIELDS):
            raise ValueError("Reference cohort, aliases, eligibility, or original HMM hits changed")
        seen.add(pair)
        if any(row["category"] not in {"TP", "FN"} for row in (a, b)):
            raise ValueError("Unexpected asserted-pair category")
        counts[a["category"] + "->" + b["category"]] += 1
        stages[a["stage"] + "->" + b["stage"]] += 1
        if a["reference_is_cross_species"] not in {"True", "False"}:
            raise ValueError("Malformed reference eligibility")
        if a["reference_is_cross_species"] == "False" and any(
                row["category"] != "FN" or row["stage"] != "reference_within_species" for row in (a, b)):
            raise ValueError("Within-species reference assertions must remain explicit native misses")
    return {"asserted_pairs": len(left), "categories": dict(counts), "stages": dict(stages),
            "rescued_asserted_pairs": counts["FN->TP"], "lost_asserted_pairs": counts["TP->FN"],
            "net_recalled_pairs": counts["FN->TP"] - counts["TP->FN"]}


def prepare(repo, source, receipt_path, previous_plan_path, previous_trace_path, output):
    if output.exists():
        raise ValueError("Require a new diagnostic namespace")
    if Path(__file__).resolve() != source / "benchmark_tools/trace_deferred_reference_pairs.py":
        raise ValueError("Prepare from frozen wrapper source")
    receipt = json.loads(receipt_path.read_text())
    check(receipt["full_report"])
    full = json.loads(Path(receipt["full_report"]["path"]).read_text())
    cells = full["cells"]
    if (tuple(c["method_key"] for c in cells) != tuple("deferred_inferred:" + label for label in LABELS)
            or any(c["status"] != "independent_content_readback_verified" or c["dataset"] != "QfO"
                   or c["coverage"]["genes"] != 984137 or c["coverage"]["species"] != 78 for c in cells)):
        raise ValueError("Require both independently admitted full-universe own-tree QfO cells")
    previous = json.loads(previous_plan_path.read_text())
    trace = json.loads(previous_trace_path.read_text())
    if (trace["status"] != "native_reference_pair_stage_trace_verified" or trace["reference_pairs"] != 23934
            or trace["within_species_reference_pairs_preserved"] != 2
            or trace["checked"][0] != record(previous_plan_path)):
        raise ValueError("Require the original complete, alias-corrected reference trace")
    checked = [record(receipt_path), receipt["full_report"], record(previous_plan_path), record(previous_trace_path),
               trace["table"], *trace["checked"], record(source / "benchmark_tools" / PROTOCOL)]
    captures = [c for c in full["captures"] if c["source_path"].endswith("execution_plan.json")]
    if len(captures) != 1:
        raise ValueError("Ambiguous admitted inferred-plan snapshot")
    check(captures[0]["snapshot"])
    checked.append(captures[0]["snapshot"])
    parent = json.loads(Path(captures[0]["snapshot"]["path"]).read_text())
    if Path(parent["repo"]).resolve() != repo or tuple(c["label"] for c in parent["cases"]) != LABELS:
        raise ValueError("Wrong inferred-parent repository or ordered cases")
    sources = []
    for old in previous["sources"]:
        relative = Path(old["path"]).relative_to(repo)
        current = record(source / relative)
        if current["sha256"] != old["sha256"]:
            raise ValueError("Existing reference-trace helper changed: " + str(relative))
        sources.append(current)
    sources.extend([record(Path(__file__)), record(source / "benchmark_tools/trace_deferred_reference_pairs.slurm")])
    cases = []
    for cell, case in zip(cells, parent["cases"]):
        if (argument(case["command"], "--species-tree-mode") != "infer"
                or argument(case["command"], "--root-rule") != "species_overlap"
                or argument(case["command"], "--pair-rule") != "positive_paralogy"):
            raise ValueError("Require matched own-tree duplication-aware inference")
        folder = Path(case["output"])
        if any(Path(cell["products"][name]["path"]).parent != folder / "output/orthohmm_phylogeny" for name in PRODUCTS):
            raise ValueError("Admitted native output differs from the planned arm")
        for item in [*cell["products"].values(), *cell["checked"]]:
            check(item)
        checked.extend(cell["products"].values())
        inputs = deepcopy(previous["inputs"])
        raw = list((folder / "qfo_scoring/results/VGNC").glob("*raw.txt.gz"))
        databases = list((folder / "qfo_scoring/other").glob("*.db"))
        if len(raw) != 1 or len(databases) != 1:
            raise ValueError("Ambiguous participant database or raw VGNC table")
        inputs.update(comparison=receipt["full_report"], candidates=record(Path(argument(case["command"], "--candidate-clusters"))),
            roots=cell["products"][PRODUCTS[0]], pairs=cell["products"][PRODUCTS[1]], raw=record(raw[0]),
            database=record(databases[0]), protocol=record(source / "benchmark_tools" / PROTOCOL))
        if case["label"] == "control" and any(inputs[key]["sha256"] != previous["inputs"][key]["sha256"]
                                               for key in ("candidates", "roots", "pairs")):
            raise ValueError("Matched control no longer reproduces original native inputs/products")
        checked.extend(inputs.values())
        cases.append({"label": case["label"], "trace_plan": {
            "status": "prepared_unrun", "inputs": inputs, "sources": sources,
            "expected_universe": previous["expected_universe"], "expected_pairs": cell["coverage"]["pairs"],
            "expected_metrics": cell["details"]["VGNC F"], "expected_f1": cell["scores"]["QfO: VGNC F"]}})
    checked = list({item["path"]: item for item in checked}.values())
    for item in [*checked, *sources]:
        check(item)
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    for case in cases:
        case["trace_plan"]["source_revision"] = revision
    plan = {"status": "prepared_unrun", "repo": str(repo), "source_root": str(source), "source_revision": revision,
        "output": str(output), "inputs": checked, "sources": sources, "cases": cases,
        "original_control_table": trace["table"], "cpus": 1, "memory_gib": 8,
        "scope": "Development-exposed diagnostic reusing completed native inference; not new search, accuracy gain, or inference timing."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", plan)
    return plan


def run(plan_path):
    plan = json.loads(plan_path.read_text())
    output, source = Path(plan["output"]), Path(plan["source_root"])
    if (Path(__file__).resolve() != source / "benchmark_tools/trace_deferred_reference_pairs.py"
            or (output / "results.json").exists() or tuple(c["label"] for c in plan["cases"]) != LABELS):
        raise ValueError("Require frozen source, untouched output, and two ordered cases")
    checked = [record(plan_path), *plan["inputs"], *plan["sources"]]
    report = {"status": "running", "job_id": os.environ.get("SLURM_JOB_ID"), "execution_plan": checked[0],
        "source_revision": plan["source_revision"], "scope": plan["scope"], "arms": [],
        "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", report)
    tables = []
    try:
        for case in plan["cases"]:
            for item in checked:
                check(item)
            child = output / (case["label"] + "_plan.json")
            write_json(child, case["trace_plan"])
            row = {"label": case["label"], "status": "trace_running"}
            report["arms"].append(row)
            write_json(output / "results.json", report)
            result = tracing.run(child, output / case["label"])
            if (result["reference_pairs"] != 23934 or result["within_species_reference_pairs_preserved"] != 2
                    or case["label"] == "control" and result["table"]["sha256"] != plan["original_control_table"]["sha256"]):
                raise ValueError("Reference universe or exact control table reproduction failed")
            checked.extend([record(child), result["table"], record(output / case["label"] / "report.json")])
            with Path(result["table"]["path"]).open() as stream:
                tables.append(list(csv.DictReader(stream, delimiter="\t")))
            if len(tables[-1]) != result["reference_pairs"]:
                raise ValueError("Physical reference table row count differs from the verified cohort")
            row.update(status="complete_reference_trace_verified", table=result["table"], summary=result["summary"],
                raw_categories=result["raw_categories"], resources=result["resources"])
            write_json(output / "results.json", report)
        report["paired_transitions"] = paired_transitions(*tables)
        for item in checked:
            check(item)
        report.update(status="complete_paired_reference_localization_verified", checked=checked)
    except Exception as error:
        report.update(status="failed_with_explicit_gap", error_type=type(error).__name__, error=str(error))
        if report["arms"] and report["arms"][-1]["status"] != "complete_reference_trace_verified":
            report["arms"][-1].update(status="failed", error_type=type(error).__name__, error=str(error))
    write_json(output / "results.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    prep = modes.add_parser("prepare")
    for name in ("repo", "source-root", "receipt", "previous-plan", "previous-trace", "output"):
        prep.add_argument("--" + name, type=Path, required=True)
    worker = modes.add_parser("run")
    worker.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source_root.resolve(), args.receipt.resolve(), args.previous_plan.resolve(),
                      args.previous_trace.resolve(), args.output.resolve()) if args.mode == "prepare" else run(args.execution_plan.resolve()))
    print(json.dumps({"status": result["status"]}))
    raise SystemExit(0 if result["status"] in {"prepared_unrun", "complete_paired_reference_localization_verified"} else 1)
