"""Measure a protected union of retained unfiltered HMM candidate arrays."""

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_search_recovery import audit_direction, read_raw
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_search_recovery import host_activity, link_changes, runtime_records
from benchmark_tools.search_decision_trace import classify_search_result
from benchmark_tools.trace_ob_search_decisions import load_watched

PROTOCOL = "hillclimb_protected_candidate_union_diagnostic_20261016.json"
SOURCES = ("diagnose_protected_candidate_union.py", "diagnose_protected_candidate_union.slurm", PROTOCOL,
    "audit_search_recovery.py", "prepare_ob_candidate_neighborhood.py", "run_membership_hillclimb.py",
    "run_search_recovery.py", "search_decision_trace.py", "trace_ob_search_decisions.py")


def protected_union(left, right, queries, targets):
    for result in (left, right):
        classify_search_result(result, queries, targets, [], 0.0001)
        if any(getattr(result, field).dtype != np.dtype("float64") for field in ("scores", "evalues")):
            raise ValueError("Require retained float64 score/E-value arrays")
        if np.any(np.bincount(result.query_indices, minlength=len(queries)) > 100):
            raise ValueError("Original arm exceeds its 100-candidate per-query cap")
    keys = [r.query_indices.astype(np.int64) * len(targets) + r.target_indices for r in (left, right)]
    common, a, b = np.intersect1d(*keys, assume_unique=True, return_indices=True)
    for field in ("scores", "evalues"):
        if not np.array_equal(getattr(left, field)[a].view(np.uint64),
                              getattr(right, field)[b].view(np.uint64)):
            raise ValueError("Common candidate numerical values are not bit-for-bit identical: " + field)
    extra = ~np.isin(keys[1], common, assume_unique=True)
    combined_keys = np.concatenate((keys[0], keys[1][extra]))
    order = np.argsort(combined_keys, kind="stable")
    result = SimpleNamespace(candidate_count=len(order), **{
        field: np.concatenate((getattr(left, field), getattr(right, field)[extra]))[order]
        for field in ("query_indices", "target_indices", "scores", "evalues")})
    classify_search_result(result, queries, targets, [], 0.0001)
    counts = np.bincount(result.query_indices, minlength=len(queries))
    if np.any(counts > 200) or len(np.unique(combined_keys)) != result.candidate_count:
        raise ValueError("Union exceeds its bound or contains duplicate candidates")
    return result, {"baseline": left.candidate_count, "relaxed": right.candidate_count,
        "common": len(common), "baseline_only": left.candidate_count - len(common),
        "relaxed_only": right.candidate_count - len(common), "union": result.candidate_count,
        "per_query_union_histogram": {str(k): v for k, v in sorted(Counter(map(int, counts)).items())},
        "common_values_bitwise_equal": True}


def advance_gate(candidates, thresholds):
    return (candidates["union"] <= 1.25 * candidates["baseline"]
            and all(t["versus_control"]["lost"] == 0 for t in thresholds)
            and any(t["versus_control"]["gained"] > 0 for t in thresholds))


def context(repo, source):
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    source_pins = [record(source / "benchmark_tools" / name) for name in SOURCES]
    for pin in source_pins:
        relative = Path(pin["path"]).relative_to(source)
        blob = subprocess.check_output(["git", "-C", str(source), "show", f"{revision}:{relative}"])
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Diagnostic differs from committed frozen source")
    protocol = json.loads((source / "benchmark_tools" / PROTOCOL).read_text())
    pins = list(source_pins)
    evidence = []
    for item in protocol["prior_evidence"]:
        pin = {**item, "path": str(repo / item["path"])}
        check(pin)
        pins.append(pin)
        evidence.append(json.loads(Path(pin["path"]).read_text()))
    audited = evidence[0]
    spec = protocol["inputs"]
    execution_pin = record(repo / spec["execution_plan"])
    if (execution_pin["sha256"] != spec["execution_plan_sha256"]
            or audited["execution_plan"] != execution_pin
            or audited["status"] != "all_search_arms_independently_verified"):
        raise ValueError("Original execution or independent audit differs")
    execution = json.loads(Path(execution_pin["path"]).read_text())
    if execution["source_revision"] != spec["scientific_revision"] or runtime_records() != execution["runtime"]:
        raise ValueError("Original diagnostic scientific source or runtime differs")
    pins += [execution_pin, *execution["inputs"], *execution["sources"], *execution["runtime"]["files"]]
    reports = []
    for label, expected_hash, minimum in zip(spec["arms"],
            (spec["control_report_sha256"], spec["relaxed_report_sha256"]), (4, 1)):
        path = Path(execution["output"]) / label / "report.json"
        pin = record(path)
        rows = [a for a in audited["arms"] if a["label"] == label]
        report = json.loads(path.read_text())
        if (len(rows) != 1 or rows[0]["status"] != "independently_verified"
                or rows[0]["retained_report"] != pin or pin["sha256"] != expected_hash
                or report["execution_plan"] != execution_pin or report["source_revision"] != spec["scientific_revision"]
                or report["status"] != "fresh_search_completed_pending_independent_audit"
                or report["label"] != label or report["explicit_control_fixture_exact"] is not True
                or report["settings"] != {**protocol["fixed_settings"], "min_total_hits": minimum}):
            raise ValueError("Original arm differs from authenticated audited search")
        pins.append(pin)
        reports.append(report)
    for pin in pins:
        check(pin)
    return protocol, execution, reports, pins, revision


def analyze(protocol, execution, reports):
    watched = load_watched(execution["pair_trace"]["path"])
    species, owners = {}, {}
    for pin in execution["fastas"]:
        name = Path(pin["path"]).name
        if name in species:
            raise ValueError("Duplicate FASTA species")
        species[name] = [s.id for s in SeqIO.parse(pin["path"], "fasta")]
        for gene in species[name]:
            if gene in owners:
                raise ValueError("Duplicate FASTA gene")
            owners[gene] = name
    if len(species) != 12 or len(owners) != 251378 or len(watched) != protocol["known_results_not_new_predictions"]["watched_pairs"]:
        raise ValueError("Wrong complete diagnostic FASTA/watched universe")
    directions = {}
    for q, t in watched:
        directions.setdefault((owners[q], owners[t]), []).append((q, t))
    expected = sorted(directions)
    if len(expected) != protocol["known_results_not_new_predictions"]["directions"] or any(
            [(r["query_species"], r["target_species"]) for r in report["directions"]] != expected for report in reports):
        raise ValueError("Missing, duplicated or reordered diagnostic directions")
    pins, output, totals = [], [], {t: {"decisions": Counter(), "versus_control": Counter(),
        "historically_absent_accepted": 0, "all_significant_hits": 0} for t in protocol["thresholds"]}
    candidates, histogram = Counter(), Counter()
    for i, key in enumerate(expected):
        pairs = sorted(directions[key])
        wanted = {q for q, _ in pairs}
        queries, targets = [q for q in species[key[0]] if q in wanted], species[key[1]]
        items = [r["directions"][i] for r in reports]
        audits = [audit_direction(item, queries, targets, pairs, watched, protocol["thresholds"]) for item in items]
        raw = [read_raw(item["raw"]["path"])[0] for item in items]
        union, stats = protected_union(*raw, queries, targets)
        histogram.update(stats["per_query_union_histogram"])
        candidates.update({k: stats[k] for k in ("baseline", "relaxed", "common", "baseline_only", "relaxed_only", "union")})
        threshold_rows = []
        for t in protocol["thresholds"]:
            rows = classify_search_result(union, queries, targets, pairs, t)
            row = {"threshold": t, "decisions": dict(Counter(r["decision"] for r in rows)),
                "versus_control": link_changes(audits[0][t]["rows"], rows),
                "historically_absent_accepted": sum(r["decision"] == "accepted" and not watched[r["query"], r["target"]]["historical_present"] for r in rows),
                "all_significant_hits": int(np.count_nonzero(union.evalues < t))}
            for name in ("decisions", "versus_control"):
                totals[t][name].update(row[name])
            for name in ("historically_absent_accepted", "all_significant_hits"):
                totals[t][name] += row[name]
            threshold_rows.append(row)
        pins += [pin for item in items for pin in (item["raw"], *[r["table"] for r in item["thresholds"]])]
        output.append({"query_species": key[0], "target_species": key[1], "candidates": stats, "thresholds": threshold_rows})
    thresholds = [{"threshold": t, **{k: dict(v) if isinstance(v, Counter) else v for k, v in values.items()}}
                  for t, values in totals.items()]
    return {"directions": output, "candidates": dict(candidates), "thresholds": thresholds,
        "per_query_union_histogram": dict(sorted(histogram.items(), key=lambda x: int(x[0]))),
        "advance_gate_passed": advance_gate(candidates, thresholds)}, pins


def run(repo, output):
    source = Path(__file__).resolve().parent.parent
    if source == repo or output.exists():
        raise ValueError("Require a separate frozen executor and fresh report namespace")
    if (os.environ.get("SLURM_CPUS_PER_TASK"), os.environ.get("SLURM_MEM_PER_NODE"),
            os.environ.get("SLURMD_NODENAME"), os.environ.get("CUDA_VISIBLE_DEVICES")) != ("1", "16384", "bizon", ""):
        raise ValueError("Require CPU-only 1-CPU 16-GiB local allocation")
    output.mkdir(parents=True)
    started, usage = time.monotonic(), resource.getrusage(resource.RUSAGE_SELF)
    report = {"status": "validating_cached_protected_union_inputs", "job_id": os.environ.get("SLURM_JOB_ID"),
        "native_inference_rerun": False, "accuracy_evaluated": False, "defaults_promoted": False, "goal_completed": False,
        "host_before": host_activity()}
    write_json(output / "results.json", report)
    try:
        protocol, execution, reports, pins, revision = context(repo, source)
        report.update(source_revision=revision, checked=pins)
        details, raw_pins = analyze(protocol, execution, reports)
        pins += raw_pins
        for pin in pins:
            check(pin)
        if runtime_records() != execution["runtime"]:
            raise ValueError("Runtime changed during cached diagnostic")
        report.update(status="cached_protected_union_verified_not_fresh_search_or_accuracy", **details,
            postflight_pins_verified=True, checked=list({p["path"]: p for p in pins}.values()),
            scope=protocol["scope"], next_action="prespecify_fresh_policy_pilot" if details["advance_gate_passed"] else "stop_this_branch_without_changing_gate")
    except Exception as error:
        report.update(status="failed_cached_protected_union_diagnostic", error_type=type(error).__name__, error=str(error))
    after = resource.getrusage(resource.RUSAGE_SELF)
    report.update(host_after=host_activity(), resources={"diagnostic_only_wall_s": time.monotonic() - started,
        "self_cpu_s": after.ru_utime + after.ru_stime - usage.ru_utime - usage.ru_stime, "self_max_rss_kib": after.ru_maxrss})
    write_json(output / "results.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = run(args.repo.resolve(), args.output.resolve())
    print(result["status"])
    raise SystemExit(0 if result.get("postflight_pins_verified") else 1)
