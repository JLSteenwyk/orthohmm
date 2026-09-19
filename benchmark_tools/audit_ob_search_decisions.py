"""Independently recount completed search-decision artifacts, without rescoring."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.trace_ob_search_decisions import load_watched, SETTINGS, TRACE_SHA, RUNTIME_SHA
from benchmark_tools.build_publication_runtime import verify_runtime
from benchmark_tools.verify_ygob_validation import require_completed_job


def reconstruct(raw, query_ids, target_ids, watched):
    if raw["query_ids"].tolist() != query_ids or raw["target_ids"].tolist() != target_ids:
        raise ValueError("Raw identities differ from input order/universe")
    q, t, scores, evalues = [raw[k] for k in ("query_indices", "target_indices", "scores", "evalues")]
    if any(a.ndim != 1 for a in (q, t, scores, evalues)) or len({len(a) for a in (q, t, scores, evalues)}) != 1:
        raise ValueError("Inconsistent raw candidate dimensions")
    count = raw["candidate_count"]
    if count.ndim != 0 or count.dtype.kind not in "iu" or int(count) != len(q):
        raise ValueError("Incomplete raw candidate output")
    if q.dtype.kind not in "iu" or t.dtype.kind not in "iu":
        raise ValueError("Noninteger candidate indices")
    if (np.any(q < 0) or np.any(q >= len(query_ids)) or np.any(t < 0) or np.any(t >= len(target_ids))
            or not np.isfinite(scores).all() or not np.isfinite(evalues).all() or np.any(evalues < 0)):
        raise ValueError("Invalid raw candidate values")
    lookup = {(query_ids[int(qi)], target_ids[int(ti)]): (float(s), float(e))
              for qi, ti, s, e in zip(q, t, scores, evalues)}
    if len(lookup) != len(q):
        raise ValueError("Duplicate raw directed candidate")
    rows = []
    for (left, right), prior in sorted(watched.items()):
        if (left, right) in lookup:
            score, evalue = lookup[left, right]
            decision = "accepted" if evalue < 0.0001 else "scored_not_significant"
        else:
            score = evalue = None
            decision = "not_selected_by_prefilter"
        rows.append(dict(query=left, target=right, decision=decision, score=score,
                         evalue=evalue, historical_present=prior["historical_present"],
                         families=",".join(sorted(prior["families"]))))
    return rows


def audit(root, report_path, report_sha, job, output):
    if output.exists():
        raise FileExistsError(output)
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, job)
    if (scheduler["State"], scheduler["ExitCode"], scheduler["NodeList"], scheduler["AllocCPUS"]) != ("COMPLETED", "0:0", "bizon", "4"):
        raise ValueError("Wrong diagnostic completion or allocation")
    report = read_frozen(report_path, report_sha)
    if (report["status"] != "search_decisions_observed_pending_independent_audit"
            or report["job_id"] != str(job) or report["settings"] != SETTINGS
            or report["accuracy_evaluated"] is not False or report["publication_ready"] is not False
            or report["threads"] != 4):
        raise ValueError("Wrong diagnostic identity/settings")
    executor = root / "benchmarks/work/publication_ob_search_decisions_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != "2a6513d04585decad1afe53dac2e74bf0c7b3968":
        raise ValueError("Wrong diagnostic executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    if report["source"] != record(executor / "benchmark_tools/trace_ob_search_decisions.py"):
        raise ValueError("Wrong driver source")
    trace_path = root / "benchmark_tools/results/ob_family_trace_verified_20260916.json"
    trace = read_frozen(trace_path, TRACE_SHA)
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = read_frozen(runtime_path, RUNTIME_SHA)
    if report["runtime"] != record(runtime_path) or report["core_commit"] != runtime["commit"]:
        raise ValueError("Wrong runtime binding")
    verify_runtime(runtime_path, Path(runtime["root"]))
    checked = [record(report_path), record(trace_path), trace["pair_trace"], *report["checked_records"]]
    for item in checked:
        check(item)
    ids, owners = {}, {}
    for item in trace["inputs"]:
        if Path(item["path"]).suffix != ".fa":
            continue
        check(item)
        name = Path(item["path"]).name
        ids[name] = [r.id for r in SeqIO.parse(item["path"], "fasta")]
        for gene in ids[name]:
            if gene in owners:
                raise ValueError("Duplicate FASTA identity")
            owners[gene] = name
    watched = load_watched(trace["pair_trace"]["path"])
    expected = {}
    for pair, prior in watched.items():
        expected.setdefault((owners[pair[0]], owners[pair[1]]), {})[pair] = prior
    observed_keys = [(d["query_species"], d["target_species"]) for d in report["directions"]]
    if observed_keys != sorted(expected) or len(expected) != 144:
        raise ValueError("Incomplete or reordered direction inventory")
    counts, transitions, total = Counter(), Counter(), 0
    for index, direction in enumerate(report["directions"]):
        qsp, tsp = direction["query_species"], direction["target_species"]
        pairs = expected[qsp, tsp]
        queries = {q for q, _ in pairs}
        query_ids = [q for q in ids[qsp] if q in queries]
        for label, suffix in (("raw", "npz"), ("table", "tsv")):
            item = direction[label]
            if Path(item["path"]) != report_path.parent / f"direction_{index:03d}.{suffix}":
                raise ValueError("Unexpected artifact location")
            check(item)
            checked.append(item)
        with np.load(direction["raw"]["path"], allow_pickle=False) as raw:
            rows = reconstruct(raw, query_ids, ids[tsp], pairs)
            if direction["candidate_count"] != int(raw["candidate_count"]):
                raise ValueError("Candidate count differs")
        with Path(direction["table"]["path"]).open(newline="") as stream:
            retained = list(csv.DictReader(stream, delimiter="\t"))
        for row in retained:
            for name in ("score", "evalue"):
                row[name] = float(row[name]) if row[name] else None
            if row["historical_present"] not in ("True", "False"):
                raise ValueError("Invalid historical presence flag")
            row["historical_present"] = row["historical_present"] == "True"
        if retained != rows:
            raise ValueError("Watched-pair table differs from independent raw reconstruction")
        local = Counter(row["decision"] for row in rows)
        mismatch = sum((row["decision"] == "accepted") != row["historical_present"] for row in rows)
        if (direction["decisions"] != dict(local) or direction["historical_presence_mismatches"] != mismatch
                or direction["watched_pairs"] != len(rows) or direction["query_count"] != len(query_ids)
                or direction["target_count"] != len(ids[tsp])):
            raise ValueError("Direction summary differs from raw recount")
        counts.update(local)
        transitions.update(("historical_present" if r["historical_present"] else "historical_absent")
                           + ":" + r["decision"] for r in rows)
        total += len(rows)
    if total != len(watched):
        raise ValueError("Incomplete watched-pair coverage")
    for item in checked:
        check(item)
    result = dict(status="search_decisions_independently_recounted", source=record(__file__),
                  scheduler=scheduler, accounting=accounting, source_report=record(report_path),
                  checked_records=checked, directions=144, directed_pairs=total,
                  decisions=dict(counts), historical_transitions=dict(transitions),
                  accuracy_evaluated=False, publication_ready=False,
                  limitations=["Recounts persisted raw output; does not independently rescore sequences.",
                      "Historical disagreements are retained; no historical execution equivalence claim.",
                      "No counterfactual scoring, causal accuracy conclusion or comparative timing."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "report", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("report-sha256", "job"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    audit(args.root.resolve(), args.report.resolve(), args.report_sha256, args.job, args.output.resolve())
