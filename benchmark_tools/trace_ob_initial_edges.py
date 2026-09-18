"""Reconstruct frozen initial RBNH edges and trace admitted reference pairs."""

import argparse
from collections import Counter
import csv
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import pickle
import subprocess
import sys

import numpy as np

CORE = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
CORE_SHA = "1a35944ab7fea859143f599b1262aa111272787f6f735b9acc37d299427f2ad6"
TRACE_SHA = "bda00fb593b357bc8f07e43544feae598150ae891ed82c988fb339a92349de0c"
CACHE_SHA = "78a5af40ea2683a69e1baefefb0966c3549cffe71b24bda03b9ba4a6e29e65e1"


def record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def classify(forward, reverse, left_threshold, right_threshold, edge):
    if any(math.isnan(value) or value <= 0 for value in (left_threshold, right_threshold)):
        raise ValueError("Invalid endpoint threshold")
    scores = [score for score in (forward, reverse) if score is not None]
    if any(not math.isfinite(s) or s <= 0 for s in scores):
        raise ValueError("Invalid cached score")
    threshold = min(left_threshold, right_threshold)
    expected = any(score >= threshold for score in scores)
    if bool(edge) != expected:
        raise ValueError("Threshold decision differs from native edge")
    if not scores:
        return "no_direct_hit"
    if edge:
        return "initial_edge"
    return "no_finite_endpoint_threshold" if math.isinf(threshold) else "below_endpoint_threshold"


def capture_edges(function, names, species, q, t, scores):
    captured = []

    def tracer(frame, event, arg):
        if frame.f_code is function.__code__:
            if event == "return" and "thresholds" in frame.f_locals:
                captured.append(frame.f_locals["thresholds"].copy())
            return tracer
        return None

    previous = sys.gettrace()
    try:
        sys.settrace(tracer)
        edges = function(names, species, q, t, scores)
    finally:
        sys.settrace(previous)
    if len(captured) != 1:
        raise ValueError("Native threshold snapshot not reached exactly once")
    return edges, captured[0]


def run(repo, core, output):
    if output.exists():
        raise FileExistsError(output)
    if subprocess.check_output(["git", "-C", str(core), "rev-parse", "HEAD"], text=True).strip() != CORE:
        raise ValueError("Wrong frozen core revision")
    subprocess.run(["git", "-C", str(core), "diff", "--exit-code", "HEAD", "--", "orthohmm"], check=True)
    core_record = record(core / "orthohmm/accuracy.py")
    if core_record["sha256"] != CORE_SHA:
        raise ValueError("Wrong frozen graph source")
    if any(name == "orthohmm" or name.startswith("orthohmm.") for name in sys.modules):
        raise ValueError("Require fresh process without an imported OrthoHMM")
    sys.path.insert(0, str(core))
    accuracy = importlib.import_module("orthohmm.accuracy")
    if Path(accuracy.__file__).resolve() != Path(core_record["path"]):
        raise ValueError("Wrong imported core")
    trace_path = repo / "benchmark_tools/results/ob_family_trace_verified_20260916.json"
    trace_record = record(trace_path)
    if trace_record["sha256"] != TRACE_SHA:
        raise ValueError("Wrong admitted reference trace")
    trace = json.loads(trace_path.read_text())
    cache = repo / "benchmarks/results/hits_BLOSUM62_mc100.pkl"
    cache_record = record(cache)
    if cache_record["sha256"] != CACHE_SHA:
        raise ValueError("Untrusted hit pickle")
    pair_record = trace["pair_trace"]
    if record(pair_record["path"]) != pair_record:
        raise ValueError("Trace pair table changed")
    # Only this previously admitted, checksum-pinned local cache is unpickled.
    with cache.open("rb") as stream:
        payload = pickle.load(stream)
    names = sorted({str(gene) for gene in payload["all_gene_ids"]})
    ids = {gene: i for i, gene in enumerate(names)}
    species_ids = {name: i for i, name in enumerate(sorted({str(payload["gene_to_species"][g]) for g in names}))}
    species = np.array([species_ids[str(payload["gene_to_species"][g])] for g in names], dtype=np.int32)
    hits = payload["all_hits"]
    q = np.fromiter((ids[str(a)] for a, _ in hits), dtype=np.int32, count=len(hits))
    t = np.fromiter((ids[str(b)] for _, b in hits), dtype=np.int32, count=len(hits))
    scores = np.fromiter(hits.values(), dtype=np.float64, count=len(hits))
    if len(names) != 251378 or not np.isfinite(scores).all() or (scores <= 0).any():
        raise ValueError("Invalid retained hit universe")
    edges, thresholds = capture_edges(accuracy.build_rbnh_edges, names, species, q, t, scores)
    # Reference labels are used only after the complete native graph is built.
    with Path(pair_record["path"]).open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    wanted = {tuple(sorted((ids[row["left"]], ids[row["right"]]))) for row in rows}
    ref_ids = {idx for pair in wanted for idx in pair}
    mask = np.zeros(len(names), dtype=bool)
    mask[list(ref_ids)] = True
    selected = mask[edges.sources] & mask[edges.targets]
    present = {tuple(sorted((int(a), int(b)))) for a, b in zip(edges.sources[selected], edges.targets[selected])}
    output.mkdir(parents=True, exist_ok=False)
    summary, families = Counter(), {}
    table = output / "reference_pair_edges.tsv"
    fields = ["refog", "left", "right", "forward_score", "reverse_score", "left_threshold", "right_threshold", "decision", "multipass_refined", "root_hogs"]
    with table.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            a, b = row["left"], row["right"]
            forward, reverse = hits.get((a, b)), hits.get((b, a))
            for value, key in ((forward, "forward_normalized_hit"), (reverse, "reverse_normalized_hit")):
                if (None if row[key] == "NA" else float(row[key])) != value:
                    raise ValueError("Admitted trace differs from retained search hit")
            left, right = float(thresholds[ids[a]]), float(thresholds[ids[b]])
            decision = classify(forward, reverse, left, right, tuple(sorted((ids[a], ids[b]))) in present)
            key = decision + "/root_" + row["root_hogs"]
            summary[key] += 1
            families.setdefault(row["refog"], Counter())[key] += 1
            writer.writerow({"refog": row["refog"], "left": a, "right": b,
                             "forward_score": "NA" if forward is None else forward,
                             "reverse_score": "NA" if reverse is None else reverse,
                             "left_threshold": left, "right_threshold": right, "decision": decision,
                             "multipass_refined": row["multipass_refined"], "root_hogs": row["root_hogs"]})
    evidence = [core_record, trace_record, cache_record, pair_record]
    for item in evidence:
        if record(item["path"]) != item:
            raise ValueError("Source changed during graph trace")
    result = {"status": "frozen_initial_rbnh_reference_pairs_traced", "source": record(__file__),
              "job_id": os.environ.get("SLURM_JOB_ID"), "numpy": np.__version__, "python": sys.version,
              "checked_records": evidence, "native_edges": len(edges), "hits": len(hits),
              "pair_memberships": len(rows), "families": families, "summary": summary, "table": record(table),
              "limitations": ["Graph reconstruction, not a new search or clustering run.",
                  "Initial edges exclude singleton assignment and subsequent profile-added edges.",
                  "Exact ties follow retained cache insertion order and frozen native code.",
                  "No direct hit does not distinguish prefilter versus scoring rejection.",
                  "Raw family pair memberships include low-certainty and within-species pairs, not official recall.",
                  "Edge presence is not biological correctness or proof of causal final grouping effects."]}
    (output / "report.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "core", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.core.resolve(), args.output.resolve())
