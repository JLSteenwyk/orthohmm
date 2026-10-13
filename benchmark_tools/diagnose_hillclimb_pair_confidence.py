"""Describe GO/EC pair-set changes by native event confidence, without tuning."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.compare_qfo_scored_pairs import read_scores
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json
from qfo_benchmark.og_to_pairwise import _strip_to_uniprot

HEADER = ["gene_a", "species_a", "gene_b", "species_b", "confidence"]


def read_confidence(path, watched, expected_count):
    found, counts, normalized = {}, Counter(), {}
    def accession(gene):
        if gene not in normalized:
            normalized[gene] = _strip_to_uniprot(gene)
        return normalized[gene]
    with path.open() as stream:
        rows = csv.reader(stream, delimiter="\t")
        if next(rows, None) != HEADER:
            raise ValueError("Unexpected native pair-confidence header")
        for row in rows:
            if (len(row) != 5 or not all(row) or row[0] >= row[2] or row[1] == row[3]
                    or row[4] not in {"high", "medium"}):
                raise ValueError("Malformed native confidence row")
            counts[row[4]] += 1
            pair = tuple(sorted((accession(row[0]), accession(row[2]))))
            if pair[0] == pair[1]:
                raise ValueError("Collapsed native accession pair")
            if pair in watched:
                if pair in found:
                    raise ValueError("Duplicate normalized watched pair")
                found[pair] = row[4]
    if sum(counts.values()) != expected_count or set(found) != watched:
        raise ValueError("Native count or scored-pair confidence coverage differs")
    if len(set(normalized.values())) != len(normalized):
        raise ValueError("Non-injective native accession normalization")
    return found, dict(counts)


def summarize(left, right, left_confidence, right_confidence):
    buckets = defaultdict(lambda: [0, 0])
    for pair in left.keys() | right.keys():
        if pair in left and pair in right:
            if left[pair] != right[pair]:
                raise ValueError("Shared serialized annotation score changed")
            category, confidence, score = "shared", left_confidence[pair] + "->" + right_confidence[pair], left[pair]
        elif pair in left:
            category, confidence, score = "removed", left_confidence[pair], left[pair]
        else:
            category, confidence, score = "added", right_confidence[pair], right[pair]
        buckets[category, confidence][0] += 1
        buckets[category, confidence][1] += score
    return [{"category": category, "confidence": confidence, "scored_relations": count,
             "score_sum_millionths": total, "serialized_mean": total / count / 1e6}
            for (category, confidence), (count, total) in sorted(buckets.items())]


def run(receipt, output):
    if output.exists():
        raise ValueError("Preserve existing diagnostic")
    summary = json.loads(receipt.read_text())
    check(summary["full_report"])
    full = json.loads(Path(summary["full_report"]["path"]).read_text())
    cells = full["cells"]
    if (tuple(c["method_key"] for c in cells) != ("deferred_inferred:control", "deferred_inferred:deferred500")
            or any(c["status"] != "independent_content_readback_verified" or c["dataset"] != "QfO" for c in cells)):
        raise ValueError("Require two independently admitted own-tree QfO cells")
    evidence = [record(receipt), summary["full_report"]]
    captures = [c for c in full["captures"] if c["source_path"].endswith("execution_plan.json")]
    if len(captures) != 1:
        raise ValueError("Ambiguous producer plan")
    check(captures[0]["snapshot"])
    evidence.append(captures[0]["snapshot"])
    plan = json.loads(Path(captures[0]["snapshot"]["path"]).read_text())
    if tuple(c["label"] for c in plan["cases"]) != ("control", "deferred500"):
        raise ValueError("Changed planned cell order")
    scores = {metric: [] for metric in ("GO", "EC")}
    for cell, case in zip(cells, plan["cases"]):
        for metric in scores:
            paths = list((Path(case["output"]) / "qfo_scoring/results" / metric).glob("*raw.txt.gz"))
            if len(paths) != 1:
                raise ValueError("Ambiguous raw functional endpoint")
            raw = read_scores(paths[0], metric)
            expected = cell["scores"]["QfO: " + metric]
            if (len(raw) != cell["details"][metric]["assessed_relations"]
                    or abs(sum(raw.values()) / len(raw) / 1e6 - expected) > 5.051e-7):
                raise ValueError("Raw count or mean differs from independently admitted endpoint")
            scores[metric].append(raw)
            evidence.append(record(paths[0]))
    confidence, native_counts = [], []
    for index, cell in enumerate(cells):
        identity = cell["products"][PRODUCTS[2]]
        check(identity)
        evidence.append(identity)
        watched = set().union(*(set(values[index]) for values in scores.values()))
        found, counts = read_confidence(Path(identity["path"]), watched, cell["coverage"]["pairs"])
        confidence.append(found)
        native_counts.append(counts)
    comparisons = {metric: summarize(*values, *confidence) for metric, values in scores.items()}
    for item in evidence:
        check(item)
    result = {"status": "descriptive_functional_pair_confidence_localization", "evidence": evidence,
        "sources": [record(__file__), record(Path(__file__).with_name("compare_qfo_scored_pairs.py")),
                    record(Path(__file__).parent.parent / "qfo_benchmark/og_to_pairwise.py")],
        "native_confidence_counts": native_counts, "comparisons": comparisons,
        "limitations": ["GO/EC are eligible annotated subsets, not orthology truth or all predicted pairs.",
            "Native high/medium labels describe reconciliation evidence, not calibrated posterior probabilities or bootstrap support.",
            "Shared/added/removed conditioning is descriptive; original benchmark denominators remain unchanged.",
            "Reference annotations do not control inference; no new search, inference, scoring, tuning, or uncertainty estimation."],
        "defaults_promoted": False, "goal_completed": False}
    write_json(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.receipt, args.output)["comparisons"]))
