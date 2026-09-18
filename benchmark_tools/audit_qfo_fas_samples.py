"""Recompute retained FAS sample statistics without assuming independent pairs."""

import argparse
from collections import Counter
import csv
import gzip
import json
import math
from pathlib import Path
import statistics
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_swiss_comparators import COMPARISON_SHA, METHODS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def read_sample(path):
    pairs = {}
    with gzip.open(path, "rt") as stream:
        reader = csv.reader(stream, delimiter="\t")
        if next(reader, None) != ["Acc1", "Acc2", "FAS"]:
            raise ValueError("Unexpected FAS header")
        for row in reader:
            if len(row) != 3:
                raise ValueError("Malformed FAS row")
            a, b, value = row
            score = float(value)
            if not a or a >= b or "_" in a or "_" in b:
                raise ValueError("Noncanonical accession pair")
            if not math.isfinite(score) or not 0 <= score <= 1:
                raise ValueError("Invalid FAS score")
            if (a, b) in pairs:
                raise ValueError("Duplicate FAS pair")
            pairs[a, b] = score
    if len(pairs) < 2:
        raise ValueError("Insufficient sample")
    return pairs


def summarize(pairs, participant):
    values = list(pairs.values())
    mean = statistics.mean(values)
    sem = statistics.stdev(values) / math.sqrt(len(values))
    for actual, key in ((mean, "metric_y"), (sem, "stderr_y")):
        if not math.isclose(actual, participant[key], rel_tol=0, abs_tol=1e-12):
            raise ValueError("Saved raw FAS sample does not reproduce aggregate")
    eligible = participant["metric_x"]
    if not math.isfinite(eligible) or eligible != int(eligible) or eligible < len(values):
        raise ValueError("Invalid eligible pair count")
    degrees = Counter(accession for pair in pairs for accession in pair)
    return {"mean": mean, "native_pair_iid_sem": sem, "sample_pairs": len(values),
            "reported_eligible_pairs": int(eligible), "sample_fraction": len(values) / eligible,
            "sample_proteins": len(degrees), "proteins_in_multiple_sample_pairs": sum(v > 1 for v in degrees.values()),
            "maximum_sample_protein_degree": max(degrees.values())}


def audit(repo):
    comparison_path = repo / "benchmark_tools/results/publication_comparison_orthomcl_complete_20260916.json"
    identity = record(comparison_path)
    if identity["sha256"] != COMPARISON_SHA:
        raise ValueError("Changed frozen comparison")
    comparison = json.loads(comparison_path.read_text())
    if tuple(row["key"] for row in comparison["methods"]) != METHODS:
        raise ValueError("Changed method inventory")
    sources = [(row["key"], row["qfo"]["metric_details"]["FAS"]["source"]) for row in comparison["methods"]]
    sources.extend((f"checked_v2_{i}", record(repo / f"qfo_benchmark/scoring/checked_v2_{i}/results/FAS/FAS.json")) for i in range(4))
    rows, samples = [], {}
    for name, source in sources:
        check(source)
        aggregate = json.loads(Path(source["path"]).read_text())["datalink"]["inline_data"]
        if aggregate["visualization"]["x_axis"] != "NR_ORTHOLOGS" or aggregate["visualization"]["y_axis"] != "FAS":
            raise ValueError("Unexpected FAS axes")
        participants = aggregate["challenge_participants"]
        raw_paths = list(Path(source["path"]).parent.glob("*raw.txt.gz"))
        if len(participants) != 1 or len(raw_paths) != 1:
            raise ValueError("Ambiguous participant or raw file")
        before = record(raw_paths[0])
        samples[name] = read_sample(raw_paths[0])
        rows.append({"method": name, "aggregate": source, "raw": before,
                     **summarize(samples[name], participants[0])})
        check(before)
    overlaps = [{"a": a, "b": b, "shared_sample_pairs": len(samples[a].keys() & samples[b].keys())}
                for i, a in enumerate(samples) for b in list(samples)[i + 1:]]
    return {"status": "saved_sample_means_and_native_sem_verified", "source": record(__file__),
            "comparison": identity, "native_scorer": record(repo / "qfo_benchmark/benchmark-webservice/fas_benchmark.py"),
            "methods": rows, "sample_overlaps": overlaps,
            "limitations": ["NR_ORTHOLOGS counts precomputed plus annotation-eligible missing pairs, not the scored sample or all predictions.",
                            "Native scorer shuffles without setting a seed and caps newly computed pairs at 9000; historical RNG state is not recovered.",
                            "Native SEM assumes independent pair values; independence cannot be assumed for repeated proteins and homologous families.",
                            "Pair overlap is descriptive, not an independent unit definition or a paired comparison confidence interval.",
                            "This verifies saved-sample arithmetic, not prediction membership, annotation completeness, FAS score correctness, or sampling representativeness."]}


def render_table(result):
    lines = ["# Retained FAS Samples", "", "Generated by `audit_qfo_fas_samples.py` from retained raw scores.", "",
             "| Method | Scored pairs | Reported eligible pairs | Sample % | Mean FAS |",
             "|---|---:|---:|---:|---:|"]
    for row in result["methods"]:
        lines.append(f"| {row['method']} | {row['sample_pairs']:,} | {row['reported_eligible_pairs']:,} | "
                     f"{100 * row['sample_fraction']:.4f} | {row['mean']:.9f} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--table", type=Path)
    args = parser.parse_args()
    result = audit(args.repo)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    if args.table:
        with args.table.open("x") as stream:
            stream.write(render_table(result))
