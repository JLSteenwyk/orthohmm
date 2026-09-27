"""Descriptive GO/EC scored-pair overlap, not paired uncertainty estimation."""

import argparse
import csv
import gzip
import json
from pathlib import Path
import re

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def read_scores(path, metric):
    if metric not in {"GO", "EC"}:
        raise ValueError("Require GO or EC")
    scores = {}
    with gzip.open(path, "rt") as stream:
        headers = [next(stream, "").rstrip("\n") for _ in range(3)]
        if (not headers[0].startswith(f"# {metric} Similarities between orthologs from ")
                or not headers[1].startswith("# Computing timestamp: ")
                or headers[2] != f"# Protein ID 1<tab>Protein ID 2<tab>{metric} Similarity"):
            raise ValueError("Unexpected raw header")
        for fields in csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE):
            if len(fields) != 3:
                raise ValueError("Malformed scored row")
            a, b, value = fields
            if not a or not b or a == b or not re.fullmatch(r"(?:0\.\d{6}|1\.000000)", value):
                raise ValueError("Invalid pair or six-decimal score")
            pair = tuple(sorted((a, b)))
            if pair in scores:
                raise ValueError("Duplicate undirected pair")
            # Integer millionths preserve the serialized sums exactly.
            scores[pair] = int(value[0]) * 1000000 + int(value[2:])
    if not scores:
        raise ValueError("Empty scored set")
    return scores


def compare(left, right):
    shared = left.keys() & right.keys()
    left_sum, right_sum = sum(left.values()), sum(right.values())
    shared_left = sum(left[p] for p in shared)
    shared_right = sum(right[p] for p in shared)
    nl, nr, ns = len(left), len(right), len(shared)
    if not nl or not nr:
        raise ValueError("Both scored sets must be nonempty")
    components = dict(
        shared_with_original_denominators=shared_left / (nl * 1e6) - shared_right / (nr * 1e6),
        left_only=(left_sum - shared_left) / (nl * 1e6),
        negative_right_only=-(right_sum - shared_right) / (nr * 1e6))
    return dict(left_pairs=nl, right_pairs=nr, shared_pairs=ns,
        left_only_pairs=nl-ns, right_only_pairs=nr-ns,
        left_score_sum_millionths=left_sum, right_score_sum_millionths=right_sum,
        shared_left_sum_millionths=shared_left, shared_right_sum_millionths=shared_right,
        shared_pairs_with_different_serialized_scores=sum(left[p] != right[p] for p in shared),
        maximum_shared_absolute_difference_millionths=max((abs(left[p]-right[p]) for p in shared), default=None),
        shared_conditional_mean_difference=(shared_left-shared_right)/(ns*1e6) if ns else None,
        left_mean=left_sum/(nl*1e6), right_mean=right_sum/(nr*1e6),
        original_mean_difference=left_sum/(nl*1e6)-right_sum/(nr*1e6),
        original_mean_difference_components=components)


def audit(left, right, metric, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    inputs = [record(left), record(right)]
    result = compare(read_scores(left, metric), read_scores(right, metric))
    for item in inputs:
        check(item)
    report = dict(status="descriptive_scored_pair_comparison", metric=metric,
        inputs=inputs, source=record(__file__), result=result, uncertainty_admitted=False,
        limitations=["Rounded eligible-pair scores only, not all predicted relations.",
            "Shared-pair conditioning changes the endpoint and is not a replacement benchmark score.",
            "Decomposition is arithmetic, not causal attribution or independent-unit inference.",
            "No underlying annotation score recomputation or biological confidence interval."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("left", "right", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--metric", required=True, choices=("GO", "EC"))
    args = parser.parse_args()
    audit(args.left, args.right, args.metric, args.output)
