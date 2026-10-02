"""Render all prespecified matched-recall graph contrasts from scored records."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.reproduce_matched_graph_statistics import reproduce
from benchmark_tools.score_matched_graph import CONDITIONS, summarize


LABELS = ["Baseline", "Divergent", "Divergent + turnover", "Missing 20%", "Taxon-count control", "Turnover", "Uneven taxa", "Overall"]


def _replay_differences(reported, recomputed, path=()):
    """Explain exact JSON-value mismatches without changing replay admission."""
    if isinstance(reported, dict) and isinstance(recomputed, dict):
        differences = []
        for key in sorted(reported.keys() | recomputed.keys()):
            if key not in reported or key not in recomputed:
                differences.append(dict(path=[*path, key],
                    reported_present=key in reported, recomputed_present=key in recomputed,
                    reported=reported.get(key), recomputed=recomputed.get(key)))
            else:
                differences.extend(_replay_differences(reported[key], recomputed[key], (*path, key)))
        return differences
    if isinstance(reported, list) and isinstance(recomputed, list) and len(reported) == len(recomputed):
        return [difference for index, (a, b) in enumerate(zip(reported, recomputed))
                for difference in _replay_differences(a, b, (*path, index))]
    if reported != recomputed:
        return [dict(path=list(path), reported=reported, recomputed=recomputed)]
    return []


def validate_replay(result, replay_policy="exact"):
    if replay_policy not in ("exact", "count-level"):
        raise ValueError("Unknown replay policy")
    recomputed = summarize(result["records"])
    differences = _replay_differences(
        {key: result[key] for key in ("contrasts", "bootstrap")},
        {key: recomputed[key] for key in ("contrasts", "bootstrap")})
    count_validation = reproduce(result) if replay_policy == "count-level" else None
    tolerance = count_validation["absolute_tolerance"] if count_validation else 0
    rejected = differences
    if count_validation:
        rejected = [row for row in differences if not (
            row["path"][0] == "contrasts"
            and type(row["reported"]) is float and type(row["recomputed"]) is float
            and np.isfinite(row["reported"]) and np.isfinite(row["recomputed"])
            and abs(row["reported"] - row["recomputed"]) <= tolerance)]
    if rejected:
        raise ValueError("Reported effects differ from paired score records: " + json.dumps(differences, sort_keys=True))
    return dict(policy=replay_policy, absolute_tolerance=tolerance, relative_tolerance=0,
                exact_differences=differences, count_validation=count_validation,
                scorer=record(Path(summarize.__code__.co_filename)),
                compared_sections=["contrasts", "bootstrap"])


def plot(result, *, replay_policy="exact"):
    validate_replay(result, replay_policy)
    return _draw(result)


def _draw(result):
    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 6.8), gridspec_kw={"width_ratios": [1.05, 1]})
    fig.subplots_adjust(left=.20, right=.97, bottom=.23, top=.79, wspace=.20)
    labels = [*CONDITIONS, "overall"]
    y = np.arange(len(labels))
    for arm, color, offset, label in (("hmm", "#187f70", -.12, "HMM"), ("diamond", "#666666", .12, "DIAMOND")):
        left.scatter([100 * result["contrasts"][c]["f1"][arm + "_mean"] for c in labels], y + offset,
                     s=40, color=color, label=label, zorder=3)
    left.set_yticks(y, LABELS)
    left.invert_yaxis()
    left.set_xlim(0, 100)
    left.set_xlabel("Mean pair F1 (%)")
    left.legend(loc="lower left", frameon=False, ncol=2, fontsize=10)
    left.set_title("A  Final graph partitions", loc="left", fontsize=12)
    for i, condition in enumerate(labels):
        row = result["contrasts"][condition]["f1"]
        lo, hi = row["bonferroni_8_ci"]
        color = "#187f70" if condition != "overall" else "#1d1d1d"
        right.plot([lo, hi], [i, i], color=color, linewidth=2)
        right.scatter(row["difference_percentage_points"], i, color=color, s=40, zorder=3)
        right.scatter(row["seed_differences_percentage_points"], np.full(5, i + .17), color="#aeb8b5", s=12)
    right.axvline(0, color="#888888", linestyle="--", linewidth=.8)
    right.set_yticks(y, [""] * len(y))
    right.invert_yaxis()
    right.set_xlabel("HMM minus DIAMOND F1 (percentage points)")
    right.set_title("B  Paired differences", loc="left", fontsize=12)
    for ax in (left, right):
        ax.set_ylim(len(labels) - .5, -.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color="#e6e6e6", linewidth=.6)
        ax.set_axisbelow(True)
    fig.suptitle("Matched-recall initial-search control", fontsize=16, x=.20, ha="left", y=.97)
    fig.text(.20, .88, "35 reporting datasets; five seed blocks across seven conditions\nSame production graph/refinement; profile expansion and phylogeny off", fontsize=10)
    fig.text(.20, .09, "Intervals: paired seed-block percentile bootstrap, Bonferroni-adjusted for eight F1 contrasts.\n"
             "Small dots: five seed effects. Five blocks give limited tail resolution and approximate coverage.\n"
             "Development-exposed simulations; not real-data equivalence or full-pipeline superiority.", fontsize=9)
    return fig


def run(path, output, *, replay_policy="exact"):
    if output.exists():
        raise FileExistsError(output)
    identity = record(path)
    result = json.loads(path.read_text())
    validation = validate_replay(result, replay_policy)
    figure = _draw(result)
    files = []
    try:
        output.mkdir(parents=True, exist_ok=False)
        for extension in ("png", "pdf", "svg"):
            destination = output / ("matched_graph." + extension)
            figure.savefig(destination, dpi=180)
            files.append(record(destination))
    finally:
        plt.close(figure)
    if record(path) != identity:
        raise ValueError("Input changed during rendering")
    with (output / "manifest.json").open("x") as stream:
        json.dump(dict(source=record(__file__), input=identity, outputs=files,
                       recomputed_from_paired_score_records=True, replay_validation=validation,
                       visual_review_complete=False), stream, indent=2, sort_keys=True)
        stream.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replay-policy", choices=("exact", "count-level"), default="exact")
    args = parser.parse_args()
    run(args.results.resolve(), args.output.absolute(), replay_policy=args.replay_policy)
