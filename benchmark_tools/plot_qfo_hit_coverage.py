"""Plot source-bound corrected QfO search coverage, not orthology accuracy."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

SUMMARY_SHA = "ba2a963b49146d723038756e723bddd3a6230172e98f4355c824a8b4d1f1038f"
VARIANTS = ("hmm", "all_hits", "top100")
LABELS = ("HMM initial search", "DIAMOND all-hits", "DIAMOND top100")


def plot(report):
    if (report["status"] != "source_bound_label_free_coverage_summary"
            or report["accuracy_evaluated"] is not False
            or report["publication_ready"] is not False
            or tuple(row["variant"] for row in report["rows"]) != VARIANTS):
        raise ValueError("Wrong search diagnostic identity")
    hits, coverage = [], []
    for row in report["rows"]:
        if (row["genes"] != 984137 or type(row["directed_hits"]) is not int
                or row["directed_hits"] < 0
                or type(row["queries_without_cross_species_hits"]) is not int
                or not 0 <= row["queries_without_cross_species_hits"] <= row["genes"]):
            raise ValueError("Invalid plot counts")
        hits.append(row["directed_hits"] / 1e6)
        coverage.append(100 * (row["genes"] - row["queries_without_cross_species_hits"]) / row["genes"])
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    fig.subplots_adjust(left=.20, right=.97, bottom=.29, top=.72, wspace=.25)
    colors = ("#007d83", "#b64b42", "#747474")
    for i, (ax, values) in enumerate(zip(axes, (hits, coverage))):
        ax.barh(range(3), values, color=colors, height=.52, zorder=3)
        ax.set_yticks(range(3), LABELS if i == 0 else ["", "", ""])
        ax.set_ylim(2.6, -.6)
        ax.set_xlim(0, max(hits) * 1.22 if i == 0 else 100)
        ax.set_xlabel("Directed hits (millions)" if i == 0 else "Queries with cross-species hits (%)")
        ax.set_title("A  Search-hit counts" if i == 0 else "B  Query coverage", loc="left", fontsize=12)
        ax.grid(axis="x", alpha=.2, zorder=0)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for y, value in enumerate(values):
            ax.text(value + (8 if i == 0 else 1.5), y, f"{value:.2f}", va="center", fontsize=10)
    fig.suptitle("Corrected QfO: initial search diagnostics", x=.04, y=.97, ha="left", fontsize=17)
    fig.text(.04, .86, "984,137 proteins | 78 species | no reference labels used", fontsize=11)
    fig.text(.04, .18, "HMM checkpoint precedes multi-sequence profile expansion; this is not full-pipeline recovery.", fontsize=10)
    fig.text(.04, .12, "Top100 is a post-search cap per query and target species, not a global 100-hit limit.", fontsize=10)
    fig.text(.04, .06, "Coverage is not recall. Equal E-value cutoffs do not establish matched sensitivity or calibration.", fontsize=10)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = read_frozen(args.results, SUMMARY_SHA)
    for item in [*report["checked_records"], report["source"], report["table"]]:
        check(item)
    figure = plot(report)
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / ("search_coverage." + extension)
        figure.savefig(path, dpi=180)
        outputs.append(record(path))
    plt.close(figure)
    (args.output / "manifest.json").write_text(json.dumps({
        "results": record(args.results), "source": record(__file__), "outputs": outputs,
        "matplotlib": matplotlib.__version__, "accuracy_evaluated": False,
        "limits": report["limitations"], "publication_ready": False,
    }, indent=2, sort_keys=True) + "\n")
