"""Plot independently recounted, reference-conditioned search transitions."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_simulation_methods import read_frozen

AUDIT_SHA = "0cda6c99b7bab61c86f2319653828fe6ada74c118166408dbab3535bf003d4bd"
STATES = ("accepted", "scored_not_significant", "not_selected_by_prefilter")


def counts(report):
    if (report["status"] != "forced_candidate_scores_independently_recounted"
            or report["candidate_only_interpretation_authorized"] is not True
            or report["numerical_disagreements"] != []
            or report["directions"] != 144 or report["directed_pairs"] != 81466):
        raise ValueError("Require complete, numerically unchanged diagnostic")
    transitions = report["forced_transitions"]
    allowed = {a + ":" + b for a in STATES for b in STATES[:2]}
    if set(transitions) - allowed or any(type(v) is not int or v < 0 for v in transitions.values()):
        raise ValueError("Invalid transition counts")
    matrix = [[transitions.get(a + ":" + b, 0) for b in STATES[:2]] for a in STATES]
    if (sum(map(sum, matrix)) != report["directed_pairs"]
            or sum(map(sum, matrix[:2])) != report["previously_scored"]
            or matrix[0][1] or matrix[1][0]
            or any(sum(row[j] for row in matrix) != report["decisions"][s]
                   for j, s in enumerate(STATES[:2]))):
        raise ValueError("Transition marginals disagree")
    return matrix


def plot(report):
    matrix = counts(report)
    fig, ax = plt.subplots(figsize=(10, 5.2))
    fig.subplots_adjust(left=.27, right=.96, top=.70, bottom=.30)
    colors = ("#007d83", "#b64b42")
    for y, row in enumerate(matrix):
        left = 0
        for j, value in enumerate(row):
            ax.barh(y, value, left=left, color=colors[j], height=.55,
                    label=("Pass E < 1e-4", "Do not pass")[j] if y == 0 else None)
            if value:
                inside = value >= 5000
                ax.text(left + value / 2 if inside else left + value + 600, y,
                        f"{value:,}", ha="center" if inside else "left", va="center",
                        color="white" if inside else "#222222", fontsize=11)
            left += value
    ax.set_yticks(range(3), ["Previously accepted", "Previously nonsignificant", "Previously prefilter-excluded"])
    ax.invert_yaxis()
    ax.set_xlim(0, 51000)
    ax.set_xticks([0, 10000, 20000, 30000, 40000, 50000], ["0", "10,000", "20,000", "30,000", "40,000", "50,000"])
    ax.set_xlabel("Directed reference-family pairs after forced scoring", fontsize=10)
    ax.tick_params(length=0)
    ax.set_axisbelow(True)
    ax.grid(axis="x", alpha=.18)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.suptitle("OrthoBench: prefilter exclusions hide significant HMM hits", x=.04,
                 ha="left", y=.97, fontsize=15)
    fig.text(.04, .88, "81,466 watched pairs | 144 species directions | frozen HMM scoring and full target databases", fontsize=10)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(.26, .84), ncol=2, frameon=False)
    fig.text(.04, .15, "All 33,098 previously scored pairs retain exactly identical scores, E-values and decisions.", fontsize=10)
    fig.text(.04, .08, "Reference-conditioned diagnostic; not unbiased recall, true-ortholog counts or a whole-pipeline F1 gain.", fontsize=10)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    report = read_frozen(args.audit, AUDIT_SHA)
    figure = plot(report)
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "svg", "pdf"):
        path = args.output / ("forced_candidates." + extension)
        figure.savefig(path, dpi=180)
        outputs.append(record(path))
    plt.close(figure)
    (args.output / "manifest.json").write_text(json.dumps(dict(
        audit=record(args.audit), source=record(__file__), outputs=outputs,
        matplotlib=matplotlib.__version__, counts=counts(report),
        limitations=report["limitations"]), indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
