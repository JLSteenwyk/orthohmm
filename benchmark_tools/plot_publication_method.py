"""Draw the frozen publication method, distinguishing groups from ortholog pairs."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from benchmark_tools.prepare_ob_candidate_neighborhood import record

NODES = (
    ("input", .38, .925, "Protein FASTA files by species", "#eeeeee"),
    ("search", .38, .840, "k-mer candidate filter + HMM sequence search", "#dcefee"),
    ("graph", .38, .750, "Normalize hits; reciprocal graph\nLeiden CPM + singleton expansion", "#eeeeee"),
    ("profiles", .38, .655, "Build cluster profiles; HMM search\nAdd edges, recluster and refine", "#dcefee"),
    ("groups", .38, .560, "High-sensitivity orthogroups", "#eeeeee"),
    ("candidates", .38, .465, "satellite_v2 candidate-family expansion\nRetain seed-membership constraints", "#f7e7e1"),
    ("gene_trees", .38, .365, "Candidate-family alignments and gene trees\nMAFFT + FastTree", "#e8e5f0"),
    ("reconcile", .38, .255, "Root gene trees; reconcile species histories\nPositive-paralogy pair inference", "#e8e5f0"),
    ("constraints", .38, .160, "Apply high-confidence-pair membership policy\nDetach unsupported satellite memberships", "#f7e7e1"),
    ("outputs", .38, .065, "Root hierarchical orthogroups (HOGs)\nPhylogenetically inferred ortholog pairs", "#eeeeee"),
)
EDGES = tuple(zip([n[0] for n in NODES[:-1]], [n[0] for n in NODES[1:]]))


def plot():
    figure, ax = plt.subplots(figsize=(11, 13))
    figure.subplots_adjust(left=.04, right=.98, top=.88, bottom=.10)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    nodes = {name: (x, y) for name, x, y, _, _ in NODES}
    width, height = .46, .058
    for name, x, y, label, color in NODES:
        ax.add_patch(Rectangle((x - width / 2, y - height / 2), width, height,
                               facecolor=color, edgecolor="#777777", linewidth=.8))
        ax.text(x, y, label, ha="center", va="center", fontsize=10)
    for source, target in EDGES:
        x, y = nodes[source]
        tx, ty = nodes[target]
        ax.annotate("", xy=(tx, ty + height / 2), xytext=(x, y - height / 2),
                    arrowprops={"arrowstyle": "->", "color": "#555555", "lw": 1.1})
    ax.add_patch(Rectangle((.66, .315), .32, .095, facecolor="#e8e5f0", edgecolor="#777777", linewidth=.8))
    ax.text(.82, .3625, "Species-tree inference\nfrom selected candidate families\nMinimum-variance rooting",
            ha="center", va="center", fontsize=10)
    ax.plot([.61, .82], [.465, .465], color="#555555", lw=1.1)
    ax.annotate("", xy=(.82, .410), xytext=(.82, .465),
                arrowprops={"arrowstyle": "->", "color": "#555555", "lw": 1.1})
    ax.plot([.82, .82], [.315, .255], color="#555555", lw=1.1)
    ax.annotate("", xy=(.61, .255), xytext=(.82, .255),
                arrowprops={"arrowstyle": "->", "color": "#555555", "lw": 1.1})
    ax.plot([.15, .07, .07], [.465, .465, .16], linestyle="--", color="#ad604b", lw=1.1)
    ax.annotate("", xy=(.15, .16), xytext=(.07, .16),
                arrowprops={"arrowstyle": "->", "color": "#ad604b", "lw": 1.1, "linestyle": "--"})
    ax.text(.045, .31, "Constraint trace", rotation=90, ha="center", va="center", fontsize=10, color="#854632")
    figure.suptitle("OrthoHMM: frozen publication workflow", x=.04, ha="left", y=.97, fontsize=18)
    figure.text(.04, .925, "HMM-centered grouping followed by candidate expansion and inferred-tree refinement", fontsize=11)
    figure.text(.04, .070, "Schematic of the full satellite_v2 configuration; high-sensitivity output stops at orthogroups.", fontsize=10)
    figure.text(.04, .043, "Small families may bypass tree inference. Exact-input checkpoints change cost, not the intended analysis.", fontsize=10)
    figure.text(.04, .016, "Output levels are distinct: group co-membership is not itself a phylogenetic orthology assignment.", fontsize=10)
    return figure


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    root = args.root.resolve()
    evidence = [record(root / "benchmark_tools/results" / name) for name in (
        "orthobench_factorial_prepared_20260916.json", "orthobench_factorial_results_20260916.json",
        "ob_family_trace_verified_20260916.json", "ob_reconciliation_trace_20260916.json")]
    evidence.append(record(root / "benchmarks/work/publication_method_native_v2/benchmark_tools/replay_high_sensitivity.py"))
    figure = plot()
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / ("publication_method." + extension)
        figure.savefig(path, dpi=180)
        outputs.append(record(path))
    plt.close(figure)
    (args.output / "manifest.json").write_text(json.dumps({"source": record(__file__), "evidence": evidence,
        "outputs": outputs, "matplotlib": matplotlib.__version__, "scope": "Conceptual frozen-method diagram, not execution or accuracy validation"},
        indent=2, sort_keys=True) + "\n")
