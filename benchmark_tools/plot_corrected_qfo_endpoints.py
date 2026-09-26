"""Plot complete corrected-QfO endpoint trade-offs from admitted retained scores."""

import argparse
import json
import math
from pathlib import Path

from benchmark_tools.audit_failed_recovery_refinement import record
from benchmark_tools import plot_publication_accuracy as helper

INPUT_SHA = "042aa221d01554114a1b0b413ca3e2ff56dad5f8c06362c69a1bf49f2de642fc"
HELPER_SHA = "92322aa7f10c3c3e96af7813b1111ae60e5903b87df04345e7cc822107533380"


def plotting_data(report):
    rows = report["methods"]
    keys = [row[0] for row in helper.METHODS]
    if (report.get("status") != "corrected_qfo_publication_comparison"
            or report.get("publication_ready") is not False or report.get("admitted_methods") != 8
            or [row["key"] for row in rows] != keys):
        raise ValueError("Require the complete corrected method inventory")
    data = []
    for row, label, color, marker in zip(rows, helper.LABELS, helper.COLORS, helper.MARKERS):
        if row["status"] != "admitted":
            raise ValueError("Unadmitted method")
        points = {}
        for metric in helper.QFO_AXES:
            name = metric.removesuffix(" F")
            detail = row["details"][name]
            score = helper.bounded(row["scores"][name])
            if metric.endswith(" F"):
                x, y = helper.bounded(detail["recall"]), helper.bounded(detail["precision"])
                expected = 2 * x * y / (x + y) if x + y else 0
                if detail["statistic"] != "F1" or not math.isclose(score, expected, abs_tol=1e-12):
                    raise ValueError("F1 disagrees with precision/recall")
            else:
                x, y = helper.bounded(detail["assessed_relations"], float("inf")), score
                if x != int(x) or detail["statistic"] != ("FAS" if name == "FAS" else "avg Schlicker"):
                    raise ValueError("Unexpected functional endpoint")
            points[metric] = dict(x=x, y=y)
        data.append(dict(key=row["key"], label=label, color=color, marker=marker, qfo=points,
                         prediction_semantics=row["prediction_semantics"]))
    return data


def render(source, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    inputs = [record(source), record(Path(helper.__file__).resolve())]
    if [item["sha256"] for item in inputs] != [INPUT_SHA, HELPER_SHA]:
        raise ValueError("Changed endpoint source or plotting helper")
    data = plotting_data(json.loads(source.read_text()))
    output.mkdir(parents=True)
    helper.plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                               "pdf.fonttype": 42, "svg.hashsalt": "corrected-qfo-complete-v1"})
    fig = helper.qfo_endpoints(data)
    fig.suptitle("Corrected QfO | Eight methods | Development-exposed", fontsize=14)
    fig.axes[-1].set_xlabel("Eligible relations (log scale)")
    fig.texts[-1].set_text("No cross-endpoint error bars: uncertainty definitions differ.\n"
                          "FAS x-axis is the reported eligible count, not the sampled count.")
    outputs = []
    for suffix in ("png", "pdf", "svg"):
        path = output / ("corrected_qfo_endpoints." + suffix)
        metadata = {"CreationDate": None, "ModDate": None} if suffix == "pdf" else {"Date": None} if suffix == "svg" else None
        fig.savefig(path, dpi=200, facecolor="white", metadata=metadata)
        if suffix == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
        outputs.append(record(path.resolve()))
    helper.plt.close(fig)
    if [record(item["path"]) for item in inputs] != inputs:
        raise ValueError("Plot inputs changed")
    manifest = dict(status="corrected_qfo_endpoints_rendered", publication_ready=False,
                    inputs=inputs, source=record(Path(__file__).resolve()), plotted_data=data,
                    outputs=outputs, matplotlib_version=helper.matplotlib.__version__,
                    limitations=["Retained-score visualization; no native scoring rerun.",
                        "QfO endpoints assess differing relation sets and prediction semantics.",
                        "FAS uses an unseeded sample; its eligible count is not its sample size.",
                        "No independence-aware uncertainty implied by scatter points."])
    with (output / "manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.source.resolve(), args.output.absolute())
