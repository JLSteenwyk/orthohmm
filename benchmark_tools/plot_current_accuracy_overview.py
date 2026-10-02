"""Reconcile a current overview with retained tables, not new native scores."""

import argparse
import csv
import json
import math
from pathlib import Path

from benchmark_tools import plot_publication_accuracy as helper
from benchmark_tools import plot_corrected_qfo_endpoints as corrected
from benchmark_tools.export_current_benchmark_scores import METRICS, index
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

SOURCES = {
    "table": ("current_benchmark_scores_20260926_v2/manifest.json", "8d56dae74721c4530ab2b86c5008b2d0ebac07878e574fd85b1ad84e4b908f07"),
    "qfo": ("qfo_corrected_comparison_20260926_v7/manifest.json", "042aa221d01554114a1b0b413ca3e2ff56dad5f8c06362c69a1bf49f2de642fc"),
    "kingdoms": ("three_kingdoms_comparison_matched_20260918/comparison.json", "ba1ed664ad8ee89ba65b72e7d4d4341eed864957b0720ab1ff4fbe91a398945a"),
    "register": ("orthobench_provenance_register_20260927/register.json", "21b1a736521ddfd87e3def36387258b8f509a73ea7ab58af8678bb101aca490d"),
    "historical": ("figures_accuracy_orthomcl_complete_20260916/manifest.json", "c9fb8d06ee37f43855b61d39963af24ff227339ff178bc45144ce951d9d89cdc"),
    "qfo_figure": ("figures_corrected_qfo_endpoints_20260926/manifest.json", "c9e0e8e93253efb6a56aff9d2ae68db22a5c834575c974a081576d5f7fde94b0"),
}


def equal(a, b):
    if (type(a) not in (int, float) or type(b) not in (int, float)
            or not math.isfinite(a) or not math.isfinite(b)
            or not math.isclose(a, b, rel_tol=0, abs_tol=1e-12)):
        raise ValueError("Retained score or coordinate disagreement")


def validate(reports):
    table, qfo, kingdoms, register, historical, qfo_figure = [reports[k] for k in SOURCES]
    if (table["status"] != "current_retained_scores_consolidated"
            or kingdoms["status"] != "three_kingdoms_comparison_recounted"
            or kingdoms["uniform_historical_input_consumption_proven"] is not False
            or register["status"] != "orthobench_provenance_register_partial"
            or register["complete_transitive_provenance"] is not False
            or register["controlled_comparative_resources"] is not False
            or any(r["publication_ready"] is not False for r in (table, qfo, kingdoms, historical, qfo_figure))):
        raise ValueError("Changed retained-evidence scope")
    inventory = [m[0] for m in helper.METHODS]
    current, qfo_rows, ob_rows = [index(r["rows"] if "rows" in r else r["methods"])
                               for r in (table, qfo, register)]
    tk_rows = index([r for r in kingdoms["rows"] if r["use"] == "comparison"])
    old = index(historical["plotted_data"])
    if any(set(r) != set(inventory) for r in (current, qfo_rows, ob_rows, tk_rows, old)):
        raise ValueError("Different method inventories")
    obsolete = next(r for r in kingdoms["rows"] if r["key"] == "sonicparanoid_2_0_9_historical")
    if (obsolete["use"] != "historical diagnostic only"
            or obsolete["input_status"] != "known Danio input mismatch"):
        raise ValueError("Changed excluded historical SonicParanoid evidence")
    equal(obsolete["counts"]["f_score"], old["sonicparanoid_2_0_9"]["three_kingdoms_f1"])
    expected_points = corrected.plotting_data(qfo)
    if qfo_figure["status"] != "corrected_qfo_endpoints_rendered" or qfo_figure["plotted_data"] != expected_points:
        raise ValueError("Corrected QfO figure differs from its complete source")
    data, plotted, changes = [], [], []
    for key, label, color, marker in zip(inventory, helper.LABELS, helper.COLORS, helper.MARKERS):
        row, native, tk = current[key], qfo_rows[key], tk_rows[key]
        if row["qfo_prediction_semantics"] != native["prediction_semantics"]:
            raise ValueError("Changed QfO output semantics")
        if set(row["scores"]) != {"OrthoBench", "ThreeKingdoms", "QfO_secondary_mean", *METRICS}:
            raise ValueError("Missing score-table endpoint")
        for value in row["scores"].values():
            equal(value, value)
            if not 0 <= value <= 1:
                raise ValueError("Changed unit scale")
        for metric in METRICS:
            equal(row["scores"][metric], native["scores"][metric])
        equal(row["scores"]["QfO_secondary_mean"], native["secondary_mean"])
        equal(native["secondary_mean"], sum(native["scores"][m] for m in METRICS) / 6)
        ob = row["orthobench_retained_evidence"]
        equal(row["scores"]["OrthoBench"], ob["f_score_percent"] / 100)
        equal(row["scores"]["OrthoBench"], ob_rows[key]["weighted_refog_f1"])
        values = {m: ob[m + "_percent"] for m in ("f_score", "precision", "recall")}
        for metric, value in values.items():
            equal(value, old[key]["orthobench"][metric])
            if not 0 <= value <= 100:
                raise ValueError("Invalid OrthoBench percent")
        counts = tk["counts"]
        for field in ("true_positive_gene_pairs", "false_positive_gene_pairs", "false_negative_gene_pairs"):
            if type(counts[field]) is not int or counts[field] < 0:
                raise ValueError("Invalid BUSCO-reference pair count")
        tp, fp, fn = [counts[k] for k in ("true_positive_gene_pairs", "false_positive_gene_pairs", "false_negative_gene_pairs")]
        if (counts["reference_genes"] != 2035 or counts["reference_orthogroups"] != 255
                or tp + fn != 7352 or tk["semantics"] != "reference-gene group co-membership, including within-species pairs"
                or row["three_kingdoms"] != {k: tk[k] for k in ("run", "input_status", "semantics", "groups")}):
            raise ValueError("Changed Three Kingdoms reference or provenance")
        equal(counts["precision"], tp / (tp + fp) if tp + fp else 0)
        equal(counts["recall"], tp / (tp + fn) if tp + fn else 0)
        equal(counts["f_score"], 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0)
        equal(row["scores"]["ThreeKingdoms"], counts["f_score"])
        if not math.isclose(old[key]["three_kingdoms_f1"], counts["f_score"], rel_tol=0, abs_tol=1e-12):
            changes.append(dict(key=key, historical_f1=old[key]["three_kingdoms_f1"],
                current_f1=counts["f_score"], difference_pp=100*(counts["f_score"]-old[key]["three_kingdoms_f1"]),
                historical_input_status="known Danio input mismatch", current_input_status=tk["input_status"]))
        data.append(dict(key=key, label=label, color=color, marker=marker, orthobench=values,
                         three_kingdoms_f1=counts["f_score"]))
        for metric in ("precision", "recall"):
            plotted.append(dict(key=key, benchmark="OrthoBench", metric=metric, value=values[metric], units="percent"))
        plotted.append(dict(key=key, benchmark="ThreeKingdoms", metric="pair_F1", value=100*counts["f_score"], units="percent"))
    if [r["key"] for r in changes] != ["sonicparanoid_2_0_9"]:
        raise ValueError("Unexpected historical-to-current overview change inventory")
    return data, plotted, changes


def plot(data):
    figure = helper.overview(data)
    figure.suptitle("Current retained group recovery | Eight methods", fontsize=14)
    figure.texts[1].set_text("BUSCO-only pairs include within-species pairs and exclude non-reference-gene errors.")
    figure.text(.5, .105, "Descriptive points only; historical input-consumption gaps remain. OrthoFinder sequence is an MCL checkpoint.",
                ha="center", fontsize=9)
    return figure


def export(results, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    reports, inputs = {}, []
    for key, (name, sha) in SOURCES.items():
        path = results / name
        pin = record(path.resolve())
        if pin["sha256"] != sha:
            raise ValueError("Changed overview source: " + key)
        reports[key] = json.loads(path.read_text())
        inputs.append(pin)
    for module, digest in ((helper, corrected.HELPER_SHA), (corrected, "db12a834f814288b590cd7a0a6f6c325ae5dd6f6e16ae532bed896fbbb90783a")):
        pin = record(Path(module.__file__).resolve())
        if pin["sha256"] != digest:
            raise ValueError("Changed retained plotting implementation")
        inputs.append(pin)
    inputs.append(record(Path(__file__).resolve()))
    data, plotted, changes = validate(reports)
    # Recheck the existing corrected figure, without rerendering or accepting its historical absolute paths.
    for saved in reports["qfo_figure"]["outputs"]:
        path = results / "figures_corrected_qfo_endpoints_20260926" / Path(saved["path"]).name
        pin = record(path.resolve())
        if (pin["bytes"], pin["sha256"]) != (saved["bytes"], saved["sha256"]):
            raise ValueError("Changed existing corrected-QfO artifact")
        inputs.append(pin)
    output.mkdir(parents=True, exist_ok=False)
    with (output / "plotted_values.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(plotted[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(plotted)
    figure = plot(data)
    try:
        figure.savefig(output / "current_accuracy_overview.png", dpi=180)
        figure.savefig(output / "current_accuracy_overview.pdf", metadata={"CreationDate": None, "ModDate": None})
    finally:
        helper.plt.close(figure)
    for pin in inputs:
        check(pin)
    manifest = dict(status="current_accuracy_overview_reconciled", checked_inputs=inputs,
        plotted_data=data, plotted_values=plotted, historical_changes=changes,
        outputs=[record(p.resolve()) for p in sorted(output.iterdir())],
        cross_dataset_score_cells_verified=72, corrected_qfo_points_verified=48,
        overview_plotted_values=24, native_predictions_or_scores_changed=False,
        uncertainty_recomputed=False, historical_input_consumption_proven=False,
        native_inference_or_scoring_rerun=False, controlled_timing=False,
        existing_corrected_qfo_figure_rerendered=False, figure_visually_reviewed=False,
        publication_ready=False, matplotlib_version=helper.matplotlib.__version__,
        limitations=["Current score table is reused unchanged; no cross-dataset mean or ranking is defined.",
            "OrthoBench is weighted group recovery; Three Kingdoms is projected BUSCO-reference pair recovery.",
            "The contemporary matched SonicParanoid row replaces only its old diagnostic in this new overview.",
            "Other Three Kingdoms historical input-consumption gaps remain; not a newly matched eight-tool experiment.",
            "Corrected QfO score/coordinate and artifact verification is not raw admission, rescoring or uncertainty.",
            "FastOMA uses supplied-tree configurations; OrthoFinder's sequence checkpoint is diagnostic."])
    with (output / "manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.results.resolve(), args.output.absolute())
