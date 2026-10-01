"""Export admitted seven-arm QfO scores and reproduced SwissTrees uncertainty."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmark_tools import audit_qfo_parameter_swiss as legacy
from benchmark_tools import run_private_cpm_parameter_uncertainty as recovered
from benchmark_tools.bootstrap_qfo_parameter_neighborhood import ARMS, validated_values
from benchmark_tools.bootstrap_qfo_swiss_stages import aggregate
from benchmark_tools.export_qfo_threshold_endpoints import ENDPOINTS, scores
from benchmark_tools.plot_ob_parameter_neighborhood import LABELS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

METRICS = ("F1", "PPV", "TPR")
REPRODUCER_SHA = "6b337c936bf1cf6896fea06de6e1924092c1dc6e5b3871ad621edbdc2138ddf8"
SOURCES = {"run_qfo_parameter_uncertainty.py": "2c551e6a961ce130f7973c064096056beef159ca2b1212bd0a414169b962c719",
           "run_private_cpm_parameter_uncertainty.py": "608ac7512993280190350d15c77aab56b26770e885dcae0bba62c88354d0a6ca"}


def validate(report, reproduction):
    controls = dict(status="corrected_qfo_parameter_uncertainty_audited", replicates=100000,
        seed=20260925, alpha=.05, multiplicity_endpoints=18, quantile_method="linear",
        units="raw 0-to-1 metric units")
    if (any(report.get(key) != value for key, value in controls.items())
            or report["scientific_inputs_admitted"] is not True or report["protocol_controls_match"] is not True
            or report["publication_ready"] is not False or report["protocol"]["sha256"] != recovered.PROTOCOL_SHA
            or report["plan"]["sha256"] != recovered.PLAN_SHA
            or report["source"]["sha256"] not in SOURCES.values()):
        raise ValueError("Require admitted frozen parameter uncertainty")
    if (reproduction["status"] != "qfo_parameter_uncertainty_numerically_reproduced"
            or reproduction["source"]["sha256"] != REPRODUCER_SHA
            or reproduction["planned_endpoints"] != 18 or reproduction["absolute_tolerance"] != 1e-12
            or reproduction["publication_ready"] is not False):
        raise ValueError("Require unchanged independent numerical reproduction")
    values = validated_values(report["reconstructed_counts"])
    expected_arms = [{key: row[key] for key in ("arm", "status", "reason") if key in row}
                     for row in report["reconstructed_counts"]["arms"]]
    if (report["arms"] != expected_arms or report["families"] != report["reconstructed_counts"]["families"]
            or set(report["point_estimates"]) != set(ARMS) or len(report["comparisons"]) != 6):
        raise ValueError("Changed parameter arm or family inventory")
    points = {}
    for arm in ARMS:
        saved = report["point_estimates"][arm]
        if arm not in values:
            if saved is not None:
                raise ValueError("Unavailable parameter score imputed")
            points[arm] = None
            continue
        expected = aggregate(values[arm].mean(axis=0))
        if (set(saved) != set(METRICS) or any(type(saved[key]) not in (int, float)
                or not math.isfinite(saved[key]) or not math.isclose(saved[key], expected[index], abs_tol=1e-12, rel_tol=0)
                for index, key in enumerate(METRICS))):
            raise ValueError("Displayed points differ from reconstructed counts")
        points[arm] = expected
    estimated = 0
    for row, arm in zip(report["comparisons"], ARMS[1:]):
        if (row["candidate"], row["reference"]) != (arm, "control"):
            raise ValueError("Changed parameter contrast identity")
        if arm not in values or "control" not in values:
            reason = "baseline_not_admitted" if "control" not in values else "variant_not_admitted"
            if (row["status"] != "not_estimable" or row["reason"] != reason
                    or row["metrics"] is not None or row["family_differences"] is not None):
                raise ValueError("Unavailable contrast must remain explicit and null")
            continue
        estimated += 1
        if row["status"] != "estimated" or set(row["metrics"]) != set(METRICS):
            raise ValueError("Available contrast missing")
        for index, key in enumerate(METRICS):
            item = row["metrics"][key]
            nominal, adjusted = item["paired_percentile_ci"], item["bonferroni_percentile_ci"]
            if (len(nominal) != 2 or len(adjusted) != 2
                    or any(type(value) not in (int, float) or not math.isfinite(value)
                           for value in (item["difference"], *nominal, *adjusted))
                    or not -1 <= adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1] <= 1
                    or not math.isclose(item["difference"], points[arm][index] - points["control"][index], abs_tol=1e-12, rel_tol=0)):
                raise ValueError("Invalid intervals or contrast arithmetic")
            directions = [item[key] for key in ("family_wins", "family_ties", "family_losses")]
            if any(type(value) is not int or value < 0 for value in directions) or sum(directions) != 18:
                raise ValueError("Invalid family direction counts")
    if (type(report["estimated_contrasts"]) is not int or report["estimated_contrasts"] != estimated
            or report["complete_panel"] is not (estimated == 6) or report["uncertainty_admitted"] is not (estimated > 0)
            or type(reproduction["endpoints"]) is not int or reproduction["endpoints"] != estimated * 3):
        raise ValueError("Incorrect panel completeness or reproduced endpoint count")
    return points


def plot(report, reproduction):
    points = validate(report, reproduction)
    figure = plt.figure(figsize=(12, 9))
    grid = figure.add_gridspec(2, 3, left=.21, right=.97, top=.83, bottom=.21,
        height_ratios=[1, 1.2], hspace=.5, wspace=.32)
    table_ax = figure.add_subplot(grid[0, :])
    table_ax.axis("off")
    rows = [[LABELS[arm], *(["NA"] * 3 if points[arm] is None else [f"{value * 100:.3f}" for value in points[arm]])] for arm in ARMS]
    table = table_ax.table(cellText=rows, colLabels=["Configuration", "F1 (%)", "Precision (%)", "Recall (%)"],
        colWidths=[.43, .19, .19, .19], cellLoc="center", bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    for (row, column), cell in table.get_celld().items():
        cell.set_edgecolor("#dddddd")
        cell.set_linewidth(.5)
        cell.set_facecolor("#eeeeee" if row == 0 else "white")
        if row == 0:
            cell.set_text_props(weight="bold")
    extent = max([1.] + [abs(value * 100) for row in report["comparisons"] if row["status"] == "estimated"
        for key in METRICS for value in row["metrics"][key]["bonferroni_percentile_ci"]])
    limit = float(np.ceil(extent * 1.1))
    for index, key in enumerate(METRICS):
        ax = figure.add_subplot(grid[1, index])
        ax.axvline(0, color="#888888", linewidth=1)
        for y, row in enumerate(report["comparisons"]):
            color = ("#007d83", "#b64b42", "#6b5b95")[y // 2]
            if row["status"] == "not_estimable":
                ax.text(0, y, "Not estimable", ha="center", va="center", fontsize=8, color="#666666",
                    bbox=dict(facecolor="white", edgecolor="none", pad=2))
                continue
            item = row["metrics"][key]
            ax.plot(np.asarray(item["bonferroni_percentile_ci"]) * 100, [y, y], color=color, linewidth=1.3)
            ax.plot(np.asarray(item["paired_percentile_ci"]) * 100, [y, y], color=color, linewidth=4, solid_capstyle="butt")
            ax.plot(item["difference"] * 100, y, "o", color=color, markersize=4)
        ax.set_xlim(-limit, limit)
        ax.set_ylim(5.6, -.6)
        ax.set_yticks(range(6), [LABELS[arm] for arm in ARMS[1:]] if index == 0 else [""] * 6)
        ax.set_title(("F1", "Precision", "Recall")[index], loc="left", fontsize=12)
        ax.set_xlabel("Variant minus control (pp)", fontsize=10)
        ax.grid(axis="x", alpha=.2)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    figure.suptitle("QfO SwissTrees: limited parameter sensitivity", x=.04, ha="left", y=.97, fontsize=17)
    figure.text(.04, .92, "Six one-at-a-time changes of 20% | 18 reference families | development-exposed", fontsize=11)
    figure.text(.04, .875, f"{report['estimated_contrasts']}/6 contrasts estimated; unavailable arms are not zero. No default selection.", fontsize=10)
    notes = ["Thick: paired 95% CI. Thin: Bonferroni CI over all 18 planned endpoints; 100,000 draws, seed 20260925.",
        "F1 is recomputed from macro precision and recall per draw; missing contrasts retain their places.",
        "No independent confirmation, selection adjustment, OrthoFinder superiority test or controlled timing is shown.",
        "Intervals including zero do not establish equivalence. No other-endpoint or secondary-mean uncertainty is shown."]
    for y, note in zip((.145, .105, .065, .025), notes):
        figure.text(.04, y, note, fontsize=9)
    return figure


def native_rows(report, checked):
    pin = report["admission_inventory"]
    inventory = read_frozen(Path(pin["path"]), pin["sha256"])
    if inventory["status"] != "qfo_parameter_score_admission_inventory" or [row["arm"] for row in inventory["arms"]] != list(ARMS):
        raise ValueError("Wrong native parameter admission inventory")
    checked.append(pin)
    if report["source"]["sha256"] == SOURCES["run_private_cpm_parameter_uncertainty.py"]:
        prior_pin = next(ref for ref in report["checked_inputs"] if ref["path"].endswith("/" + recovered.PRIOR))
        recovered.validate_inventory(inventory, read_frozen(Path(prior_pin["path"]), prior_pin["sha256"]))
        checked.append(prior_pin)
    rows = []
    for entry, availability in zip(inventory["arms"], report["arms"]):
        arm = entry["arm"]
        row = dict(arm=arm, status=entry["status"], reason=entry.get("reason", ""))
        if entry["status"] == "not_admitted":
            if availability != entry:
                raise ValueError("Unavailable arm differs from audited counts")
            row.update({key: None for key in (*ENDPOINTS, "secondary_mean")})
        elif entry["status"] == "admitted":
            if availability != {"arm": arm, "status": "counts_verified"}:
                raise ValueError("Native arm lacks verified family counts")
            checked.append(entry["admission"])
            admitted = read_frozen(Path(entry["admission"]["path"]), entry["admission"]["sha256"])
            pair = admitted["pairs_manifest"]
            conversion = read_frozen(Path(pair["path"]), pair["sha256"])
            checked.extend([pair, *legacy.file_records(admitted)])
            if admitted["status"] == "private_recovered_cpm_assessment_admitted":
                if arm != "cpm_high" or report["source"]["sha256"] != SOURCES["run_private_cpm_parameter_uncertainty.py"]:
                    raise ValueError("Private admission cannot enter a legacy or different arm")
                candidate_pin = conversion["candidate_admission"]
                candidate = read_frozen(Path(candidate_pin["path"]), candidate_pin["sha256"])
                recovered.validate_private(admitted, conversion, candidate)
                checked.extend([candidate_pin, entry["admission_submission"]])
            else:
                legacy.validate_admission(arm, admitted, conversion)
            row.update(scores(admitted["assessment"]))
            if not math.isclose(row["SwissTrees"], report["point_estimates"][arm]["F1"], rel_tol=0, abs_tol=1e-6):
                raise ValueError("Native SwissTrees table differs from reconstructed counts")
        else:
            raise ValueError("Unknown native parameter availability")
        rows.append(row)
    return rows


def write_tables(output, report, rows):
    interval_rows = []
    for row in report["comparisons"]:
        for key in METRICS:
            item = row["metrics"][key] if row["status"] == "estimated" else None
            interval_rows.append(dict(candidate=row["candidate"], reference="control", metric=key, status=row["status"],
                difference_pp=None if item is None else 100 * item["difference"],
                nominal_low_pp=None if item is None else 100 * item["paired_percentile_ci"][0],
                nominal_high_pp=None if item is None else 100 * item["paired_percentile_ci"][1],
                adjusted_low_pp=None if item is None else 100 * item["bonferroni_percentile_ci"][0],
                adjusted_high_pp=None if item is None else 100 * item["bonferroni_percentile_ci"][1],
                family_wins=None if item is None else item["family_wins"],
                family_ties=None if item is None else item["family_ties"],
                family_losses=None if item is None else item["family_losses"], reason=row.get("reason", "")))
    for name, data in (("scores.tsv", rows), ("intervals.tsv", interval_rows)):
        with (output / name).open("x") as stream:
            writer = csv.DictWriter(stream, list(data[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(data)
    lines = ["# Corrected QfO Parameter Neighborhood", "",
        "| Arm | Status | GO similarity | EC similarity | VGNC F1 | SwissTrees F1 | TreeFam-A F1 | FAS | Secondary mean |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in rows:
        lines.append("| " + " | ".join([row["arm"], row["status"], *["NA" if row[key] is None else f"{row[key]:.6f}" for key in (*ENDPOINTS, "secondary_mean")]]) + " |")
    lines.extend(["", "GO/EC similarity and FAS are not F1. The six-score mean is a project-defined secondary summary.",
        "", "SwissTrees table values use rounded native aggregates; paired intervals use reconstructed family counts.",
        "", "Unavailable arms are not zero. All 18 planned interval endpoints are retained.",
        "", "FAS sampling is unseeded; small FAS differences cannot be attributed to parameters alone.",
        "", "Development-exposed sensitivity analysis, not independent validation, equivalence or general superiority."])
    with (output / "scores.md").open("x") as stream:
        stream.write("\n".join(lines) + "\n")
    return interval_rows


def export(path, sha, reproduction_path, reproduction_sha, output):
    if not output.is_absolute() or output.resolve() != output:
        raise ValueError("Require direct absolute parameter figure destination")
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    report, reproduction = read_frozen(path, sha), read_frozen(reproduction_path, reproduction_sha)
    input_pin = record(path)
    if reproduction["input"] != input_pin:
        raise ValueError("Numerical reproduction belongs to a different analysis")
    validate(report, reproduction)
    checked = [record(__file__), input_pin, record(reproduction_path), reproduction["source"],
        report["source"], *report["helpers"], *report["checked_inputs"],
        *[record(module.__file__) for name, module in sorted(sys.modules.items())
          if name.startswith("benchmark_tools.") and getattr(module, "__file__", None)]]
    for ref in checked:
        check(ref)
    rows = native_rows(report, checked)
    for ref in checked:
        check(ref)
    figure = plot(report, reproduction)
    output.mkdir(parents=True, exist_ok=False)
    try:
        interval_rows = write_tables(output, report, rows)
        for extension in ("png", "pdf", "svg"):
            figure.savefig(output / ("qfo_parameter_neighborhood." + extension), dpi=180)
    finally:
        plt.close(figure)
    for ref in checked:
        check(ref)
    outputs = [record(output / name) for name in ("scores.tsv", "scores.md", "intervals.tsv",
        "qfo_parameter_neighborhood.png", "qfo_parameter_neighborhood.pdf", "qfo_parameter_neighborhood.svg")]
    result = dict(status="qfo_parameter_results_exported", source=record(__file__), input=input_pin,
        reproduction=record(reproduction_path), rows=rows, intervals=interval_rows, outputs=outputs,
        checked_records=checked, planned_endpoints=18, estimated_endpoints=report["estimated_contrasts"] * 3,
        complete_panel=report["complete_panel"], matplotlib=matplotlib.__version__, publication_ready=False,
        limitations=["Visualization/export of admitted results, not a new inference, score or uncertainty admission.",
            "Reuse of a hash-bound numerical reproduction; no bootstrap rerun.",
            "Only SwissTrees has family-level intervals here; development exposure and family dependence remain.",
            "No default selection, equivalence, superiority or controlled-resource claim."])
    with (output / "manifest.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("results", "reproduction", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    for name in ("sha256", "reproduction-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    export(args.results.resolve(), args.sha256, args.reproduction.resolve(), args.reproduction_sha256, args.output)
