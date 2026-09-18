"""Export admitted sequence-search controls without claiming HMM superiority."""

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.export_qfo_corrected_comparison import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_qfo_native_assessment import validate_directory

VARIANTS = ("all_hits", "top100")


def binding(report):
    variant = report["variant"]
    if (variant not in VARIANTS or report["status"] != "corrected_sequence_assessment_admitted"
            or report["accuracy_admitted"] is not True or report["publication_ready"] is not False):
        raise ValueError("Require admitted sequence-control assessment")
    conversion = report["conversion"]
    participant = "ohmm_qfo_corrected_sequence_" + variant
    if (conversion["variant"] != variant or conversion["participant"] != participant
            or conversion["status"] != "corrected_sequence_group_pairs_prepared_unscored"
            or conversion["semantics"] != "cross-species group-derived clique pairs"
            or report["assessment"]["participant"] != participant):
        raise ValueError("Sequence identity or prediction semantics differs")
    counts = [conversion[key] for key in ("total_pairs", "retained_pairs", "expected_pairs", "removed_mapping_pairs")]
    if any(type(value) is not int or value < 0 for value in counts) or not counts[0] == counts[1] == counts[2] or counts[3]:
        raise ValueError("Invalid pair accounting")
    return variant, participant


def export(sources, output):
    if output.exists():
        raise FileExistsError(output)
    rows, checked = {}, []
    for path, sha in sources:
        path = Path(path).resolve()
        report = read_frozen(path, sha)
        variant, participant = binding(report)
        if variant in rows:
            raise ValueError("Duplicate sequence variant")
        evidence = [record(path), report["pairs_manifest"], report["execution_report"],
                    report["environment_manifest"], report["native_trace"], *report["metric_files"]]
        for item in evidence:
            check(item)
        conversion = json.loads(Path(report["pairs_manifest"]["path"]).read_text())
        if conversion != report["conversion"]:
            raise ValueError("Conversion manifest differs from admission")
        execution = json.loads(Path(report["execution_report"]["path"]).read_text())
        environment = json.loads(Path(report["environment_manifest"]["path"]).read_text())
        assessment, paths = validate_directory(Path(execution["results"]), participant,
                                               Path(environment["pipeline"]) / "reference_data")
        if assessment != report["assessment"]:
            raise ValueError("Native metrics differ from admitted assessment")
        evidence.extend(record(path) for path in paths)
        row = dict(variant=variant, status="admitted", admission=record(path),
                   submitted_pairs=conversion["total_pairs"], removed_mapping_pairs=0,
                   semantics=conversion["semantics"],
                   scores={key: assessment["endpoints"][key]["score"] for key in ENDPOINTS},
                   endpoint_details=assessment["endpoints"],
                   secondary_mean=assessment["secondary_six_metric_mean"])
        rows[variant] = row
        checked.extend(evidence)
    table = [rows.get(variant, dict(variant=variant, status="not_admitted",
                scores={key: None for key in ENDPOINTS}, secondary_mean=None, submitted_pairs=None))
             for variant in VARIANTS]
    for item in checked:
        check(item)
    limits = ["DIAMOND sequence-search controls with frozen OrthoHMM downstream grouping, not standalone competing tools.",
        "GO/EC similarity and FAS are not F1; the six-metric mean is project-defined and secondary.",
        "Native metric standard errors are not paired method-difference confidence intervals.",
        "Equal search cutoffs do not establish matched sensitivity or computational effort.",
        "The HMM baseline and paired uncertainty are not supplied by this table; no HMM advantage is established.",
        "Pending means no admitted score, not zero and not a scheduler-state assertion.",
        "Rechecks admitted hashes and native metric content, not the entire inference workflow."]
    output.mkdir(parents=True, exist_ok=False)
    header = ["Variant", "Status", "GO similarity", "EC similarity", "VGNC F1", "SwissTrees F1",
              "TreeFam-A F1", "FAS", "Secondary mean", "Submitted pairs"]
    values = [[row["variant"], row["status"], *[row["scores"][key] for key in ENDPOINTS],
               row["secondary_mean"], row["submitted_pairs"]] for row in table]
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(["pending" if value is None else value for value in row] for row in values)
    lines = ["# Corrected QfO Sequence-Search Controls", "", "| " + " | ".join(header) + " |",
             "| " + " | ".join(["---"] * len(header)) + " |"]
    for row in values:
        lines.append("| " + " | ".join("pending" if value is None else f"{value:.6f}"
                      if isinstance(value, float) else str(value) for value in row) + " |")
    lines.extend(["", *["- " + text for text in limits]])
    (output / "scores.md").write_text("\n".join(lines) + "\n")
    result = dict(status="admitted_sequence_control_scores", rows=table, checked_records=checked,
        source=record(__file__), helpers=[record(Path(__file__).with_name(name)) for name in (
            "validate_qfo_native_assessment.py", "qfo_summarize_scores.py", "export_qfo_corrected_comparison.py")],
        limitations=limits, publication_ready=False,
        outputs=[record(output / name) for name in ("scores.tsv", "scores.md")])
    (output / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assessment", nargs=2, action="append", required=True, metavar=("PATH", "SHA256"))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    export(args.assessment, args.output.resolve())
