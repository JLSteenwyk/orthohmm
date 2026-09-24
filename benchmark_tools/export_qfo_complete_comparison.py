"""Combine corrected comparator and publication-cell admissions without retuning."""

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.export_qfo_corrected_comparison import ENDPOINTS, RECOVERY_NOTE, extract as comparator
from benchmark_tools.export_qfo_corrected_factorial import extract as factorial
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.publication_comparison import METHODS
from benchmark_tools.run_simulation_methods import read_frozen

REPLAY_SHA = "1eec5ffb675fff234e5ce0db7e65abfa9bec09bdea8f5bd13ca72ad72d7ec40e"
PUBLICATION_CELLS = {4: "orthohmm_high_sensitivity", 7: "orthohmm_phylogeny_satellite_v2"}


def extract(report, conversion):
    if report.get("status") != "corrected_factorial_assessment_admitted":
        return comparator(report, conversion)
    row = factorial(report, conversion)
    if row["index"] not in PUBLICATION_CELLS:
        raise ValueError("Ablation cell is not a frozen publication method")
    return {"key": PUBLICATION_CELLS[row["index"]], "status": row["status"],
            "factorial_cell": row["cell"], "participant": row["participant"],
            "scores": row["scores"], "details": row["details"],
            "secondary_mean": row["secondary_mean"],
            "submitted_pairs": row["total_pairs"], "retained_pairs": row["retained_pairs"],
            "removed_mapping_pairs": row["removed_mapping_pairs"],
            "prediction_semantics": row["prediction_semantics"]}


def export(sources, replay_path, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    replay = read_frozen(replay_path, REPLAY_SHA)
    if replay["status"] != "corrected_checked_replay_admitted":
        raise ValueError("Require independently admitted corrected replay")
    replay_record = record(replay_path)
    checked, rows = [replay_record], {}
    for path, digest in sources:
        path = Path(path).resolve()
        report = read_frozen(path, digest)
        admission = record(path)
        pair_record = report["pairs_manifest"]
        check(pair_record)
        conversion = read_frozen(Path(pair_record["path"]), pair_record["sha256"])
        row = extract(report, conversion)
        if row["key"] in rows:
            raise ValueError("Duplicate corrected publication method")
        rows[row["key"]] = {**row, "admission": admission, "conversion": pair_record}
        checked.extend([admission, pair_record])
    table = []
    for key, label, _, semantics in METHODS:
        row = rows.get(key, {"key": key, "status": "not_admitted",
            "scores": {e: None for e in ENDPOINTS}, "secondary_mean": None,
            "submitted_pairs": None, "retained_pairs": None,
            "removed_mapping_pairs": None, "prediction_semantics": semantics})
        table.append({**row, "label": label})
    if set(rows) - {r["key"] for r in table}:
        raise ValueError("Method is absent from publication registry")
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "export_qfo_corrected_comparison.py", "export_qfo_corrected_factorial.py",
        "run_qfo_corrected_factorial_assessment.py", "publication_comparison.py")]
    for item in checked + helpers:
        check(item)
    output.mkdir(parents=True, exist_ok=False)
    headers = ["Method", "Status", "GO similarity", "EC similarity", "VGNC F1",
        "SwissTrees F1", "TreeFam-A F1", "FAS", "Secondary mean", "Submitted pairs",
        "Retained pairs", "Mapping losses", "Prediction semantics"]
    values = [[r["label"], r["status"], *[r["scores"][e] for e in ENDPOINTS],
        r["secondary_mean"], r["submitted_pairs"], r["retained_pairs"],
        r["removed_mapping_pairs"], r["prediction_semantics"]] for r in table]
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(values)
    limitations = [
        "Only supplied corrected-release admissions are included; missing is not zero or a scheduler-state claim.",
        "OrthoHMM high sensitivity is p1_c0_r0; phylogeny satellite_v2 is p1_c1_r1. Other factorial cells are ablations, not interchangeable defaults.",
        "The pinned independent replay admission establishes equality with the corrected native high-sensitivity partition, not historical-input outputs.",
        "GO/EC similarity and FAS are not F1; the six-metric mean is a project-defined secondary summary.",
        "Prediction semantics differ by method and are shown explicitly. Pair volume is not protein coverage.",
        "This table establishes no paired significance, ranking, independent generalization or matched efficiency.",
        "Checks cover report hashes, conversion binding and native score arithmetic; upstream admission audits are not rerun."]
    if any(r.get("search_recovery") for r in table):
        limitations.append(RECOVERY_NOTE)
    lines = ["# Corrected-Release QfO Publication Comparison", "",
        "| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in values:
        lines.append("| " + " | ".join("not admitted" if v is None else
            f"{v:.6f}" if isinstance(v, float) else str(v) for v in row) + " |")
    (output / "scores.md").write_text("\n".join(lines) + "\n\n" + "\n\n".join(limitations) + "\n")
    result = {"status": "corrected_qfo_publication_comparison", "methods": table,
        "admitted_methods": len(rows), "replay_admission": replay_record,
        "checked_records": checked, "source": record(__file__), "helpers": helpers,
        "publication_ready": False, "limitations": limitations,
        "outputs": [record(output / name) for name in ("scores.tsv", "scores.md")]}
    (output / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assessment", nargs=2, action="append", required=True)
    parser.add_argument("--replay-admission", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.assessment, args.replay_admission.resolve(), args.output.absolute())
