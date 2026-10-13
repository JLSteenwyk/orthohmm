"""Localize development-exposed OrthoBench changes from independently read-back cells."""

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json

COUNTS = ("true_positive", "false_positive", "false_negative")


def reconstruct(details):
    families, totals = {}, {key: 0.0 for key in COUNTS}
    for row in details["refog_records"]:
        name, genes = row["refog"], row["genes"]
        if name in families or type(genes) is not int or genes < 2:
            raise ValueError("Duplicate family or invalid scored gene count")
        if any(type(row[k]) not in (int, float) or not math.isfinite(row[k]) or row[k] < 0 for k in COUNTS):
            raise ValueError("Invalid family pair counts")
        families[name] = row
        for key in COUNTS:
            totals[key] += row[key] / (genes - 1)
    if not families:
        raise ValueError("Empty reference family inventory")
    for key in COUNTS:
        if not math.isclose(totals[key], details["weighted_" + key], rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("Family counts do not reproduce weighted totals")
    tp, fp, fn = (totals[key] for key in COUNTS)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0
    for key, value in (("precision", precision), ("recall", recall), ("f_score", f1)):
        if not math.isclose(value * 100, details[key], rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("Reconstructed OrthoBench endpoint differs")
    return families, {"F1": f1, "precision": precision, "recall": recall, "weighted_counts": totals}


def localize(summary):
    cells = summary["cells"]
    if (len(cells) != 2 or any(c.get("dataset") != "OrthoBench"
            or c.get("status") != "independent_content_readback_verified" for c in cells)):
        raise ValueError("Require two independently verified OrthoBench cells in control/variant order")
    arms = [reconstruct(c["details"]["family_details"]) for c in cells]
    left, right = arms[0][0], arms[1][0]
    if left.keys() != right.keys() or any(left[name]["genes"] != right[name]["genes"] for name in left):
        raise ValueError("Reference families or scored gene counts changed")
    for cell, (_, scores) in zip(cells, arms):
        if not math.isclose(cell["scores"]["OrthoBench: OrthoBench F1"], scores["F1"], rel_tol=1e-12):
            raise ValueError("Claimed cell F1 differs from family reconstruction")
    changes = []
    for name in sorted(left):
        a, b = left[name], right[name]
        delta = {key: b[key] - a[key] for key in (*COUNTS, "splits")}
        if any(delta.values()):
            changes.append({"refog": name, "genes": a["genes"], "control": a, "variant": b,
                "delta": delta, "weighted_delta": {key: delta[key] / (a["genes"] - 1) for key in COUNTS}})
    return {"status": "verified_family_count_localization_not_independent_validation", "families": len(left),
        "changed_families": len(changes), "control": arms[0][1], "variant": arms[1][1],
        "F1_delta": arms[1][1]["F1"] - arms[0][1]["F1"],
        "weighted_delta": {key: arms[1][1]["weighted_counts"][key] - arms[0][1]["weighted_counts"][key] for key in COUNTS},
        "changes": sorted(changes, key=lambda row: (-row["weighted_delta"]["false_positive"], row["refog"])),
        "limitations": ["Previously examined references are development-exposed; this is mechanistic diagnosis, not confirmation.",
            "Per-family weighted count changes are additive; F1 is nonlinear and has no claimed additive family attribution.",
            "Do not introduce reference-family-specific inference rules or selectively exclude declining families."],
        "defaults_promoted": False, "goal_completed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Retain prior localization evidence")
    identity = record(args.summary)
    summary = json.loads(args.summary.read_text())
    full_identity = summary["report"]
    check(full_identity)
    full_report = json.loads(Path(full_identity["path"]).read_text())
    checked = [identity, full_identity, *(item for cell in full_report["cells"] for item in cell["checked"])]
    for item in checked:
        check(item)
    result = localize(summary)
    if result != localize(full_report):
        raise ValueError("Compact summary differs from the pinned full readback")
    for item in checked:
        check(item)
    result.update(input_summary=identity, full_readback=full_identity,
                  source=record(Path(__file__)), input_pins_reverified=True)
    write_json(args.output, result)
    print(json.dumps({key: result[key] for key in ("status", "changed_families", "F1_delta", "weighted_delta")}))


if __name__ == "__main__":
    main()
