"""Independently verify corrected VGNC deletions using exact rational ratios."""

import argparse
import csv
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CATEGORIES = ("TP", "FP", "FN")
MAPPING_SHA = "01a51a62a536bcf759f2a0ef83c56196a92022b0bb8290448eacbfefa47d6005"


def ratios(counts):
    tp, fp, fn = counts
    if any(type(v) is not int or v < 0 for v in counts) or not (tp + fp and tp + fn):
        raise ValueError("Invalid or undefined count ratios")
    return (Fraction(tp, tp + fp), Fraction(tp, tp + fn), Fraction(2 * tp, 2 * tp + fp + fn))


def close(observed, expected):
    if type(observed) not in (int, float) or not math.isfinite(observed):
        raise ValueError("Nonfinite or nonnumeric reported value")
    error = abs(observed - float(expected))
    if error > 1e-12:
        raise ValueError("Reported arithmetic differs")
    return error


def reconstruct(path, blocks):
    incident = {b: [0, 0, 0] for b in blocks}
    total, seen = [0, 0, 0], set()
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["block_left", "block_right", *CATEGORIES]:
            raise ValueError("Wrong sparse table columns")
        for row in reader:
            left, right = row["block_left"], row["block_right"]
            counts = [int(row[c]) for c in CATEGORIES]
            if (left not in blocks or right not in blocks or left > right
                    or (left, right) in seen or min(counts) < 0 or not sum(counts)
                    or (left != right and (counts[0] or counts[2]))):
                raise ValueError("Invalid sparse cell")
            seen.add((left, right))
            for i, n in enumerate(counts):
                total[i] += n
                incident[left][i] += n
                if right != left:
                    incident[right][i] += n
    deleted = {b: (removed, ratios([n - r for n, r in zip(total, removed)]))
               for b, removed in incident.items()}
    return total, ratios(total), deleted


def verify(mapping_path, report_path, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    mapping_ref, report_ref = record(mapping_path), record(report_path)
    if mapping_ref["sha256"] != MAPPING_SHA:
        raise ValueError("Wrong corrected mapping")
    mapping, report = json.loads(mapping_path.read_text()), json.loads(report_path.read_text())
    refs = [mapping_ref, report_ref, report["table"], *report["checked_records"],
            mapping["reference_table"], *[m["table"] for m in mapping["methods"]], record(Path(__file__))]
    for ref in refs:
        check(ref)
    with Path(mapping["reference_table"]["path"]).open() as handle:
        labels = [r["block"] for r in csv.DictReader(handle, delimiter="\t")]
    blocks = set(labels)
    keys = [m["key"] for m in mapping["methods"]]
    baseline = "orthofinder_3_1_5_full"
    if (len(labels) != 16844 or len(blocks) != len(labels) or len(keys) != 8 or len(set(keys)) != 8
            or report["baseline"] != baseline or set(report["methods"]) != set(keys)
            or set(report["comparisons"]) != set(keys) - {baseline}
            or report["rows"] != 134752):
        raise ValueError("Wrong result panel")
    expected, full = {}, {}
    error = 0.
    for method in mapping["methods"]:
        key = method["key"]
        total, full[key], expected[key] = reconstruct(Path(method["table"]["path"]), blocks)
        if dict(zip(CATEGORIES, total)) != report["methods"][key]["counts"]:
            raise ValueError("Wrong full counts")
        for metric, value in zip(("precision", "recall", "f1"), full[key]):
            error = max(error, close(report["methods"][key]["full_metrics"][metric], value))
    seen = set()
    with gzip.open(report["table"]["path"], "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key, block = row["method"], row["block"]
            if key not in expected or block not in blocks or (key, block) in seen:
                raise ValueError("Unknown or duplicate deletion row")
            seen.add((key, block))
            removed, values = expected[key][block]
            if [int(row["removed_" + c]) for c in CATEGORIES] != removed:
                raise ValueError("Incident counts differ")
            for metric, value in zip(("precision", "recall", "f1"), values):
                error = max(error, close(float(row[metric]), value))
            error = max(error, close(float(row["f1_change"]), values[2] - full[key][2]))
    if seen != {(key, block) for key in keys for block in blocks}:
        raise ValueError("Incomplete deletion table")
    for key, observed in report["comparisons"].items():
        differences = {b: expected[key][b][1][2] - expected[baseline][b][1][2] for b in blocks}
        for field, value in (("full_difference", full[key][2] - full[baseline][2]),
                ("minimum_deleted_difference", min(differences.values())),
                ("maximum_deleted_difference", max(differences.values()))):
            error = max(error, close(observed[field], value))
        for field, count in (("positive", sum(v > 0 for v in differences.values())),
                ("negative", sum(v < 0 for v in differences.values())),
                ("zero", sum(v == 0 for v in differences.values()))):
            if type(observed[field]) is not int or observed[field] != count:
                raise ValueError("Wrong contrast sign counts")
        for direction in ("minimum", "maximum"):
            error = max(error, close(observed[direction + "_deleted_difference"],
                                    differences[observed[direction + "_block"]]))
    for ref in refs:
        check(ref)
    result = dict(status="corrected_vgnc_rational_replay_passed", rows_checked=len(seen),
                  comparisons_checked=7, maximum_absolute_error=error, checked_records=refs,
                  uncertainty_admitted=False, publication_ready=False,
                  scope="Full counts/ratios, all deletion counts/ratios/F1 changes, paired ranges, signs and extremal block values; not top-ten influence ranking or statistical validity.")
    with output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("mapping", "report", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    verify(args.mapping, args.report, args.output)
