"""Independent rational-arithmetic check and full F1 display for the OB export."""

import argparse
import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_counts(records, selected):
    return [sum((Fraction(str(records[n][key])) / (records[n]["genes"] - 1)
                 for n in selected), Fraction())
            for key in ("true_positive", "false_positive", "false_negative")]


def verify(base, export):
    report = json.loads((export / "report.json").read_text())
    frozen = json.loads((base / "ob_error_strata_prepared_20260916.json").read_text())
    main = json.loads((base / "orthobench_paired_uncertainty_20260916.json").read_text())
    others = json.loads((base / "retained_ob_comparator_readback_20260926.json").read_text())
    scores = dict(main["scores"])
    scores.update({r["key"]: r["score"] for r in others["rows"]})
    methods = list(scores)
    expected = {(label, method) for label in frozen["strata"] for method in methods}
    rows = report["rows"]
    assert len(rows) == len(expected) == 112
    assert {(r["stratum"], r["method"]) for r in rows} == expected
    for ref in report["checked_records"]:
        assert digest(Path(ref["path"])) == ref["sha256"]
    indexed = {m: {r["refog"]: r for r in s["refog_records"]} for m, s in scores.items()}
    with (export / "scores.tsv").open() as handle:
        table = list(csv.DictReader(handle, delimiter="\t"))
    assert len(table) == len(rows)
    for row, tsv in zip(rows, table):
        selected = frozen["strata"][row["stratum"]]
        assert row["families"] == sorted(selected)
        assert row["family_count"] == len(selected) == int(tsv["families"])
        assert all(row[k] == tsv[k] for k in ("stratum", "method", "status"))
        if not selected:
            assert row["metrics_percent"] is None and row["weighted_counts"] is None
            assert all(tsv[k] == "NA" for k in ("f_score", "precision", "recall"))
            continue
        tp, fp, fn = exact_counts(indexed[row["method"]], selected)
        for key, value in zip(("tp", "fp", "fn"), (tp, fp, fn)):
            assert math.isclose(row["weighted_counts"][key], float(value), rel_tol=1e-12, abs_tol=1e-10)
        for key, numerator, denominator in (("f_score", 2 * tp, 2 * tp + fp + fn),
                ("precision", tp, tp + fp), ("recall", tp, tp + fn)):
            value = float(100 * numerator / denominator) if denominator else 0.
            assert abs(value - row["metrics_percent"][key]) < 1e-10
            assert float(tsv[key]) == row["metrics_percent"][key]
    receipt = dict(status="independent_rational_arithmetic_verified", rows=len(rows),
        empty_rows=sum(not r["families"] for r in rows),
        checked_source_records=len(report["checked_records"]),
        files={str(p): digest(p) for p in (export / "report.json", export / "scores.tsv", Path(__file__))},
        limitations="Verifies exported arithmetic, not biological truth or historical scorer execution.")
    lookup = {(r["stratum"], r["method"]): r for r in rows}
    lines = ["# OrthoBench All-Method Strata", "", "Weighted F1 percentages; development-exposed descriptive results. No new confidence intervals.", "",
             "Method columns: " + "; ".join(f"M{i + 1} = {m}" for i, m in enumerate(methods)), "",
             "| Stratum | Families | " + " | ".join(f"M{i + 1}" for i in range(8)) + " |",
             "|---|---:|" + "---:|" * 8]
    for label in sorted(frozen["strata"]):
        values = [lookup[label, m]["metrics_percent"] for m in methods]
        lines.append(f"| {label} | {len(frozen['strata'][label])} | " + " | ".join(
            "NA" if v is None else f"{v['f_score']:.3f}" for v in values) + " |")
    lines.extend(["", "Precision, recall, full-precision F1 and empty-bin status are in scores.tsv.",
                  "Family descriptors are proxies, not validated domain, fragment or duplication-history labels."])
    (export / "crosscheck.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    (export / "TABLE.md").write_text("\n".join(lines) + "\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--export", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(verify(args.base.resolve(), args.export.resolve())))
