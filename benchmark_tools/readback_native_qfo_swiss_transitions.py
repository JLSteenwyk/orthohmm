"""Independent stdlib/SQLite readback of native SwissTrees transition evidence."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3

LABELS = ("TP", "FP", "FN", "TN")
HEADER = "# Dataset<tab>Protein ID 1<tab>Protein ID 2<tab>Correctness (TP:True positive, FP: False positive, TN: True Negative. FN: False Negative)"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).absolute()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def verify(report_path, report_sha):
    report_ref = record(report_path)
    require(report_ref["sha256"] == report_sha, "Changed diagnostic report")
    report = json.loads(Path(report_path).read_text())
    require(report["schema"] == "native_qfo_swiss_pair_transitions_v1"
            and report["source"] == record(Path(__file__).with_name("trace_native_qfo_swiss_transitions.py"))
            and all(report[k] is False for k in ("new_scoring_or_admission", "uncertainty_admitted",
                "scientific_timings_admitted", "independent_confirmation", "publication_ready"))
            and type(report["new_bootstrap_draws"]) is int and report["new_bootstrap_draws"] == 0,
            "Changed diagnostic source/scope")
    checked = [*report["checked_inputs"], report["changed_relations_ledger"], report_ref]
    for ref in checked:
        require(record(ref["path"]) == ref, "Changed direct diagnostic evidence")
    require([row["cell"] for row in report["cells"]] == ["p0_c0_r0", "p0_c0_r1"]
            and report["cells"][1]["timing_eligible"] is False
            and report["cells"][1]["timing_admitted"] is False
            and report["cells"][1]["resources"] is None, "Wrong cells or repaired timing")
    families = [row["family"] for row in report["comparison"]["families"]]
    require(len(families) == 18 and len(set(families)) == 18, "Changed family inventory")
    with sqlite3.connect(":memory:") as db:
        db.execute("CREATE TABLE decisions(cell INTEGER, family TEXT, a TEXT, b TEXT, label TEXT, "
                   "PRIMARY KEY(cell,family,a,b))")
        for cell, row in enumerate(report["cells"]):
            require(row["raw"] in report["checked_inputs"], "Unchecked raw evidence")
            with gzip.open(row["raw"]["path"], "rt", newline="") as stream:
                require(stream.readline().rstrip("\r\n") == HEADER, "Changed raw header")
                for item in csv.reader(stream, delimiter="\t"):
                    require(len(item) == 4, "Invalid raw width")
                    f, a, b, label = item
                    require(f in families and a and b and a != b and label in LABELS, "Invalid native relation")
                    db.execute("INSERT INTO decisions VALUES(?,?,?,?,?)", (cell, f, min(a, b), max(a, b), label))
        db.execute("CREATE VIEW paired AS SELECT l.family,l.a,l.b,l.label AS before,r.label AS after "
                   "FROM decisions l JOIN decisions r ON l.family=r.family AND l.a=r.a AND l.b=r.b "
                   "WHERE l.cell=0 AND r.cell=1")
        n = db.execute("SELECT COUNT(*) FROM paired").fetchone()[0]
        require(n == report["comparison"]["reference_relations"] and
                db.execute("SELECT COUNT(*) FROM decisions").fetchone()[0] == 2 * n,
                "Reference relation universes differ")
        require(db.execute("SELECT COUNT(*) FROM paired WHERE "
            "(before IN ('TP','FN')) != (after IN ('TP','FN'))").fetchone()[0] == 0, "Changed reference truth")
        transitions = {(f, a, b): count for f, a, b, count in db.execute(
            "SELECT family,before,after,COUNT(*) FROM paired GROUP BY family,before,after")}
        marginals = {(cell, f, label): count for cell, f, label, count in db.execute(
            "SELECT cell,family,label,COUNT(*) FROM decisions GROUP BY cell,family,label")}
        totals = {a + "->" + b: 0 for a in LABELS for b in LABELS}
        for row in report["comparison"]["families"]:
            f = row["family"]
            table = {a + "->" + b: transitions.get((f, a, b), 0) for a in LABELS for b in LABELS}
            require(table == row["transitions"] and sum(table.values()) == row["relations"],
                    "Family transition table differs")
            for key, transition in (("removed_true_positives", "TP->FN"), ("removed_false_positives", "FP->TN"),
                                   ("added_true_positives", "FN->TP"), ("added_false_positives", "TN->FP")):
                require(row[key] == table[transition], "Family transition summary differs")
            for cell, key in enumerate(("before_counts", "after_counts")):
                require(row[key] == {label: marginals.get((cell, f, label), 0) for label in LABELS},
                        "Family marginal differs")
            for key, value in table.items():
                totals[key] += value
        require(totals == report["comparison"]["transitions"], "Pooled transitions differ")
        for key, transition in (("removed_true_positives", "TP->FN"), ("removed_false_positives", "FP->TN"),
                               ("added_true_positives", "FN->TP"), ("added_false_positives", "TN->FP")):
            require(report["comparison"][key] == totals[transition], "Pooled transition summary differs")
        require(report["comparison"]["after_predictions_subset_on_reference"] is
                (totals["FN->TP"] + totals["TN->FP"] == 0), "Changed subset statement")
        expected_changes = list(db.execute("SELECT family,a,b,before,after FROM paired "
            "WHERE before!=after ORDER BY family,a,b"))
        with open(report["changed_relations_ledger"]["path"], newline="") as stream:
            reader = csv.reader(stream, delimiter="\t")
            require(next(reader) == ["family", "protein_a", "protein_b", "before", "after"], "Changed ledger header")
            actual_changes = [tuple(row) for row in reader]
        require(actual_changes == expected_changes and len(expected_changes) == report["comparison"]["changed_relations"],
                "Changed-relation ledger differs from SQL join")
        points = []
        for cell, row in enumerate(report["cells"]):
            precision, recall = [], []
            for f in families:
                tp, fp, fn = (marginals.get((cell, f, label), 0) / 2 + 1 for label in ("TP", "FP", "FN"))
                precision.append(tp / (tp + fp))
                recall.append(tp / (tp + fn))
            p, r = sum(precision) / 18, sum(recall) / 18
            point = dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)
            require(all(abs(point[k] - row["macro_statistics"][k]) <= 1e-12 for k in point),
                    "Independently reconstructed macro endpoint differs")
            require(abs(point["F1"] - row["native_endpoint_f1"]) <= 5e-8, "Native decimal F1 differs")
            points.append(point)
    for ref in checked:
        require(record(ref["path"]) == ref, "Direct evidence changed during readback")
    return dict(schema="native_qfo_swiss_transition_sql_readback_v1", report=report_ref,
        source=record(__file__), checked_inputs=checked, families_checked=18,
        raw_rows_checked=2 * n, paired_relations_checked=n, changed_relations_checked=len(expected_changes),
        transition_cells_checked=18 * 16, macro_points=points, new_scoring_or_admission=False,
        uncertainty_admitted=False, scientific_timings_admitted=False, publication_ready=False,
        limitations=["Independent SQL arithmetic/parser readback, not new scientific admission or transitive execution proof.",
                     "Same retained source data; no new biological evidence, uncertainty or inference timing."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = verify(args.report, args.report_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("families_checked", "raw_rows_checked", "paired_relations_checked",
                                            "changed_relations_checked", "transition_cells_checked")}))
