import csv
import json
from fractions import Fraction

import pytest

from benchmark_tools.check_ob_complete_strata import exact_counts, verify


def test_rational_weights():
    records = {"a": dict(genes=4, true_positive=2, false_positive=1, false_negative=4)}
    assert exact_counts(records, ["a"]) == [Fraction(2, 3), Fraction(1, 3), Fraction(4, 3)]
    assert exact_counts(records, []) == [0, 0, 0]


@pytest.mark.parametrize("tamper", [False, True])
def test_full_check_and_metric_tampering(tmp_path, tamper):
    methods = [f"m{i}" for i in range(8)]
    score = dict(refog_records=[dict(refog="a", genes=2, true_positive=1, false_positive=0, false_negative=0)])
    strata = {f"s{i}": ["a"] if i else [] for i in range(14)}
    sources = {
        "ob_error_strata_prepared_20260916.json": dict(strata=strata),
        "orthobench_paired_uncertainty_20260916.json": dict(scores={m: score for m in methods[:3]}),
        "retained_ob_comparator_readback_20260926.json": dict(rows=[dict(key=m, score=score) for m in methods[3:]])}
    rows = [dict(stratum=s, method=m, families=f, family_count=len(f),
                 status="descriptive" if f else "empty_nonestimable",
                 weighted_counts=dict(tp=1, fp=0, fn=0) if f else None,
                 metrics_percent=dict(f_score=100, precision=100, recall=100) if f else None)
            for s, f in strata.items() for m in methods]
    sources["report.json"] = dict(rows=rows, checked_records=[])
    if tamper:
        rows[-1]["metrics_percent"]["recall"] = 99
    for name, content in sources.items():
        (tmp_path / name).write_text(json.dumps(content))
    with (tmp_path / "scores.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["stratum", "method", "families", "f_score", "precision", "recall", "status"])
        for r in rows:
            writer.writerow([r["stratum"], r["method"], r["family_count"],
                             *[r["metrics_percent"][k] if r["families"] else "NA"
                               for k in ("f_score", "precision", "recall")], r["status"]])
    if tamper:
        with pytest.raises(AssertionError):
            verify(tmp_path, tmp_path)
    else:
        assert verify(tmp_path, tmp_path)["rows"] == 112
