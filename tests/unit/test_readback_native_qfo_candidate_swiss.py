import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import readback_native_qfo_candidate_swiss as module
from benchmark_tools import bind_native_qfo_swiss_uncertainty as binder
from tests.unit.test_bind_native_qfo_swiss_uncertainty import native_rows, retained_bootstrap


def raw_file(tmp_path, rows, header=module.HEADER):
    path = tmp_path / "raw.txt.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(header + "\n")
        for row in rows:
            stream.write("\t".join(row) + "\n")
    return path


def test_raw_parser_preserves_all_labels_and_orientation(tmp_path):
    rows = [("A", "b", "a", "TP"), ("A", "a", "c", "FP"),
            ("B", "a", "b", "FN"), ("B", "b", "c", "TN")]
    counts, members, decisions = module.raw(raw_file(tmp_path, rows), ["A", "B"])
    assert counts["A"] == dict(TP=1, FP=1, FN=0, TN=0)
    assert counts["B"] == dict(TP=0, FP=0, FN=1, TN=1)
    assert members == dict(A={"a", "b", "c"}, B={"a", "b", "c"})
    assert decisions["A", "a", "b"] == "TP"


@pytest.mark.parametrize("rows", [[], [("A", "a", "a", "TP")], [("X", "a", "b", "TP")],
    [("A", "", "b", "TP")], [("A", "a", "b", "XX")], [("A", "a", "b")],
    [("A", "a", "b", "TP"), ("A", "b", "a", "TP")]])
def test_bad_or_incomplete_raw_refused(tmp_path, rows):
    with pytest.raises(ValueError):
        module.raw(raw_file(tmp_path, rows), ["A"])


def test_raw_header_refused(tmp_path):
    with pytest.raises(ValueError, match="header"):
        module.raw(raw_file(tmp_path, [("A", "a", "b", "TP")], header="changed"), ["A"])


def test_rational_prior_and_nonlinear_macro_statistic():
    counts = dict(A=dict(TP=2, FP=1, FN=1, TN=0), B=dict(TP=0, FP=0, FN=8, TN=2))
    families, macro = module.statistics(counts, ["A", "B"])
    assert families["A"]["PPV"] == families["A"]["TPR"] == Fraction(4, 7)
    assert families["B"]["PPV"] == Fraction(1, 2)
    assert families["B"]["TPR"] == Fraction(1, 6)
    p, r = Fraction(15, 28), Fraction(31, 84)
    assert macro == dict(PPV=p, TPR=r, F1=2*p*r/(p+r))
    assert macro["F1"] != sum(row["F1"] for row in families.values()) / 2


@pytest.mark.parametrize("change", ["negative", "boolean", "float", "missing_label", "empty", "duplicate_family"])
def test_invalid_count_or_family_types_refused(change):
    counts, names = dict(A=dict(TP=2, FP=1, FN=1, TN=0)), ["A"]
    if change in ("negative", "boolean", "float"):
        counts["A"]["TP"] = {"negative": -1, "boolean": True, "float": 2.0}[change]
    elif change == "missing_label":
        counts["A"].pop("TN")
    elif change == "empty":
        counts["A"] = dict.fromkeys(module.LABELS, 0)
    else:
        names.append("A")
    with pytest.raises(ValueError):
        module.statistics(counts, names)


def test_three_matched_cells_add_only_candidate_and_reconciliation_contrasts(retained_bootstrap):
    counts, frozen = retained_bootstrap
    _, contrasts = binder.project(native_rows(counts, [0, 1, 2]), counts, frozen)
    assert {r["name"] for r in contrasts if r["status"] == "native_records_matched"} == {
        "C_at_P0_R0", "R_at_P0_C0"}
    for row in contrasts:
        if row["status"] == "native_records_matched":
            original = next(r for r in frozen["comparisons"] if r["name"] == row["name"])
            assert row["metrics"] == original["metrics"]
        else:
            assert row["metrics"] is row["family_differences"] is None


def test_unmatched_candidate_does_not_disable_already_matched_reconciliation(retained_bootstrap):
    counts, frozen = retained_bootstrap
    rows = native_rows(counts, [0, 1, 2])
    rows[module.CELLS[1]]["families"][0]["represented_genes"][0] = "different_gene"
    rows[module.CELLS[1]]["retained_family_records_identical"] = False
    _, contrasts = binder.project(rows, counts, frozen)
    assert [r["name"] for r in contrasts if r["status"] == "native_records_matched"] == ["R_at_P0_C0"]


@pytest.mark.parametrize("change", ["counts", "members", "family_statistic", "macro_statistic", "family_order"])
def test_cell_readback_refuses_changes(tmp_path, change):
    counts, members, _ = module.raw(raw_file(tmp_path, [("A", "a", "b", "TP")]), ["A"])
    points, macro = module.statistics(counts, ["A"])
    cell = dict(families=[dict(family="A", counts_without_prior=dict(counts["A"]), represented_genes=["a", "b"],
        statistics_with_prior={k: float(v) for k, v in points["A"].items()})],
        aggregate={k: float(v) for k, v in macro.items()})
    module.check_cell(cell, counts, members, ["A"])
    cell = copy.deepcopy(cell)
    if change == "counts":
        cell["families"][0]["counts_without_prior"]["TP"] += 1
    elif change == "members":
        cell["families"][0]["represented_genes"].append("c")
    elif change == "family_statistic":
        cell["families"][0]["statistics_with_prior"]["F1"] = .123
    elif change == "macro_statistic":
        cell["aggregate"]["F1"] = .123
    else:
        cell["families"][0]["family"] = "B"
    with pytest.raises(ValueError):
        module.check_cell(cell, counts, members, ["A"])


def test_no_overwrite_before_readback(tmp_path, monkeypatch):
    output = tmp_path / "existing.json"
    output.write_text("retain me\n")
    monkeypatch.setattr(sys, "argv", ["readback", "--audit", "absent", "--audit-sha256", "unused",
        "--binding", "absent", "--binding-sha256", "unused", "--output", str(output)])
    with pytest.raises(ValueError, match="Output already exists"):
        module.main()
    assert output.read_text() == "retain me\n"


def test_identical_count_copy_may_have_different_path(tmp_path):
    parsed = dict(counts=[1, 2, 3])
    first, second = tmp_path / "original.json", tmp_path / "copy.json"
    for path in (first, second):
        path.write_text(json.dumps(parsed))
    assert module.record(first) != module.record(second)
    module.check_count_copy(module.record(first), module.record(second), parsed)


@pytest.mark.parametrize("change", ["different_bytes", "changed_original", "changed_supplied", "parsed_content"])
def test_count_copy_requires_both_records_and_exact_content(tmp_path, change):
    first, second = tmp_path / "original.json", tmp_path / "copy.json"
    for path in (first, second):
        path.write_text('{"count": 1}\n')
    if change == "different_bytes":
        second.write_text('{"count": 2}\n')
    a, b = module.record(first), module.record(second)
    if change == "changed_original":
        first.write_text('{"count": 2}\n')
    elif change == "changed_supplied":
        second.write_text('{"count": 2}\n')
    parsed = dict(count=2 if change == "parsed_content" else 1)
    with pytest.raises(ValueError, match="count copy"):
        module.check_count_copy(a, b, parsed)


def test_manuscript_candidate_endpoint_table_matches_admitted_snapshot():
    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    snapshot = json.loads((root / "native_qfo_scientific_scores_20261006_v2/report.json").read_text())
    baseline, candidate = (next(r for r in snapshot["rows"] if r["cell"] == c) for c in module.CELLS)
    section = (root / "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md").read_text().split(
        "### Third Native QfO Cell And Candidate-Expansion Trade-Off\n", 1)[1].split(
        "### Native Functional-Pair Composition", 1)[0]
    for endpoint in ("VGNC", "SwissTrees", "TreeFam-A", "GO", "EC", "FAS"):
        row = next(line for line in section.splitlines() if line.startswith("| " + endpoint + " |"))
        fields = [f.strip() for f in row.split("|")[1:-1]]
        places = 8 if endpoint in ("GO", "EC") else 10
        assert fields[2] == f'{baseline["scores"][endpoint]:.{places}f}'
        assert fields[3] == f'{candidate["scores"][endpoint]:.{places}f}'
        assert fields[4] == f'{100*(candidate["scores"][endpoint]-baseline["scores"][endpoint]):+.6f}'
    assert "three admitted cells and four unavailable score rows" in section
    assert "No paired\ninterval is supplied for the other QfO endpoints or the secondary mean" in section


def test_manuscript_and_result_intervals_match_guarded_binding():
    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    binding = json.loads((root / "native_qfo_candidate_swiss_uncertainty_20261006_v1.json").read_text())
    effect = next(r for r in binding["contrasts"] if r["name"] == "C_at_P0_R0")
    for filename in ("PUBLICATION_MANUSCRIPT_DRAFT_20260916.md", "NATIVE_QFO_CANDIDATE_UNCERTAINTY_RESULT_20261006.md"):
        text = (root / filename).read_text()
        if filename.startswith("PUBLICATION_"):
            text = text.split("### Third Native QfO Cell And Candidate-Expansion Trade-Off\n", 1)[1]
        for metric, name in (("F1", "F1"), ("PPV", "Precision"), ("TPR", "Recall")):
            value = effect["metrics"][metric]
            label = name if filename.startswith("PUBLICATION_") else "SwissTrees " + name.lower()
            if metric == "F1" and not filename.startswith("PUBLICATION_"):
                label = "SwissTrees F1"
            row = next(line for line in text.splitlines() if line.startswith("| " + label + " |"))
            fields = [f.strip() for f in row.split("|")[1:-1]]
            assert fields[1] == f'{100*value["difference"]:+.6f}'
            for field, key in zip(fields[2:4], ("paired_percentile_ci", "bonferroni_percentile_ci")):
                low, high = value[key]
                assert field == f"[{100*low:.6f}, {100*high:.6f}]"
            assert fields[4] == "/".join(str(value[k]) for k in ("family_wins", "family_ties", "family_losses"))


def test_current_claims_preserve_incomplete_and_conditional_scope():
    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    text = (root / "PUBLICATION_CLAIMS_20260916.md").read_text().split("## Current Requirement Status", 1)[0]
    for phrase in ("three admitted fresh cells and four", "Native candidate expansion improves SwissTrees F1 conclusively",
                   "Not established: all 18 complete records", "not equivalence", "including retained failure 9"):
        assert phrase in text
