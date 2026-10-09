"""Check integration against already executed, independently validated records."""

import copy
import csv
import io
import json
from pathlib import Path

import pytest

from benchmark_tools import integrate_native_qfo_four_cell_strata as current


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"


def retained():
    docs = {}
    for key, (name, _) in current.PINS.items():
        path = RESULTS / name
        docs[key] = path.read_text() if key == "parent" else json.loads(path.read_text())
    return docs


def test_all_rows_values_and_exact_parent_restoration():
    docs = retained()
    before = copy.deepcopy(docs)
    revised, methods, results, rows, changes = current.manuscript(docs["parent"], docs["result"], docs["reader"])
    assert docs == before
    assert revised.replace(methods, "", 1).replace(results, "", 1) == docs["parent"]
    assert len(rows) == 23 and sum(r["families"] == 0 for r in rows) == 5
    assert changes == [dict(family="CASP", TP=0, FP=-3, FN=0, TN=3), dict(family="GH14", TP=-1, FP=0, FN=1, TN=0)]
    table = [line for line in results.splitlines() if line.startswith("| ")
             and line.split("|")[1].strip() in current.SUITES]
    assert len(table) == 23
    for row, line in zip(rows, table):
        cells = line.split("|")[1:-1]
        assert [c.strip() for c in cells[:3]] == [row["suite"], row["stratum"], str(row["families"])]
        assert [c.strip() for c in cells[3:]] == ["NA" if row[m] is None else f"{100*row[m]:+.3f}" for m in current.METRICS]
    assert "-0.320215" in results and "+0.160465" in results and "-0.794278" in results
    assert "No subgroup intervals were computed" in results


@pytest.mark.parametrize("change", ["schema", "reader_schema", "rows", "order", "reference", "empty", "nan", "scope",
                                   "reader_count", "draws", "localized_count", "localized_bin", "complement"])
def test_changed_evidence_not_integrated(change):
    docs = retained()
    result, reader = docs["result"], docs["reader"]
    row = next(r for r in result["differences"] if r["contrast"] == "P_at_C0_R1")
    if change == "schema": result["schema"] = "old"
    elif change == "reader_schema": reader["schema"] = "old"
    elif change == "rows": result["differences"].pop()
    elif change == "order": result["differences"] = list(reversed(result["differences"]))
    elif change == "reference": row["reference"] = "p0_c0_r0"
    elif change == "empty": next(r for r in result["differences"] if r["contrast"] == "P_at_C0_R1" and not r["families"])["F1"] = 0
    elif change == "nan": row["F1"] = float("nan")
    elif change == "scope": reader["publication_ready"] = True
    elif change == "reader_count": reader["score_rows_checked"] -= 1
    elif change == "draws": result["new_bootstrap_draws"] = True
    elif change == "localized_count": next(r for r in result["family_rows"] if r["cell"] == "p1_c0_r1" and r["family"] == "CASP")["counts_without_prior"]["TP"] += 1
    elif change == "localized_bin": result["bins"]["sequence"]["lower_entropy"].remove("CASP")
    elif change == "complement": next(r for r in result["differences"] if r["contrast"] == "P_at_C0_R1" and r["stratum"] == "short_relative")["F1"] = .01
    with pytest.raises(ValueError): current.sections(result, reader)


@pytest.mark.parametrize("parent", ["missing", current.METHOD_ANCHOR*2 + current.RESULT_ANCHOR])
def test_ambiguous_anchors_rejected(parent):
    docs = retained()
    with pytest.raises(ValueError): current.manuscript(parent, docs["result"], docs["reader"])


def temporary_inputs(tmp_path, monkeypatch):
    docs = retained()
    directory = tmp_path / "benchmark_tools/results"
    directory.mkdir(parents=True)
    pins = {}
    for key in ("result", "reader", "parent"):
        path = directory / (key + (".md" if key == "parent" else ".json"))
        if key == "reader": docs[key]["report"] = current.record(directory / "result.json")
        path.write_text(docs[key] if key == "parent" else json.dumps(docs[key], sort_keys=True, allow_nan=False))
        pins[key] = path.name, current.record(path)["sha256"]
    monkeypatch.setattr(current, "PINS", pins)
    return directory


def test_temporary_generation_and_serialized_values(tmp_path, monkeypatch):
    directory = temporary_inputs(tmp_path, monkeypatch)
    output, manuscript = directory / "test_integration", directory / "test_manuscript.md"
    manifest = current.run(tmp_path, output, manuscript)
    methods = (output / "methods_section.md").read_text()
    results = (output / "results_section.md").read_text()
    assert manuscript.read_text().replace(methods, "", 1).replace(results, "", 1) == (directory / "parent.md").read_text()
    parser = csv.DictReader(io.StringIO((output / "profile_bins.tsv").read_text()), delimiter="\t")
    rows = list(parser)
    assert parser.fieldnames == list(current.FIELDS) and len(rows) == 23
    assert sum(r["F1"] == "NA" for r in rows) == 5
    assert manifest["profile_value_cells"] == 69 and manifest["publication_ready"] is False
    assert manifest["parent_body_unchanged_except_insertions"] is True
    assert manifest["new_bootstrap_draws"] == 0 and manifest["manuscript_rendered"] is False


def test_mixed_reader_binding_and_changed_pin_rejected(tmp_path, monkeypatch):
    directory = temporary_inputs(tmp_path, monkeypatch)
    reader = json.loads((directory / "reader.json").read_text())
    reader["report"]["sha256"] = "0"*64
    (directory / "reader.json").write_text(json.dumps(reader))
    with pytest.raises(ValueError): current.run(tmp_path, directory / "first", directory / "first.md")
    current.PINS["reader"] = "reader.json", current.record(directory / "reader.json")["sha256"]
    with pytest.raises(ValueError, match="Mixed"): current.run(tmp_path, directory / "second", directory / "second.md")
    assert not (directory / "first").exists() and not (directory / "second").exists()


def test_existing_output_guards_precede_loading(tmp_path):
    with pytest.raises(FileExistsError): current.run(tmp_path, tmp_path, tmp_path / "unused.md")
