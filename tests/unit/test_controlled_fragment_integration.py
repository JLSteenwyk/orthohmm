"""Test bounded manuscript integration of the independently verified fragment panel."""

import copy
import csv
import io
import json
from pathlib import Path

import pytest

from benchmark_tools import integrate_controlled_fragment_results as current


ROOT = Path(__file__).resolve().parents[2]


def retained():
    return {key: (ROOT / name).read_text() if key == "parent" else json.loads((ROOT / name).read_text())
            for key, (name, _) in current.PINS.items()}


def test_complete_tables_restore_exact_parent_and_leave_inputs_unchanged():
    docs = retained()
    before = copy.deepcopy(docs)
    revised, inserts, tables = current.manuscript(docs["parent"], docs)
    restored = revised
    for section in inserts:
        restored = restored.replace(section, "", 1)
    assert restored == docs["parent"] and docs == before
    assert [len(t) for t in tables] == [4, 15, 12]
    lines = [line for line in inserts[1].splitlines() if line.startswith("| ")]
    assert len(lines) == 37
    means, differences, counts = tables
    for row in means:
        expected = "| " + " | ".join([current.LABELS[row["method"]],
            *[f"{100*row[k]:.4f}" for k in current.MEAN_FIELDS[1:-1]]]) + " |"
        assert expected in lines
    for row in differences:
        for kind in ("nominal", "adjusted"):
            assert f"[{100*row[kind+'_low']:+.4f}, {100*row[kind+'_high']:+.4f}]" in inserts[1]
    for row in counts:
        for arm in current.ARMS:
            assert " / ".join(str(row[arm+"_"+k]) for k in ("tp", "fp", "fn")) in inserts[1]
    assert "not pooled-pair F1" in inserts[0]
    assert "not four simultaneous workers" in inserts[0]
    assert "neither equivalence nor absence" in inserts[1]
    assert "not estimates of isolated performance" in inserts[2]
    assert "Representative search/candidate/group/reconciliation tracing" in inserts[1]


@pytest.mark.parametrize("change", ["report_scope", "reader_exit", "reader_count", "reader_binding",
    "panel_seed", "runtime_equivalence", "scheduler", "record_duplicate", "record_missing", "record_status",
    "nan", "counts", "stratum_order", "means", "mean_eligibility", "comparison_order", "comparison_scope",
    "paired_values", "intervals", "f1_interpretation"])
def test_changed_or_incomplete_evidence_rejected(change):
    docs = retained()
    r, e, p, v = (docs[k] for k in ("result", "reader", "panel", "runtime"))
    if change == "report_scope": r["publication_ready"] = True
    elif change == "reader_exit": e["independent_result_readback"]["terminal"]["exit_code"] = 1
    elif change in ("reader_count", "reader_binding"):
        observed = json.loads(e["independent_result_readback"]["terminal"]["output"])
        observed["score_records" if change == "reader_count" else "report_sha256"] = 79 if change == "reader_count" else "0"*64
        e["independent_result_readback"]["terminal"]["output"] = json.dumps(observed)
    elif change == "panel_seed": p["datasets"][0]["seed"] += 1
    elif change == "runtime_equivalence": v["inventory_amendment"]["historical_output_equivalence_established"] = True
    elif change == "scheduler": r["scheduler"]["ExitCode"] = "1:0"
    elif change == "record_duplicate": r["records"].append(r["records"][0])
    elif change == "record_missing": r["records"].pop()
    elif change == "record_status": r["records"][0]["status"] = "unavailable"
    elif change == "nan": r["records"][0]["score"]["f1"] = float("nan")
    elif change == "counts": r["records"][0]["strata"][0]["score"]["tp"] += 1
    elif change == "stratum_order": r["records"][0]["strata"].reverse()
    elif change == "means": r["metric_means"][0]["metrics"]["f1"]["mean"] += .01
    elif change == "mean_eligibility": r["metric_means"][0]["metrics"]["f1"]["eligible_seeds"].pop()
    elif change == "comparison_order": r["comparisons"].reverse()
    elif change == "comparison_scope": r["comparisons"][0]["replicates"] = 19999
    elif change == "paired_values": r["comparisons"][0]["paired_differences"][0] += .01
    elif change == "intervals": r["comparisons"][0]["adjusted_interval"][0] = .1
    elif change == "f1_interpretation": r["comparisons"][0]["adjusted_interval"][1] = -.00001
    with pytest.raises(ValueError): current.sections(docs)


@pytest.mark.parametrize("parent", ["missing", "".join(current.ANCHORS)*2])
def test_ambiguous_parent_rejected(parent):
    with pytest.raises(ValueError): current.manuscript(parent, retained())


def temporary_inputs(tmp_path, monkeypatch):
    docs = retained()
    directory = tmp_path / "benchmark_tools/results"
    directory.mkdir(parents=True)
    pins = {}
    for key in ("panel", "runtime", "result", "reader", "parent"):
        path = directory / (key + (".md" if key == "parent" else ".json"))
        if key == "result":
            for name in ("panel", "runtime"):
                ref = current.record(directory / (name + ".json"))
                docs[key]["manifest" if name == "panel" else name] = dict(
                    absolute_path=ref["path"], bytes=ref["bytes"], sha256=ref["sha256"])
        if key == "reader":
            for field, name in (("outputs", "result"), ("prepared_manifest", "panel"), ("runtime_manifest", "runtime")):
                docs[key][field]["sha256"] = current.record(directory / (name + ".json"))["sha256"]
            terminal = docs[key]["independent_result_readback"]["terminal"]
            output = json.loads(terminal["output"])
            output["report_sha256"] = docs[key]["outputs"]["sha256"]
            terminal["output"] = json.dumps(output)
        path.write_text(docs[key] if key == "parent" else json.dumps(docs[key], sort_keys=True, allow_nan=False))
        pins[key] = str(path.relative_to(tmp_path)), current.record(path)["sha256"]
    monkeypatch.setattr(current, "PINS", pins)
    return directory


def test_temporary_generation_all_serialized_tables_and_parent(tmp_path, monkeypatch):
    directory = temporary_inputs(tmp_path, monkeypatch)
    output, revised = directory / "integration", directory / "manuscript.md"
    manifest = current.run(tmp_path, output, revised)
    restored = revised.read_text()
    for section in ("methods", "results", "limitations"):
        restored = restored.replace((output / (section + "_section.md")).read_text(), "", 1)
    assert restored == (directory / "parent.md").read_text()
    docs = {k: json.loads((directory / (k + ".json")).read_text()) for k in ("result", "reader", "panel", "runtime")}
    expected = current.validated_tables(docs)
    for name, fields, table in zip(("means", "comparisons", "endpoint_counts"),
                                 (current.MEAN_FIELDS, current.DIFF_FIELDS, current.COUNT_FIELDS), expected):
        parser = csv.DictReader(io.StringIO((output / (name + ".tsv")).read_text()), delimiter="\t")
        rows = list(parser)
        assert parser.fieldnames == list(fields)
        assert rows == [{k: str(r[k]) for k in fields} for r in table]
    assert manifest["table_rows"] == dict(means=4, comparisons=15, endpoint_counts=12)
    assert manifest["new_bootstrap_draws"] == 0 and manifest["new_inference_or_scoring"] is False
    assert manifest["publication_ready"] is False and manifest["visual_reviewed"] is False
    assert manifest["parent_body_unchanged_except_insertions"] is True


def test_mixed_binding_refused_even_when_pin_updated(tmp_path, monkeypatch):
    directory = temporary_inputs(tmp_path, monkeypatch)
    reader_path = directory / "reader.json"
    reader = json.loads(reader_path.read_text())
    reader["runtime_manifest"]["sha256"] = "0"*64
    reader_path.write_text(json.dumps(reader))
    current.PINS["reader"] = str(reader_path.relative_to(tmp_path)), current.record(reader_path)["sha256"]
    with pytest.raises(ValueError, match="Mixed"):
        current.run(tmp_path, directory / "output", directory / "revised.md")
    assert not (directory / "output").exists()


def test_existing_destination_refused_before_input_loading(tmp_path):
    with pytest.raises(FileExistsError): current.run(tmp_path, tmp_path, tmp_path / "new.md")
