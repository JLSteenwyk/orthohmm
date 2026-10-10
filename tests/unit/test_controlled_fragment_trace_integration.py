import copy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import integrate_controlled_fragment_trace as current


ROOT = Path(__file__).resolve().parents[2]


def retained():
    return {k: (ROOT / name).read_text() if k == "parent" else json.loads((ROOT / name).read_text())
            for k, (name, _) in current.PINS.items()}


def test_exact_parent_restoration_with_explicit_status_update():
    docs = retained()
    before = copy.deepcopy(docs)
    revised, inserts, rows = current.manuscript(docs["parent"], docs)
    restored = revised.replace(current.NEW_STATUS, current.OLD_STATUS, 1)
    for section in inserts:
        restored = restored.replace(section, "", 1)
    assert restored == docs["parent"] and docs == before
    assert len(rows) == 11
    assert sum(r["baseline_representatives"] for r in rows) == sum(r["fragment_representatives"] for r in rows) == 53
    assert "not rates, prevalence estimates" in inserts[1]
    assert "tree-level cause is unknown" in inserts[1]
    assert "publication readiness" in inserts[2]
    assert current.OLD_STATUS not in revised


def test_all_location_cells_and_method_totals_match_native_report():
    docs = retained()
    inserts, rows = current.sections(docs)
    for row in rows:
        for arm in ("baseline", "fragment"):
            assert row[arm + "_representatives"] == sum(
                c["stages"][arm]["observed_location"] == row["observed_location"]
                for c in docs["stages"]["cases"] if c["method"] == row["method"])
        expected = f"| {current.LABELS[row['method']]} | {row['observed_location'].replace('_', ' ')} | {row['baseline_representatives']} | {row['fragment_representatives']} |"
        assert expected in inserts[1]
    assert [sum(r["baseline_representatives"] for r in rows if r["method"] == m) for m in current.METHODS] == [12, 18, 11, 12]


@pytest.mark.parametrize("change", ["bins", "selection", "report", "readback_status", "summary", "terminal", "receipt", "identity", "scope", "example"])
def test_changed_or_unverified_inputs_are_rejected(change):
    docs = retained()
    if change == "bins": docs["selection"]["planned_bins"] = 83
    elif change == "selection": docs["stages"]["selection"]["sha256"] = "0" * 64
    elif change == "report": docs["readback"]["report"]["sha256"] = "0" * 64
    elif change == "readback_status": docs["readback"]["status"] = "unverified"
    elif change == "summary": docs["readback"]["summary"].pop()
    elif change == "terminal": docs["execution"]["independent_stage_readback_recovery"]["terminal"]["exit_code"] = 1
    elif change == "receipt": docs["execution"]["independent_stage_readback_recovery"]["readback_sha256"] = "0" * 64
    elif change == "identity": docs["stages"]["cases"][0]["truth"] = False
    elif change == "scope": docs["stages"]["publication_ready"] = True
    elif change == "example": docs["stages"]["cases"][0]["stages"]["fragment"]["graph_connected"] = False
    with pytest.raises(ValueError):
        current.tables(docs)


def test_existing_output_refused_before_input_reads(tmp_path):
    with pytest.raises(FileExistsError):
        current.run(tmp_path, tmp_path, tmp_path / "manuscript")


def test_temporary_generation_preserves_inputs_and_restores_parent(tmp_path, monkeypatch):
    directory = tmp_path / "benchmark_tools/results"
    directory.mkdir(parents=True)
    docs, pins = retained(), {}
    # Logical cross-bindings stay frozen; temporary documentary paths get new pins.
    for key, (name, sha) in current.PINS.items():
        path = directory / (key + (".md" if key == "parent" else ".json"))
        path.write_text(docs[key] if key == "parent" else json.dumps(docs[key], sort_keys=True))
        pins[key] = str(path.relative_to(tmp_path)), current.record(path)["sha256"]
    original_tables = current.tables
    original_pins = current.PINS
    def tables(fixture):
        monkeypatch.setattr(current, "PINS", original_pins)
        try:
            return original_tables(fixture)
        finally:
            monkeypatch.setattr(current, "PINS", pins)
    monkeypatch.setattr(current, "tables", tables)
    monkeypatch.setattr(current, "PINS", pins)
    output, revised = directory / "integration", directory / "revised.md"
    result = current.run(tmp_path, output, revised)
    restored = revised.read_text().replace(current.NEW_STATUS, current.OLD_STATUS, 1)
    for name in ("methods", "results", "limitations"):
        restored = restored.replace((output / (name + "_section.md")).read_text(), "", 1)
    assert restored == docs["parent"]
    with (output / "representative_locations.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == result["location_rows"] == 11
    assert result["publication_ready"] is result["manuscript_rendered"] is result["visual_reviewed"] is False
    assert result["parent_body_unchanged_except_insertions_and_recorded_status"]
