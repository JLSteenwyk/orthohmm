import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import integrate_publication_supplements as current


ROOT = Path(__file__).resolve().parents[2]


def retained():
    directory = ROOT / "benchmark_tools/results"
    return {key: (directory / name).read_text() if key == "parent" else json.loads((directory / name).read_text())
            for key, (name, _) in current.PINS.items()}


def test_exact_parent_restoration_and_unchanged_inputs():
    docs = retained()
    before = copy.deepcopy(docs)
    revised, sections = current.manuscript(docs["parent"], docs)
    for section in sections:
        revised = revised.replace(section, "", 1)
    assert revised == docs["parent"] and docs == before
    assert "zero coverage" in sections[0] and "No native VGNC" in sections[0]
    assert "0.999559608605" in sections[0] and "5.999645225074e-13" in sections[0]
    assert "not a repeat of" in sections[1] and "component remains local" in sections[1]
    assert "prior manuscript, not this extended text" in sections[1]


@pytest.mark.parametrize("change", ["status", "admission", "cells", "duplicate", "coverage", "tail", "controls", "replay_status", "kernel", "native_rows"])
def test_altered_evidence_refused(change):
    docs = retained()
    p, r = docs["poisson"], docs["relocation"]
    if change == "status": p["status"] = "pending"
    elif change == "admission": p["native_intervals_admitted"] = True
    elif change == "cells": p["rows"].pop()
    elif change == "duplicate": p["rows"][1] = copy.deepcopy(p["rows"][0])
    elif change == "coverage": p["rows"][0]["covered_enumerated_mass"] = float("nan")
    elif change == "tail": p["rows"][0]["omitted_tail_mass"] = 1e-3
    elif change == "controls": p["model_violation_controls"][0]["coverage"] = .99
    elif change == "replay_status": r["status"] = "pending"
    elif change == "kernel": r["original_reader_kernels_changed"] = True
    elif change == "native_rows": r["native_readback"]["stage_rows"] = 105
    with pytest.raises(ValueError):
        current.sections(docs)


def test_ambiguous_anchor_and_already_integrated_refused():
    docs = retained()
    with pytest.raises(ValueError):
        current.manuscript(docs["parent"] + current.ANCHORS[0], docs)
    revised, _ = current.manuscript(docs["parent"], docs)
    with pytest.raises(ValueError):
        current.manuscript(revised, docs)


def test_occupied_output_refused_before_input_reads(tmp_path):
    for destination in (tmp_path, tmp_path / "receipt.json"):
        with pytest.raises(FileExistsError):
            current.run(tmp_path, destination, tmp_path)


def test_actual_generation_with_copied_reporting_inputs(tmp_path, monkeypatch):
    directory = tmp_path / "benchmark_tools/results"
    directory.mkdir(parents=True)
    docs, pins = retained(), {}
    for key in current.PINS:
        path = directory / (key + (".md" if key == "parent" else ".json"))
        path.write_text(docs[key] if key == "parent" else json.dumps(docs[key]))
        pins[key] = path.name, current.record(path)["sha256"]
    monkeypatch.setattr(current, "PINS", pins)
    output, receipt = directory / "revised.md", directory / "generation.json"
    result = current.run(tmp_path, output, receipt)
    assert json.loads(receipt.read_text()) == result
    restored = output.read_text()
    for section in result["inserted_sections"]:
        restored = restored.replace(section, "", 1)
    assert restored == docs["parent"]
    assert result["parent_unchanged_except_insertions"]
    assert result["new_bootstrap_draws"] == 0
    assert all(result[k] is False for k in ("publication_ready", "native_intervals_admitted", "manuscript_rendered", "visual_reviewed"))
    with pytest.raises(FileExistsError):
        current.run(tmp_path, output, receipt)


def test_changed_pin_refused_without_output(tmp_path, monkeypatch):
    directory = tmp_path / "benchmark_tools/results"
    directory.mkdir(parents=True)
    parent = directory / "parent.md"
    parent.write_text("changed")
    monkeypatch.setattr(current, "PINS", {"parent": (parent.name, "0" * 64)})
    output, receipt = directory / "revised.md", directory / "generation.json"
    with pytest.raises(ValueError, match="Changed input parent"):
        current.run(tmp_path, output, receipt)
    assert not output.exists() and not receipt.exists()
