"""Regression for configured ownership FASTAs absent from selection pins."""

import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools import readback_controlled_fragment_stage_table as reader
from benchmark_tools import relocate_controlled_fragment_stages_v2 as corrected
from tests.unit.test_relocate_controlled_fragment_stages import dump, panel, put


@pytest.fixture
def recovery_panel(panel):
    root, report_path, readback_path, _ = panel
    put(root / "benchmark_tools" / Path(corrected.__file__).name, Path(corrected.__file__).read_bytes())
    input_refs = [corrected.previous.pin(p) for p in sorted((root / "input").glob("*.fasta"))]
    input_paths = {r["path"] for r in input_refs}
    execution_path = root / "execution.json"
    execution = json.loads(execution_path.read_text())
    execution["verified_inputs"] = dict(inputs=input_refs)
    execution_ref = dump(execution_path, execution)
    selection_path = root / "selection.json"
    selection = json.loads(selection_path.read_text())
    for binding in selection["bindings"]:
        for arms in binding["arms"].values():
            for row in arms.values():
                row["execution"] = execution_ref
    selection["checked_inputs"] = [execution_ref]
    selection_ref = dump(selection_path, selection)
    report = json.loads(report_path.read_text())
    report["selection"] = selection_ref
    report["checked_inputs"] = [execution_ref if corrected.previous.logical(r) == str(execution_path) else r
                                for r in report["checked_inputs"] if corrected.previous.logical(r) not in input_paths]
    report_ref = dump(report_path, report)
    original = json.loads(readback_path.read_text())
    expected = reader.verify(report_path, report_ref["sha256"])
    expected.update(source=original["source"], previous_reader=original["previous_reader"])
    readback_ref = dump(readback_path, expected)
    return root, report_path, readback_path, readback_ref["sha256"]


def test_all_context_ownership_is_added_from_retained_inventories(recovery_panel):
    old = corrected.previous.plan(*recovery_panel)
    new = corrected.plan(*recovery_panel)
    assert new["configured_native_contexts"] == 8
    assert new["configured_fasta_reads"] == 16
    assert new["added_configured_fastas"] == 2
    assert new["added_checkpoint_fastas"] == old["added_checkpoint_fastas"] == 2
    assert new["preparation_revision"] == corrected.REVISION
    assert new["original_reader_kernels_changed"] is False


def test_corrected_standalone_replay_after_synthetic_original_removed(recovery_panel, tmp_path):
    original_plan = corrected.previous.plan
    prepared = corrected.prepare(*recovery_panel, tmp_path / "prepared")
    assert corrected.previous.plan is original_plan
    component = tmp_path / "relocated"
    shutil.move(tmp_path / "prepared", component)
    shutil.rmtree(recovery_panel[0])
    result_path = tmp_path / "readback.json"
    subprocess.run([sys.executable, "-I", "-B", str(component / "runner/benchmark_tools" / Path(corrected.__file__).name),
                    "replay", "--component", str(component), "--manifest-sha256", prepared["manifest"]["sha256"],
                    "--output", str(result_path)], cwd=component, check=True, capture_output=True, text=True)
    result = json.loads(result_path.read_text())
    assert result["native_readback"]["stage_rows"] == 8
    assert result["native_readback"]["native_contexts"] == 8
    assert result["added_configured_fastas"] == 2
    assert result["original_reader_kernels_changed"] is False


def test_missing_execution_inventory_is_not_replaced_with_current_file_hashes(recovery_panel, tmp_path):
    root, report_path, readback_path, expected_sha = recovery_panel
    path = root / "execution.json"
    execution = json.loads(path.read_text())
    execution["verified_inputs"]["inputs"].pop()
    changed_ref = dump(path, execution)
    selection_path = root / "selection.json"
    selection = json.loads(selection_path.read_text())
    for binding in selection["bindings"]:
        for arms in binding["arms"].values():
            for row in arms.values():
                row["execution"] = changed_ref
    selection["checked_inputs"] = [changed_ref]
    selection_ref = dump(selection_path, selection)
    report = json.loads(report_path.read_text())
    report["selection"] = selection_ref
    report["checked_inputs"] = [changed_ref if corrected.previous.logical(r) == str(path) else r for r in report["checked_inputs"]]
    report_ref = dump(report_path, report)
    expected = json.loads(readback_path.read_text())
    expected["report"] = report_ref
    expected_ref = dump(readback_path, expected)
    with pytest.raises(ValueError, match="retained execution inventory"):
        corrected.prepare(root, report_path, readback_path, expected_ref["sha256"], tmp_path / "component")
    assert not (tmp_path / "component").exists()


def test_frozen_writer_plan_is_restored_after_copy_failure(recovery_panel, tmp_path, monkeypatch):
    original_plan = corrected.previous.plan
    def fail(*args):
        assert corrected.previous.plan is not original_plan
        raise RuntimeError("synthetic copy failure")
    monkeypatch.setattr(corrected.previous, "prepare", fail)
    with pytest.raises(RuntimeError, match="synthetic copy failure"):
        corrected.prepare(*recovery_panel, tmp_path / "component")
    assert corrected.previous.plan is original_plan


def test_frozen_first_driver_and_native_sources_stay_unchanged():
    assert corrected.previous.identity(corrected.previous.__file__)["sha256"] == "9b00de075bf638a0ffd47a4d93ac5a7a97fa27e5afd430b66e4b5f317bbafe25"
    assert reader.previous.pin(reader.previous.__file__)["sha256"] == reader.PREVIOUS_SHA
