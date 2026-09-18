import json
from pathlib import Path

import pytest

from benchmark_tools.audit_publication_figures import inspect_manifest, records, PANELS
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture(tmp_path):
    directory = tmp_path / "figures"
    directory.mkdir()
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = directory / ("figure." + extension)
        path.write_bytes(b"synthetic format placeholder")
        outputs.append(record(path))
    source = tmp_path / "source.py"
    source.write_text("pass\n")
    data = {"outputs": outputs, "source": record(source)}
    path = directory / "manifest.json"
    path.write_text(json.dumps(data))
    return path, data


def test_checks_identity_and_separates_tracking(tmp_path):
    path, data = fixture(tmp_path)
    result = inspect_manifest(path, tmp_path, {"source.py"})
    assert result["status"] == "all_recorded_bytes_match"
    assert result["output_count"] == 3
    assert sum(row["tracked_in_main_repository"] for row in result["files"]) == 1
    assert len(result["files"]) == 4


@pytest.mark.parametrize("mutation", ["changed", "missing"])
def test_preserves_integrity_failures(tmp_path, mutation):
    path, data = fixture(tmp_path)
    source = Path(data["source"]["path"])
    if mutation == "changed":
        source.write_text("changed\n")
    else:
        source.unlink()
    result = inspect_manifest(path, tmp_path, set())
    assert result["status"] == "integrity_failure"
    assert result["files"][-1]["status"] == mutation


@pytest.mark.parametrize("mutation", ["duplicate", "format", "outside", "conflict", "schema"])
def test_rejects_invalid_inventory(tmp_path, mutation):
    path, data = fixture(tmp_path)
    if mutation == "duplicate":
        data["outputs"].append(data["outputs"][0])
    elif mutation == "format":
        data["outputs"].pop()
    elif mutation == "outside":
        data["outputs"][0]["path"] = str(tmp_path / "elsewhere.png")
    elif mutation == "conflict":
        data["extra"] = {**data["source"], "bytes": 999}
    else:
        data["source"]["unexpected"] = True
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        inspect_manifest(path, tmp_path, set())


def test_nested_record_inventory():
    item = {"path": "/tmp/a", "bytes": 1, "sha256": "a"}
    assert list(records({"nested": [item, {"deeper": item}], "text": "ignored"})) == [item, item]


def test_original_qfo_factorial_is_explicitly_retained():
    assert "qfo_factorial_swiss_figure_20260918" in PANELS
    assert len(set(PANELS)) == len(PANELS)
    assert not any("corrected_factorial" in panel for panel in PANELS)
