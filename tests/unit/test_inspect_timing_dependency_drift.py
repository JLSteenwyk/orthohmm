import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import inspect_timing_dependency_drift as audit
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def test_no_overwrite(tmp_path):
    with pytest.raises(FileExistsError):
        audit.inspect(tmp_path / "missing", "0" * 64, tmp_path / "missing", tmp_path)


def test_wrong_pin(tmp_path):
    source = tmp_path / "lookup.json"
    source.write_text("{}")
    with pytest.raises(ValueError):
        audit.inspect(source, "0" * 64, tmp_path / "missing", tmp_path / "out")


@pytest.mark.parametrize("ambiguous", [False, True])
def test_drift_report(tmp_path, monkeypatch, ambiguous):
    package = tmp_path / "pyparsing"
    package.mkdir()
    source = package / "__init__.py"
    source.write_text("old")
    old = record(source)
    source.write_text("new")
    prior = tmp_path / "prior.json"
    prior.write_text(json.dumps(dict(modules={"pyparsing": str(source)}, files=[old])))
    lookup = tmp_path / "lookup.json"
    lookup.write_text(json.dumps(dict(interpreters={"orthohmm": dict(reports=[record(prior)])})))
    diff = tmp_path / "diff.json"
    diff.write_text(json.dumps(dict(added=["new"], removed=["old"])))
    metadata = tmp_path / "METADATA"
    metadata.write_text("Name: example\nVersion: 1\n")
    distribution = SimpleNamespace(version="1", files=[] if ambiguous else [Path("example.dist-info/METADATA")],
                                   locate_file=lambda _: metadata)
    monkeypatch.setattr(audit.importlib.metadata, "distribution", lambda _: distribution)
    output = tmp_path / "result.json"
    args = (lookup, record(lookup)["sha256"], diff, output)
    if ambiguous:
        with pytest.raises(ValueError, match="metadata"):
            audit.inspect(*args)
        assert not output.exists()
    else:
        result = audit.inspect(*args)
        assert result["changed_imported_files"] == 1
        assert result["added_entry_count"] == result["removed_entry_count"] == 1
        assert result["scientific_execution_authorized"] is False
        assert source.read_text() == "new"
