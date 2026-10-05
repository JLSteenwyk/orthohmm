import copy
import hashlib
import json
import sys

import pytest

from benchmark_tools import recheck_benchmark_files as audit


@pytest.fixture
def sample(tmp_path):
    path = tmp_path / "input.bin"
    path.write_bytes(b"retained input")
    ref = {"path": str(path), "bytes": path.stat().st_size,
           "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    rows = [{"dataset": dataset, "key": f"method{method}", "input_records": [ref],
             "output_records": [ref]} for dataset in ("OrthoBench", "QfO", "ThreeKingdoms")
            for method in range(8)]
    rows[15]["input_records"] = {"bpo": ref}
    return path, ref, {"rows": rows}


def test_shared_files_and_named_downstream_inputs_are_not_independent_repeats(sample):
    _, ref, register = sample
    pins, rows = audit.selections(register)
    assert pins == {ref["path"]: ref} and len(rows) == 24
    assert rows[15]["files"][0]["entry"] == "bpo"
    assert audit.inspect(ref)["status"] == "matches"


@pytest.mark.parametrize("kind", ["missing", "mismatch", "directory", "symlink"])
def test_file_outcomes_are_retained(sample, kind):
    path, ref, _ = sample
    if kind == "missing":
        path.unlink()
    elif kind == "mismatch":
        path.write_bytes(b"different")
    elif kind == "directory":
        path.unlink()
        path.mkdir()
    else:
        link = path.with_name("alias")
        link.symlink_to(path)
        ref = dict(ref, path=str(link))
    expected = "nonregular" if kind in ("directory", "symlink") else kind
    assert audit.inspect(ref)["status"] == expected


def test_file_mutation_during_read_is_not_a_match(sample, monkeypatch):
    _, ref, _ = sample
    original = audit.identity
    calls = 0
    def unstable(info):
        nonlocal calls
        calls += 1
        value = original(info)
        if calls == 2:
            value[-1] += 1
        return value
    monkeypatch.setattr(audit, "identity", unstable)
    assert audit.inspect(ref)["status"] == "changed_during_check"


@pytest.mark.parametrize("kind", ["missing_row", "duplicate_method", "different_panel", "conflicting_pin", "bad_pin", "bad_container"])
def test_invalid_register_refused(sample, kind):
    _, _, original = sample
    register = copy.deepcopy(original)
    row = register["rows"][0]
    if kind == "missing_row":
        register["rows"].pop()
    elif kind == "duplicate_method":
        row["key"] = "method1"
    elif kind == "different_panel":
        row["key"] = "other"
    elif kind == "conflicting_pin":
        row["input_records"][0] = dict(row["input_records"][0], sha256="0" * 64)
    elif kind == "bad_pin":
        row["input_records"][0]["bytes"] = True
    else:
        row["input_records"] = "unknown"
    with pytest.raises(ValueError):
        audit.selections(register)


def test_complete_run_does_not_upgrade_historical_or_scientific_provenance(sample, tmp_path):
    _, _, data = sample
    path = tmp_path / "register.json"
    content = json.dumps(data).encode()
    path.write_bytes(content)
    result = audit.run(path, hashlib.sha256(content).hexdigest())
    assert result["status_counts"] == {"matches": 1}
    assert result["all_selected_files_match"] is True
    assert result["historical_input_consumption_proven"] is False
    assert result["scores_recomputed"] is False and result["publication_ready"] is False
    with pytest.raises(ValueError, match="checksum"):
        audit.run(path, "0" * 64)


def test_unreadable_file_is_not_a_missing_or_matching_file(sample, monkeypatch):
    _, ref, _ = sample
    def denied(*args):
        raise PermissionError(13, "permission denied")
    monkeypatch.setattr(audit.os, "open", denied)
    result = audit.inspect(ref)
    assert result["status"] == "unreadable" and result["errno"] == 13


def test_missing_raw_file_remains_explicit(sample, tmp_path):
    raw, _, data = sample
    raw.unlink()
    path = tmp_path / "register.json"
    content = json.dumps(data).encode()
    path.write_bytes(content)
    result = audit.run(path, hashlib.sha256(content).hexdigest())
    assert result["status_counts"] == {"missing": 1}
    assert not result["all_selected_files_match"]


def test_main_preserves_results_and_refuses_raw_output_alias(sample, tmp_path, monkeypatch):
    raw, _, data = sample
    path = tmp_path / "register.json"
    content = json.dumps(data).encode()
    path.write_bytes(content)
    output = tmp_path / "result.json"
    argv = ["recheck", "--register", str(path), "--sha256", hashlib.sha256(content).hexdigest(),
            "--output", str(output)]
    monkeypatch.setattr(sys, "argv", argv)
    audit.main()
    saved = output.read_bytes()
    with pytest.raises(ValueError, match="already exists"):
        audit.main()
    assert output.read_bytes() == saved
    raw.unlink()
    argv[-1] = str(raw)
    with pytest.raises(ValueError, match="aliases selected"):
        audit.main()
    assert not raw.exists()
    child = tmp_path / "child"
    child.mkdir()
    argv[-1] = str(child / ".." / raw.name)
    with pytest.raises(ValueError, match="aliases selected"):
        audit.main()
    assert not raw.exists()
