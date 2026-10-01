import hashlib
from pathlib import Path

import pytest

from benchmark_tools import qfo_parameter_source_lineage as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture(tmp_path, monkeypatch):
    old, current = b"historical report exporter\n", b"current recovered reporting exporter\n"
    for prefix, value in (("OLD", old), ("CURRENT", current)):
        monkeypatch.setattr(module, prefix + "_BYTES", len(value))
        monkeypatch.setattr(module, prefix + "_SHA", hashlib.sha256(value).hexdigest())
    for name, value in ((module.SOURCE, current), (module.ARCHIVE, old)):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(value)
    calls = []
    def git(argv):
        calls.append(argv)
        revision = argv[-1].split(":", 1)[0]
        assert argv[:-1] == ["git", "-C", str(tmp_path), "show"]
        assert argv[-1].endswith(":" + module.SOURCE)
        return {module.OLD_COMMIT: old, module.CURRENT_COMMIT: current}[revision]
    monkeypatch.setattr(module.subprocess, "check_output", git)
    original = module.identity(tmp_path, module.SOURCE, len(old), module.OLD_SHA)
    data = tmp_path / "scientific_data.tsv"
    data.write_text("unchanged data\n")
    return [original, record(data), record(tmp_path / module.SOURCE)], calls


def test_explicit_binding_preserves_original_and_checks_both_git_blobs(tmp_path, monkeypatch):
    refs, calls = fixture(tmp_path, monkeypatch)
    original = dict(refs[0])
    checked, lineage = module.check_records(refs + [refs[0]], tmp_path, True)
    assert refs[0] == original
    assert original not in checked
    assert record(tmp_path / module.ARCHIVE) in checked
    assert record(tmp_path / module.SOURCE) in checked
    assert refs[1] in checked
    assert lineage["original_record"] == original
    assert lineage["historical_copy"] == record(tmp_path / module.ARCHIVE)
    assert lineage["current_record"] == refs[2]
    assert lineage["scientific_source_substitution"] is False
    assert len(calls) == 2


def test_default_route_does_not_silently_bind_stale_source(tmp_path, monkeypatch):
    refs, calls = fixture(tmp_path, monkeypatch)
    with pytest.raises(ValueError):
        module.check_records(refs, tmp_path)
    assert not calls


@pytest.mark.parametrize("problem", ["missing_original", "unknown_original", "current_changed",
    "archive_changed", "old_git_changed", "current_git_changed", "archive_symlink",
    "current_symlink", "data_changed", "unknown_same_path", "other_source_changed"])
def test_binding_never_weakens_other_evidence_or_known_source_identities(tmp_path, monkeypatch, problem):
    refs, _ = fixture(tmp_path, monkeypatch)
    if problem == "missing_original":
        refs.pop(0)
    elif problem == "unknown_original":
        refs[0] = {**refs[0], "sha256": "wrong"}
    elif problem in ("current_changed", "archive_changed"):
        (tmp_path / (module.SOURCE if problem == "current_changed" else module.ARCHIVE)).write_bytes(b"changed")
    elif problem in ("old_git_changed", "current_git_changed"):
        original = module.subprocess.check_output
        bad = module.OLD_COMMIT if problem == "old_git_changed" else module.CURRENT_COMMIT
        monkeypatch.setattr(module.subprocess, "check_output",
            lambda argv: b"changed Git blob" if argv[-1].startswith(bad + ":") else original(argv))
    elif problem in ("archive_symlink", "current_symlink"):
        path = tmp_path / (module.ARCHIVE if problem == "archive_symlink" else module.SOURCE)
        target = tmp_path / "symlink_target.py"
        target.write_bytes(path.read_bytes())
        path.unlink()
        path.symlink_to(target)
    elif problem == "data_changed":
        Path(refs[1]["path"]).write_text("changed scientific data\n")
    elif problem == "unknown_same_path":
        refs.append({**refs[2], "sha256": "unknown"})
    else:
        path = tmp_path / "scientific_helper.py"
        path.write_text("before\n")
        refs.append(record(path))
        path.write_text("after\n")
    with pytest.raises(ValueError):
        module.check_records(refs, tmp_path, True)


def test_nonbinding_records_deduplicate_only_identical_evidence(tmp_path):
    path = tmp_path / "data.tsv"
    path.write_text("data\n")
    ref = record(path)
    checked, lineage = module.check_records([ref, ref], tmp_path)
    assert checked == [ref] and lineage is None
    with pytest.raises(ValueError, match="Conflicting"):
        module.check_records([ref, {**ref, "sha256": "different"}], tmp_path)
