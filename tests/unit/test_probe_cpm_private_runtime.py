import json
import os
from pathlib import Path
import signal
import sys

import pytest

from benchmark_tools import probe_cpm_private_runtime as control


def test_environment_matches_allocator_control_without_forced_gc(monkeypatch):
    for key in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH"):
        monkeypatch.setenv(key, "unwanted")
    env, overrides, removed = control.environment(Path("/frozen"))
    assert all(key not in env for key in removed)
    assert overrides["PYTHONMALLOC"] == "debug" and overrides["PYTHONFAULTHANDLER"] == "1"
    assert overrides["PYTHONPATH"] == "/frozen"
    assert all(overrides[key] == "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"))


def package(tmp_path):
    site = tmp_path / "site"
    site.mkdir()
    path = site / "example.py"
    path.write_bytes(b"pass\n")
    ref = control.record(path)
    row = dict(member="example.py", bytes=ref["bytes"], sha256=ref["sha256"])
    return site, dict(packages=[dict(matched_files=1, matched=[row])])


def test_package_pins_reuse_retained_payload_identities(tmp_path):
    site, audit = package(tmp_path)
    rows = control.package_records(audit, site)
    control.check(rows)
    (site / "example.py").write_bytes(b"changed")
    with pytest.raises(ValueError, match="Changed"):
        control.check(rows)


def test_base_alias_normalization_keeps_expected_hashes(tmp_path):
    base = tmp_path / "base"
    base.mkdir()
    binary = base / "python3.10"
    binary.write_bytes(b"original")
    alias = base / "python"
    alias.symlink_to(binary)
    expected = control.record(binary)
    plan = dict(checked_records=[expected, dict(expected, path=str(alias)),
                                dict(expected, path="/unrelated/runtime")])
    assert control.bound_base_records(plan, base) == [expected]
    binary.write_bytes(b"changed")
    assert control.bound_base_records(plan, base) == [expected]
    with pytest.raises(ValueError, match="Changed"):
        control.check(control.bound_base_records(plan, base))


def test_base_alias_cannot_escape_restored_runtime(tmp_path):
    base = tmp_path / "base"
    base.mkdir()
    outside = tmp_path / "outside"
    outside.write_bytes(b"outside")
    alias = base / "python"
    alias.symlink_to(outside)
    plan = dict(checked_records=[dict(control.record(outside), path=str(alias))])
    with pytest.raises(ValueError, match="escapes"):
        control.bound_base_records(plan, base)


@pytest.mark.parametrize("change", ["count", "duplicate", "traversal", "absolute", "backslash", "symlink"])
def test_invalid_package_binding_rejected(tmp_path, change):
    site, audit = package(tmp_path)
    group = audit["packages"][0]
    if change == "count": group["matched_files"] = 2
    elif change == "duplicate": audit["packages"].append(group.copy())
    elif change in {"traversal", "absolute", "backslash"}:
        group["matched"][0]["member"] = {"traversal": "../outside", "absolute": "/outside", "backslash": "a\\b"}[change]
    else:
        path = site / "example.py"
        outside = tmp_path / "outside"
        path.rename(outside)
        path.symlink_to(outside)
    with pytest.raises(ValueError):
        control.package_records(audit, site)


def test_conflicting_source_pins_are_not_silently_deduplicated():
    row = dict(path="/same", bytes=1, sha256="a")
    assert control.unique_records([row, row.copy()]) == [row]
    with pytest.raises(ValueError, match="Conflicting"):
        control.unique_records([row, dict(row, sha256="b")])


def test_exact_child_metadata_and_partition_required(tmp_path):
    path = tmp_path / "partition"
    path.write_bytes(b"a b\n")
    ref = control.record(path)
    original = dict(genes=2, groups=1, output=dict(ref, path="/original"))
    result = dict(genes=2, groups=1, output=ref)
    assert control.validate_result(result, original, path) == ref
    with pytest.raises(ValueError, match="metadata"):
        control.validate_result(dict(result, groups=2), original, path)
    with pytest.raises(ValueError, match="partition"):
        control.validate_result(result, dict(original, output=dict(ref, sha256="0" * 64)), path)


def test_native_child_completion_and_signal_failure_are_retained(tmp_path):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    code, timeout = control.execute([sys.executable, "-B", "-c", "print(42)"],
        tmp_path, env, tmp_path / "out", tmp_path / "err", 10)
    assert (code, timeout) == (0, False) and (tmp_path / "out").read_text().strip() == "42"
    code, timeout = control.execute([sys.executable, "-B", "-c", "import os,signal;os.kill(os.getpid(),signal.SIGTERM)"],
        tmp_path, env, tmp_path / "signal.out", tmp_path / "signal.err", 10)
    assert (code, timeout) == (-signal.SIGTERM, False)


def test_timeout_kills_only_owned_child_process_group(tmp_path):
    code, timeout = control.execute([sys.executable, "-B", "-c", "import time;time.sleep(10)"],
        tmp_path, os.environ.copy(), tmp_path / "out", tmp_path / "err", .05)
    assert code == -signal.SIGKILL and timeout is True


def test_preflight_failure_retained_without_refinement_attempt(tmp_path, monkeypatch):
    monkeypatch.setattr(control.os, "sched_setaffinity", lambda *_: None)
    def fail(*_): raise ValueError("bad prerequisite")
    monkeypatch.setattr(control, "_run", fail)
    output = tmp_path / "fresh"
    with pytest.raises(ValueError, match="bad prerequisite"):
        control.run(tmp_path, output, "unused")
    report = json.loads((output / "report.json").read_bytes())
    assert report["status"] == "private_runtime_preflight_failed"
    assert report["refinement_attempts"] == 0 and report["seed_admitted"] is False
    before = (output / "report.json").read_bytes()
    with pytest.raises(FileExistsError):
        control.run(tmp_path, output, "unused")
    assert (output / "report.json").read_bytes() == before
