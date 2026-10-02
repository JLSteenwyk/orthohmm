import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import run_ci_cpu_wheel as module


@pytest.fixture
def source(tmp_path, monkeypatch):
    root = tmp_path / "repo"
    names = sorted(module.REQUIRED | {"orthohmm/search/csrc/pair_align.c"})
    content = {name: (name + "\n").encode() for name in names}
    for name, raw in content.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    for name in ("requirements.txt", "tests/requirements.txt"):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# synthetic test requirements\n")

    def git(argv, **kwargs):
        if "rev-parse" in argv:
            return "a" * 40 + "\n"
        if "ls-files" in argv:
            return "\0".join(names) + "\0"
        if "show" in argv:
            return content[argv[-1].split(":", 1)[1]]
        pytest.fail("Unexpected source command: " + repr(argv))

    monkeypatch.setattr(module.subprocess, "check_output", git)
    return root, names


def test_stage_copies_committed_bytes_and_excludes_untracked_binaries(source, tmp_path):
    root, names = source
    (root / "orthohmm/inherited.so").write_bytes(b"Do not copy")
    result = module.stage(root, tmp_path / "staged")
    assert result["commit"] == "a" * 40
    assert {row["relative_path"] for row in result["files"]} == set(names)
    assert not (tmp_path / "staged/orthohmm/inherited.so").exists()
    for row in result["files"]:
        module.check(row["staged"])
        assert row["original"]["sha256"] == row["staged"]["sha256"]


@pytest.mark.parametrize("problem", ["dirty", "symlink", "incomplete", "traversal", "tracked_binary", "duplicate"])
def test_stage_rejects_unbound_or_unsafe_source(source, tmp_path, problem):
    root, names = source
    if problem == "dirty":
        (root / "setup.py").write_text("changed")
    elif problem == "symlink":
        (root / "setup.py").unlink()
        (root / "setup.py").symlink_to(root / "README.md")
    elif problem == "incomplete":
        names.remove("setup.py")
    elif problem == "traversal":
        names.append("orthohmm/../../outside")
    elif problem == "tracked_binary":
        names.append("orthohmm/inherited.so")
    else:
        names.append("setup.py")
    with pytest.raises(ValueError):
        module.stage(root, tmp_path / "staged")
    assert not (tmp_path / "staged").exists()


@pytest.mark.parametrize("problem", [None, "pip_failure", "timeout", "multiple_wheels", "verification", "source_drift"])
def test_isolated_workflow_and_failure_receipt(source, tmp_path, monkeypatch, problem):
    root, _ = source
    output = tmp_path / "run"
    monkeypatch.setattr(module.sys, "platform", "linux")
    commands = []

    def execute(argv, log, *, cwd, env):
        assert cwd == output
        commands.append((argv, env))
        log.write_text("synthetic command log\n")
        if problem == "pip_failure" and "install_dependencies" in log.name:
            raise subprocess.CalledProcessError(1, argv)
        if problem == "timeout" and "install_dependencies" in log.name:
            raise subprocess.TimeoutExpired(argv, 900)
        if "wheel" in argv:
            assert env["PATH"] == "/usr/bin:/bin" and env["ORTHOHMM_CPU_TARGET"] == "baseline"
            wheels = Path(argv[argv.index("--wheel-dir") + 1])
            wheels.mkdir()
            (wheels / "orthohmm.whl").write_bytes(b"synthetic wheel")
            if problem == "multiple_wheels":
                (wheels / "extra.whl").write_bytes(b"unexpected")

    def verify(repo, python, wheel, directory):
        assert repo == root and python == output / "venv/bin/python"
        assert wheel == output / "wheels/orthohmm.whl"
        assert directory == output / "verification"
        if problem == "verification":
            raise ValueError("synthetic verification failure")
        if problem == "source_drift":
            (root / "setup.py").write_text("changed")
        return {"status": "synthetic_test_only"}

    monkeypatch.setattr(module, "execute", execute)
    monkeypatch.setattr(module, "verify", verify)
    if problem:
        with pytest.raises((ValueError, subprocess.CalledProcessError, subprocess.TimeoutExpired)):
            module.run(root, output)
    else:
        result = module.run(root, output)
        assert result["status"] == "linux_cpu_wheel_installation_verified"
        assert len(commands) == 5
        assert not result["accuracy_evaluated"] and not result["controlled_timing"]
    receipt = json.loads((output / "result.json").read_text())
    assert receipt["status"] == ("cpu_wheel_attempt_failed" if problem else "linux_cpu_wheel_installation_verified")
    assert not receipt["publication_ready"]
    with pytest.raises(FileExistsError):
        module.run(root, output)


def test_non_linux_fails_before_any_attempt(tmp_path, monkeypatch):
    monkeypatch.setattr(module.sys, "platform", "darwin")
    with pytest.raises(NotImplementedError):
        module.run(tmp_path, tmp_path / "output")
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("problem", ["relative", "traversal", "indirect", "file", "directory", "dangling"])
def test_output_guards_precede_source_or_installation(tmp_path, monkeypatch, problem):
    monkeypatch.setattr(module.sys, "platform", "linux")
    monkeypatch.setattr(module, "stage", lambda *args: pytest.fail("Staged source before output validation"))
    output = tmp_path / "output"
    if problem == "relative":
        output = Path("relative-output")
    elif problem == "traversal":
        output = tmp_path / "escape" / ".." / "output"
    elif problem == "indirect":
        alias = tmp_path / "alias"
        alias.symlink_to(tmp_path, target_is_directory=True)
        output = alias / "output"
    elif problem == "file":
        output.write_text("preserve")
    elif problem == "directory":
        output.mkdir()
    else:
        output.symlink_to(tmp_path / "missing")
    with pytest.raises((ValueError, FileExistsError)):
        module.run(tmp_path, output)
    if problem == "file":
        assert output.read_text() == "preserve"


def test_failed_subprocess_log_is_retained(tmp_path):
    log = tmp_path / "failed.log"
    with pytest.raises(subprocess.CalledProcessError) as error:
        module.execute([sys.executable, "-I", "-c", "import sys; print('test-only failure'); sys.exit(3)"],
                       log, cwd=tmp_path, env=module.os.environ.copy())
    assert error.value.returncode == 3
    assert "test-only failure" in log.read_text()
