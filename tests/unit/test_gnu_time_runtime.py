import os
from pathlib import Path
import subprocess
import types

import pytest

from tests import gnu_time_runtime as runtime


def executable(tmp_path, name):
    path = tmp_path / name
    path.write_text("fixture, not executed\n")
    path.chmod(0o700)
    return str(path)


def probe(monkeypatch, callback):
    monkeypatch.setattr(runtime, "subprocess", types.SimpleNamespace(
        run=callback, TimeoutExpired=subprocess.TimeoutExpired))


def test_explicit_gnu_binary_is_verified_without_fallback(tmp_path, monkeypatch):
    binary = executable(tmp_path, "GNU time with spaces")
    monkeypatch.setenv("ORTHOHMM_TEST_GNU_TIME", binary)
    monkeypatch.setattr(runtime, "shutil", types.SimpleNamespace(
        which=lambda name: pytest.fail("Explicit selection must not search PATH")))
    calls = []

    def version(argv, **kwargs):
        calls.append((argv, kwargs))
        return subprocess.CompletedProcess(argv, 0, "GNU Time fixture\n", "")

    probe(monkeypatch, version)
    assert runtime.select_gnu_time() == binary
    assert calls == [([binary, "--version"], dict(capture_output=True, text=True,
                                               check=False, timeout=5))]


@pytest.mark.parametrize("invalid", ["empty", "relative", "missing", "nonexecutable",
                                    "bsd_error", "bsd_success", "os_error", "timeout"])
def test_invalid_explicit_binary_fails_closed(tmp_path, monkeypatch, invalid):
    binary = executable(tmp_path, "time")
    if invalid == "empty":
        binary = ""
    elif invalid == "relative":
        binary = "time"
    elif invalid == "missing":
        Path(binary).unlink()
    elif invalid == "nonexecutable":
        Path(binary).chmod(0o600)
    monkeypatch.setenv("ORTHOHMM_TEST_GNU_TIME", binary)
    monkeypatch.setattr(runtime, "shutil", types.SimpleNamespace(
        which=lambda name: pytest.fail("Invalid explicit selection must not fall back")))
    calls = []

    def version(argv, **kwargs):
        calls.append(argv)
        if invalid == "os_error":
            raise PermissionError("fixture")
        if invalid == "timeout":
            raise subprocess.TimeoutExpired(argv, 5)
        return subprocess.CompletedProcess(argv, int(invalid == "bsd_error"),
                                           "BSD time", "illegal option")

    probe(monkeypatch, version)
    with pytest.raises(ValueError, match="ORTHOHMM_TEST_GNU_TIME"):
        runtime.select_gnu_time()
    assert len(calls) == int(invalid in {"bsd_error", "bsd_success", "os_error", "timeout"})


def test_discovery_prefers_gtime(tmp_path, monkeypatch):
    binary = executable(tmp_path, "gtime")
    monkeypatch.delenv("ORTHOHMM_TEST_GNU_TIME", raising=False)
    monkeypatch.setattr(runtime, "shutil", types.SimpleNamespace(
        which=lambda name: binary if name == "gtime" else "/usr/bin/time"))
    calls = []

    def version(argv, **kwargs):
        calls.append(argv)
        return subprocess.CompletedProcess(argv, 0, "GNU Time fixture", "")

    probe(monkeypatch, version)
    assert runtime.select_gnu_time() == binary
    assert calls == [[binary, "--version"]]


def test_bsd_version_failure_is_unavailable_not_called_process_error(monkeypatch):
    monkeypatch.delenv("ORTHOHMM_TEST_GNU_TIME", raising=False)
    monkeypatch.setattr(runtime, "shutil", types.SimpleNamespace(which=lambda name: "/usr/bin/time"))
    calls = []

    def version(argv, **kwargs):
        calls.append(argv)
        return subprocess.CompletedProcess(argv, 1, "", "illegal option")

    probe(monkeypatch, version)
    with pytest.raises(FileNotFoundError, match="GNU Time"):
        runtime.select_gnu_time()
    assert calls == [["/usr/bin/time", "--version"]]


@pytest.mark.parametrize("method", ["run", "Popen"])
@pytest.mark.parametrize("argv", [["/usr/bin/time", "-v", "-o", "log with spaces", "/bin/true"],
                                  ("/usr/bin/time", "--version"),
                                  ["/bin/echo", "/usr/bin/time"], [], "/usr/bin/time --version"])
def test_proxy_changes_only_leading_literal_executor(method, argv):
    calls = []
    sentinel = object()

    def launch(args, *positional, **keywords):
        calls.append((args, positional, keywords))
        return sentinel

    original = types.SimpleNamespace(run=launch, Popen=launch, PIPE=subprocess.PIPE)
    proxy = runtime.GnuTimeSubprocess("/selected/gtime", original)
    before = list(argv) if isinstance(argv, (list, tuple)) else argv
    assert getattr(proxy, method)(argv, "positional", env={"PATH": "/frozen"}) is sentinel
    expected = (["/selected/gtime", *argv[1:]]
                if isinstance(argv, (list, tuple)) and argv and argv[0] == "/usr/bin/time" else argv)
    assert calls == [(expected, ("positional",), {"env": {"PATH": "/frozen"}})]
    assert (list(argv) if isinstance(argv, (list, tuple)) else argv) == before
    assert proxy.PIPE == subprocess.PIPE
    assert original.run is launch and original.Popen is launch


def test_explicit_binding_is_module_local_and_restored(bind_test_gnu_time, monkeypatch):
    original_run, original_popen = subprocess.run, subprocess.Popen
    module = types.SimpleNamespace(subprocess=subprocess)
    proxy = bind_test_gnu_time(module)
    assert module.subprocess is proxy
    assert proxy.subprocess_module is subprocess
    assert subprocess.run is original_run and subprocess.Popen is original_popen
    monkeypatch.undo()
    assert module.subprocess is subprocess


@pytest.mark.parametrize("exit_code", [0, 7])
def test_bound_proxy_runs_real_child_and_keeps_native_log(tmp_path, gnu_time_binary, exit_code):
    import sys

    proxy = runtime.GnuTimeSubprocess(gnu_time_binary, subprocess)
    output = tmp_path / "native time with spaces.log"
    command = ["/usr/bin/time", "-v", "-o", str(output), sys.executable,
               "-c", f"print('native evidence'); raise SystemExit({exit_code})"]
    result = proxy.run(command, capture_output=True, text=True, env=os.environ.copy())
    assert result.returncode == exit_code
    assert result.stdout == "native evidence\n"
    assert f"Exit status: {exit_code}" in output.read_text()
    assert command[0] == "/usr/bin/time"
