import json
import signal
import subprocess
import sys

import pytest

from benchmark_tools.run_wgd_application import METHODS, check_copies, execute, pinned, select
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools import run_wgd_application as application


def test_pinned_rejects_mutation(tmp_path):
    path = tmp_path / "manifest.json"
    path.write_text('{"a": 1}')
    identity = record(path)
    assert pinned(identity) == {"a": 1}
    path.write_text('{"a": 2}')
    with pytest.raises(ValueError, match="Changed pinned"):
        pinned(identity)


def test_select_requires_authorization_and_full_inventory(tmp_path):
    path = tmp_path / "plan.json"
    path.write_text(json.dumps({"runs": [{"method": m} for m in METHODS], "output_root": str(tmp_path / "out")}))
    spec = {"command_plan": record(path), "execution_authorized": True, "purpose": "wgd_application"}
    assert select(spec, 3)[2] == tmp_path / "out/sonicparanoid"
    with pytest.raises(ValueError, match="index"):
        select(spec, -1)
    spec["execution_authorized"] = False
    with pytest.raises(ValueError, match="authorized"):
        select(spec, 0)


def test_copy_check_rejects_changes_and_extra_fasta(tmp_path):
    path = tmp_path / "a.fasta"
    path.write_text(">a\nMA\n")
    original = record(path)
    assert len(check_copies(tmp_path, [original])) == 1
    (tmp_path / "b.fasta").write_text(">b\nMA\n")
    with pytest.raises(ValueError, match="copies"):
        check_copies(tmp_path, [original])


@pytest.mark.parametrize("code", [0, 3])
def test_execute_preserves_exit_and_log(tmp_path, code):
    log = tmp_path / "native.log"
    result = execute([sys.executable, "-c", f"print('native'); raise SystemExit({code})"], tmp_path, log, 5)
    assert result["exit_code"] == code and not result["timed_out"]
    assert result["cleanup_source"] == record(application.owned_process_group.__file__)
    assert log.read_text() == "native\n"
    with pytest.raises(FileExistsError):
        execute([sys.executable, "-c", "pass"], tmp_path, log, 5)


def test_timeout_is_failure(tmp_path):
    result = execute([sys.executable, "-c", "import time; time.sleep(30)"], tmp_path, tmp_path / "native.log", 0.05)
    assert result["timed_out"] and result["exit_code"] != 0


@pytest.mark.parametrize("gone", [False, True])
def test_timeout_permission_error_is_not_silently_ignored(tmp_path, monkeypatch, gone):
    signals = []

    class Command:
        pid = 999999
        returncode = None
        polls = 0

        def wait(self, timeout=None):
            if timeout is not None:
                raise subprocess.TimeoutExpired("fixture", timeout)
            return self.returncode

        def poll(self):
            self.polls += 1
            self.returncode = -signal.SIGTERM
            return self.returncode

    command = Command()

    def popen(*args, **kwargs):
        assert kwargs["start_new_session"] is True
        return command

    def killpg(pid, signum):
        assert pid == command.pid
        signals.append(signum)
        if command.polls and gone:
            raise ProcessLookupError("gone after reap")
        raise PermissionError("fixture permission denied")

    monkeypatch.setattr(application.subprocess, "Popen", popen)
    monkeypatch.setattr(application.os, "killpg", killpg)
    log = tmp_path / "native.log"
    if gone:
        result = execute(["fixture"], tmp_path, log, .05)
        assert result["timed_out"] and result["exit_code"] == -signal.SIGTERM
        assert signals == [signal.SIGTERM, 0]
    else:
        with pytest.raises(PermissionError, match="fixture permission denied"):
            execute(["fixture"], tmp_path, log, .05)
    assert log.exists()
