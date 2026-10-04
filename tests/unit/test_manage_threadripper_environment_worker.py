import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools.manage_threadripper_environment_worker import EnvironmentWorker
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def owner(tmp_path, code, monkeypatch=None):
    request = tmp_path / "request.json"
    policy = tmp_path / "policy.json"
    request.write_text('{"synthetic_test_only":true}')
    policy.write_text('{"synthetic_test_only":true}')
    calls = []
    def launch(command, **kwargs):
        calls.append((command, kwargs))
        assert kwargs["cwd"] == tmp_path
        assert kwargs["env"]["PYTHONPATH"] == str(tmp_path)
        assert kwargs["env"]["PYTHONDONTWRITEBYTECODE"] == "1"
        assert kwargs["env"]["PYTHONNOUSERSITE"] == "1"
        assert kwargs["env"]["PYTHONHASHSEED"] == "0"
        assert not any(k in kwargs["env"] for k in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"))
        # A real owned subprocess, but explicitly synthetic work, not the responder.
        return subprocess.Popen([sys.executable, "-I", "-B", "-c", code], **kwargs)
    return EnvironmentWorker(tmp_path, record(request), record(policy), tmp_path, popen=launch), calls


def receipt(tmp_path):
    result = json.loads((tmp_path / "environment_worker_lifecycle.json").read_text())
    check(result["log"])
    assert result["terminal"] and not result["automatic_retry"]
    return result


def test_owned_worker_response_join_and_environment(tmp_path, monkeypatch):
    for key in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        monkeypatch.setenv(key, "/synthetic-invalid-inheritance")
    worker, calls = owner(tmp_path, "import pathlib; print('synthetic log'); pathlib.Path('response.json').write_text('{\"synthetic\":true}')")
    with worker:
        assert worker.wait_response(tmp_path / "response.json", 5) == {"synthetic": True}
        worker.finish()
    result = receipt(tmp_path)
    assert result["status"] == "completed" and result["exit_code"] == 0
    assert "synthetic log" in (tmp_path / "environment_worker.log").read_text()
    assert len(calls) == 1
    with pytest.raises(RuntimeError, match="cannot be reused"):
        worker.__enter__()
    assert len(calls) == 1


def test_child_failure_without_response(tmp_path):
    worker, calls = owner(tmp_path, "raise SystemExit(7)")
    with pytest.raises(RuntimeError, match="without a response"):
        with worker:
            worker.wait_response(tmp_path / "missing.json", 5)
    result = receipt(tmp_path)
    assert result["exit_code"] == 7 and result["status"] == "failed_or_cancelled"
    assert len(calls) == 1


def test_nonzero_child_cannot_pass_join(tmp_path):
    worker, _ = owner(tmp_path, "raise SystemExit(3)")
    with pytest.raises(RuntimeError, match="status 3"):
        with worker: worker.finish()
    assert receipt(tmp_path)["exit_code"] == 3


@pytest.mark.parametrize("kind", [RuntimeError, KeyboardInterrupt])
def test_preparation_failure_reaps_only_owned_waiting_child(tmp_path, kind):
    worker, _ = owner(tmp_path, "import time; time.sleep(60)")
    with pytest.raises(kind, match="synthetic parent failure"):
        with worker: raise kind("synthetic parent failure")
    result = receipt(tmp_path)
    assert result["cleanup_requested"] and worker.process.poll() is not None
    assert result["parent_error_type"] == kind.__name__


def test_response_timeout_retains_failure_and_reaps(tmp_path):
    worker, _ = owner(tmp_path, "import time; time.sleep(60)")
    with pytest.raises(TimeoutError):
        with worker: worker.wait_response(tmp_path / "missing.json", .05)
    assert receipt(tmp_path)["cleanup_requested"]


def test_return_without_successful_join_is_not_completed(tmp_path):
    worker, _ = owner(tmp_path, "import time; time.sleep(60)")
    with pytest.raises(RuntimeError, match="not successfully joined"):
        with worker: pass
    assert receipt(tmp_path)["status"] == "failed_or_cancelled"


def test_spawn_failure_closes_log_and_retains_failure(tmp_path):
    worker, _ = owner(tmp_path, "pass")
    def fail(*args, **kwargs): raise OSError("synthetic spawn failure")
    worker.popen = fail
    with pytest.raises(OSError, match="spawn failure"):
        with worker: pass
    assert worker.log.closed and receipt(tmp_path)["parent_error_type"] == "OSError"


def test_term_ignoring_child_is_killed_and_reaped(tmp_path):
    code = "import signal,time,pathlib; signal.signal(signal.SIGTERM,signal.SIG_IGN); pathlib.Path('ready.json').write_text('{}'); time.sleep(60)"
    worker, _ = owner(tmp_path, code)
    with pytest.raises(RuntimeError, match="synthetic failure"):
        with worker:
            worker.wait_response(tmp_path / "ready.json", 5)
            raise RuntimeError("synthetic failure")
    result = receipt(tmp_path)
    assert result["kill_requested"] and worker.process.poll() is not None


def test_join_timeout_still_cleans_up_owned_child(tmp_path):
    worker, _ = owner(tmp_path, "import time; time.sleep(60)")
    with pytest.raises(subprocess.TimeoutExpired):
        with worker:
            original = worker.process.wait
            first = True
            def wait(timeout=None):
                nonlocal first
                if first:
                    first = False
                    raise subprocess.TimeoutExpired("synthetic join timeout", timeout)
                return original(timeout=timeout)
            worker.process.wait = wait
            worker.finish()
    assert receipt(tmp_path)["cleanup_requested"] and worker.process.poll() is not None


def preparing_owner(tmp_path, monkeypatch, *, delay=.1):
    worker, calls = owner(tmp_path, '')
    request_path = Path(worker.request_ref['path'])
    request_path.write_text(json.dumps(dict(job_id=42, index=21, synthetic_test_only=True)))
    worker.request_ref = record(request_path)
    payload = dict(schema='threadripper_environment_worker_prepared_v1',
        status='prepared_waiting_for_release_request', request=worker.request_ref,
        policy=worker.policy_ref, job_id=42, index=21, job_scope='/synthetic/job_42',
        native_release_authorized=False, scientific_timings_admitted=False)
    code = f"""
import json, os, pathlib, time
time.sleep({delay!r})
data = json.loads({json.dumps(payload)!r})
data.update(pid=os.getpid(), prepared_unix_ns=time.time_ns(),
            boot_id=pathlib.Path('/proc/sys/kernel/random/boot_id').read_text().strip())
path = pathlib.Path('environment_worker_prepared.json')
pending = path.with_suffix('.pending')
pending.write_text(json.dumps(data))
pending.rename(path)
while not pathlib.Path('synthetic_release.json').exists():
    time.sleep(.01)
pathlib.Path('response.json').write_text('{{"synthetic":true}}')
"""
    def launch(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.Popen([sys.executable, '-I', '-B', '-c', code], **kwargs)
    worker.popen = launch
    original_read = Path.read_text
    def read(path, *args, **kwargs):
        if worker.process is not None and path == Path(f'/proc/{worker.process.pid}/cgroup'):
            return '0::/synthetic/job_42/step_batch/task_0\n'
        return original_read(path, *args, **kwargs)
    monkeypatch.setattr(Path, 'read_text', read)
    return worker, calls


def test_preparation_barrier_precedes_release_and_retains_bound_receipt(tmp_path, monkeypatch):
    worker, calls = preparing_owner(tmp_path, monkeypatch)
    with worker:
        prepared = worker.wait_prepared(seconds=5)
        check(prepared)
        value = json.loads(Path(prepared['path']).read_text())
        assert value['pid'] == worker.process.pid and not value['native_release_authorized']
        assert not (tmp_path / 'synthetic_release.json').exists()
        with pytest.raises(RuntimeError, match='newly owned'):
            worker.wait_prepared(seconds=5)
        (tmp_path / 'synthetic_release.json').write_text('{}')
        assert worker.wait_response(tmp_path / 'response.json', seconds=5) == {'synthetic': True}
        worker.finish()
    assert receipt(tmp_path)['prepared'] == prepared
    assert receipt(tmp_path)['status'] == 'completed' and len(calls) == 1


@pytest.mark.parametrize('key,value', [
    ('schema', 'wrong'), ('status', 'passed'), ('job_id', 43), ('index', 20),
    ('pid', 1), ('request', {}), ('policy', {}), ('job_scope', '/synthetic/job_43'),
    ('boot_id', 'wrong-boot'), ('prepared_unix_ns', 0), ('prepared_unix_ns', True),
    ('native_release_authorized', True), ('scientific_timings_admitted', True),
])
def test_preparation_cannot_borrow_or_admit_another_worker(tmp_path, monkeypatch, key, value):
    worker, _ = preparing_owner(tmp_path, monkeypatch, delay=0)
    with pytest.raises(ValueError):
        with worker:
            path = tmp_path / 'environment_worker_prepared.json'
            prepared = worker.wait_response(path, seconds=5)
            prepared[key] = value
            path.write_text(json.dumps(prepared))
            worker.wait_prepared(seconds=5)
    assert receipt(tmp_path)['status'] == 'failed_or_cancelled'
    assert not (tmp_path / 'synthetic_release.json').exists()


def test_preparation_timeout_reaps_child_without_starting_native_work(tmp_path, monkeypatch):
    worker, _ = preparing_owner(tmp_path, monkeypatch, delay=60)
    with pytest.raises(TimeoutError):
        with worker:
            worker.wait_prepared(seconds=.05)
    assert receipt(tmp_path)['cleanup_requested']
    assert not (tmp_path / 'synthetic_release.json').exists()


def test_dead_prepared_worker_cannot_release_native_work(tmp_path, monkeypatch):
    worker, _ = preparing_owner(tmp_path, monkeypatch, delay=0)
    with pytest.raises((ValueError, FileNotFoundError)):
        with worker:
            worker.wait_response(tmp_path / 'environment_worker_prepared.json', seconds=5)
            worker.process.terminate()
            worker.process.wait(timeout=5)
            worker.wait_prepared(seconds=5)
    assert receipt(tmp_path)['status'] == 'failed_or_cancelled'
