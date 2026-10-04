import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from benchmark_tools.results import capture_shared_terminal_controller_20261004 as collector


@pytest.mark.parametrize('index', [True, 0, 20, 27, -1, None])
def test_invalid_identity_never_queries(index, monkeypatch):
    monkeypatch.setattr(collector.subprocess, 'run', lambda *a, **k: pytest.fail('Must not query'))
    with pytest.raises(ValueError):
        collector.capture(index)


def fixture(tmp_path, monkeypatch, *, state='COMPLETED', wrong_digest=False, failed_query=False):
    monkeypatch.setattr(collector.panel, 'WORK', tmp_path)
    request = collector.panel.save(tmp_path / 'request.json', dict(index=24, job_id=999,
        allocation_cwd=str(collector.ROOT), scheduler_command=str(collector.panel.SCRIPT)))
    collector.panel.save(tmp_path / 'launch_24.json', dict(index=24, job_id=999,
        request=request, status='repaired_shared_job_bound_and_released',
        execution_scope=collector.panel.executor.SHARED_SCOPE, automatic_retry=False,
        uncontended_timing=False))
    raw = (f'JobId=999 JobState={state} ExitCode={"0:0" if state=="COMPLETED" else "1:0"} '
        f'Partition=gpu NodeList=bizon NumNodes=1 NumCPUs=64 NumTasks=1 CPUs/Task=64 '
        f'OverSubscribe=OK MinMemoryNode=128G Requeue=0 Restarts=0 '
        f'Command={collector.panel.SCRIPT} WorkDir={collector.ROOT} '
        f'TimeLimit={collector.panel.executor.TIME_LIMIT} '
        f'Comment={"wrong" if wrong_digest else request["sha256"]} '
        'StartTime=2026-10-04T00:00:00 EndTime=2026-10-04T01:00:00')
    calls = []
    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        if failed_query:
            raise subprocess.CalledProcessError(1, argv, stderr='Invalid job id')
        return SimpleNamespace(returncode=0, stdout=raw, stderr='')
    monkeypatch.setattr(collector.subprocess, 'run', run)
    return request, raw, calls


@pytest.mark.parametrize('state', ['COMPLETED', 'FAILED', 'TIMEOUT', 'CANCELLED'])
def test_terminal_capture_preserves_raw_failure_and_never_admits(tmp_path, monkeypatch, state):
    request, raw, calls = fixture(tmp_path, monkeypatch, state=state)
    pin = collector.capture(24)
    record = json.loads(Path(pin['path']).read_text())
    assert record['stdout'] == raw and record['request'] == request
    assert record['native_review_completed'] is False
    assert record['next_submission_authorized'] is False
    assert record['scientific_timings_admitted'] is False
    assert calls == [(['scontrol', 'show', 'job', '999', '--oneliner'],
        dict(capture_output=True, text=True, timeout=5, check=True))]
    assert record['started_unix_ns'] <= record['finished_unix_ns']


@pytest.mark.parametrize('mode', ['live', 'wrong_digest', 'failed_query', 'existing'])
def test_invalid_or_existing_evidence_is_not_written(tmp_path, monkeypatch, mode):
    _, _, calls = fixture(tmp_path, monkeypatch,
        state='RUNNING' if mode=='live' else 'COMPLETED',
        wrong_digest=mode=='wrong_digest', failed_query=mode=='failed_query')
    output = tmp_path / 'terminal_controller_999.json'
    if mode == 'existing': output.write_text('preserved')
    error = FileExistsError if mode=='existing' else subprocess.CalledProcessError if mode=='failed_query' else ValueError
    with pytest.raises(error): collector.capture(24)
    assert output.read_text() == 'preserved' if mode=='existing' else not output.exists()
    if mode=='existing': assert not calls
