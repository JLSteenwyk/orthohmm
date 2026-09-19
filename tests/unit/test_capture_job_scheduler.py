import pytest

from benchmark_tools.capture_job_scheduler import terminal_record
from benchmark_tools.capture_array_scheduler import REQUIRED
from benchmark_tools import capture_job_scheduler as module
from types import SimpleNamespace
import json


def raw(**changes):
    fields = {key: "1" for key in REQUIRED - {"ArrayJobId", "ArrayTaskId"}}
    fields.update(JobId="21912", JobState="COMPLETED", ExitCode="0:0", **changes)
    return " ".join(f"{k}={v}" for k,v in fields.items())


def test_exact_terminal_job_only():
    text = raw()
    assert terminal_record(text, 21912) == text + "\n"
    assert terminal_record(text, 21913) is None
    assert terminal_record(text.replace("COMPLETED", "RUNNING"), 21912) is None


@pytest.mark.parametrize("fault", ["duplicate", "array", "missing", "exit"])
def test_malformed_record(fault):
    text = raw()
    if fault == "duplicate":
        text += " JobId=21912"
    elif fault == "array":
        text += " ArrayJobId=21912"
    elif fault == "missing":
        text = text.replace("Command=1", "Other=1")
    else:
        text = text.replace("ExitCode=0:0", "ExitCode=bad")
    with pytest.raises(ValueError):
        terminal_record(text, 21912)


def test_capture_retains_all_exact_terminal_records(tmp_path, monkeypatch):
    def run(argv, **kwargs):
        job = argv[3]
        return SimpleNamespace(returncode=0, stderr="", stdout=raw().replace("JobId=21912", f"JobId={job}"))
    monkeypatch.setattr(module.subprocess, "run", run)
    result = module.capture([21912, 21913, 21914], tmp_path / "capture")
    assert result["status"] == "complete"
    assert result["polls"] == 1
    assert result["missing"] == []
    assert not result["scientific_timings_admitted"]
    for job in result["jobs"]:
        assert f"JobId={job}" in (tmp_path / "capture" / f"scheduler_{job}.txt").read_text()
    assert json.loads((tmp_path / "capture/capture.json").read_text())["status"] == "complete"
