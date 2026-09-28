import json
import os

import pytest

from benchmark_tools import measure_threadripper_run as module


@pytest.mark.parametrize("failure", [None, "before", "prepare", "native", "after", "nonzero", "timeout"])
def test_composed_boundary_and_failure_retention(tmp_path, monkeypatch, failure):
    root, target = tmp_path / "run", tmp_path / "inputs"
    run = dict(measurement_directory=str(root / "measurement"), native_argv=["/native"],
               dataset=dict(inputs=[]))
    events = []
    monkeypatch.setattr(module, "paths", lambda r: (root, target))
    monkeypatch.setattr(module, "assert_environment", lambda *args: None)
    def verify(*args):
        events.append("baseline")
        assert "run_verification_caches" in os.environ["NUMBA_CACHE_DIR"]
    monkeypatch.setattr(module, "verify_environment", verify)
    calls = 0
    def checker(specs):
        nonlocal calls
        calls += 1
        events.append("before" if calls == 1 else "after")
        if failure == ("before" if calls == 1 else "after"):
            raise ValueError("runtime changed")
        return {"verified": True}
    def prepare(*args):
        events.append("prepare")
        assert root.is_dir() and not (root / "measurement").exists()
        if failure == "prepare":
            raise ValueError("copy failed")
        target.mkdir()
        return {"inference_started": False}
    monkeypatch.setattr(module, "prepare", prepare)
    monkeypatch.setattr(module, "check_prepared", lambda *args: events.append("input_check"))
    def collector(command, directory, job, cpus, memory, timeout, cadence, **kwargs):
        events.append("native")
        assert command == ["/native"] and job == 1
        assert (cpus, memory, timeout, cadence) == (32, 128*1024**3, 85800, 1.)
        assert os.environ["NUMBA_CACHE_DIR"] == str(root / "native_numba_cache")
        assert all(os.environ[key] == str(root / "native_tmp") for key in ("TMPDIR", "TMP", "TEMP"))
        assert not list((root / "native_numba_cache").iterdir())
        if failure == "native":
            raise RuntimeError("collector failed")
        if failure in ("nonzero", "timeout"):
            return dict(status="command_failed", native=dict(exit_code=124 if failure == "timeout" else 7,
                                                             timed_out=failure == "timeout"))
        return {"status": "command_exited_zero"}
    result = module.measure_run(run, {}, [("manifest", "hash")], collector, 1, runtime_checker=checker)
    assert json.loads((root / "verification.json").read_text()) == result
    assert result["scientific_results_admitted"] is False
    if failure is None:
        assert events == ["before", "baseline", "prepare", "input_check", "native", "after", "baseline", "input_check"]
        assert result["status"] == "command_exited_zero"
    elif failure in ("nonzero", "timeout"):
        assert result["status"] == "command_failed"
        assert result["measurement"]["native"]["timed_out"] == (failure == "timeout")
        assert events.count("native") == 1 and "after" in events
    elif failure == "after":
        assert result["status"] == "runtime_changed_or_unverifiable"
        assert result["measurement"]["status"] == "command_exited_zero"
    else:
        assert result["status"] == "verified_wrapper_failed"
        if failure == "before":
            assert "prepare" not in events and "native" not in events
        else:
            assert "after" in events


def test_cwd_rejected_before_environment_lookup(tmp_path):
    with pytest.raises(ValueError, match="cwd"):
        module.assert_environment({"cwd": str(tmp_path / "other")}, {})


def test_environment_rejected(monkeypatch):
    monkeypatch.setenv("OMP_NUM_THREADS", "8")
    with pytest.raises(ValueError, match="environment"):
        module.assert_environment({"cwd": str(module.Path.cwd())},
                                  {"environment_overrides": {"OMP_NUM_THREADS": "1"}})
