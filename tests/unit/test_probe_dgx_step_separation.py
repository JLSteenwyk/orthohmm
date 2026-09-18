import json

import pytest

from benchmark_tools import probe_dgx_step_separation as module


def sample(scope):
    return {"raw": {"cgroup_membership": "0::" + scope + "\n"}}


def test_separate_step_scopes():
    parent = sample("/system.slice/slurm/job_123/step_batch/user/task_0")
    native = sample("/system.slice/slurm/job_123/step_0/user/task_0")
    assert module.validate_scopes(parent, native, 123)["native_step"].endswith("/step_0")


@pytest.mark.parametrize("scope", ["/job_124/step_0", "/job_123/step_batch", "/job_123/step_extern",
                                    "/job_123/user", "/job_123", "/../job_123/step_0"])
def test_invalid_scopes(scope):
    with pytest.raises(ValueError):
        module.validate_scopes(sample("/job_123/step_batch"), sample(scope), 123)


def test_atomic_publication_and_no_overwrite(tmp_path):
    path = tmp_path / "ready.json"
    module.save(path, {"ready": True})
    assert module.wait_file(path) == {"ready": True}
    assert not path.with_name("ready.json.pending").exists()
    with pytest.raises(FileExistsError):
        module.save(path, {"changed": True})
    assert json.loads(path.read_text()) == {"ready": True}


def test_wait_timeout(tmp_path):
    with pytest.raises(TimeoutError):
        module.wait_file(tmp_path / "missing", seconds=.01)


def test_unscheduled_rejected(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(KeyError):
        module.run(tmp_path / "out")
    assert not (tmp_path / "out").exists()
