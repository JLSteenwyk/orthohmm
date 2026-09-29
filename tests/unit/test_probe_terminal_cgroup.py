from pathlib import Path

import pytest

from benchmark_tools.probe_terminal_cgroup import job_scope


def test_exact_job_ancestor():
    assert job_scope("0::/system/job_123/step_batch/user/task_0\n", 123) == Path("/sys/fs/cgroup/system/job_123")


@pytest.mark.parametrize("membership", ["0::/job_124/step_batch", "0::/job_123/../step_batch", "0::relative/job_123/step_batch", "", "0::/job_123/step_batch\n0::/other"])
def test_invalid_scope_rejected(membership):
    with pytest.raises(ValueError):
        job_scope(membership, 123)
