import pytest

from benchmark_tools.probe_threadripper_allocation import validate


def sample():
    return dict(host="bizon", affinity=list(range(32)),
                topology=[dict(cpu=i, package=0, core=i) for i in range(32)],
                cgroup="0::/slurm/job_1/step_0\n",
                slurm=dict(SLURM_JOB_ID="1", SLURM_CPUS_PER_TASK="32", SLURM_MEM_PER_NODE="131072"),
                ancestors=[{"memory.max": "max"}, {"memory.max": str(128 * 1024**3)}])


def test_valid():
    validate(sample(), sample())


@pytest.mark.parametrize("change", ["siblings", "ram", "unlimited", "cores", "host", "child", "request"])
def test_rejects_wrong_allocation(change):
    parent, child = sample(), sample()
    if change == "siblings":
        child["topology"][1]["core"] = 0
    elif change == "ram":
        child["ancestors"].append({"memory.max": str(64 * 1024**3)})
    elif change == "unlimited":
        child["ancestors"] = [{"memory.max": "max"}]
    elif change == "cores":
        child["affinity"] = list(range(192))
    elif change == "host":
        child["host"] = "other"
    elif change == "child":
        child["cgroup"] = "0::/other\n"
    else:
        child["slurm"]["SLURM_CPUS_PER_TASK"] = "20"
    with pytest.raises(ValueError):
        validate(parent, child)
