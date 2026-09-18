import copy

import pytest

from benchmark_tools import probe_cgroup_frontier as module


def fixture(root):
    for scope in ("", "system.slice", "system.slice/slurm", "system.slice/slurm/job_1",
                  "system.slice/slurm/job_2", "system.slice/daemon", "user.slice"):
        directory = root / scope
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "cpu.stat").write_text("usage_usec 1000000\n")
        (directory / "cgroup.procs").write_text("")


def points(tmp_path):
    fixture(tmp_path)
    left = module.snapshot(tmp_path, "/system.slice/slurm/job_1")
    increments = {"": 500000, "system.slice/slurm/job_1": 100000,
                  "system.slice/slurm/job_2": 200000, "system.slice/daemon": 100000,
                  "user.slice": 200000}
    for scope, value in increments.items():
        (tmp_path / scope / "cpu.stat").write_text(f"usage_usec {1000000+value}\n")
    return left, module.snapshot(tmp_path, "/system.slice/slurm/job_1")


def test_disjoint_accounting_and_negative_residual(tmp_path):
    left, right = points(tmp_path)
    result = module.compare(left, right)
    assert set(result["scope_cpu_s"]) == {"/system.slice/slurm/job_1", "/system.slice/slurm/job_2",
                                          "/system.slice/daemon", "/user.slice"}
    assert result["target_cpu_s"] == pytest.approx(.1)
    assert result["outside_target_frontier_cpu_s"] == pytest.approx(.5)
    assert result["root_minus_frontier_cpu_s"] == pytest.approx(-.1)
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("scope", ["/", "relative", "/a/../b", "/a/", "/a//b", "//a"])
def test_scope_normalization(scope):
    with pytest.raises(ValueError):
        module.validate_scope(scope)


@pytest.mark.parametrize("problem", ["boot", "inode", "counter", "overlap", "order", "duplicate"])
def test_reject_invalid_observations(tmp_path, problem):
    left, right = points(tmp_path)
    if problem == "boot":
        right["boot_id"] = "changed"
    elif problem == "inode":
        for key in ("inventory_before", "inventory_after"):
            right[key]["identities"]["/user.slice"] = [1, 99]
    elif problem == "counter":
        right["rows"][0]["raw"] = "usage_usec 0\n"
    elif problem == "overlap":
        right = copy.deepcopy(left)
    elif problem == "order":
        right["rows"][0]["started_ns"] = 0
    else:
        right["rows"].append(copy.deepcopy(right["rows"][-1]))
    with pytest.raises(ValueError):
        module.compare(left, right)


def test_new_scope_rejected_between_observations(tmp_path):
    fixture(tmp_path)
    left = module.snapshot(tmp_path, "/system.slice/slurm/job_1")
    new = tmp_path / "new.slice"
    new.mkdir()
    (new / "cpu.stat").write_text("usage_usec 0\n")
    right = module.snapshot(tmp_path, "/system.slice/slurm/job_1")
    with pytest.raises(ValueError, match="identity changed"):
        module.compare(left, right)


def test_symlink_rejected(tmp_path):
    fixture(tmp_path)
    (tmp_path / "alias").symlink_to(tmp_path / "user.slice", target_is_directory=True)
    with pytest.raises(ValueError, match="Symlink"):
        module.inventory(tmp_path, "/system.slice/slurm/job_1")


def test_ancestor_direct_processes_not_claimed_as_frontier(tmp_path):
    fixture(tmp_path)
    (tmp_path / "system.slice/cgroup.procs").write_text("123\n456\n")
    value = module.inventory(tmp_path, "/system.slice/slurm/job_1")
    assert value["ancestor_direct_process_counts"]["/system.slice"] == 2
    assert "/system.slice" not in value["identities"]
