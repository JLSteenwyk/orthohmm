import pytest

from benchmark_tools.observe_thread_affinity import observe


def fixture(tmp_path):
    root = tmp_path / "cgroup"
    scope = root / "job_1" / "step_0"
    scope.mkdir(parents=True)
    (scope / "cgroup.threads").write_text("11\n")
    child = scope / "child"
    child.mkdir()
    (child / "cgroup.threads").write_text("12\n")
    proc = tmp_path / "proc"
    for tid, membership in ((11, "/job_1/step_0"), (12, "/job_1/step_0/child")):
        path = proc / str(tid)
        path.mkdir(parents=True)
        (path / "stat").write_text(f"{tid} (odd ) name) S " + "0 " * 18 + "123 0\n")
        (path / "cgroup").write_text(f"0::{membership}\n")
    return root, scope, proc


def test_descendant_threads_and_subset_affinity(tmp_path):
    root, scope, proc = fixture(tmp_path)
    row = observe(scope, range(32), cgroup_root=root, proc_root=proc,
                  get_affinity=lambda tid: {tid - 11})
    assert row["status"] == "observed_within_affinity"
    assert [t["start_ticks"] for t in row["threads"]] == [123, 123]
    assert row["initial_tids"] == [11, 12]
    assert row["full_run_affinity_verified"] is False


@pytest.mark.parametrize("case", ["escape", "exit", "migration", "reuse", "new_thread", "denied"])
def test_retains_violations_and_gaps(tmp_path, case):
    root, scope, proc = fixture(tmp_path)
    def affinity(tid):
        if tid == 12:
            if case == "escape":
                return {0, 96}
            if case == "exit":
                raise ProcessLookupError("exited")
            if case == "denied":
                raise PermissionError("denied")
            if case == "migration":
                (proc / "12" / "cgroup").write_text("0::/job_1/step_01\n")
            if case == "reuse":
                path = proc / "12" / "stat"
                path.write_text(path.read_text().replace("123", "124"))
            if case == "new_thread":
                (scope / "cgroup.threads").write_text("11\n13\n")
        return {0}
    row = observe(scope, range(32), cgroup_root=root, proc_root=proc, get_affinity=affinity)
    assert row["status"] == ("violation" if case == "escape" else "incomplete")
    if case == "escape":
        assert row["violating_tids"] == [12]
    else:
        assert row["errors"]


def test_empty_inventory_is_not_compliance(tmp_path):
    root, scope, proc = fixture(tmp_path)
    for path in scope.rglob("cgroup.threads"):
        path.write_text("")
    assert observe(scope, [0], cgroup_root=root, proc_root=proc)["status"] == "incomplete"


def test_reject_whole_host(tmp_path):
    root, _, _ = fixture(tmp_path)
    with pytest.raises(ValueError, match="whole-host"):
        observe(root, [0], cgroup_root=root)
