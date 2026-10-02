import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools import replay_lineage_read_crossing_control as module
from benchmark_tools.replay_lineage_read_crossing_control import replay, replay_trial

REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "benchmark_tools/results"
HISTORICAL_COMMAND = (
    "/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/"
    "benchmark_tools/run_dgx_lineage_read_crossing.sh"
)

CHILD = r'''
import json, os, pathlib, sys
repo = pathlib.Path(sys.argv[1]).resolve()
original = pathlib.Path(sys.argv[4]).resolve()
blocked = []
def audit(event, args):
    if event == "open" and isinstance(args[0], (str, bytes, os.PathLike)):
        path = pathlib.Path(os.path.abspath(os.fsdecode(args[0])))
        if path == original or original in path.parents:
            blocked.append(str(path))
            raise PermissionError("Original checkout Python reads forbidden")
    if event == "subprocess.Popen":
        argv = args[1]
        if (len(argv) != 3 or argv[:2] != ["git", "show"]
                or not argv[2].startswith("6599c6e:benchmark_tools/")):
            raise PermissionError("Only frozen-source Git reads allowed")
    if event in ("os.system", "os.exec", "os.posix_spawn"):
        raise PermissionError("Other child execution forbidden")
sys.addaudithook(audit)
try:
    (original / "tests/conftest.py").read_bytes()
except PermissionError:
    pass
assert len(blocked) == 1
blocked.clear()
sys.path.insert(0, str(repo))
from benchmark_tools import replay_lineage_read_crossing_control as module
result = module.replay(pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3]), repo)
origins = {name: str(pathlib.Path(value.__file__).resolve())
           for name, value in sys.modules.items()
           if getattr(value, "__file__", None)
           and (name.startswith("benchmark_tools") or name == "probe_host_counters")}
assert origins and all(pathlib.Path(p).is_relative_to(repo / "benchmark_tools")
                       for p in origins.values()), origins
assert not blocked, blocked
result["test_import_evidence"] = dict(origins=origins, canary_blocked=True,
    later_original_python_open_events=len(blocked), only_git_show_children=True,
    python=sys.version, git_objects_shared_with_original=True)
print(json.dumps(result, sort_keys=True, allow_nan=False))
'''


@pytest.fixture
def frozen_repo(tmp_path):
    target = tmp_path / "frozen_repo"
    subprocess.run(["git", "clone", "--shared", "--no-checkout", str(REPO), str(target)],
                   check=True, capture_output=True, timeout=30)
    names = module.SOURCES | {"__init__.py", "prepare_ob_candidate_neighborhood.py",
                              "run_dgx_lineage_read_crossing.sh"}
    for name in sorted(names):
        data = subprocess.check_output(["git", "show", module.COMMIT + ":benchmark_tools/" + name],
                                       cwd=target, timeout=10)
        path = target / "benchmark_tools" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    shutil.copyfile(module.__file__, target / "benchmark_tools/replay_lineage_read_crossing_control.py")
    return target


@pytest.fixture
def scheduler(tmp_path, frozen_repo):
    text = (RESULTS / "lineage_read_crossing_scheduler_22018.txt").read_text()
    old = "   Command=" + HISTORICAL_COMMAND + "\n"
    assert text.count(old) == 1
    target = tmp_path / "scheduler.txt"
    target.write_text(text.replace(old, "   Command=" + str(
        frozen_repo / "benchmark_tools/run_dgx_lineage_read_crossing.sh") + "\n", 1))
    return target


def child_replay(archive, scheduler, frozen_repo, *, python=sys.executable):
    return subprocess.run([python, "-I", "-B", "-c", CHILD, str(frozen_repo),
                           str(archive), str(scheduler), str(REPO)], cwd=frozen_repo,
                          capture_output=True, text=True, timeout=30)


@pytest.fixture
def archive(tmp_path):
    target = tmp_path / "controls"
    shutil.copytree(RESULTS / "lineage_read_crossing_22018", target)
    return target


def mutate(path, change):
    data = json.loads(path.read_text())
    change(data)
    path.write_text(json.dumps(data))


@pytest.mark.parametrize("index", [0, 1, 2])
def test_real_trial_raw_replay(index, archive):
    identity = json.loads((archive / "identity.json").read_text())
    result = replay_trial(archive / f"trial_{index}", index, identity)
    assert result["result"]["control_met"] is True
    assert result["result"]["spans"]["before_to_crossing"]["root_minus_target_cpu_usec"] < 0
    assert result["result"]["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["command", "service", "removal", "nested_scope", "nested_result",
                                   "event", "outer_scope", "result"])
def test_tampered_trial_rejected(archive, fault):
    trial = archive / "trial_0"
    cases = {
        "command": ("lifecycle/command.json", lambda d: d["command"].append("unexpected")),
        "service": ("lifecycle/service.json", lambda d: d.update(returncode=1)),
        "removal": ("lifecycle/removal.json", lambda d: d["observations"][-1].update(exists=True)),
        "nested_scope": ("lifecycle/before.json", lambda d: d["target"].update(target="/other")),
        "nested_result": ("lifecycle_result.json", lambda d: d["result"].update(manager_cpu_s=0)),
        "event": ("event.json", lambda d: d.update(finished_ns=d["started_ns"])),
        "outer_scope": ("before.json", lambda d: d.update(target="/other")),
        "result": ("result.json", lambda d: d.update(control_met=False)),
    }
    filename, change = cases[fault]
    mutate(trial / filename, change)
    identity = json.loads((archive / "identity.json").read_text())
    with pytest.raises(ValueError):
        replay_trial(trial, 0, identity)


def test_complete_raw_replay_with_pinned_sources(archive, scheduler, frozen_repo):
    child = child_replay(archive, scheduler, frozen_repo)
    assert child.returncode == 0, child.stderr
    result = json.loads(child.stdout)
    assert len(result["trials"]) == 3
    assert [row["result"] for row in result["trials"]] == json.loads(
        (archive / "report.json").read_text())["trials"]
    assert result["all_controls_met"] is True
    assert result["environmental_validity_established"] is False
    assert result["scientific_timings_admitted"] is False
    evidence = result["test_import_evidence"]
    assert evidence["canary_blocked"] and evidence["only_git_show_children"]
    assert evidence["later_original_python_open_events"] == 0
    assert {Path(p).name for p in evidence["origins"].values()} >= {
        Path(name).name for name in module.SOURCES if name.endswith(".py")}


def test_missing_trial_evidence_cannot_be_skipped(archive):
    (archive / "trial_2/result.json").unlink()
    with pytest.raises(ValueError, match="35-file"):
        replay(archive, RESULTS / "lineage_read_crossing_scheduler_22018.txt", REPO)


def test_wrong_source_hash_rejected(archive, scheduler, frozen_repo):
    mutate(archive / "identity.json", lambda d: d["sources"].update(
        {next(iter(d["sources"])): "0" * 64}))
    child = child_replay(archive, scheduler, frozen_repo)
    assert child.returncode != 0
    assert "Deployed source differs from frozen commit" in child.stderr


def test_wrong_scheduler_allocation_rejected(archive, scheduler, frozen_repo):
    scheduler.write_text(scheduler.read_text().replace(
        "NumCPUs=20", "NumCPUs=2"))
    child = child_replay(archive, scheduler, frozen_repo)
    assert child.returncode != 0
    assert "Terminal scheduler identity or allocation differs" in child.stderr


def test_wrong_scheduler_command_rejected(archive, scheduler, frozen_repo):
    scheduler.write_text(scheduler.read_text().replace(
        "run_dgx_lineage_read_crossing.sh", "other.sh"))
    child = child_replay(archive, scheduler, frozen_repo)
    assert child.returncode != 0
    assert "Terminal scheduler identity or allocation differs" in child.stderr


@pytest.mark.parametrize("fault", ["same_size", "truncated", "missing"])
def test_changed_local_frozen_source_rejected(archive, scheduler, frozen_repo, fault):
    path = frozen_repo / "benchmark_tools/results/LINEAGE_READ_CROSSING_PROTOCOL_20260919.md"
    if fault == "same_size":
        path.write_bytes(b"X" + path.read_bytes()[1:])
    elif fault == "truncated":
        path.write_bytes(path.read_bytes()[:-1])
    else:
        path.unlink()
    child = child_replay(archive, scheduler, frozen_repo)
    assert child.returncode != 0
    assert ("Local replay dependency differs from frozen source" if fault != "missing"
            else "FileNotFoundError") in child.stderr


@pytest.mark.parametrize("fault", ["same_size", "current_helper"])
def test_changed_imported_helper_rejected(archive, scheduler, frozen_repo, fault):
    path = frozen_repo / "benchmark_tools/probe_cgroup_frontier.py"
    original = path.read_bytes()
    if fault == "same_size":
        assert original.count(b"Read a disjoint") == 1
        changed = original.replace(b"Read a disjoint", b"Reed a disjoint", 1)
        assert len(changed) == len(original)
    else:
        changed = (REPO / "benchmark_tools/probe_cgroup_frontier.py").read_bytes()
    assert changed != original
    path.write_bytes(changed)
    child = child_replay(archive, scheduler, frozen_repo)
    assert child.returncode != 0
    assert "Local replay dependency differs from frozen source" in child.stderr
