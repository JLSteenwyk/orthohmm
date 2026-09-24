import copy
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from benchmark_tools import admit_cpm_checkpoint_recovery as module
from benchmark_tools import verify_qfo_replay_launcher as runtime_module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.mark.parametrize("problem", [None, "preflight", "fresh_runtime", "optimizer", "profile_copy", "reference",
    "checkpoint", "comparison", "child_failure", "child_metadata", "child_partition", "mutation", "runtime"])
def test_admission_flow(tmp_path, monkeypatch, problem):
    # Test orchestration and real hashes/partition comparisons, with scientific
    # execution mocked. Separate tests exercise actual frozen native functions.
    root = tmp_path
    executor = root / "benchmarks/work/cpm_checkpoint_recovery_v1_20260923"
    source = executor / "benchmark_tools/run_cpm_checkpoint_recovery.py"
    source.parent.mkdir(parents=True)
    source.write_text("import json\ndef optimizer_evidence(root, directory, preflight, environment):\n"
                      "    return json.loads((directory / 'fixture_optimizer.json').read_text())\n")
    monkeypatch.setattr(module, "RUNNER_SHA", record(source)["sha256"])
    directory = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
    directory.mkdir(parents=True)
    (directory / "payload").mkdir()
    (directory / "payload/gene_names.txt").write_text("a\nb\n")
    stages = []
    for name in ("multipass", "multipass_refined", "profiles", "profiles_refined"):
        path = directory / f"orthogroups_{name}.txt"
        path.write_text("a b\n")
        stages.append(dict(output=record(path)))
    (directory / "refinement_repeat.txt").write_text("a b\n")
    optimizer = dict(partition=record(directory / "orthogroups_profiles.txt"), checked_records=[])
    write(directory / "fixture_optimizer.json", optimizer if problem != "optimizer" else {})
    reference_path = root / "benchmarks/work/qfo_cpm_refinement_check_22153/worker.json"
    write(reference_path, dict(numeric_checkpoint={"verified": True}, modules=[]))
    reference = record(reference_path)
    child = dict(genes=984137, groups=1, refinement_directed_hits=0, numeric_checkpoint={"verified": True},
                 modules=[], accuracy_evaluated=False)
    if problem == "checkpoint":
        child["numeric_checkpoint"] = {}
    for label, filename in (("refinement", "orthogroups_profiles_refined.txt"), ("refinement_repeat", "refinement_repeat.txt")):
        write(directory / f"{label}.json", {**child, "output": record(directory / filename)})
    runtime = {"verified": True}
    preflight = dict(status="cpm_checkpoint_preflight_verified_unscored", preflight_passed=True,
        checked_records=[reference], context={"cwd": str(root)}, runtime=runtime, saved_graph={}, saved_payload="saved",
        original_scheduler={"State": "FAILED"})
    preflight_path = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_preflight_v1/status.json"
    write(preflight_path, preflight)
    parent = dict(preflight=record(preflight_path), checked_records=[], phases=[], runtime_after=runtime, optimizer=optimizer,
        refinement_reports=[record(directory / name) for name in ("refinement.json", "refinement_repeat.json")],
        refinement_comparison={"genes": 2, "groups": 1, "partition_equal": True})
    if problem == "profile_copy":
        parent["optimizer"]["partition"] = {**optimizer["partition"], "sha256": "changed"}
        write(directory / "fixture_optimizer.json", parent["optimizer"])
    elif problem == "comparison":
        parent["refinement_comparison"]["partition_equal"] = False
    write(directory / "status.json", parent)
    monkeypatch.setattr(module, "parent_contract", lambda *args: stages)
    monkeypatch.setattr(module.subprocess, "check_output", lambda command, **kwargs:
        "JobID|State|ExitCode|AllocCPUS|ReqMem|NodeList\n123|COMPLETED|0:0|1|64G|bizon\n"
        if command[0] == "sacct" else "commit")
    def prepare(root, output):
        if problem == "preflight":
            raise ValueError("fresh preflight failed")
        result = copy.deepcopy(preflight)
        if problem == "fresh_runtime":
            result["runtime"] = {}
        elif problem == "reference":
            result["checked_records"] = []
        write(output / "status.json", result)
        return result
    monkeypatch.setattr(module, "prepare", prepare)
    calls = []
    def run(command, **kwargs):
        if command[0] == "git":
            return SimpleNamespace(returncode=0)
        calls.append(command)
        if problem == "child_failure":
            raise subprocess.CalledProcessError(1, command)
        output = Path(command[command.index("--output") + 1])
        path = output / "refinement_repeat.txt"
        path.write_text("a\nb\n" if problem == "child_partition" else "a b\n")
        result = {**child, "output": record(path)}
        if problem == "child_metadata":
            result["modules"] = ["wrong source"]
        write(output / "refinement_repeat.json", result)
        if problem == "mutation":
            (directory / "orthogroups_profiles_refined.txt").write_text("changed")
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(module.subprocess, "run", run)
    monkeypatch.setattr(runtime_module, "verify", lambda *args: {} if problem == "runtime" else runtime)
    output = root / "admission"
    if problem:
        with pytest.raises((ValueError, subprocess.CalledProcessError)):
            module.admit(root, "123", "commit", output)
        report = json.loads((output / "status.json").read_text())
        assert report["status"] == "recovery_admission_failed"
        assert report["seed_admitted"] is False
        assert report["accuracy_evaluated"] is False
    else:
        report = module.admit(root, "123", "commit", output)
        assert report["status"] == "cpm_checkpoint_recovered_seed_admitted_unscored"
        assert report["seed_partition"] == stages[-1]["output"]
        assert report["seed_admitted"] is True
        assert report["accuracy_evaluated"] is False
        assert len(calls) == 1 and calls[0][-1] == "repeat-refinement"
