import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import run_cpm_checkpoint_recovery as module
from benchmark_tools import prepare_cpm_checkpoint_recovery as preflight_module
from benchmark_tools import verify_qfo_replay_launcher as runtime_module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.mark.parametrize("problem", [None, "allocation", "revision", "preflight", "optimizer", "observation",
    "refine", "repeat-refinement", "mismatch", "checkpoint", "module", "runtime", "changed_input", "existing"])
def test_orchestration(tmp_path, monkeypatch, problem):
    # Exercise parent handoffs and real file hashes. Scientific subprocesses are
    # mocked here; frozen native execution is tested separately.
    root = tmp_path
    output = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
    for name, value in dict(SLURM_CPUS_PER_TASK="1", SLURM_MEM_PER_NODE="65536",
                            SLURMD_NODENAME="bizon", SLURM_JOB_ID="fixture").items():
        monkeypatch.setenv(name, value)
    if problem == "allocation":
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "wrong" if problem == "revision" else "commit")
    original = root / "original"
    (original / "replay").mkdir(parents=True)
    for name in ("orthogroups_multipass.txt", "orthogroups_multipass_refined.txt"):
        (original / "replay" / name).write_text("a b\n")
    saved = root / "saved"
    saved.mkdir()
    for name in ("sources.npy", "targets.npy", "weights.npy"):
        (saved / name).write_text("fixture")
    (saved / "gene_names.txt").write_text("a\nb\n")
    reference_path = root / "benchmarks/work/qfo_cpm_refinement_check_22153/worker.json"
    write(reference_path, dict(numeric_checkpoint={"verified": True}, modules=[]))
    reference = record(reference_path)
    expected_runtime = {"verified": True}
    def prepare(root, destination):
        if problem == "preflight":
            raise ValueError("preflight failed")
        destination.mkdir(parents=True)
        result = dict(status="cpm_checkpoint_preflight_verified_unscored", preflight_passed=True,
            checked_records=[reference, record(saved / "weights.npy")], saved_payload=str(saved), saved_graph={},
            context={"cwd": str(root), "output_root": str(original)}, runtime=expected_runtime)
        write(destination / "status.json", result)
        return result
    monkeypatch.setattr(preflight_module, "prepare", prepare)
    calls = []
    def run(command, **kwargs):
        if command[0] == "git":
            return SimpleNamespace(returncode=0)
        mode = command[-1]
        calls.append(mode)
        if problem == mode or (problem == "optimizer" and mode == "optimize"):
            return SimpleNamespace(returncode=1)
        if mode == "optimize":
            partition = output / "orthohmm_working_res/orthohmm_edges_clustered.txt"
            partition.parent.mkdir()
            partition.write_text("a b\n")
        else:
            filename = "orthogroups_profiles_refined.txt" if mode == "refine" else "refinement_repeat.txt"
            path = output / filename
            path.write_text("a\nb\n" if problem == "mismatch" and mode == "repeat-refinement" else "a b\n")
            result = dict(genes=984137, groups=1, refinement_directed_hits=0, accuracy_evaluated=False,
                output=record(path), numeric_checkpoint={"verified": True}, modules=[])
            if problem == "checkpoint":
                result["numeric_checkpoint"] = {}
            elif problem == "module":
                result["modules"] = [{"unexpected": True}]
            write(output / ("refinement.json" if mode == "refine" else "refinement_repeat.json"), result)
        if problem == "changed_input":
            (saved / "weights.npy").write_text("changed")
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(module.subprocess, "run", run)
    def optimizer(*args):
        if problem == "observation":
            raise ValueError("missing graph evidence")
        return dict(partition=record(output / "orthohmm_working_res/orthohmm_edges_clustered.txt"), checked_records=[])
    monkeypatch.setattr(module, "optimizer_evidence", optimizer)
    monkeypatch.setattr(runtime_module, "verify", lambda *args: {} if problem == "runtime" else expected_runtime)
    if problem == "existing":
        output.mkdir(parents=True)
    if problem:
        with pytest.raises((ValueError, RuntimeError, FileExistsError)):
            module.run(root, "commit")
        if (output / "status.json").exists():
            report = json.loads((output / "status.json").read_text())
            assert report["status"] == "checkpoint_recovery_failed"
            assert report["downstream_admitted"] is False
            assert report["accuracy_evaluated"] is False
        assert calls.count("optimize") <= 1
        if problem in ("allocation", "revision", "preflight", "existing"):
            assert not calls
        elif problem in ("optimizer", "observation"):
            assert calls == ["optimize"]
    else:
        report = module.run(root, "commit")
        assert calls == ["optimize", "refine", "repeat-refinement"]
        assert report["status"] == "cpm_checkpoint_recovered_pending_independent_admission"
        assert [s["origin"] for s in report["stages"]] == ["reused", "reused", "recovered", "recovered"]
        assert report["refinement_comparison"]["partition_equal"] is True
        assert report["downstream_admitted"] is False


def test_isolated_cli_help(tmp_path):
    done = subprocess.run([sys.executable, "-I", str(Path(module.__file__).resolve()), "--help"],
                          cwd=tmp_path, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert "--mode" in done.stdout
