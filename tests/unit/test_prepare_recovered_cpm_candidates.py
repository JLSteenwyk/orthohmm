import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import prepare_recovered_cpm_candidates as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("problem", [None, "runtime_before", "runtime_after", "numeric_after",
                                    "build", "changed_input", "wrong_import", "short_universe"])
def test_candidate_orchestration(tmp_path, monkeypatch, problem):
    from benchmark_tools import prepare_qfo_cpm_candidates as builder
    from benchmark_tools import checked_replay_payload_worker as corrected
    from benchmark_tools import cpm_replay_context as context
    from benchmark_tools import verify_qfo_replay_launcher as runtime
    from benchmark_tools import run_simulation_methods as frozen
    from benchmark_tools import audit_accuracy_checkpoint as numeric

    for key, value in dict(SLURM_JOB_ID="fixture", SLURM_CPUS_PER_TASK="2", SLURM_MEM_PER_NODE="65536",
            SLURM_JOB_NODELIST="bizon", PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1").items():
        monkeypatch.setenv(key, value)
    def file(relative, content="fixture"):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return record(path)
    seed = file("seed.txt", "a b\n")
    plan_record = file("plan.json")
    checkpoint = file("checkpoint/manifest.json")
    file("benchmarks/results/qfo_corrected_factorial_v1/manifest.json")
    file("benchmark_tools/results/publication_native_runtime_20260916.json")
    admission = dict(seed_partition=seed, checked_records=[seed])
    monkeypatch.setattr(module, "evidence", lambda root: admission)
    plan = dict(runtime={"fixed": True}, checkpoint_manifest=checkpoint)
    monkeypatch.setattr(corrected, "corrected_evidence", lambda *a: (plan, plan_record, None, None))
    monkeypatch.setattr(context, "evidence", lambda *a: {"checked_records": []})
    baseline = dict(runtime_before=plan["runtime"], runtime_after=plan["runtime"],
        numeric_checkpoint={"summary": {"species": 78}},
        candidate_arms={"p1_c1": {"expansion": {"parameters": {"fixed": "parameters"}}}})
    monkeypatch.setattr(frozen, "read_frozen", lambda *a: baseline)
    runtime_calls = []
    def verify(*args):
        runtime_calls.append(args)
        return {} if problem == "runtime_before" or (problem == "runtime_after" and len(runtime_calls) == 2) else plan["runtime"]
    monkeypatch.setattr(runtime, "verify", verify)
    audit_calls = []
    def audit(*args):
        audit_calls.append(args)
        return {"summary": {"species": 77 if problem == "numeric_after" and len(audit_calls) == 2 else 78}}
    monkeypatch.setattr(numeric, "audit", audit)
    launcher = tmp_path / "benchmarks/work/publication_qfo_replay_native_v1"
    engine = SimpleNamespace(__file__=str(launcher / "orthohmm/orthohmm.py"))
    accuracy = SimpleNamespace(__file__=str(launcher / "orthohmm/accuracy.py"),
        load_accuracy_checkpoint=lambda *a, **k: (range(3 if problem == "short_universe" else 984137), [], [], [], []))
    if problem == "wrong_import":
        accuracy.__file__ = str(tmp_path / "accuracy.py")
    monkeypatch.setattr(module.importlib, "import_module", lambda name: engine if name.endswith(".orthohmm") else accuracy)
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "fixture-commit")
    build_calls = []
    def build(actual_engine, parameters, actual_seed, names, species, hits, directory, auditor):
        build_calls.append(directory)
        assert actual_engine is engine
        assert parameters == {"fixed": "parameters"}
        assert actual_seed == seed
        assert len(names) == 984137
        if problem == "build":
            raise ValueError("build failed")
        directory.mkdir()
        output = directory / "partition.txt"
        output.write_text("a b\n")
        if problem == "changed_input":
            Path(seed["path"]).write_text("changed")
        return dict(candidate_arm={"output_files": [record(output)]}, content_audit={"fixture": True})
    monkeypatch.setattr(builder, "build", build)
    manifest = tmp_path / "benchmarks/results/qfo_cpm_recovered_candidates_v1/manifest.json"
    if problem:
        with pytest.raises(ValueError):
            module.prepare(tmp_path)
        if problem in ("runtime_before", "wrong_import", "short_universe"):
            assert not manifest.exists()
            assert not build_calls
        else:
            report = json.loads(manifest.read_text())
            assert report["status"] == "recovered_candidate_preparation_failed"
            assert report["accuracy_evaluated"] is False
            assert report["publication_ready"] is False
    else:
        report = module.prepare(tmp_path)
        assert json.loads(manifest.read_text()) == report
        assert report["status"] == "recovered_cpm_candidates_prepared_pending_admission"
        assert report["recovery_admission"] == admission
        assert report["accuracy_evaluated"] is False
        assert report["publication_ready"] is False
        with pytest.raises(FileExistsError):
            module.prepare(tmp_path)
