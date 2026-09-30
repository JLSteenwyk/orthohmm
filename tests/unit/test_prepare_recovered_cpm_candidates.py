import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import textwrap
from types import SimpleNamespace

import pytest

from benchmark_tools import prepare_recovered_cpm_candidates as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("problem", [None, "runtime_before", "runtime_after", "numeric_after",
                                    "build", "changed_input", "wrong_import", "short_universe"])
@pytest.mark.parametrize("amended", [False, True])
def test_candidate_orchestration(tmp_path, monkeypatch, problem, amended):
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
    options = {}
    if amended:
        from benchmark_tools import helper_recovered_cpm_seed_evidence as handoff
        def forbidden(root):
            raise AssertionError("Explicit amendment cannot use historical gate")
        monkeypatch.setattr(module, "evidence", forbidden)
        readback = tmp_path / "readback.json"
        def selected(root, path, digest, protocol):
            assert (root, path, digest, protocol) == (tmp_path, readback, "fixed-readback", "fixed-protocol")
            return admission
        monkeypatch.setattr(handoff, "evidence", selected)
        options = dict(recovery_readback=readback, readback_sha="fixed-readback", protocol_sha="fixed-protocol")
    else:
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
    output = "qfo_cpm_helper_recovered_candidates_v1" if amended else "qfo_cpm_recovered_candidates_v1"
    manifest = tmp_path / "benchmarks/results" / output / "manifest.json"
    if problem:
        with pytest.raises(ValueError):
            module.prepare(tmp_path, **options)
        if problem in ("runtime_before", "wrong_import", "short_universe"):
            assert not manifest.exists()
            assert not build_calls
        else:
            report = json.loads(manifest.read_text())
            assert report["status"] == "recovered_candidate_preparation_failed"
            assert report["accuracy_evaluated"] is False
            assert report["publication_ready"] is False
    else:
        report = module.prepare(tmp_path, **options)
        assert json.loads(manifest.read_text()) == report
        assert report["status"] == "recovered_cpm_candidates_prepared_pending_admission"
        assert report["recovery_admission"] == admission
        assert report["accuracy_evaluated"] is False
        assert report["publication_ready"] is False
        assert report["seed_handoff"] == ("explicit_helper_runtime_seed_amendment" if amended else "historical_admission_22155")
        with pytest.raises(FileExistsError):
            module.prepare(tmp_path, **options)


@pytest.mark.parametrize("mask", range(1, 7))
def test_partial_amendment_rejected_before_environment_or_output(tmp_path, monkeypatch, mask):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    fields = ("recovery_readback", "readback_sha", "protocol_sha")
    options = {field: "fixture" for index, field in enumerate(fields) if mask & (1 << index)}
    with pytest.raises(ValueError, match="together"):
        module.prepare(tmp_path, **options)
    assert not (tmp_path / "benchmarks").exists()


@pytest.mark.parametrize("preloaded_wrong_package", [False, True])
def test_real_transitive_status_import_respects_selected_scientific_package(tmp_path, preloaded_wrong_package):
    root = Path(module.__file__).resolve().parent.parent
    launcher = tmp_path / "benchmarks/work/publication_qfo_replay_native_v1/orthohmm"
    launcher.parent.mkdir(parents=True)
    shutil.copytree(root / "orthohmm", launcher,
                    ignore=shutil.ignore_patterns("__pycache__", "*.so", "*.nbc", "*.nbi"))
    for path in (tmp_path / "seed.txt", tmp_path / "plan.json", tmp_path / "checkpoint/manifest.json",
                 tmp_path / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json",
                 tmp_path / "benchmark_tools/results/publication_native_runtime_20260916.json"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture")
    script = textwrap.dedent('''
        import json, sys
        from pathlib import Path
        from unittest.mock import patch
        source, root, contaminated = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3] == "yes"
        sys.path.insert(0, str(source))
        from benchmark_tools import prepare_recovered_cpm_candidates as driver
        from benchmark_tools import prepare_qfo_cpm_candidates as builder
        from benchmark_tools import checked_replay_payload_worker as corrected
        from benchmark_tools import cpm_replay_context as context
        from benchmark_tools import verify_qfo_replay_launcher as runtime
        from benchmark_tools import run_simulation_methods as frozen
        from benchmark_tools.prepare_ob_candidate_neighborhood import record
        if contaminated:
            import orthohmm
        seed, plan_record = record(root / "seed.txt"), record(root / "plan.json")
        plan = dict(runtime={"fixed": True}, checkpoint_manifest=record(root / "checkpoint/manifest.json"))
        baseline = dict(runtime_before=plan["runtime"], runtime_after=plan["runtime"],
            numeric_checkpoint={"summary": {"species": 78}},
            candidate_arms={"p1_c1": {"expansion": {"parameters": {"fixed": "parameters"}}}})
        real_import = driver.importlib.import_module
        def scientific(name):
            imported = real_import(name)
            if name == "orthohmm.accuracy":
                imported.load_accuracy_checkpoint = lambda *a, **k: (range(984137), [], [], [], [])
                from benchmark_tools import audit_accuracy_checkpoint as numeric
                numeric.audit = lambda *a, **k: {"summary": {"species": 78}}
            return imported
        with patch.object(driver, "evidence", return_value=dict(seed_partition=seed, checked_records=[seed])), \
             patch.object(corrected, "corrected_evidence", return_value=(plan, plan_record, None, None)), \
             patch.object(context, "evidence", return_value={"checked_records": []}), \
             patch.object(runtime, "verify", return_value=plan["runtime"]), \
             patch.object(frozen, "read_frozen", return_value=baseline), \
             patch.object(builder, "build", return_value={"candidate_arm": {"output_files": []}}), \
             patch.object(driver.importlib, "import_module", side_effect=scientific):
            if contaminated:
                try:
                    driver.prepare(root)
                except ValueError as error:
                    assert str(error) == "Wrong frozen candidate scientific import"
                else:
                    raise AssertionError("Preloaded wrong package must remain rejected")
                assert not (root / "benchmarks/results/qfo_cpm_recovered_candidates_v1").exists()
            else:
                report = driver.prepare(root)
                assert report["status"] == "recovered_cpm_candidates_prepared_pending_admission"
                expected = root / "benchmarks/work/publication_qfo_replay_native_v1/orthohmm"
                assert Path(sys.modules["orthohmm.orthohmm"].__file__).resolve().parent == expected
                assert Path(sys.modules["orthohmm.phylogeny_pipeline"].__file__).resolve().parent == expected
                assert report["accuracy_evaluated"] is False
                assert report["publication_ready"] is False
        print(json.dumps(dict(status="real_import_order_regression_passed", contaminated=contaminated)))
    ''')
    environment = dict(os.environ, SLURM_JOB_ID="fixture", SLURM_CPUS_PER_TASK="2", SLURM_MEM_PER_NODE="65536",
        SLURM_JOB_NODELIST="bizon", PYTHONHASHSEED="0", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    done = subprocess.run([sys.executable, "-I", "-B", "-c", script, str(root), str(tmp_path),
                           "yes" if preloaded_wrong_package else "no"],
                          cwd=tmp_path, env=environment, capture_output=True, text=True, timeout=60)
    assert done.returncode == 0, done.stdout + done.stderr
    assert json.loads(done.stdout)["status"] == "real_import_order_regression_passed"
