from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import recover_cap400_scores as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_membership_hillclimb import write_json


def admitted():
    native = {"genes": 984137, "species": 78, "groups": 370680, "pairs": 5836092,
              "phylogeny": "/original/native", "products": {}, "summary": {}, "checked": []}
    evidence = {"status": runner.READY, "scores": {}, "native_inference_rerun": False,
        "producer": {"sha256": "original"}, "scheduler": {"State": "FAILED"},
        "scientific_revision": "original", "native": native, "checked": []}
    admission = {k: deepcopy(evidence[k]) for k in ("status", "scores", "native_inference_rerun",
        "producer", "scheduler", "scientific_revision")}
    admission["coverage"] = {k: native[k] for k in ("genes", "species", "groups", "pairs")}
    return admission, evidence


def test_exact_unscored_native_admission_required():
    runner.require_admission(*admitted())


@pytest.mark.parametrize("index,field,value", [(0, "status", "running"), (1, "status", "running"),
    (0, "scores", {"GO": 1}), (1, "scores", {"GO": 1}),
    (0, "native_inference_rerun", True), (1, "native_inference_rerun", True),
    (0, "producer", {}), (0, "scheduler", {}), (0, "scientific_revision", "changed"),
    (0, "coverage", {})])
def test_changed_admission_refused(index, field, value):
    pair = list(admitted())
    pair[index][field] = value
    with pytest.raises(ValueError, match="exact unscored admission"):
        runner.require_admission(*pair)


@pytest.mark.parametrize("field,value", [("genes", 984136), ("species", 77)])
def test_incomplete_native_universe_refused_even_when_claims_agree(field, value):
    admission, evidence = admitted()
    evidence["native"][field] = admission["coverage"][field] = value
    with pytest.raises(ValueError, match="exact unscored admission"):
        runner.require_admission(admission, evidence)


def test_real_converter_and_helper_import_are_bound_to_this_checkout():
    source = Path(runner.__file__).resolve().parent.parent
    assert runner.verify_imported_converter(source) == record(source / "qfo_benchmark/og_to_pairwise.py")


@pytest.mark.parametrize("changed", ["converter_path", "helper_path", "binding"])
def test_different_imported_converter_or_binding_refused(tmp_path, monkeypatch, changed):
    function = lambda value: value
    converter = SimpleNamespace(__file__=str(tmp_path / "qfo_benchmark/og_to_pairwise.py"),
                                _strip_to_uniprot=function)
    helper = SimpleNamespace(__file__=str(tmp_path / "benchmark_tools/prepare_qfo_factorial_pairs.py"),
                             _strip_to_uniprot=function)
    if changed == "converter_path":
        converter.__file__ = str(tmp_path / "untrusted.py")
    elif changed == "helper_path":
        helper.__file__ = str(tmp_path / "untrusted.py")
    else:
        helper._strip_to_uniprot = lambda value: value
    monkeypatch.setattr(runner.importlib, "import_module", lambda name:
        converter if name == "qfo_benchmark.og_to_pairwise" else helper)
    with pytest.raises(ValueError, match="different converter"):
        runner.verify_imported_converter(tmp_path)


@pytest.fixture
def execution(tmp_path, monkeypatch):
    source, output = tmp_path / "executor", tmp_path / "recovery"
    output.mkdir()
    artifact = tmp_path / "immutable-native"
    artifact.write_text("original native fixture\n")
    pin = record(artifact)
    raw = {name: .81 for name in ("VGNC F", "SwissTrees F", "TreeFam-A F", "EC", "GO", "FAS")}
    parent_path = tmp_path / "parent_plan.json"
    parent = {"repo": str(tmp_path), "output": str(tmp_path / "parent"),
        "scientific_revision": "original", "runtime": {"files": []},
        "control": {"scores": {key: .8 for key in raw}}}
    write_json(parent_path, parent)
    plan = {"repo": str(tmp_path), "source_root": str(source), "source_revision": "executor",
        "scientific_revision": "original", "output": str(output), "admission": pin,
        "parent_plan": record(parent_path), "parent_report": pin, "parent_root": parent["output"],
        "sources": [pin], "runtime": parent["runtime"], "scoring_id": "fixed-short-id",
        "native_inference_rerun": False}
    path = output / "execution_plan.json"
    write_json(path, plan)
    admission, prior = admitted()
    prior.update(execution_plan=plan["parent_plan"], producer=plan["parent_report"])
    admission.update(producer=plan["parent_report"], audit=pin)
    monkeypatch.setattr(runner, "__file__", str(source / "benchmark_tools/recover_cap400_scores.py"))
    monkeypatch.setattr(runner, "read_admission", lambda *args: (deepcopy(admission), deepcopy(prior)))
    monkeypatch.setattr(runner, "audit", lambda *args: deepcopy(prior))
    monkeypatch.setattr(runner, "runtime", lambda: {"files": []})
    monkeypatch.setattr(runner, "host_activity", lambda: {"scope": "fixture"})
    monkeypatch.setattr(runner, "short_id", lambda *args: "fixed-short-id")
    monkeypatch.setattr(runner, "verify_imported_converter", lambda *args: pin)
    monkeypatch.setattr(runner, "verify_helpers", lambda *args: [pin])
    calls = []
    def score(repo, phylo, folder, identifier):
        calls.append((repo, phylo, folder, identifier))
        assert phylo == Path(prior["native"]["phylogeny"])
        return {"scores": deepcopy(raw)}
    monkeypatch.setattr(runner, "score_qfo", score)
    monkeypatch.setattr(runner, "readback_scores", lambda *args: (deepcopy(raw), {}, [pin]))
    for key, value in {"SLURM_CPUS_PER_TASK": "8", "SLURM_MEM_PER_NODE": "65536",
        "SLURMD_NODENAME": "bizon", "CUDA_VISIBLE_DEVICES": "", "SLURM_JOB_ID": "new-job"}.items():
        monkeypatch.setenv(key, value)
    return path, artifact, calls, prior, raw


def test_scoring_only_success_retains_original_native_content(execution):
    path, artifact, calls, _, raw = execution
    before = record(artifact)
    result = runner.run(path)
    assert result["status"] == runner.SUCCESS and result["postflight_pins_verified"]
    assert result["scores"] == raw and len(calls) == 1
    assert all(delta == .81 - .8 for delta in result["deltas"].values())
    assert record(artifact) == before
    assert result["native_inference_rerun"] is result["defaults_promoted"] is result["goal_completed"] is False


def test_changed_native_content_refused_before_scoring(execution, monkeypatch):
    path, _, calls, prior, _ = execution
    changed = deepcopy(prior)
    changed["native"]["pairs"] += 1
    monkeypatch.setattr(runner, "audit", lambda *args: changed)
    result = runner.run(path)
    assert result["status"] == "failed_cap400_scoring_only_recovery" and not calls
    assert result["scores"] == {}


def test_scoring_error_is_retained_without_retry(execution, monkeypatch):
    path, _, _, _, _ = execution
    def fail(*args):
        raise RuntimeError("scorer failed")
    monkeypatch.setattr(runner, "score_qfo", fail)
    result = runner.run(path)
    assert result["error"] == "scorer failed" and result["scores"] == {}
    with pytest.raises(ValueError, match="attempted recovery namespace"):
        runner.run(path)


def test_postflight_changed_pin_prevents_success(execution, monkeypatch):
    path, artifact, _, _, raw = execution
    def mutate(*args):
        artifact.write_text("changed fixture\n")
        return deepcopy(raw), {}, []
    monkeypatch.setattr(runner, "readback_scores", mutate)
    result = runner.run(path)
    assert result["status"] == "failed_cap400_scoring_only_recovery"
    assert "postflight_error" in result and not result.get("postflight_pins_verified", False)


def test_wrong_allocation_refused_before_output(execution, monkeypatch):
    path, _, calls, _, _ = execution
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    with pytest.raises(ValueError, match="8-CPU 64-GiB"):
        runner.run(path)
    assert not calls and not (path.parent / "results.json").exists()


def test_changed_scoring_namespace_refused_before_scoring(execution):
    path, _, calls, _, _ = execution
    plan = json.loads(path.read_text())
    plan["scoring_id"] = "different"
    write_json(path, plan)
    result = runner.run(path)
    assert result["status"] == "failed_cap400_scoring_only_recovery" and not calls


def test_missing_raw_endpoint_refused(execution, monkeypatch):
    path, _, _, _, raw = execution
    incomplete = {key: value for key, value in raw.items() if key != "GO"}
    monkeypatch.setattr(runner, "readback_scores", lambda *args: (incomplete, {}, []))
    result = runner.run(path)
    assert result["status"] == "failed_cap400_scoring_only_recovery" and result["scores"] == {}


def test_changed_executor_inventory_refused_before_scoring(execution, monkeypatch):
    path, _, calls, _, _ = execution
    monkeypatch.setattr(runner, "verify_helpers", lambda *args: [])
    result = runner.run(path)
    assert result["status"] == "failed_cap400_scoring_only_recovery" and not calls


def test_inference_rerun_plan_refused(execution):
    path, _, calls, _, _ = execution
    plan = json.loads(path.read_text())
    plan["native_inference_rerun"] = True
    write_json(path, plan)
    result = runner.run(path)
    assert result["status"] == "failed_cap400_scoring_only_recovery" and not calls


@pytest.fixture
def frozen_helpers(tmp_path, monkeypatch):
    source, original, science = (tmp_path / name for name in ("source", "original", "science"))
    files = {"benchmark_tools/" + name: "# original helper\n" for name in runner.SCORING_HELPERS}
    files.update({"orthohmm/__init__.py": "# package fixture\n",
        "qfo_benchmark/og_to_pairwise.py": "# original converter\n",
        "benchmark_tools/" + runner.PROTOCOL: "{}\n",
        "benchmark_tools/recover_cap400_scores.slurm": "#!/bin/bash\n"})
    for name, content in files.items():
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    for name in runner.SCORING_HELPERS:
        path = original / "benchmark_tools" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(files["benchmark_tools/" + name])
    parent = {"source_root": str(original), "scientific_source": str(science),
              "scientific_revision": "science"}
    def git(command, **kwargs):
        if "rev-parse" in command:
            return "executor\n"
        _, relative = command[-1].split(":", 1)
        return files[relative].encode()
    monkeypatch.setattr(runner.subprocess, "check_output", git)
    return source, parent, files


def test_frozen_helpers_include_unchanged_original_converter(frozen_helpers):
    source, parent, _ = frozen_helpers
    pins = runner.verify_helpers(source, parent, "executor")
    assert record(source / "qfo_benchmark/og_to_pairwise.py") in pins


def test_missing_converter_detected_before_long_native_audit(frozen_helpers):
    source, parent, _ = frozen_helpers
    (source / "qfo_benchmark/og_to_pairwise.py").unlink()
    with pytest.raises(FileNotFoundError):
        runner.verify_helpers(source, parent, "executor")


def test_changed_frozen_helper_rejected(frozen_helpers):
    source, parent, _ = frozen_helpers
    (source / "benchmark_tools/run_membership_hillclimb.py").write_text("# changed\n")
    with pytest.raises(ValueError, match="committed content"):
        runner.verify_helpers(source, parent, "executor")


def test_changed_original_helper_rejected(frozen_helpers):
    source, parent, _ = frozen_helpers
    (Path(parent["source_root"]) / "benchmark_tools/run_membership_hillclimb.py").write_text("# changed\n")
    with pytest.raises(ValueError, match="original scoring helper"):
        runner.verify_helpers(source, parent, "executor")


def test_changed_executor_revision_rejected(frozen_helpers):
    source, parent, _ = frozen_helpers
    with pytest.raises(ValueError, match="executor revision"):
        runner.verify_helpers(source, parent, "different")


def test_custom_recovery_protocol_and_slurm_are_pinned(frozen_helpers):
    source, parent, files = frozen_helpers
    for relative in ("benchmark_tools/shared_protocol.json", "benchmark_tools/shared.slurm"):
        files[relative] = "# shared recovery fixture\n"
        (source / relative).write_text(files[relative])
    pins = runner.verify_helpers(source, parent, "executor",
                                 protocol="shared_protocol.json", slurm="shared.slurm")
    assert record(source / "benchmark_tools/shared_protocol.json") in pins
    assert record(source / "benchmark_tools/shared.slurm") in pins


def test_missing_custom_protocol_refused_before_execution(frozen_helpers):
    source, parent, _ = frozen_helpers
    with pytest.raises(FileNotFoundError):
        runner.verify_helpers(source, parent, "executor", protocol="missing_protocol.json")
