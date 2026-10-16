from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import recover_shared_qfo_export as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.report_hillclimb import endpoint_vector
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json
from benchmark_tools.shared_qfo_export import HEADER, CONFIDENCE_HEADER, verify_projection

ROW = {"JobIDRaw": "24247", "State": "FAILED", "ExitCode": "1:0",
       "AllocCPUS": "16", "NodeList": "bizon"}


@pytest.fixture
def failure():
    evidence = {"status": "shared_qfo_terminal_duplicate_export_failure_reproduced_with_original_provenance",
        "job_id": "24247", "scientific_revision": runner.REVISION,
        "pair_and_confidence_stream_audit": {"rows_each": runner.ROWS, "distinct_pairs": runner.DISTINCT,
            "adjacent_duplicate_rows_each": 1, "duplicate_confidence_conflicts": 0}}
    report = {"job_id": "24247", "status": "failed_or_unverified_shared_qfo_continuation",
        "error_type": "ValueError", "error": "Native pairs are not canonical, unique and sorted",
        "scores": {}, "execution": {"returncode": 0}}
    return evidence, report


def test_exact_failed_parent_required(failure):
    assert runner.require_failure(*failure, [ROW]) == ROW


@pytest.mark.parametrize("rows", [[], [ROW, ROW], [{**ROW, "State": "RUNNING"}],
    [{**ROW, "ExitCode": "0:0"}], [{**ROW, "AllocCPUS": "8"}], [{**ROW, "NodeList": "other"}]])
def test_wrong_scheduler_refused(failure, rows):
    with pytest.raises(ValueError):
        runner.require_failure(*failure, rows)


@pytest.mark.parametrize("field,value", [("job_id", "other"), ("status", "running"),
    ("error_type", "RuntimeError"), ("error", "different"), ("scores", {"GO": .5}),
    ("execution", {"returncode": 1})])
def test_other_failure_is_not_recoverable(failure, field, value):
    evidence, report = failure
    report[field] = value
    with pytest.raises(ValueError):
        runner.require_failure(evidence, report, [ROW])


@pytest.mark.parametrize("field,value", [("rows_each", 3), ("distinct_pairs", 2),
    ("adjacent_duplicate_rows_each", 2), ("duplicate_confidence_conflicts", 1)])
def test_unplanned_export_shape_refused(failure, field, value):
    evidence, report = failure
    evidence["pair_and_confidence_stream_audit"][field] = value
    with pytest.raises(ValueError):
        runner.require_failure(evidence, report, [ROW])


def test_derived_view_keeps_literal_records_and_all_other_products(tmp_path, monkeypatch):
    folder, output = tmp_path / "native", tmp_path / "recovery"
    phylo = folder / "output/orthohmm_phylogeny"
    phylo.mkdir(parents=True)
    output.mkdir()
    pair, other = b"a\tA\tb\tB\n", b"a\tA\tc\tC\n"
    annotated, other_annotated = pair[:-1] + b"\tmedium\n", other[:-1] + b"\thigh\n"
    (phylo / PRODUCTS[1]).write_bytes(HEADER + pair * 2 + other)
    (phylo / PRODUCTS[2]).write_bytes(CONFIDENCE_HEADER + annotated * 2 + other_annotated)
    for name in (PRODUCTS[0], "provenance_manifest.json", "reconciliation_summary.json", "species_tree.rooted.nwk"):
        (phylo / name).write_text("original immutable fixture\n")
    (folder / "metrics.json").write_text("{}\n")
    before = {p: record(p) for p in phylo.iterdir()}
    monkeypatch.setattr(runner, "DUPLICATE", ("a", "A", "b", "B"))
    monkeypatch.setattr(runner, "ROWS", 3)
    monkeypatch.setattr(runner, "DISTINCT", 2)
    monkeypatch.setattr(runner, "DUPLICATE_LINE", 3)
    def validate(case, parent, pair_export):
        assert verify_projection(pair_export)["distinct_rows"] == 2
        return {"status": "derived_native_content_verified", "genes": 984137,
                "species": 78, "groups": 1, "pairs": 2}
    monkeypatch.setattr(runner, "shared_native", validate)
    native, manifest = runner.derive({"case": {"output": str(folder)}}, {}, output)
    assert native["status"] == "derived_native_content_verified"
    assert manifest["normalization"]["removed_rows"] == 1
    assert {p: record(p) for p in phylo.iterdir()} == before
    view = Path(native["phylogeny"])
    assert (view / PRODUCTS[1]).read_bytes() == HEADER + pair + other
    assert (view / PRODUCTS[2]).read_bytes() == CONFIDENCE_HEADER + annotated + other_annotated
    assert (view / PRODUCTS[0]).read_bytes() == (phylo / PRODUCTS[0]).read_bytes()
    assert not manifest["native_inference_rerun"]
    with pytest.raises(FileExistsError):
        runner.derive({"case": {"output": str(folder)}}, {}, output)


@pytest.fixture
def execution(tmp_path, monkeypatch):
    source, output = tmp_path / "executor", tmp_path / "recovery"
    output.mkdir()
    artifact = tmp_path / "immutable"
    artifact.write_text("original native fixture\n")
    pin = record(artifact)
    raw = {name: .81 for name in ("VGNC F", "SwissTrees F", "TreeFam-A F", "EC", "GO", "FAS")}
    original = {"output": str(tmp_path / "parent"), "runtime": {"files": []}}
    failed = {"control": {"report": pin, "cell": {"scores": endpoint_vector("QfO", {k: .8 for k in raw})}}}
    parent = Path(original["output"])
    parent.mkdir()
    write_json(parent / "results.json", {"status": "retained_failed_fixture"})
    plan = {"repo": str(tmp_path), "source_root": str(source), "source_revision": "executor",
        "scientific_revision": runner.REVISION, "output": str(output), "failure": pin,
        "sources": [pin], "runtime": original["runtime"], "scoring_id": "fixed-id",
        "native_inference_rerun": False, "scope": "explicit derived export fixture"}
    path = output / "execution_plan.json"
    write_json(path, plan)
    native = {"genes": 984137, "species": 78, "groups": 370016, "pairs": 5561114,
              "phylogeny": str(output / "canonical_export"), "products": {}, "summary": {}}
    manifest = {"unchanged_copied_products": [], "normalization": {"derived": {"pairs": pin}}}
    def derive(*args):
        write_json(output / "normalization.json", manifest)
        return deepcopy(native), deepcopy(manifest)
    monkeypatch.setattr(runner, "__file__", str(source / "benchmark_tools/recover_shared_qfo_export.py"))
    monkeypatch.setattr(runner, "parent_context", lambda *args: (original, {}, failed, [pin]))
    monkeypatch.setattr(runner, "runtime", lambda: {"files": []})
    monkeypatch.setattr(runner, "host_activity", lambda: {})
    monkeypatch.setattr(runner, "verify_helpers", lambda *args, **kwargs: [pin])
    monkeypatch.setattr(runner, "verify_imported_converter", lambda *args: pin)
    monkeypatch.setattr(runner, "short_id", lambda *args: "fixed-id")
    monkeypatch.setattr(runner, "derive", derive)
    calls = []
    def score(repo, phylo, folder, identifier):
        assert phylo == Path(native["phylogeny"]) and identifier == "fixed-id"
        calls.append(identifier)
        return {"scores": deepcopy(raw)}
    monkeypatch.setattr(runner, "score_qfo", score)
    monkeypatch.setattr(runner, "readback_scores", lambda *args: (deepcopy(raw), {}, [pin]))
    for key, value in {"SLURM_CPUS_PER_TASK": "8", "SLURM_MEM_PER_NODE": "65536",
        "SLURMD_NODENAME": "bizon", "CUDA_VISIBLE_DEVICES": "", "SLURM_JOB_ID": "new-job"}.items():
        monkeypatch.setenv(key, value)
    return path, artifact, calls, raw


def test_derived_export_scoring_success_preserves_parent_and_complete_vector(execution):
    path, artifact, calls, raw = execution
    before = record(artifact)
    result = runner.run(path)
    assert result["status"] == runner.SUCCESS and result["postflight_pins_verified"]
    assert result["scores"] == raw and calls == ["fixed-id"]
    assert all(value == .81 - .8 for value in result["deltas"].values())
    assert record(artifact) == before
    assert result["native_inference_rerun"] is result["defaults_promoted"] is result["goal_completed"] is False


def test_scoring_not_called_when_normalization_or_native_validation_fails(execution, monkeypatch):
    path, _, calls, _ = execution
    def fail(*args):
        raise ValueError("lossless gate failed")
    monkeypatch.setattr(runner, "derive", fail)
    result = runner.run(path)
    assert result["status"] == "failed_shared_qfo_lossless_export_recovery" and not calls
    assert result["scores"] == {} and result["error"] == "lossless gate failed"


def test_missing_endpoint_not_admitted(execution, monkeypatch):
    path, _, _, raw = execution
    incomplete = {k: v for k, v in raw.items() if k != "GO"}
    monkeypatch.setattr(runner, "readback_scores", lambda *args: (incomplete, {}, []))
    result = runner.run(path)
    assert result["status"] == "failed_shared_qfo_lossless_export_recovery" and result["scores"] == {}


def test_postflight_changed_original_invalidates_recovery(execution, monkeypatch):
    path, artifact, _, raw = execution
    def mutate(*args):
        artifact.write_text("changed fixture\n")
        return raw, {}, []
    monkeypatch.setattr(runner, "readback_scores", mutate)
    result = runner.run(path)
    assert result["status"] == "failed_shared_qfo_lossless_export_recovery" and "postflight_error" in result


def test_attempted_namespace_not_reused(execution):
    path, _, _, _ = execution
    runner.run(path)
    with pytest.raises(ValueError, match="attempted recovery"):
        runner.run(path)


def test_wrong_resources_refused(execution, monkeypatch):
    path, _, calls, _ = execution
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    with pytest.raises(ValueError, match="8-CPU 64-GiB"):
        runner.run(path)
    assert not calls


def test_relabeling_as_inference_rerun_refused(execution):
    path, _, calls, _ = execution
    plan = json.loads(path.read_text())
    plan["native_inference_rerun"] = True
    write_json(path, plan)
    result = runner.run(path)
    assert result["status"] == "failed_shared_qfo_lossless_export_recovery" and not calls


@pytest.fixture
def matched_control(tmp_path, monkeypatch):
    root = tmp_path / "control"
    root.mkdir()
    path = root / "execution_plan.json"
    write_json(path, {})
    write_json(root / "results.json", {"job_id": "24246", "execution_plan": record(path)})
    pin = record(root / "results.json")
    cell = {"scores": {"QfO: GO": .49}, "products": {}, "coverage": {"genes": 984137},
            "source_revision": "canonical", "checked": [pin]}
    original = {"recovery_root": str(root), "inputs": [record(path)]}
    report = {"control": {"report": pin, "cell": deepcopy(cell)}}
    monkeypatch.setattr(runner, "check_plan", lambda *args: ({}, []))
    monkeypatch.setattr(runner, "canonical_recovery_cell", lambda *args: deepcopy(cell))
    return original, report, cell


def test_original_control_independently_read_back(tmp_path, matched_control):
    original, report, _ = matched_control
    assert runner.verify_control(tmp_path, original, report)


@pytest.mark.parametrize("field,value", [("scores", {"QfO: GO": 1}), ("products", {"changed": True}),
    ("coverage", {"genes": 1}), ("source_revision", "changed")])
def test_changed_control_refused(tmp_path, matched_control, field, value):
    original, report, cell = matched_control
    cell[field] = value
    with pytest.raises(ValueError, match="matched-control readback differs"):
        runner.verify_control(tmp_path, original, report)
