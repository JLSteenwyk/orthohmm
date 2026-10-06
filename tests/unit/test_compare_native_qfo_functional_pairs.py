import copy
import gzip
import json
from pathlib import Path

import pytest

from benchmark_tools import compare_native_qfo_functional_pairs as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_compare_qfo_scored_pairs import raw
from tests.unit.test_export_native_qfo_factorial_scores import write

ROOT = Path(__file__).resolve().parents[2]


def fixture(tmp_path, monkeypatch):
    snapshot = json.loads((ROOT / "benchmark_tools/results/native_qfo_scientific_scores_20261006_v1/report.json").read_text())
    for row in snapshot["rows"][:2]:
        cell = row["cell"]
        references = []
        for metric in ("GO", "EC", "FAS"):
            directory = tmp_path / cell / "results" / metric
            directory.mkdir(parents=True)
            path = directory / "fixture_raw.txt.gz"
            if metric == "FAS":
                with gzip.open(path, "wt") as stream:
                    stream.write("Acc1\tAcc2\tFAS\nA\tB\t0.5\nA\tC\t0.9\n")
            else:
                raw(path, "B\tA\t0.500000\nA\tC\t0.900000\n", metric)
                row["endpoint_details"][metric]["assessed_relations"] = 2
            row["scores"][metric] = .7
            references.append(record(path))
        execution_ref = write(tmp_path / cell / "execution.json", dict(
            exit_code=0, status="process_succeeded_pending_independent_admission", outputs=references))
        admission = dict(cell=cell, native_index=row["index"], native_job_id=row["native_job_id"],
            accuracy_admitted=True, execution_report=execution_ref, checked_records=[execution_ref, *references],
            fas_sample=dict(raw=references[-1], sample_pairs=2))
        row["admission"] = write(tmp_path / cell / "admission.json", admission)
    ref = write(tmp_path / "snapshot.json", snapshot)
    monkeypatch.setattr(module.reporter, "collect", lambda *args: copy.deepcopy(snapshot))
    return snapshot, ref


def test_fixture_binds_six_raw_tables_without_uncertainty_or_timing(tmp_path, monkeypatch):
    _, ref = fixture(tmp_path, monkeypatch)
    result = module.run(Path(ref["path"]), ref["sha256"], tmp_path / "result.json")
    assert len(result["endpoints"]) == 6 and len(result["comparisons"]) == 3
    assert result["cells"] == ["p0_c0_r1", "p0_c0_r0"]
    for key in ("uncertainty_admitted", "scientific_timings_admitted", "publication_ready", "new_scoring_or_admission"):
        assert result[key] is False
    assert result["new_bootstrap_draws"] == 0
    assert all(row["proteins_in_multiple_pairs"] == 1 for row in result["endpoints"])
    assert result["comparisons"][0]["result"]["original_mean_difference"] == 0
    assert result["comparisons"][2]["result"]["shared_sample_pairs"] == 2


@pytest.mark.parametrize("change", ["not_admitted", "duplicate_cell", "resources", "timing_admitted",
    "timing_eligible", "score", "count", "replay", "source", "publication_ready"])
def test_snapshot_and_native_failure_contracts(tmp_path, monkeypatch, change):
    snapshot, _ = fixture(tmp_path, monkeypatch)
    if change == "not_admitted": snapshot["rows"][1]["accuracy_admitted"] = False
    elif change == "duplicate_cell": snapshot["rows"][2]["cell"] = "p0_c0_r1"
    elif change in ("resources", "timing_admitted", "timing_eligible"):
        snapshot["rows"][1][change] = {} if change == "resources" else True
    elif change == "score": snapshot["rows"][0]["scores"]["GO"] = .9
    elif change == "count": snapshot["rows"][0]["endpoint_details"]["GO"]["assessed_relations"] = 3
    elif change == "source": snapshot["source"]["sha256"] = "0" * 64
    elif change == "publication_ready": snapshot[change] = True
    else: monkeypatch.setattr(module.reporter, "collect", lambda *args: {"supplied_admissions": 0})
    ref = write(tmp_path / "changed_snapshot.json", snapshot)
    output = tmp_path / "result.json"
    with pytest.raises(ValueError): module.run(Path(ref["path"]), ref["sha256"], output)
    assert not output.exists()


@pytest.mark.parametrize("change", ["admission_identity", "execution_unpinned", "execution_failed",
    "execution_duplicate_raw", "raw_unadmitted", "fas_binding", "raw_mutated"])
def test_raw_binding_requires_admission_and_original_execution_pins(tmp_path, monkeypatch, change):
    snapshot, _ = fixture(tmp_path, monkeypatch)
    row = snapshot["rows"][0]
    admission_path = Path(row["admission"]["path"])
    admission = json.loads(admission_path.read_text())
    execution_path = Path(admission["execution_report"]["path"])
    execution = json.loads(execution_path.read_text())
    if change == "admission_identity": admission["native_job_id"] = 99
    elif change == "execution_unpinned": admission["checked_records"].remove(admission["execution_report"])
    elif change in ("execution_failed", "execution_duplicate_raw"):
        if change == "execution_failed": execution["exit_code"] = 1
        else: execution["outputs"].append(execution["outputs"][0])
        old = admission["execution_report"]
        admission["execution_report"] = write(execution_path, execution)
        admission["checked_records"] = [admission["execution_report"] if r == old else r for r in admission["checked_records"]]
    elif change == "raw_unadmitted": admission["checked_records"].remove(execution["outputs"][0])
    elif change == "fas_binding": admission["fas_sample"]["raw"] = execution["outputs"][0]
    else: Path(execution["outputs"][0]["path"]).write_bytes(b"changed raw bytes")
    row["admission"] = write(admission_path, admission)
    with pytest.raises(ValueError): module.raw_bindings(row, [])


def test_fas_overlap_retains_original_denominators_and_shared_score_changes():
    left = {("A", "B"): .5, ("A", "C"): .9}
    right = {("A", "B"): .4, ("B", "C"): .2, ("C", "D"): .3}
    result = module.fas_overlap(left, right)
    assert result["shared_sample_pairs"] == 1
    assert result["shared_pairs_with_different_serialized_scores"] == 1
    assert result["shared_conditional_mean_difference"] == pytest.approx(.1)
    assert result["original_sample_mean_difference"] == pytest.approx(.4)
    assert sum(result["original_sample_mean_difference_components"].values()) == pytest.approx(.4)
    assert result["shared_fraction_of_left"] == .5
    assert result["shared_fraction_of_right"] == pytest.approx(1/3)


def test_fas_disjoint_samples_have_no_shared_conditional_statistic():
    result = module.fas_overlap({("A", "B"): .5, ("A", "C"): .9},
                                {("D", "E"): .5, ("D", "F"): .9})
    assert result["shared_sample_pairs"] == 0
    assert result["shared_conditional_mean_difference"] is result["maximum_shared_absolute_difference"] is None


def test_existing_output_refused_before_input_access(tmp_path):
    output = tmp_path / "retain.json"
    output.write_text("retain")
    with pytest.raises(ValueError, match="Output already exists"):
        module.run(Path("missing"), "missing", output)
    assert output.read_text() == "retain"
