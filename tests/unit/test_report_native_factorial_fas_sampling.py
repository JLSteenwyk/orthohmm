import copy
import gzip
import json
from pathlib import Path

import pytest

from benchmark_tools import report_native_factorial_fas_sampling as current
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture_sample():
    values = [.1, .3, .4, .8, .5, .7, .9]
    pairs = {(f"A{i:02d}", f"B{i:02d}"): x for i, x in enumerate(values)}
    z = sum(values) / len(values)
    row = dict(index=6, cell="p0_c0_r0", participant="fixture", native_job_id=123,
               prediction_semantics="cross-species group-derived clique pairs",
               fas_sample=dict(reported_eligible_pairs=14, sample_pairs=7), scores=dict(FAS=z))
    text = ("Namespace(participant='fixture', limited_species=False)\n"
            "4 pairs precomputed, 10 missing (will compute); 0 no feature annotations\n"
            "we will compute 10 new pairs and sample 4 precomputed pairs\n"
            "FAS score[precomputed]: 0.400000 +- 0.100000 [N=4]\n"
            "FAS score[missing]: 0.700000 +- 0.100000 [N=3]\n"
            f"FAS_mean: {z} +- 0.1; nr_orthologs: 14; sample_size: 7 vs 7\n")
    return row, pairs, text


def test_sample_stratum_order_uses_exact_raw_not_rounded_log_means():
    row, pairs, text = fixture_sample()
    result = current.sample_summary(text.replace("0.400000", "0.4000004"), pairs, row)
    assert result["precomputed_sample_mean"] == pytest.approx(.4)
    assert result["returned_sample_mean"] == pytest.approx(.7)
    assert (result["P"], result["M"], result["k"], result["c"], result["r"]) == (4, 10, 4, 10, 3)
    assert result["omitted_new_numeric_scores"] == 7


@pytest.mark.parametrize("before,after", [
    ("participant='fixture'", "participant='another'"),
    ("limited_species=False", "limited_species=True"),
    ("sample 4", "sample 5"), ("compute 10 new", "compute 9 new"),
    ("[N=4]", "[N=5]"), ("[N=3]", "[N=2]"),
    ("0.400000", "0.800000"), ("nr_orthologs: 14", "nr_orthologs: 15"),
    ("sample_size: 7 vs 7", "sample_size: 7 vs 8"),
])
def test_log_mapping_and_arithmetic_changes_rejected(before, after):
    row, pairs, text = fixture_sample()
    with pytest.raises(ValueError):
        current.sample_summary(text.replace(before, after), pairs, row)


def test_ambiguous_log_batch_failure_or_reordered_raw_rejected():
    row, pairs, text = fixture_sample()
    for changed in (text + text, text + "Computing fas.runMultiTaxa failed: None"):
        with pytest.raises(ValueError):
            current.sample_summary(changed, pairs, row)
    with pytest.raises(ValueError):
        current.sample_summary(text, dict(reversed(list(pairs.items()))), row)


def bound_fixture(tmp_path):
    row, pairs, text = fixture_sample()
    work = tmp_path / "work"
    task = work / "ab" / "cdef12345"
    task.mkdir(parents=True)
    log = task / ".command.log"
    log.write_text(text)
    raw = tmp_path / "raw.txt.gz"
    with gzip.open(raw, "wt") as stream:
        stream.write("Acc1\tAcc2\tFAS\n")
        for (a, b), v in pairs.items():
            stream.write(f"{a}\t{b}\t{v}\n")
    raw_ref = record(raw)
    row["fas_sample"]["raw"] = raw_ref
    execution = tmp_path / "execution.json"
    execution.write_text(json.dumps(dict(exit_code=0, status="process_succeeded_pending_independent_admission",
                                        work=str(work), outputs=[raw_ref])))
    exec_ref = record(execution)
    scorer = record(Path(current.__file__))
    admission = dict(accuracy_admitted=True, cell=row["cell"], native_index=row["index"],
                     native_job_id=row["native_job_id"], participant=row["participant"],
                     fas_protocol=dict(source=scorer), execution_report=exec_ref,
                     checked_records=[exec_ref, raw_ref], fas_sample=dict(raw=raw_ref),
                     native_tasks=[dict(name="fas_benchmark (1)", status="COMPLETED", exit="0", hash="ab/cdef")])
    path = tmp_path / "admission.json"
    path.write_text(json.dumps(admission))
    row["admission"] = record(path)
    return row, record(log), scorer, admission, path


def test_raw_and_task_bound_to_successful_admission_without_historical_log_claim(tmp_path):
    row, log_ref, scorer, _, _ = bound_fixture(tmp_path)
    evidence = []
    result = current.read_cell(row, log_ref, scorer, evidence)
    assert result["native_log_historically_hash_bound"] is False
    assert len(evidence) == 4 and result["observed_native_mean"] == row["scores"]["FAS"]


@pytest.mark.parametrize("change", ["admission", "execution", "raw", "task", "path", "scorer", "participant"])
def test_changed_native_bindings_refused(tmp_path, change):
    row, log_ref, scorer, admission, path = bound_fixture(tmp_path)
    if change == "admission":
        admission["accuracy_admitted"] = False
    elif change == "execution":
        admission["checked_records"].remove(admission["execution_report"])
    elif change == "raw":
        admission["checked_records"].remove(row["fas_sample"]["raw"])
    elif change == "task":
        admission["native_tasks"][0]["exit"] = "1"
    elif change == "path":
        admission["native_tasks"][0]["hash"] = "ef/another"
    elif change == "scorer":
        scorer = dict(scorer, sha256="0" * 64)
    elif change == "participant":
        admission["participant"] = "wrong"
    path.write_text(json.dumps(admission))
    row["admission"] = record(path)
    with pytest.raises(ValueError):
        current.read_cell(row, log_ref, scorer, [])


def panel_rows():
    row, pairs, text = fixture_sample()
    counts = current.sample_summary(text, pairs, row)
    return [dict(index=index, cell=cell, **counts) for index, cell in current.CELLS]


def test_four_cells_all_six_contrasts_and_scope():
    result = current.panel(panel_rows())
    assert len(result["methods"]) == 4 and len(result["contrasts"]) == 6
    assert result["components"] == 12 and result["joint_error_bound"] == .05
    assert result["conditional_design_ranges_computed"] is True
    assert all(r["zero_included"] is True for r in result["contrasts"])
    for key in ("historical_scores_rerun", "observed_scores_changed", "biological_generalization_intervals",
                "unconditional_historical_interval_admission", "other_endpoint_uncertainty_admitted",
                "failed_factorial_cells_repaired", "publication_ready"):
        assert result[key] is False
    assert current.render(result).count("| p0_") == 9


def test_numerical_failure_preserved_with_only_dependent_contrasts_unavailable(monkeypatch):
    rows = panel_rows()
    rows[0]["P"] = -1
    result = current.panel(rows)
    assert result["status"] == "conditional_design_ranges_partial"
    assert result["methods"][0]["error"]["type"] == "ValueError"
    assert "traceback" in result["methods"][0]["error"]
    assert sum(r["conditional_expected_difference_bounds"] is None for r in result["contrasts"]) == 3
    assert result["methods"][0]["observed_native_mean"] == rows[0]["observed_native_mean"]


def test_reordered_or_missing_cells_not_silently_selected():
    rows = panel_rows()
    for bad in (rows[:3], list(reversed(rows)), rows + [copy.deepcopy(rows[0])]):
        with pytest.raises(ValueError):
            current.panel(bad)


def test_fresh_output_and_preserved_fatal_failure(tmp_path, monkeypatch):
    output = tmp_path / "result"
    def fail(_):
        raise ValueError("fixture input defect")
    monkeypatch.setattr(current, "build", fail)
    with pytest.raises(ValueError, match="fixture input defect"):
        current.run(tmp_path, output)
    retained = json.loads((output / "failure.json").read_text())
    assert retained["historical_scores_rerun"] is False
    assert not (output / "report.json").exists()
    with pytest.raises(ValueError, match="Output already exists"):
        current.run(tmp_path, output)


def test_success_writes_bound_result_without_new_scoring(tmp_path, monkeypatch):
    result = current.panel(panel_rows())
    monkeypatch.setattr(current, "build", lambda _: result)
    output = tmp_path / "result"
    assert current.run(tmp_path, output) == result
    assert json.loads((output / "report.json").read_text()) == result
    assert (output / "summary.md").read_text() == current.render(result)
