"""Real conversion kernels with synthetic admission context, not production scores."""

import gzip
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import prepare_native11_composed_qfo_pairs as current
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.validate_native_factorial_outputs import Evidence
from tests.unit.test_prepare_native_factorial_qfo_pairs import conversion_fixture


def test_actual_group_kernel_counts_and_full_input_coverage(conversion_fixture, tmp_path):
    inputs, owners, groups, native, mapping = conversion_fixture
    result = current.materialize(groups, inputs, owners, mapping, tmp_path)
    assert result["expected_pairs"] == result["total_pairs"] == result["retained_pairs"] == 5
    assert result["pair_coverage"]["input_accessions"] == 5
    assert result["pair_coverage"]["accessions_in_any_pair"] == 4
    assert result["pair_coverage"]["fraction_inputs_in_any_pair"] == .8
    text = (tmp_path / "pairs.qfo.partial.tsv").read_text()
    assert "A1\tA2\n" not in text
    assert record(tmp_path / "pairs.partial.tsv")["sha256"] == record(tmp_path / "pairs.qfo.partial.tsv")["sha256"]
    with pytest.raises(ValueError, match="fresh"):
        current.materialize(groups, inputs, owners, mapping, tmp_path)


def test_singletons_have_zero_pairs_without_fabricated_score(conversion_fixture, tmp_path):
    inputs, owners, groups, native, mapping = conversion_fixture
    groups.write_text("\n".join(owners) + "\n")
    result = current.materialize(groups, inputs, owners, mapping, tmp_path)
    assert result["empty_predictions"] is True
    assert result["expected_pairs"] == 0
    assert result["pair_coverage"]["input_accessions"] == 5
    assert result["pair_coverage"]["fraction_inputs_in_any_pair"] == 0
    assert "accuracy_evaluated" not in result


def test_unmapped_input_rejected_even_if_singleton(conversion_fixture, tmp_path):
    inputs, owners, groups, native, mapping = conversion_fixture
    del mapping["U1"]
    with pytest.raises(ValueError, match="absent from frozen"):
        current.materialize(groups, inputs, owners, mapping, tmp_path)
    assert not (tmp_path / "pairs.partial.tsv").exists()


@pytest.fixture
def joined(conversion_fixture, tmp_path, monkeypatch):
    inputs, owners, groups, native, mapping = conversion_fixture
    root = tmp_path / "repo"
    root.mkdir()
    tools = root / "benchmark_tools"
    tools.mkdir()
    original = Path(current.__file__).parent
    for name in ("native11_composed_review_binding.py", "prepare_native_factorial_qfo_pairs.py",
        "prepare_qfo_corrected_group_pairs.py", "prepare_qfo_factorial_pairs.py", "qfo_filter_pairs.py",
        "simulation_method_outputs.py", "validate_native_factorial_outputs.py"):
        shutil.copyfile(original / name, tools / name)
    (root / "qfo_benchmark").mkdir()
    shutil.copyfile(original.parent / "qfo_benchmark/og_to_pairwise.py", root / "qfo_benchmark/og_to_pairwise.py")
    output = root / "native"
    shutil.copytree(inputs, output / "input")
    prediction = output / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
    prediction.parent.mkdir(parents=True)
    shutil.copyfile(groups, prediction)
    run = dict(index=11, cell="p1_c1_r0", dataset="qfo_corrected", repeat=0,
        output_root=str(output), inputs=[record(p) for p in sorted(inputs.glob("*.fasta"))],
        native_order=[p.name for p in sorted(inputs.glob("*.fasta"))], genes=5, proteomes=4)
    _, digest, counts = current.universe(dict(run, input_directory=str(output / "input")), Evidence())
    outputs = dict(gene_ownership_sha256=digest, per_species_counts=counts,
        checked_files=[record(prediction)], native_cpu_ids=[1], allocated_ready={})
    prepared = root / "prepared.json"
    prepared.write_text(json.dumps(dict(input_fastas=run["inputs"])))
    mapping_path = root / "mapping.json.gz"
    with gzip.open(mapping_path, "wt") as stream:
        json.dump(dict(mapping=mapping), stream)
    results = tools / "results"
    results.mkdir()
    environment = results / "qfo_assessment_environment_20260917.json"
    environment.write_text(json.dumps(dict(reference_files=[record(mapping_path)])))
    request_ref, review_ref = {"request": 1}, {"composed_review": 1}
    request = dict(plan={"plan": 1}, amendment={"amendment": 1}, job_id=23985)
    binding = dict(composed_schema_preserved=True, original_review_translated=False, next_identity_authorized=False)
    context = [request, {}, {}, run, {}, outputs, "group", {}, [], binding]
    monkeypatch.setattr(current, "native_binding", lambda *a: tuple(context))
    monkeypatch.setattr(current, "ROOT", root)
    monkeypatch.setattr(current, "FIXED_INPUTS", {"qfo_preparation": ("prepared.json", record(prepared)["sha256"])})
    monkeypatch.setattr(current, "ENV_SHA", record(environment)["sha256"])
    destination = root / "conversion"
    monkeypatch.setattr(current, "DESTINATION", destination)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setenv("SLURM_JOB_ID", "25001")
    return request_ref, review_ref, destination, record(current.__file__)["sha256"], context


def test_joined_real_conversion_preserves_explicit_schema_and_no_score(joined):
    request_ref, review_ref, destination, source_sha, context = joined
    ref = current.prepare(request_ref, review_ref, destination, source_sha)
    report = json.loads(Path(ref["path"]).read_text())
    assert report["schema"] == current.SCHEMA and report["status"] == current.STATUS
    assert report["terminal_review"] == review_ref
    assert report["conversion_kind"] == "group" and report["semantics"] == current.SEMANTICS
    assert report["total_pairs"] == report["retained_pairs"] == report["expected_pairs"] == 5
    assert report["participant"] == current.PARTICIPANT
    assert report["pair_coverage"]["input_accessions"] == 5
    assert report["conversion_finished_monotonic_ns"] >= report["conversion_started_monotonic_ns"]
    assert (destination / "pairs.qfo.tsv").is_file()
    assert not (destination / "pairs.qfo.partial.tsv").exists()
    for key in ("accuracy_evaluated", "publication_ready", "next_identity_authorized", "automatic_retry",
                "native_inference_reexecuted", "original_review_translated"):
        assert report[key] is False
    with pytest.raises(ValueError, match="fresh"):
        current.prepare(request_ref, review_ref, destination, source_sha)


@pytest.mark.parametrize("change", ["ownership", "prediction", "kind", "source", "unscheduled"])
def test_invalid_inputs_refuse_before_conversion(joined, change, monkeypatch):
    request_ref, review_ref, destination, source_sha, context = joined
    if change == "ownership":
        context[5]["gene_ownership_sha256"] = "wrong"
    elif change == "prediction":
        context[5]["checked_files"] = []
    elif change == "kind":
        context[6] = "native"
    elif change == "source":
        source_sha = "wrong"
    else:
        monkeypatch.delenv("SLURM_JOB_ID")
    with pytest.raises(ValueError):
        current.prepare(request_ref, review_ref, destination, source_sha)
    assert not destination.exists()


def test_failed_materialization_retains_preflight_and_failure_report(joined, monkeypatch):
    request_ref, review_ref, destination, source_sha, context = joined
    def failed(*args):
        raise ValueError("fixture conversion failure")
    monkeypatch.setattr(current, "materialize", failed)
    with pytest.raises(ValueError, match="fixture conversion failure"):
        current.prepare(request_ref, review_ref, destination, source_sha)
    assert (destination / "preflight.json").exists()
    report = json.loads((destination / "results.json").read_text())
    assert report["status"] == "native11_composed_qfo_conversion_failed_retained"
    assert report["accuracy_evaluated"] is False
