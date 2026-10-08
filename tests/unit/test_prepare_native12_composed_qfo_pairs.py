"""Real native-pair conversion under explicitly synthetic review context."""

import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from benchmark_tools import prepare_native12_composed_qfo_pairs as current
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.validate_native_factorial_outputs import Evidence
from tests.unit.test_prepare_native_factorial_qfo_pairs import conversion_fixture


def test_actual_native_kernel_does_not_expand_group_cliques(conversion_fixture, tmp_path):
    inputs, owners, groups, native, mapping = conversion_fixture
    result = current.materialize(native, inputs, owners, mapping, tmp_path, 1)
    assert result["expected_pairs"] == result["total_pairs"] == result["retained_pairs"] == 1
    assert result["command"] is None
    assert result["pair_coverage"]["input_accessions"] == 5
    assert result["pair_coverage"]["accessions_in_any_pair"] == 2
    assert result["pair_coverage"]["fraction_inputs_in_any_pair"] == .4
    assert (tmp_path / "pairs.qfo.partial.tsv").read_text() == "A1\tB1\n"
    assert record(tmp_path / "pairs.partial.tsv")["sha256"] == record(tmp_path / "pairs.qfo.partial.tsv")["sha256"]
    assert not (tmp_path / "conversion.log").exists()
    with pytest.raises(ValueError, match="fresh"):
        current.materialize(native, inputs, owners, mapping, tmp_path, 1)


def test_empty_native_predictions_retain_full_coverage_denominator(conversion_fixture, tmp_path):
    inputs, owners, groups, native, mapping = conversion_fixture
    native.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\n")
    result = current.materialize(native, inputs, owners, mapping, tmp_path, 0)
    assert result["empty_predictions"] is True and result["expected_pairs"] == 0
    assert result["pair_coverage"]["input_accessions"] == 5
    assert result["pair_coverage"]["fraction_inputs_in_any_pair"] == 0
    assert "accuracy_evaluated" not in result


@pytest.mark.parametrize("expected", [2, True, None, -1, 1.0])
def test_invalid_or_mismatched_native_count_refuses(conversion_fixture, tmp_path, expected):
    inputs, owners, groups, native, mapping = conversion_fixture
    with pytest.raises(ValueError):
        current.materialize(native, inputs, owners, mapping, tmp_path, expected)
    assert not (tmp_path / "pairs.qfo.partial.tsv").exists()


def test_unmapped_singleton_refuses_before_materialization(conversion_fixture, tmp_path):
    inputs, owners, groups, native, mapping = conversion_fixture
    del mapping["U1"]
    with pytest.raises(ValueError, match="absent from frozen"):
        current.materialize(native, inputs, owners, mapping, tmp_path, 1)
    assert not (tmp_path / "pairs.partial.tsv").exists()


def test_clique_file_cannot_substitute_for_native_pairs(conversion_fixture, tmp_path):
    inputs, owners, groups, native, mapping = conversion_fixture
    with pytest.raises(ValueError):
        current.materialize(groups, inputs, owners, mapping, tmp_path, 5)
    assert not (tmp_path / "pairs.qfo.partial.tsv").exists()


def test_dangling_partial_output_is_not_fresh(conversion_fixture, tmp_path):
    inputs, owners, groups, native, mapping = conversion_fixture
    (tmp_path / "pairs.partial.tsv").symlink_to(tmp_path / "missing.tsv")
    with pytest.raises(ValueError, match="fresh"):
        current.materialize(native, inputs, owners, mapping, tmp_path, 1)


@pytest.fixture
def joined(conversion_fixture, tmp_path, monkeypatch):
    inputs, owners, groups, native, mapping = conversion_fixture
    root = tmp_path / "repo"
    tools = root / "benchmark_tools"
    tools.mkdir(parents=True)
    original = Path(current.__file__).parent
    for name in ("native12_composed_review_binding.py", "prepare_native_factorial_qfo_pairs.py",
        "prepare_qfo_factorial_pairs.py", "qfo_filter_pairs.py", "simulation_method_outputs.py",
        "validate_native_factorial_outputs.py"):
        shutil.copyfile(original / name, tools / name)
    output = root / "native"
    shutil.copytree(inputs, output / "input")
    prediction = output / "native/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv"
    prediction.parent.mkdir(parents=True)
    shutil.copyfile(native, prediction)
    run = dict(index=12, cell="p1_c1_r1", dataset="qfo_corrected", repeat=0,
        output_root=str(output), inputs=[record(path) for path in sorted(inputs.glob("*.fasta"))],
        native_order=[path.name for path in sorted(inputs.glob("*.fasta"))], genes=5, proteomes=4)
    _, digest, counts = current.universe(dict(run, input_directory=str(output / "input")), Evidence())
    outputs = dict(gene_ownership_sha256=digest, per_species_counts=counts,
        checked_files=[record(prediction)], native_cpu_ids=[1], allocated_ready={}, phylogeny=dict(native_pair_rows=1))
    prepared = root / "prepared.json"
    prepared.write_text(json.dumps(dict(input_fastas=run["inputs"])))
    mapping_path = root / "mapping.json.gz"
    with gzip.open(mapping_path, "wt") as stream:
        json.dump(dict(mapping=mapping), stream)
    results = tools / "results"
    results.mkdir()
    environment = results / "qfo_assessment_environment_20260917.json"
    environment.write_text(json.dumps(dict(reference_files=[record(mapping_path)])))
    request_ref, review_ref, held_ref, release_ref = ({"request": 12}, {"review": 12}, {"held": 1}, {"release": 1})
    request = dict(plan={"plan": 1}, amendment={"amendment": 1}, job_id=24036)
    binding = dict(review_producer_job_id=25000, review_held=held_ref, review_release=release_ref,
        composed_schema_preserved=True, original_review_translated=False, next_identity_authorized=False)
    context = [request, {}, {}, run, {}, outputs, "native", {}, [], binding]
    monkeypatch.setattr(current, "native_binding", lambda *args: tuple(context))
    monkeypatch.setattr(current, "ROOT", root)
    monkeypatch.setattr(current, "FIXED_INPUTS", {"qfo_preparation": ("prepared.json", record(prepared)["sha256"])})
    monkeypatch.setattr(current, "ENV_SHA", record(environment)["sha256"])
    destination = root / "conversion"
    monkeypatch.setattr(current, "DESTINATION", destination)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setenv("SLURM_JOB_ID", "25001")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "32768")
    args = (request_ref, review_ref, 25000, held_ref, release_ref, destination, record(current.__file__)["sha256"])
    return dict(args=args, context=context, prediction=prediction, destination=destination)


def test_joined_conversion_retains_new_schema_and_only_native_pairs(joined):
    ref = current.prepare(*joined["args"])
    report = json.loads(Path(ref["path"]).read_text())
    assert report["schema"] == current.SCHEMA and report["status"] == current.STATUS
    assert report["terminal_review"] == joined["args"][1]
    assert report["composed_binding"]["review_producer_job_id"] == 25000
    assert report["conversion_kind"] == "native" and report["semantics"] == current.SEMANTICS
    assert report["native_job_id"] == 24036 and report["native_index"] == 12
    assert report["total_pairs"] == report["retained_pairs"] == report["expected_pairs"] == 1
    assert report["participant"] == current.PARTICIPANT
    assert report["pair_coverage"]["input_accessions"] == 5
    assert report["conversion_finished_monotonic_ns"] >= report["conversion_started_monotonic_ns"]
    assert Path(report["filtered_pairs"]["path"]).read_text() == "A1\tB1\n"
    assert not (joined["destination"] / "conversion.log").exists()
    for key in ("accuracy_evaluated", "publication_ready", "next_identity_authorized", "automatic_retry",
                "native_inference_reexecuted", "original_review_translated"):
        assert report[key] is False
    with pytest.raises(ValueError, match="fresh"):
        current.prepare(*joined["args"])


@pytest.mark.parametrize("change", ["ownership", "prediction", "changed_prediction", "kind", "cell",
    "index", "source", "unscheduled", "memory", "count_missing", "count_bool", "count_negative"])
def test_invalid_preflight_refuses_before_destination_creation(joined, change, monkeypatch):
    args = list(joined["args"])
    context = joined["context"]
    if change == "ownership":
        context[5]["gene_ownership_sha256"] = "wrong"
    elif change == "prediction":
        context[5]["checked_files"] = []
    elif change == "changed_prediction":
        joined["prediction"].write_text("changed\n")
    elif change == "kind":
        context[6] = "group"
    elif change in {"index", "cell"}:
        context[3][change] = 11 if change == "index" else "p1_c1_r0"
    elif change == "source":
        args[-1] = "wrong"
    elif change == "unscheduled":
        monkeypatch.delenv("SLURM_JOB_ID")
    elif change == "memory":
        monkeypatch.setenv("SLURM_MEM_PER_NODE", "1024")
    elif change == "count_missing":
        context[5]["phylogeny"] = {}
    else:
        context[5]["phylogeny"]["native_pair_rows"] = True if change == "count_bool" else -1
    with pytest.raises(ValueError):
        current.prepare(*args)
    assert not joined["destination"].exists()


def test_native_count_mismatch_retains_actual_partial_output_and_failure(joined):
    joined["context"][5]["phylogeny"]["native_pair_rows"] = 2
    with pytest.raises(ValueError, match="count differs"):
        current.prepare(*joined["args"])
    destination = joined["destination"]
    assert (destination / "preflight.json").is_file()
    assert (destination / "pairs.partial.tsv").read_text() == "A1\tB1\n"
    assert not (destination / "pairs.tsv").exists()
    report = json.loads((destination / "results.json").read_text())
    assert report["status"] == "native12_composed_qfo_conversion_failed_retained"
    assert report["error_type"] == "ValueError" and report["accuracy_evaluated"] is False


def test_postflight_drift_retains_failure_not_successful_conversion(joined, monkeypatch):
    original = current.materialize
    def drift(*args):
        result = original(*args)
        joined["prediction"].write_text("changed after conversion\n")
        return result
    monkeypatch.setattr(current, "materialize", drift)
    with pytest.raises(ValueError):
        current.prepare(*joined["args"])
    report = json.loads((joined["destination"] / "results.json").read_text())
    assert report["status"] == "native12_composed_qfo_conversion_failed_retained"
    assert not (joined["destination"] / "pairs.tsv").exists()


def test_batch_syntax_and_unscheduled_guard():
    batch = Path(current.__file__).parent / "results/native12_composed_conversion_20261008_v1.sh"
    subprocess.run(["bash", "-n", str(batch)], check=True)
    environment = {key: value for key, value in os.environ.items() if not key.startswith("SLURM_")}
    result = subprocess.run(["bash", str(batch)], env=environment, capture_output=True, text=True)
    assert result.returncode != 0
