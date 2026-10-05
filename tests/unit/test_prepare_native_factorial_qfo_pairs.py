"""Full-native conversion gates/semantics, not new scientific QfO inference."""

import gzip
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import prepare_native_factorial_qfo_pairs as pairs
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.validate_native_factorial_outputs import Evidence


def admission_fixture(index=7):
    request_ref = dict(path="/request", bytes=1, sha256="a" * 64)
    request = dict(plan=dict(path="/plan", bytes=1, sha256="b" * 64), job_id=999)
    dataset, cell = pairs.IDENTITIES[index]
    run = dict(index=index, dataset=dataset, cell=cell, repeat=0)
    review = dict(schema="native_factorial_terminal_review_v1", request=request_ref, plan=request["plan"],
        job_id=999, index=index, dataset=dataset, cell=cell, repeat=0, status="native_success",
        scheduler_state="COMPLETED", scheduler_exit_code="0:0", terminal_reviewed=True,
        native_outputs_validated=True, primary_resources_replayed=True, shared_host_resources_reviewed=True,
        execution_scope=pairs.SCOPE, resource_scopes=pairs.SCOPES, uncontended_timing=False, automatic_retry=False,
        reviews=dict.fromkeys(("runtime", "resources", "environment", "outputs_or_failure")))
    return review, request_ref, request, run


@pytest.mark.parametrize("index", range(6, 13))
def test_all_seven_bound_full_native_identities(index):
    review, ref, request, run = admission_fixture(index)
    kind = pairs.admit_conversion(review, ref, request, run)
    assert kind == ("native" if run["cell"].endswith("r1") else "group")


@pytest.mark.parametrize("key,value", [("dataset", "orthobench"), ("index", 5), ("index", True),
    ("index", 13), ("cell", "p1_c0_r0"), ("cell", "p0_c0_r0")])
def test_other_or_mismatched_identity_refused(key, value):
    review, ref, request, run = admission_fixture()
    run[key] = value
    with pytest.raises(ValueError):
        pairs.admit_conversion(review, ref, request, run)


@pytest.mark.parametrize("key,value", [("schema", "cached_replay_admission"), ("status", "native_failure_retained"),
    ("job_id", 1), ("index", 6), ("cell", "p0_c0_r0"), ("request", {}), ("plan", {}),
    ("scheduler_state", "RUNNING"), ("scheduler_exit_code", "1:0"), ("terminal_reviewed", False),
    ("native_outputs_validated", False), ("primary_resources_replayed", False),
    ("shared_host_resources_reviewed", False), ("execution_scope", "isolated"), ("resource_scopes", {}),
    ("uncontended_timing", True), ("automatic_retry", True), ("reviews", {})])
def test_unreviewed_or_failed_native_result_not_converted(key, value):
    review, ref, request, run = admission_fixture()
    review[key] = value
    with pytest.raises(ValueError, match="bound successful"):
        pairs.admit_conversion(review, ref, request, run)


@pytest.fixture
def conversion_fixture(tmp_path):
    inputs = tmp_path / "input"
    inputs.mkdir()
    owners = {"sp|A1|AA": "s0", "sp|A2|AB": "s0", "sp|B1|BA": "s1", "sp|C1|CA": "s2", "sp|U1|UA": "s3"}
    for species in sorted(set(owners.values())):
        (inputs / f"{species}.fasta").write_text("".join(
            f">{gene}\nAAAA\n" for gene in owners if owners[gene] == species))
    cluster = tmp_path / "clusters.txt"
    cluster.write_text(" ".join(list(owners)[:4]) + "\n" + list(owners)[4] + "\n")
    native = tmp_path / "native.tsv"
    native.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\nsp|A1|AA\ts0\tsp|B1|BA\ts1\n")
    mapping = {pairs._strip_to_uniprot(g): i for i, g in enumerate(owners)}
    return inputs, owners, cluster, native, mapping


def test_group_and_native_semantics_remain_different(conversion_fixture, tmp_path):
    inputs, owners, cluster, native, mapping = conversion_fixture
    group_output, native_output = tmp_path / "group_pairs.tsv", tmp_path / "native_pairs.tsv"
    group_count, command = pairs.convert("group", cluster, group_output, inputs, owners)
    native_count, native_command = pairs.convert("native", native, native_output, inputs, owners, 1)
    assert group_count == 5 and native_count == 1
    assert command[1].endswith("og_to_pairwise.py") and native_command is None
    assert native_output.read_text() == "A1\tB1\n"
    assert "A1\tA2\n" not in group_output.read_text()
    normalized = pairs.normalize_owners(owners, mapping)
    coverage = pairs.pair_coverage(group_output, normalized)
    assert coverage["accessions_in_any_pair"] == 4 and coverage["input_accessions"] == 5
    assert coverage["fraction_inputs_in_any_pair"] == .8
    assert pairs.pair_coverage(native_output, normalized)["fraction_inputs_in_any_pair"] == .4


def test_normalization_is_injective_for_entire_input_not_only_pairs():
    with pytest.raises(ValueError, match="Noninjective"):
        pairs.normalize_owners({"sp|A1|AA": "s0", "A1": "s1"}, {"A1": 1})
    with pytest.raises(ValueError, match="absent from frozen"):
        pairs.normalize_owners({"sp|U1|UA": "s3"}, {"A1": 1})


def test_complete_cluster_universe_required(conversion_fixture, tmp_path):
    inputs, owners, cluster, _, _ = conversion_fixture
    cluster.write_text(" ".join(list(owners)[:4]) + "\n")
    with pytest.raises(ValueError, match="Incomplete candidate"):
        pairs.convert("group", cluster, tmp_path / "out.tsv", inputs, owners)


def test_empty_predictions_are_not_imputed_or_lost(conversion_fixture, tmp_path):
    inputs, owners, cluster, native, mapping = conversion_fixture
    cluster.write_text("\n".join(owners) + "\n")
    native.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\n")
    normalized = pairs.normalize_owners(owners, mapping)
    for kind, source, count in (("group", cluster, None), ("native", native, 0)):
        output = tmp_path / f"{kind}_empty_converted.tsv"
        actual, _ = pairs.convert(kind, source, output, inputs, owners, count)
        assert actual == 0 and output.read_bytes() == b""
        assert pairs.pair_coverage(output, normalized)["fraction_inputs_in_any_pair"] == 0


@pytest.mark.parametrize("value", [True, -1, 1.0, None])
def test_bad_admitted_native_count_refused(conversion_fixture, tmp_path, value):
    inputs, owners, _, native, _ = conversion_fixture
    with pytest.raises(ValueError, match="admitted native pair count"):
        pairs.convert("native", native, tmp_path / "out.tsv", inputs, owners, value)


def test_native_count_mismatch_preserves_failed_output(conversion_fixture, tmp_path):
    inputs, owners, _, native, _ = conversion_fixture
    output = tmp_path / "out.tsv"
    with pytest.raises(ValueError, match="count differs"):
        pairs.convert("native", native, output, inputs, owners, 2)
    assert output.read_text() == "A1\tB1\n"


@pytest.mark.parametrize("line", ["B1\tA1\n", "A1\tA1\n", "A1\tZ1\n", "A1\tA2\n", "A1 B1\n"])
def test_independent_coverage_rejects_invalid_pairs(conversion_fixture, tmp_path, line):
    _, owners, _, _, mapping = conversion_fixture
    path = tmp_path / "pairs.tsv"
    path.write_text(line)
    with pytest.raises(ValueError):
        pairs.pair_coverage(path, pairs.normalize_owners(owners, mapping))


@pytest.fixture(params=[6, 7])
def joined_conversion(conversion_fixture, tmp_path, monkeypatch, request):
    inputs, owners, clusters, native, mapping = conversion_fixture
    original = Path(pairs.__file__).parent
    index = request.param
    _, _, _, run = admission_fixture(index)
    output = tmp_path / "native_attempt"
    shutil.copytree(inputs, output / "input")
    run.update(output_root=str(output), input_directory=str(inputs), genes=5, proteomes=4,
        inputs=[record(path) for path in sorted(inputs.glob("*.fasta"))],
        native_order=[path.name for path in sorted(inputs.glob("*.fasta"))])
    working = output / "native/orthohmm_working_res"
    phylogeny = output / "native/orthohmm_phylogeny"
    working.mkdir(parents=True)
    phylogeny.mkdir(parents=True)
    clustered = working / "orthohmm_edges_clustered.txt"
    native_pairs = phylogeny / "orthohmm_pairwise_orthologs.tsv"
    shutil.copyfile(clusters, clustered)
    shutil.copyfile(native, native_pairs)
    tools = tmp_path / "benchmark_tools"
    tools.mkdir()
    for name in ("prepare_qfo_corrected_group_pairs.py", "prepare_qfo_factorial_pairs.py", "qfo_filter_pairs.py",
                 "simulation_method_outputs.py", "validate_native_factorial_outputs.py", "review_native_factorial_attempt.py"):
        shutil.copyfile(original / name, tools / name)
    (tmp_path / "qfo_benchmark").mkdir()
    shutil.copyfile(original.parent / "qfo_benchmark/og_to_pairwise.py", tmp_path / "qfo_benchmark/og_to_pairwise.py")
    baseline_path = tmp_path / "baseline.json"
    baseline_path.write_text("{}")
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(dict(panel_root=str(tmp_path / "panel"), baseline=record(baseline_path),
        helper_sources=[record(tools / "qfo_filter_pairs.py")])))
    request_path = tmp_path / "request.json"
    request_path.write_text(json.dumps(dict(schema="native_factorial_cost_request_v1", execution_authorized=True,
        job_id=999, index=index, history=[{}] * index, plan=record(plan_path),
        scheduler_command=str(original / "run_native_factorial_cost.sh"), allocation_cwd=str(original.parent))))
    request_ref = record(request_path)
    _, digest, counts = pairs.universe(dict(run, input_directory=str(output / "input")), Evidence())
    output_path = tmp_path / "output_review.json"
    output_review = dict(native_outputs_validated=True, request=request_ref, plan=record(plan_path), job_id=999,
        index=index, cell=run["cell"], gene_ownership_sha256=digest, per_species_counts=counts,
        checked_files=[record(clustered), record(native_pairs)], evidence=[], phylogeny=dict(native_pair_rows=1))
    output_path.write_text(json.dumps(output_review))
    scheduler_path = tmp_path / "scheduler.json"
    scheduler_path.write_text("{}")
    review, _, _, _ = admission_fixture(index)
    review.update(request=request_ref, plan=record(plan_path), source=record(tools / "review_native_factorial_attempt.py"),
        scheduler=record(scheduler_path), reviews=dict(runtime=record(scheduler_path), resources=record(scheduler_path),
            environment=record(scheduler_path), outputs_or_failure=record(output_path)))
    review_path = tmp_path / "review.json"
    review_path.write_text(json.dumps(review))
    prepared_path = tmp_path / "prepared.json"
    prepared_path.write_text(json.dumps(dict(input_fastas=run["inputs"])))
    mapping_path = tmp_path / "mapping.json.gz"
    with gzip.open(mapping_path, "wt") as stream:
        json.dump(dict(mapping=mapping), stream)
    (tools / "results").mkdir()
    environment_path = tools / "results/qfo_assessment_environment_20260917.json"
    environment_path.write_text(json.dumps(dict(reference_files=[record(mapping_path)])))
    monkeypatch.setattr(pairs, "ROOT", tmp_path)
    monkeypatch.setattr(pairs, "FIXED_INPUTS", {"qfo_preparation": ("prepared.json", record(prepared_path)["sha256"])})
    monkeypatch.setattr(pairs, "ENV_SHA", record(environment_path)["sha256"])
    monkeypatch.setattr(pairs, "validate_plan", lambda plan: [run] * 13)
    monkeypatch.setattr(pairs, "verify_terminal", lambda job: dict(source="fresh_accounting_after_controller_expiry",
        verified=dict(State="COMPLETED", ExitCode="0:0")))
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setenv("SLURM_JOB_ID", "7777")
    return request_ref, record(review_path), tmp_path / "conversion", index


def test_joined_conversion_no_overwrite_and_no_score(joined_conversion):
    ref, review, output, index = joined_conversion
    result_ref = pairs.prepare(ref, review, output)
    result = json.loads(Path(result_ref["path"]).read_text())
    assert result["status"] == "full_native_factorial_qfo_pairs_prepared_unscored"
    assert result["native_index"] == index and result["native_job_id"] == 999 and result["job_id"] == "7777"
    assert result["total_pairs"] == result["retained_pairs"] == (5 if index == 6 else 1)
    assert result["pair_coverage"]["accessions_in_any_pair"] == (4 if index == 6 else 2)
    assert result["removed_mapping_pairs"] == 0
    assert result["accuracy_evaluated"] is result["publication_ready"] is result["next_identity_authorized"] is False
    assert result["native_inference_reexecuted"] is False
    assert not result["participant"].startswith("ohmm_qfo_corrected_factorial_")
    assert not (output / "pairs.partial.tsv").exists()
    with pytest.raises(FileExistsError):
        pairs.prepare(ref, review, output)


def test_converter_requires_own_scheduled_allocation(joined_conversion, monkeypatch):
    ref, review, output, _ = joined_conversion
    monkeypatch.delenv("SLURM_JOB_ID")
    with pytest.raises(ValueError, match="scheduled two-CPU"):
        pairs.prepare(ref, review, output)
    assert not output.exists()


def test_changed_actual_native_outcome_not_converted(joined_conversion, monkeypatch):
    ref, review, output, _ = joined_conversion
    monkeypatch.setattr(pairs, "verify_terminal", lambda job: dict(source="fresh_accounting_after_controller_expiry",
        verified=dict(State="FAILED", ExitCode="1:0")))
    with pytest.raises(ValueError, match="successful terminal scheduler"):
        pairs.prepare(ref, review, output)
    assert not output.exists()


def test_real_conversion_failure_is_retained_not_promoted(joined_conversion, monkeypatch):
    ref, review, output, _ = joined_conversion
    def fail(*args, **kwargs):
        raise RuntimeError("synthetic converter failure")
    monkeypatch.setattr(pairs, "convert", fail)
    with pytest.raises(RuntimeError):
        pairs.prepare(ref, review, output)
    result = json.loads((output / "results.json").read_text())
    assert result["status"] == "full_native_factorial_qfo_conversion_failed_retained"
    assert result["accuracy_evaluated"] is False
    assert not (output / "pairs.tsv").exists()
    assert (output / "preflight.json").exists()


def test_mapping_loss_is_failure_not_silent_filtering(joined_conversion, monkeypatch):
    ref, review, output, _ = joined_conversion
    original_filter = pairs.filter_pairs
    def lose(source, dest, mapping):
        total, retained = original_filter(source, dest, mapping)
        return total, retained - 1
    monkeypatch.setattr(pairs, "filter_pairs", lose)
    with pytest.raises(ValueError, match="mapping loses"):
        pairs.prepare(ref, review, output)
    result = json.loads((output / "results.json").read_text())
    assert result["status"] == "full_native_factorial_qfo_conversion_failed_retained"
    assert (output / "pairs.partial.tsv").exists()
    assert not (output / "pairs.tsv").exists()
