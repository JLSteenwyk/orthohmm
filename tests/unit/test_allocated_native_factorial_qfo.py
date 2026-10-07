"""Synthetic new handoffs with real conversion and independent endpoint kernels.

The amendment and native scheduler are stubs, not production readiness evidence.
"""

from copy import deepcopy
import gzip
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest

from benchmark_tools import prepare_allocated_native_factorial_qfo_pairs as pairs
from benchmark_tools import run_allocated_native_factorial_qfo_assessment as runner
from benchmark_tools import admit_allocated_native_factorial_qfo_assessment as admission
from benchmark_tools import native_factorial_allocated_execution as contract
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_native_factorial_cost import IDENTITIES
from benchmark_tools.validate_native_factorial_outputs import Evidence
from benchmark_tools.validate_qfo_native_assessment import AXES
from tests.unit.test_prepare_native_factorial_qfo_pairs import conversion_fixture
from tests.unit.test_run_native_factorial_qfo_assessment import (
    short_root, put, write_native_endpoints, stage_fixture as historical_stage)


def review_fixture(index=10):
    request_ref = dict(path="/request", bytes=1, sha256="a" * 64)
    request = dict(plan=dict(path="/plan", bytes=1, sha256="b" * 64),
        amendment=dict(path="/amendment", bytes=1, sha256="c" * 64), job_id=999)
    run = dict(index=index, dataset="qfo_corrected", cell=IDENTITIES[index][1], repeat=0)
    review = dict(schema="allocated_native_factorial_terminal_review_v1", request=request_ref,
        plan=request["plan"], amendment=request["amendment"], job_id=999, index=index, dataset=run["dataset"],
        cell=run["cell"], repeat=0, status="native_success", scheduler_state="COMPLETED", scheduler_exit_code="0:0",
        terminal_reviewed=True, native_outputs_validated=True, primary_resources_replayed=True,
        shared_host_resources_reviewed=True, execution_scope=pairs.SCOPE, resource_scopes=pairs.SCOPES,
        uncontended_timing=False, automatic_retry=False, accuracy_evaluated=False,
        scientific_timings_admitted=False, publication_ready=False,
        reviews=dict.fromkeys(("runtime", "resources", "environment", "outputs_or_failure")))
    return review, request_ref, request, run


@pytest.mark.parametrize("index", [10, 11, 12])
def test_new_conversion_route_retains_scientific_semantics(index):
    review, ref, request, run = review_fixture(index)
    assert pairs.admit_conversion(review, ref, request, run) == ("group" if index == 11 else "native")


@pytest.mark.parametrize("key,value", [("schema", "native_factorial_terminal_review_v1"),
    ("amendment", {}), ("request", {}), ("plan", {}), ("job_id", 998), ("index", 9),
    ("dataset", "orthobench"), ("cell", "p0_c0_r1"), ("repeat", True),
    ("status", "native_failure_retained"), ("scheduler_state", "RUNNING"), ("scheduler_exit_code", "1:0"),
    ("terminal_reviewed", False), ("native_outputs_validated", False), ("primary_resources_replayed", False),
    ("shared_host_resources_reviewed", False), ("execution_scope", "isolated"), ("resource_scopes", {}),
    ("uncontended_timing", True), ("automatic_retry", True), ("accuracy_evaluated", True),
    ("scientific_timings_admitted", True), ("publication_ready", True), ("reviews", {})])
def test_wrong_terminal_route_or_claim_is_refused(key, value):
    review, ref, request, run = review_fixture()
    review[key] = value
    with pytest.raises(ValueError):
        pairs.admit_conversion(review, ref, request, run)


@pytest.mark.parametrize("index", [6, 7, 8, 9])
def test_historical_native_identity_is_not_a_new_conversion(index):
    review, ref, request, run = review_fixture(index)
    with pytest.raises(ValueError):
        pairs.admit_conversion(review, ref, request, run)


@pytest.fixture(params=[10, 11, 12])
def joined(short_root, conversion_fixture, monkeypatch, request):
    root = short_root
    inputs, owners, clusters, native_pairs, mapping = conversion_fixture
    index = request.param
    original = Path(pairs.__file__).parent
    tools = root / "benchmark_tools"
    tools.mkdir()
    for ref in runner.historical_helpers():
        shutil.copyfile(ref["path"], tools / Path(ref["path"]).name)
    for name in ("prepare_native_factorial_qfo_pairs.py", "prepare_qfo_corrected_group_pairs.py",
        "prepare_qfo_factorial_pairs.py", "qfo_filter_pairs.py", "simulation_method_outputs.py",
        "validate_native_factorial_outputs.py", "review_native_factorial_attempt.py",
        "review_allocated_native_factorial_attempt.py", "validate_allocated_native_factorial_outputs.py",
        "prepare_allocated_native_factorial_qfo_pairs.py", "run_allocated_native_factorial_qfo_assessment.py",
        "admit_allocated_native_factorial_qfo_assessment.py"):
        shutil.copyfile(original / name, tools / name)
    for module in (pairs, runner, admission):
        monkeypatch.setattr(module, "__file__", str(tools / Path(module.__file__).name))
    (root / "qfo_benchmark").mkdir()
    shutil.copyfile(original.parent / "qfo_benchmark/og_to_pairwise.py", root / "qfo_benchmark/og_to_pairwise.py")
    output = root / "native"
    shutil.copytree(inputs, output / "input")
    cluster_ref = put(output / "native/orthohmm_working_res/orthohmm_edges_clustered.txt", clusters.read_text())
    native_ref = put(output / "native/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv",
        native_pairs.read_text() + "sp|A1|AA\ts0\tsp|C1|CA\ts2\n")
    ready_ref = put(output / "measurement/ready.json", dict(synthetic_fixture=True))
    run = dict(index=index, dataset="qfo_corrected", cell=IDENTITIES[index][1], repeat=0,
        output_root=str(output), genes=5, proteomes=4, inputs=[record(p) for p in sorted(inputs.glob("*.fasta"))],
        native_order=[p.name for p in sorted(inputs.glob("*.fasta"))])
    baseline_ref = put(root / "baseline.json", {})
    plan = dict(panel_root=str(root / "panel"), baseline=baseline_ref, helper_sources=[], evidence=[], runs=[run]*13)
    plan_ref = put(root / "plan.json", plan)
    execution = dict(historical_plan=plan_ref, historical_prefix=[{}]*10, new_sources=contract.sources())
    amendment_ref = put(root / "amendment.json", execution)
    request_value = dict(schema="allocated_native_factorial_request_v1", execution_authorized=True,
        job_id=999, index=index, plan=plan_ref, amendment=amendment_ref, history=[{}]*index,
        scheduler_command=str(contract.SCRIPT), allocation_cwd=str(contract.ROOT), automatic_retry=False)
    request_ref = put(root / "request.json", request_value)
    _, digest, counts = pairs.universe(dict(run, input_directory=str(output / "input")), Evidence())
    output_review = dict(schema="allocated_native_factorial_output_review_v1", native_outputs_validated=True,
        source=record(tools / "validate_allocated_native_factorial_outputs.py"),
        semantic_validator_source=record(tools / "validate_native_factorial_outputs.py"),
        request=request_ref, plan=plan_ref, amendment=amendment_ref, index=index, cell=run["cell"], job_id=999,
        gene_ownership_sha256=digest, per_species_counts=counts, checked_files=[cluster_ref, native_ref], evidence=[],
        phylogeny=dict(native_pair_rows=2), execution_scope=pairs.SCOPE, native_cpu_ids=list(range(52,84)),
        allocated_ready=ready_ref)
    output_ref = put(root / "output_review.json", output_review)
    scheduler_ref = put(root / "scheduler.json", {})
    review, _, _, _ = review_fixture(index)
    review.update(request=request_ref, plan=plan_ref, amendment=amendment_ref,
        source=record(tools / "review_allocated_native_factorial_attempt.py"),
        common_reviewer_source=record(tools / "review_native_factorial_attempt.py"), scheduler=scheduler_ref,
        native_cpu_ids=output_review["native_cpu_ids"], evidence=[],
        reviews=dict(runtime=scheduler_ref, resources=scheduler_ref, environment=scheduler_ref, outputs_or_failure=output_ref))
    review_ref = put(root / "review.json", review)
    prepared_ref = put(root / "prepared.json", dict(input_fastas=run["inputs"]))
    pipeline = root / "pipeline"
    main = put(pipeline / "main.nf", "// Synthetic pinned pipeline\n")
    fas = put(pipeline / "fas_benchmark.py", "MAX_PAIRS_COMPUTE = 9000\n")
    mapping_path = pipeline / "reference_data/2020/mapping.json.gz"
    mapping_path.parent.mkdir(parents=True)
    with gzip.open(mapping_path, "wt") as stream:
        json.dump(dict(mapping=mapping), stream)
    templates = [put(pipeline / "reference_data/data" / (challenge + ".json"),
        dict(type="aggregation", challenge_ids=[challenge], datalink=dict(inline_data=dict(
            visualization=dict(x_axis=axes[0], y_axis=axes[1], type="2D-plot"))))) for challenge, axes in AXES.items()]
    swiss = put(pipeline / "reference_data/2020/ReconciledTrees_SwissTrees.drw",
        "ReconciledTrees['X'] := RecTreeCase('X',...)\n")
    manifest = dict(status="local_qfo_assessment_environment_frozen", accuracy_evaluated=False,
        pipeline=str(pipeline), source=fas, execution_config=main, singularity_config=main,
        pipeline_files=[main, fas], reference_files=[record(mapping_path), swiss, *templates],
        java_files=[], images=[], executables=[], singularity_support=[], environment_overrides=dict(NXF_OFFLINE="true"))
    env_ref = put(tools / "results/qfo_assessment_environment_20260917.json", manifest)
    monkeypatch.setattr(pairs, "ROOT", root)
    monkeypatch.setattr(runner, "ROOT", root)
    monkeypatch.setattr(pairs, "amendment", lambda ref: (execution, plan))
    monkeypatch.setattr(pairs, "FIXED_INPUTS", {"qfo_preparation": ("prepared.json", prepared_ref["sha256"])})
    monkeypatch.setattr(pairs, "ENV_SHA", env_ref["sha256"])
    monkeypatch.setattr(runner, "ENV_SHA", env_ref["sha256"])
    monkeypatch.setattr(pairs, "verify_terminal", lambda job: dict(source="live_controller",
        verified=dict(fields=dict(JobState="COMPLETED", ExitCode="0:0", Comment=request_ref["sha256"]))))
    monkeypatch.setattr(runner, "accounting", lambda job: ("synthetic conversion accounting",
        dict(JobIDRaw="777", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2")))
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setenv("SLURM_JOB_ID", "777")
    return root, request_ref, review_ref, run, manifest


def converted(joined):
    root, request_ref, review_ref, run, manifest = joined
    ref = pairs.prepare(request_ref, review_ref, root / "conversion")
    return root, ref, json.loads(Path(ref["path"]).read_text()), manifest


def test_real_converter_retains_mapping_and_all_input_denominator(joined):
    root, ref, stage, _ = converted(joined)
    assert stage["schema"] == "allocated_native_factorial_qfo_conversion_v1"
    assert stage["status"] == "allocated_native_factorial_qfo_pairs_prepared_unscored"
    assert stage["retained_pairs"] == (5 if stage["native_index"] == 11 else 2)
    assert stage["pair_coverage"]["input_accessions"] == 5
    assert stage["pair_coverage"]["fraction_inputs_in_any_pair"] == (.8 if stage["native_index"] == 11 else .6)
    assert stage["pairs"]["sha256"] == stage["filtered_pairs"]["sha256"]
    assert stage["accuracy_evaluated"] is stage["next_identity_authorized"] is stage["publication_ready"] is False
    assert stage["amendment"] == json.loads(Path(joined[1]["path"]).read_text())["amendment"]
    with pytest.raises(FileExistsError):
        pairs.prepare(joined[1], joined[2], root / "conversion")


@pytest.mark.parametrize("key,value", [("schema", "native_factorial_output_review_v1"),
    ("source", {}), ("semantic_validator_source", {}), ("amendment", {}), ("request", {}),
    ("plan", {}), ("job_id", 998), ("native_cpu_ids", [0]), ("allocated_ready", {}),
    ("native_outputs_validated", False), ("execution_scope", "isolated")])
def test_new_native_output_gate_cannot_consume_old_or_mismatched_receipts(joined, key, value):
    root, request_ref, review_ref, _, _ = joined
    review = json.loads(Path(review_ref["path"]).read_text())
    path = Path(review["reviews"]["outputs_or_failure"]["path"])
    outputs = json.loads(path.read_text())
    outputs[key] = value
    review["reviews"]["outputs_or_failure"] = put(path, outputs)
    review_ref = put(Path(review_ref["path"]), review)
    with pytest.raises(ValueError, match="output/source/placement"):
        pairs.native_binding(request_ref, review_ref)


@pytest.mark.parametrize("problem", ["scheduler", "comment", "owner", "mapping_loss", "source", "failed_conversion"])
def test_converter_failure_boundary(joined, monkeypatch, problem):
    root, request_ref, review_ref, _, _ = joined
    if problem in {"scheduler", "comment"}:
        monkeypatch.setattr(pairs, "verify_terminal", lambda job: dict(source="live_controller",
            verified=dict(fields=dict(JobState="FAILED" if problem == "scheduler" else "COMPLETED",
                ExitCode="0:0", Comment="wrong"))))
    elif problem in {"owner", "source"}:
        review = json.loads(Path(review_ref["path"]).read_text())
        if problem == "source":
            review["source"] = review["common_reviewer_source"]
        else:
            output_path = Path(review["reviews"]["outputs_or_failure"]["path"])
            outputs = json.loads(output_path.read_text())
            outputs["gene_ownership_sha256"] = "bad"
            review["reviews"]["outputs_or_failure"] = put(output_path, outputs)
        review_ref = put(Path(review_ref["path"]), review)
    elif problem == "mapping_loss":
        original = pairs.filter_pairs
        def lose(*args):
            total, retained = original(*args)
            return total, retained-1
        monkeypatch.setattr(pairs, "filter_pairs", lose)
    else:
        def fail(*args):
            raise ValueError("synthetic converter failure")
        monkeypatch.setattr(pairs, "convert", fail)
    with pytest.raises(ValueError):
        pairs.prepare(request_ref, review_ref, root / "conversion")
    if problem in {"mapping_loss", "failed_conversion"}:
        result = json.loads((root / "conversion/results.json").read_text())
        assert result["status"] == "allocated_native_factorial_qfo_conversion_failed_retained"
        assert result["accuracy_evaluated"] is False
        assert not (root / "conversion/pairs.tsv").exists()
    else:
        assert not (root / "conversion").exists()


def stage_fixture(index=10, empty=False):
    stage, scheduler = historical_stage(index, empty)
    stage.update(schema="allocated_native_factorial_qfo_conversion_v1",
        status="allocated_native_factorial_qfo_pairs_prepared_unscored",
        conversion_started_monotonic_ns=10, conversion_finished_monotonic_ns=20)
    return stage, scheduler


@pytest.mark.parametrize("index", [10,11,12])
@pytest.mark.parametrize("empty", [False,True])
def test_new_stage_identity_and_empty_counts(index, empty):
    stage, scheduler = stage_fixture(index, empty)
    runner.validate_stage(stage, scheduler)


@pytest.mark.parametrize("key,value", [("schema", "full_native_factorial_qfo_conversion_v1"),
    ("status", "full_native_factorial_qfo_pairs_prepared_unscored"), ("native_index", 9),
    ("native_index", True), ("cell", "p0_c0_r1"), ("participant", "cached"), ("conversion_kind", "group"),
    ("accuracy_evaluated", True), ("automatic_retry", True), ("empty_predictions", True),
    ("removed_mapping_pairs", 1), ("total_pairs", True), ("retained_pairs", 3),
    ("filtered_pairs", dict(bytes=1, sha256="a"*64)), ("conversion_started_monotonic_ns", True),
    ("conversion_finished_monotonic_ns", 9)])
def test_new_stage_rejects_historical_schema_or_bad_semantics(key, value):
    stage, scheduler = stage_fixture()
    stage[key] = value
    with pytest.raises(ValueError):
        runner.validate_stage(stage, scheduler)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
    ("NodeList", "dgx"), ("AllocCPUS", "8"), ("JobIDRaw", "888")])
def test_successful_conversion_requires_actual_own_two_cpu_scheduler(key, value):
    stage, scheduler = stage_fixture()
    scheduler[key] = value
    with pytest.raises(ValueError):
        runner.validate_stage(stage, scheduler)


def test_readonly_assessment_preflight_uses_new_namespace_and_same_six_endpoint_command(joined):
    root, ref, stage, _ = converted(joined)
    report = runner.run(root, ref, "777", check_only=True)
    assert report["schema"] == "allocated_native_factorial_qfo_execution_v1"
    assert report["work"].endswith(f"/aq{stage['native_index']:02d}")
    assert report["results"].endswith(f"/allocated_native_{stage['native_index']:02d}")
    assert report["command"][report["command"].index("--challenges_ids")+1] == "GO EC VGNC SwissTrees TreeFam-A FAS"
    assert report["command"][report["command"].index("--event_year")+1] == "2020"
    assert "-resume" not in report["command"] and report["accuracy_admitted"] is False
    assert report["fas_protocol"]["newly_computed_pair_cap"] == 9000
    assert not Path(report["cwd"]).exists()


@pytest.mark.parametrize("key,value", [("amendment", {}), ("native_cpu_ids", [0]), ("allocated_ready", {}),
    ("conversion_kernel_source", {}), ("source", {}), ("gene_ownership_sha256", "bad"),
    ("input_fastas", []), ("native_job_id", 998), ("native_input", {}), ("mapping", {})])
def test_assessment_binding_refuses_resealed_stage_changes(joined, key, value):
    root, ref, stage, _ = converted(joined)
    stage[key] = value
    ref = put(Path(ref["path"]), stage)
    with pytest.raises(ValueError):
        runner.prepare(root, ref, "777")


def execute_fixture(converted_fixture, monkeypatch, problem=None):
    root, ref, stage, manifest = converted_fixture
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setenv("SLURM_JOB_ID", "778")
    def execute(command, **kwargs):
        assert "-resume" not in command and kwargs["env"]["NXF_OFFLINE"] == "true"
        # Two sample pairs must belong to the submitted fixture predictions.
        write_native_endpoints(Path(command[command.index("--results_dir")+1]), stage, manifest)
        if problem == "tamper":
            Path(stage["filtered_pairs"]["path"]).write_text("changed")
        if problem == "interrupt":
            raise KeyboardInterrupt()
        return SimpleNamespace(returncode=1 if problem == "exit" else 0)
    monkeypatch.setattr(runner.subprocess, "run", execute)
    return runner.run(root, ref, "777")


@pytest.mark.parametrize("problem", ["exit", "tamper", "interrupt"])
def test_assessment_failure_retained_without_accuracy_claim(joined, monkeypatch, problem):
    fixture = converted(joined)
    root, ref, stage, _ = fixture
    with pytest.raises(KeyboardInterrupt if problem == "interrupt" else ValueError):
        execute_fixture(fixture, monkeypatch, problem)
    path = root / "benchmarks/results/allocated_native_qfo_assessment_v1" / stage["cell"]
    report = json.loads((path / "results.json").read_text())
    assert report["status"] == "failed" and report["accuracy_admitted"] is False
    assert (path / "preflight.json").exists()


def assessment_scheduler(monkeypatch, state="COMPLETED"):
    monkeypatch.setattr(admission, "accounting", lambda job: ("synthetic assessment accounting",
        dict(JobIDRaw="778", State=state, ExitCode="0:0", NodeList="bizon", AllocCPUS="8")))


def destination(root):
    return root / "benchmarks/results/allocated_native_qfo_admission_v1/test_admission"


def test_joined_independent_endpoints_admitted_only_after_validation(joined, monkeypatch):
    fixture = converted(joined)
    root, ref, stage, _ = fixture
    execution = execute_fixture(fixture, monkeypatch)
    assessment_scheduler(monkeypatch)
    result = admission.admit(root, ref, "777", "778", destination(root))
    assert result["status"] == "allocated_native_factorial_qfo_assessment_admitted"
    assert set(result["assessment"]["endpoints"]) == set(AXES)
    assert result["assessment"]["secondary_six_metric_mean"] == .5
    assert result["fas_sample"]["sample_membership_verified"] is True
    assert result["amendment"] == stage["amendment"]
    assert result["accuracy_admitted"] is True and result["publication_ready"] is False
    assert execution["accuracy_admitted"] is False
    with pytest.raises(FileExistsError):
        admission.admit(root, ref, "777", "778", destination(root))


@pytest.mark.parametrize("problem", ["scheduler", "failed", "command", "preflight", "outputs", "trace", "metrics",
    "schema", "source", "fas_membership", "fas_mean"])
def test_independent_admission_rejects_changed_or_incomplete_execution(joined, monkeypatch, problem):
    fixture = converted(joined)
    root, ref, _, _ = fixture
    execution = execute_fixture(fixture, monkeypatch)
    assessment_scheduler(monkeypatch, "RUNNING" if problem == "scheduler" else "COMPLETED")
    output = Path(execution["cwd"])
    if problem in {"failed", "command", "outputs", "schema", "source"}:
        if problem == "failed": execution["status"] = "failed"
        elif problem == "command": execution["command"].append("-resume")
        elif problem == "outputs": execution["outputs"] = []
        elif problem == "schema": execution["schema"] = "full_native_factorial_qfo_execution_v1"
        else: execution["source"] = ref
        put(output / "results.json", execution)
    elif problem == "preflight":
        preflight = json.loads((output / "preflight.json").read_text())
        preflight["amendment"] = {}
        put(output / "preflight.json", preflight)
    elif problem != "scheduler":
        results = Path(execution["results"])
        if problem == "trace":
            path = results / "stats/trace_fixture.txt"
            path.write_text(path.read_text().replace("COMPLETED", "CACHED", 1))
        elif problem == "metrics":
            path = results / "assessment_out/Assessment_datasets.json"
            rows = json.loads(path.read_text())
            rows.pop()
            put(path, rows)
        else:
            path = results / "results/FAS/sample_raw.txt.gz"
            with gzip.open(path, "wt") as stream:
                stream.write("Acc1\tAcc2\tFAS\nA1\tB1\t0.4\nA1\t" +
                    ("D1\t0.6\n" if problem == "fas_membership" else "C1\t0.7\n"))
        execution["outputs"] = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
        put(output / "results.json", execution)
    with pytest.raises(ValueError):
        admission.admit(root, ref, "777", "778", destination(root))
    report = json.loads((destination(root) / "results.json").read_text())
    assert report["status"] == "allocated_native_factorial_qfo_admission_failed_retained"
    assert report["accuracy_admitted"] is False
