"""Recovered accuracy handoffs, never new inference or timing admission."""

from copy import deepcopy
import gzip
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import admit_measurement_failed_native_qfo_assessment as admission
from benchmark_tools import run_measurement_failed_native_qfo_assessment as runner
from benchmark_tools import run_native_factorial_qfo_assessment as original
from benchmark_tools.prepare_measurement_failed_native_qfo_pairs import admit_recovery
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_run_native_factorial_qfo_assessment import (
    AXES, joined, put, short_root, stage_fixture, write_native_endpoints,  # noqa: F401
)
from tests.unit.test_prepare_measurement_failed_native_qfo_pairs import recovered_fixture


def recovered_stage(index=7, empty=False):
    stage, scheduler = stage_fixture(index, empty)
    stage.update(schema="measurement_failed_native_qfo_conversion_v1",
        status="measurement_failed_native_qfo_pairs_prepared_unscored",
        participant="ohmm_qfo_recovered_native_" + stage["cell"], resources=None)
    for key in ("original_native_scheduler_success", "scientific_timings_admitted", "eligible_for_timing_comparison"):
        stage[key] = False
    return stage, scheduler


@pytest.mark.parametrize("index", range(6, 13))
@pytest.mark.parametrize("empty", [False, True])
def test_recovered_identities_and_empty_predictions(index, empty):
    stage, scheduler = recovered_stage(index, empty)
    runner.validate_stage(stage, scheduler)
    with pytest.raises(ValueError):
        original.validate_stage(stage, scheduler)


@pytest.mark.parametrize("key,value", [
    ("schema", "full_native_factorial_qfo_conversion_v1"), ("status", "running"), ("native_index", True),
    ("native_index", 5), ("cell", "p1_c0_r0"), ("participant", "ohmm_qfo_full_native_p0_c0_r1"),
    ("conversion_kind", "group"), ("semantics", "root HOG cliques"), ("accuracy_evaluated", True),
    ("native_inference_reexecuted", True), ("automatic_retry", True), ("next_identity_authorized", True),
    ("publication_ready", True), ("original_native_scheduler_success", True), ("scientific_timings_admitted", True),
    ("eligible_for_timing_comparison", True), ("resources", {}), ("job_id", "778"), ("total_pairs", True),
    ("retained_pairs", 1), ("expected_pairs", 3), ("removed_mapping_pairs", 1), ("empty_predictions", True),
    ("filtered_pairs", dict(bytes=9, sha256="a" * 64)),
    ("pair_coverage", dict(pair_rows=2, input_accessions=3, accessions_in_any_pair=2, fraction_inputs_in_any_pair=1.)),
])
def test_wrong_scope_or_counts_rejected(key, value):
    stage, scheduler = recovered_stage()
    stage[key] = value
    with pytest.raises(ValueError):
        runner.validate_stage(stage, scheduler)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"), ("NodeList", "dgx"),
    ("AllocCPUS", "8"), ("JobIDRaw", "778")])
def test_wrong_conversion_scheduler_rejected(key, value):
    stage, scheduler = recovered_stage()
    scheduler[key] = value
    with pytest.raises(ValueError):
        runner.validate_stage(stage, scheduler)


@pytest.fixture
def recovered_joined(joined, monkeypatch):
    root, pairs_ref, old, manifest = joined
    stage, scheduler = recovered_stage()
    stage.update({k: old[k] for k in ("pairs", "filtered_pairs", "native_input", "gene_ownership_sha256",
        "input_fastas", "plan", "request", "environment_manifest", "mapping")})
    stage["source"] = record(Path(runner.__file__).with_name("prepare_measurement_failed_native_qfo_pairs.py"))
    review, _, request, run = recovered_fixture()
    request.update(job_id=999, plan=stage["plan"], index=7)
    run.update(genes=3, inputs=stage["input_fastas"])
    review.update(request=stage["request"], plan=stage["plan"], job_id=999)
    recovery_ref = put(root / "recovery.json", review)
    stage.update(scientific_recovery=recovery_ref, checked_records=[stage["source"]])
    pairs_ref = put(Path(pairs_ref["path"]), stage)
    outputs = dict(gene_ownership_sha256=stage["gene_ownership_sha256"], checked_files=[stage["native_input"]],
        phylogeny=dict(native_pair_rows=2))
    def bind(request_ref, review_ref):
        check_request = json.loads(Path(request_ref["path"]).read_text())
        assert check_request["plan"] == request["plan"]
        reviewed = json.loads(Path(review_ref["path"]).read_text())
        admit_recovery(reviewed, request_ref, request, run)
        return request, reviewed, {}, run, "native", outputs, dict(source="fixture",
            verified=dict(State="FAILED", ExitCode="1:0")), [stage["request"], recovery_ref]
    monkeypatch.setattr(runner, "bind_recovery", bind)
    monkeypatch.setattr(runner, "ENV_SHA", stage["environment_manifest"]["sha256"])
    monkeypatch.setattr(runner, "accounting", lambda job: ("synthetic conversion accounting", deepcopy(scheduler)))
    return root, pairs_ref, stage, manifest


def execute(data, monkeypatch, exit_code=0, tamper=False):
    root, ref, stage, manifest = data
    def invoke(command, **kwargs):
        assert "-resume" not in command
        assert command[command.index("--event_year") + 1] == "2020"
        assert command[command.index("--challenges_ids") + 1] == "GO EC VGNC SwissTrees TreeFam-A FAS"
        assert kwargs["env"]["NXF_OFFLINE"] == "true"
        write_native_endpoints(Path(command[command.index("--results_dir") + 1]), stage, manifest)
        if tamper:
            Path(stage["filtered_pairs"]["path"]).write_text("changed")
        return SimpleNamespace(returncode=exit_code)
    monkeypatch.setattr(runner.subprocess, "run", invoke)
    return runner.run(root, ref, "777")


def destination(root):
    return root / "benchmarks/results/measurement_failed_native_qfo_admission_v1/test_admission"


def assessment_scheduler(monkeypatch, state="COMPLETED"):
    monkeypatch.setattr(admission, "accounting", lambda job: ("synthetic assessment accounting",
        dict(JobIDRaw="778", State=state, ExitCode="0:0", NodeList="bizon", AllocCPUS="8")))


def test_check_only_is_distinct_read_only_and_keeps_failure(recovered_joined):
    root, ref, _, _ = recovered_joined
    report = runner.run(root, ref, "777", check_only=True)
    assert report["status"] == "prepared_unrun" and report["accuracy_admitted"] is False
    assert report["work"].endswith("/mr07") and report["results"].endswith("/recovered_native_07")
    assert report["native_scheduler"]["verified"]["State"] == "FAILED" and report["resources"] is None
    assert report["fas_protocol"]["newly_computed_pair_cap"] == 9000
    assert not Path(report["cwd"]).exists()


@pytest.mark.parametrize("problem", ["source", "input", "mapping", "job", "owner", "output", "native_count", "recovery"])
def test_bound_stage_changes_refused(recovered_joined, problem):
    root, ref, stage, _ = recovered_joined
    modified = deepcopy(stage)
    if problem == "source": modified["source"] = modified["native_input"]
    elif problem == "input": modified["input_fastas"] = []
    elif problem == "mapping": modified["mapping"] = modified["native_input"]
    elif problem == "job": modified["native_job_id"] = 1
    elif problem == "owner": modified["gene_ownership_sha256"] = "bad"
    elif problem == "output": modified["native_input"] = modified["pairs"]
    elif problem == "native_count":
        modified.update(total_pairs=1, retained_pairs=1, expected_pairs=1)
        modified["pair_coverage"]["pair_rows"] = 1
    else:
        reviewed = json.loads(Path(stage["scientific_recovery"]["path"]).read_text())
        reviewed["scheduler_success"] = True
        modified["scientific_recovery"] = put(Path(stage["scientific_recovery"]["path"]), reviewed)
    with pytest.raises(ValueError):
        runner.prepare(root, put(Path(ref["path"]), modified), "777")


@pytest.mark.parametrize("allocation", [None, "2", "32"])
def test_requires_eight_cpu_allocation(recovered_joined, monkeypatch, allocation):
    root, ref, _, _ = recovered_joined
    if allocation is None: monkeypatch.delenv("SLURM_JOB_ID")
    else: monkeypatch.setenv("SLURM_CPUS_PER_TASK", allocation)
    with pytest.raises(ValueError):
        runner.run(root, ref, "777")


def test_execution_and_independent_accuracy_admission_keep_timing_failed(recovered_joined, monkeypatch):
    root, ref, _, _ = recovered_joined
    report = execute(recovered_joined, monkeypatch)
    assert report["accuracy_admitted"] is report["scientific_timings_admitted"] is False
    with pytest.raises(FileExistsError): runner.run(root, ref, "777")
    assessment_scheduler(monkeypatch)
    result = admission.admit(root, ref, "777", "778", destination(root))
    assert result["status"] == "measurement_failed_native_qfo_assessment_admitted"
    assert result["accuracy_admitted"] is True and result["resources"] is None
    assert set(result["assessment"]["endpoints"]) == set(AXES)
    assert result["assessment"]["secondary_six_metric_mean"] == .5
    assert result["fas_sample"]["sample_membership_verified"] is True
    assert result["native_scheduler"]["verified"]["State"] == "FAILED"
    for k in ("original_native_scheduler_success", "scientific_timings_admitted", "eligible_for_timing_comparison",
        "next_identity_authorized", "publication_ready", "native_inference_reexecuted", "automatic_retry"):
        assert result[k] is False
    with pytest.raises(FileExistsError): admission.admit(root, ref, "777", "778", destination(root))


@pytest.mark.parametrize("problem", ["exit", "tamper", "interrupt"])
def test_execution_failure_is_retained_unadmitted(recovered_joined, monkeypatch, problem):
    root, ref, stage, _ = recovered_joined
    if problem == "interrupt":
        def interrupt(*args, **kwargs): raise KeyboardInterrupt()
        monkeypatch.setattr(runner.subprocess, "run", interrupt)
        with pytest.raises(KeyboardInterrupt): runner.run(root, ref, "777")
    else:
        with pytest.raises(ValueError): execute(recovered_joined, monkeypatch, 1 if problem == "exit" else 0,
            tamper=problem == "tamper")
    output = root / "benchmarks/results/measurement_failed_native_qfo_assessment_v1" / stage["cell"]
    result = json.loads((output / "results.json").read_text())
    assert result["status"] == "failed" and result["accuracy_admitted"] is False
    assert result["resources"] is None and result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("problem", ["scheduler", "failed", "command", "preflight", "outputs", "trace",
    "metrics", "fas_membership", "fas_mean", "inventory", "source", "timing", "resources"])
def test_independent_validation_refuses_changes(recovered_joined, monkeypatch, problem):
    root, ref, _, _ = recovered_joined
    report = execute(recovered_joined, monkeypatch)
    assessment_scheduler(monkeypatch, "RUNNING" if problem == "scheduler" else "COMPLETED")
    output = Path(report["cwd"])
    if problem in {"failed", "command", "outputs", "source", "timing", "resources"}:
        if problem == "failed": report["status"] = "failed"
        elif problem == "command": report["command"].append("-resume")
        elif problem == "outputs": report["outputs"] = []
        elif problem == "source": report["source"] = ref
        elif problem == "timing": report["scientific_timings_admitted"] = True
        else: report["resources"] = dict(cpu_s=0)
        put(output / "results.json", report)
    elif problem == "preflight":
        preflight = json.loads((output / "preflight.json").read_text())
        preflight["native_index"] = 6
        put(output / "preflight.json", preflight)
    elif problem != "scheduler":
        results = Path(report["results"])
        if problem == "inventory": put(results / "extra.json", "{}")
        elif problem == "trace":
            path = results / "stats/trace_fixture.txt"
            path.write_text(path.read_text().replace("COMPLETED", "CACHED", 1))
        elif problem == "metrics":
            path = results / "assessment_out/Assessment_datasets.json"
            rows = json.loads(path.read_text());rows.pop();put(path, rows)
        else:
            path = results / "results/FAS/sample_raw.txt.gz"
            with gzip.open(path, "wt") as stream:
                stream.write("Acc1\tAcc2\tFAS\nA1\tB1\t0.4\nA1\t" +
                    ("D1\t0.6\n" if problem == "fas_membership" else "C1\t0.7\n"))
        if problem != "inventory":
            report["outputs"] = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
            put(output / "results.json", report)
    with pytest.raises(ValueError): admission.admit(root, ref, "777", "778", destination(root))
    result = json.loads((destination(root) / "results.json").read_text())
    assert result["status"] == "measurement_failed_native_qfo_admission_failed_retained"
    assert result["accuracy_admitted"] is result["scientific_timings_admitted"] is False


def test_no_admission_inside_execution(recovered_joined):
    root, ref, _, _ = recovered_joined
    with pytest.raises(ValueError): admission.admit(root, ref, "777", "778", root / "native/assessment")
    assert not (root / "native/assessment").exists()


def test_conflicting_evidence_pins_refused(recovered_joined):
    root, ref, stage, _ = recovered_joined
    changed = deepcopy(stage)
    changed["checked_records"].append(dict(stage["source"], sha256="0" * 64))
    with pytest.raises(ValueError, match="Conflicting"):
        runner.prepare(root, put(Path(ref["path"]), changed), "777")


def test_darwin_path_limit_is_retained(recovered_joined):
    root, ref, stage, manifest = recovered_joined
    with pytest.raises(ValueError, match="path exceeds"):
        runner.execution_spec(root / ("x" * 130), ref, stage, manifest, {})


def test_live_conversion_is_not_scored(recovered_joined, monkeypatch):
    root, ref, _, _ = recovered_joined
    monkeypatch.setattr(runner, "accounting", lambda job: ("live", dict(State="RUNNING")))
    with pytest.raises(ValueError): runner.prepare(root, ref, "777")


def test_empty_failed_assessment_is_not_zero_imputed(recovered_joined, monkeypatch):
    root, ref, stage, _ = recovered_joined
    stage.update(total_pairs=0, retained_pairs=0, expected_pairs=0, empty_predictions=True,
        pair_coverage=dict(pair_rows=0, input_accessions=3, accessions_in_any_pair=0, fraction_inputs_in_any_pair=0.))
    stage["pairs"] = put(Path(stage["pairs"]["path"]), "")
    stage["filtered_pairs"] = put(Path(stage["filtered_pairs"]["path"]), "")
    ref = put(Path(ref["path"]), stage)
    original_bind = runner.bind_recovery
    def bind(*args):
        data = list(original_bind(*args))
        data[5] = dict(data[5], phylogeny=dict(native_pair_rows=0))
        return data
    monkeypatch.setattr(runner, "bind_recovery", bind)
    monkeypatch.setattr(runner.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=1))
    with pytest.raises(ValueError, match="assessment failed"):
        runner.run(root, ref, "777")
    path = root / "benchmarks/results/measurement_failed_native_qfo_assessment_v1" / stage["cell"] / "results.json"
    report = json.loads(path.read_text())
    assert report["status"] == "failed" and report["stage"]["empty_predictions"] is True
    assert report["outputs"] == [] and "assessment" not in report and report["accuracy_admitted"] is False
