"""Synthetic full-native handoffs; no scientific QfO job or timing benchmark."""

from copy import deepcopy
import gzip
import json
from pathlib import Path
import statistics
import tempfile
from types import SimpleNamespace

import pytest

from benchmark_tools import admit_native_factorial_qfo_assessment as admission
from benchmark_tools import run_native_factorial_qfo_assessment as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_native_factorial_cost import IDENTITIES, SCOPE
from benchmark_tools.derive_threadripper_resources import SCOPES
from benchmark_tools.validate_qfo_native_assessment import AXES


def stage_fixture(index=7, empty=False):
    count = 0 if empty else 2
    kind = "native" if IDENTITIES[index][1].endswith("r1") else "group"
    stage = dict(schema="full_native_factorial_qfo_conversion_v1",
        status="full_native_factorial_qfo_pairs_prepared_unscored", native_index=index,
        cell=IDENTITIES[index][1], native_job_id=999, job_id="777",
        participant="ohmm_qfo_full_native_" + IDENTITIES[index][1], conversion_kind=kind,
        semantics="native phylogenetically inferred pairs" if kind == "native" else "cross-species group-derived clique pairs",
        accuracy_evaluated=False, native_inference_reexecuted=False, automatic_retry=False,
        next_identity_authorized=False, publication_ready=False,
        total_pairs=count, retained_pairs=count, expected_pairs=count, removed_mapping_pairs=0,
        empty_predictions=empty, pairs=dict(bytes=0 if empty else 12, sha256="a" * 64),
        filtered_pairs=dict(bytes=0 if empty else 12, sha256="a" * 64),
        pair_coverage=dict(pair_rows=count, input_accessions=3, accessions_in_any_pair=0 if empty else 3,
            fraction_inputs_in_any_pair=0 if empty else 1.))
    scheduler = dict(JobIDRaw="777", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2")
    return stage, scheduler


@pytest.mark.parametrize("index", range(6, 13))
@pytest.mark.parametrize("empty", [False, True])
def test_seven_identities_and_empty_predictions(index, empty):
    stage, scheduler = stage_fixture(index, empty)
    runner.validate_stage(stage, scheduler)


@pytest.mark.parametrize("key,value", [("schema", "cached"), ("status", "running"), ("native_index", True),
    ("native_index", 0), ("cell", "p1_c0_r0"), ("participant", "ohmm_qfo_corrected_factorial_p0_c0_r1"),
    ("conversion_kind", "group"), ("semantics", "root HOG cliques"), ("accuracy_evaluated", True),
    ("native_inference_reexecuted", True), ("automatic_retry", True), ("next_identity_authorized", True),
    ("publication_ready", True), ("job_id", "778"), ("total_pairs", True), ("retained_pairs", 1),
    ("expected_pairs", 3), ("removed_mapping_pairs", 1), ("empty_predictions", True),
    ("filtered_pairs", dict(bytes=9, sha256="a" * 64)),
    ("pair_coverage", dict(pair_rows=2, input_accessions=3, accessions_in_any_pair=2, fraction_inputs_in_any_pair=1.))])
def test_invalid_conversion_rejected(key, value):
    stage, scheduler = stage_fixture()
    stage[key] = value
    with pytest.raises(ValueError):
        runner.validate_stage(stage, scheduler)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"), ("NodeList", "dgx"),
    ("AllocCPUS", "8"), ("JobIDRaw", "778")])
def test_conversion_scheduler_rejected(key, value):
    stage, scheduler = stage_fixture()
    scheduler[key] = value
    with pytest.raises(ValueError):
        runner.validate_stage(stage, scheduler)


def put(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(content) if not isinstance(content, str) else content)
    return record(path)


@pytest.fixture
def short_root():
    # Darwin's validated path limit applies to synthetic handoffs as well.
    with tempfile.TemporaryDirectory(prefix="nq-", dir="/tmp") as directory:
        yield Path(directory).resolve()


@pytest.fixture
def joined(short_root, monkeypatch):
    root = short_root
    stage, scheduler = stage_fixture()
    stage["source"] = record(Path(runner.__file__).with_name("prepare_native_factorial_qfo_pairs.py"))
    stage["pairs"] = put(root / "conversion/pairs.tsv", "A1\tB1\nA1\tC1\n")
    stage["filtered_pairs"] = put(root / "conversion/pairs.qfo.tsv", "A1\tB1\nA1\tC1\n")
    stage["native_input"] = put(root / "native/native.tsv", "synthetic native prediction\n")
    stage["gene_ownership_sha256"] = "c" * 64
    stage["input_fastas"] = [put(root / "input/one.fasta", ">A1\nAAAA\n")]
    plan_ref = put(root / "plan.json", dict(helper_sources=[stage["source"]]))
    request_ref = put(root / "request.json", dict(job_id=999, index=7, plan=plan_ref))
    output_ref = put(root / "outputs.json", dict(native_outputs_validated=True, request=request_ref,
        gene_ownership_sha256=stage["gene_ownership_sha256"], checked_files=[stage["native_input"]], evidence=[],
        phylogeny=dict(native_pair_rows=2)))
    sched_ref = put(root / "scheduler.json", "{}")
    review_ref = put(root / "review.json", dict(schema="native_factorial_terminal_review_v1", request=request_ref,
        plan=plan_ref, job_id=999, index=7, dataset="qfo_corrected", cell=stage["cell"], repeat=0,
        status="native_success", scheduler_state="COMPLETED", scheduler_exit_code="0:0", terminal_reviewed=True,
        native_outputs_validated=True, primary_resources_replayed=True, shared_host_resources_reviewed=True,
        execution_scope=SCOPE, resource_scopes=SCOPES, uncontended_timing=False, automatic_retry=False,
        source=record(Path(runner.__file__).with_name("review_native_factorial_attempt.py")), scheduler=sched_ref,
        reviews=dict(runtime=sched_ref, resources=sched_ref, environment=sched_ref, outputs_or_failure=output_ref)))
    stage.update(request=request_ref, terminal_review=review_ref, plan=plan_ref, checked_records=[stage["source"]])
    pipeline = root / "pipeline"
    main = put(pipeline / "main.nf", "// Synthetic pinned pipeline\n")
    fas = put(pipeline / "fas_benchmark.py", "MAX_PAIRS_COMPUTE = 9000\n")
    mapping = put(pipeline / "reference_data/2020/mapping.json.gz", "synthetic mapping\n")
    template_refs = []
    for challenge, axes in AXES.items():
        template_refs.append(put(pipeline / "reference_data/data" / (challenge + ".json"),
            dict(type="aggregation", challenge_ids=[challenge], datalink=dict(inline_data=dict(
                visualization=dict(x_axis=axes[0], y_axis=axes[1], type="2D-plot"))))))
    swiss = put(pipeline / "reference_data/2020/ReconciledTrees_SwissTrees.drw",
        "ReconciledTrees['X'] := RecTreeCase('X',...)\n")
    manifest = dict(status="local_qfo_assessment_environment_frozen", accuracy_evaluated=False, pipeline=str(pipeline),
        source=fas, execution_config=main, singularity_config=main, pipeline_files=[main, fas],
        reference_files=[mapping, swiss, *template_refs], java_files=[], images=[], executables=[], singularity_support=[],
        environment_overrides=dict(NXF_OFFLINE="true"))
    env_ref = put(root / "benchmark_tools/results/qfo_assessment_environment_20260917.json", manifest)
    stage.update(environment_manifest=env_ref, mapping=mapping)
    pairs_ref = put(root / "conversion/results.json", stage)
    monkeypatch.setattr(runner, "ENV_SHA", env_ref["sha256"])
    monkeypatch.setattr(runner, "validate_request", lambda *a: None)
    run = dict(index=7, dataset="qfo_corrected", cell=stage["cell"], repeat=0, genes=3, inputs=stage["input_fastas"])
    monkeypatch.setattr(runner, "validate_plan", lambda *a: [run] * 13)
    monkeypatch.setattr(runner, "verify_terminal", lambda job: dict(source="live_controller",
        verified=dict(fields=dict(JobState="COMPLETED", ExitCode="0:0", Comment=request_ref["sha256"]))))
    monkeypatch.setattr(runner, "accounting", lambda job: ("synthetic conversion accounting", deepcopy(scheduler)))
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setenv("SLURM_JOB_ID", "778")
    return root, pairs_ref, stage, manifest


def write_native_endpoints(results, stage, manifest):
    participant = stage["participant"]
    rows = []
    sem = statistics.stdev([.4, .6]) / 2 ** .5
    for challenge, axes in {**AXES, "SwissTrees-X": ("TPR", "PPV")}.items():
        for metric in axes:
            rows.append(dict(_id=challenge + metric, type="assessment", community_id="QfO", participant_id=participant,
                challenge_id=challenge, metrics=dict(metric_id=metric, value=2 if metric == "NR_ORTHOLOGS" else .5,
                    stderr=sem if challenge == "FAS" and metric == "FAS" else .01)))
    put(results / "assessment_out/Assessment_datasets.json", rows)
    for challenge, axes in AXES.items():
        template = json.loads((Path(manifest["pipeline"]) / "reference_data/data" / (challenge + ".json")).read_text())
        template["datalink"]["inline_data"]["challenge_participants"] = [dict(participant_id=participant,
            metric_x=2 if axes[0] == "NR_ORTHOLOGS" else .5, metric_y=.5, stderr_x=.01,
            stderr_y=sem if challenge == "FAS" else .01)]
        put(results / "results" / challenge / (challenge + ".json"), template)
    fas_path = results / "results/FAS/sample_raw.txt.gz"
    with gzip.open(fas_path, "wt") as stream:
        stream.write("Acc1\tAcc2\tFAS\nA1\tB1\t0.4\nA1\tC1\t0.6\n")
    names = ["validate_input_file", "convertPredictions", "consolidate", "vgnc_benchmark (1)", "ec_benchmark (1)",
        "go_benchmark (1)", "fas_benchmark (1)", "reference_genetrees_benchmark (SwissTrees)",
        "reference_genetrees_benchmark (TreeFam-A)", *[f"scheduleMetrics ({i})" for i in range(1, 7)]]
    put(results / "stats/trace_fixture.txt", "task_id\tname\tstatus\texit\n" + "".join(
        f"{i}\t{name}\tCOMPLETED\t{'-' if name.startswith('scheduleMetrics') else '0'}\n" for i, name in enumerate(names)))


def execute_fixture(joined, monkeypatch, exit_code=0, tamper=False):
    root, pairs_ref, stage, manifest = joined
    def execute(command, **kwargs):
        assert "-resume" not in command and command[command.index("--event_year") + 1] == "2020"
        assert command[command.index("--challenges_ids") + 1] == "GO EC VGNC SwissTrees TreeFam-A FAS"
        assert kwargs["env"]["NXF_OFFLINE"] == "true"
        write_native_endpoints(Path(command[command.index("--results_dir") + 1]), stage, manifest)
        if tamper:
            Path(stage["filtered_pairs"]["path"]).write_text("changed")
        return SimpleNamespace(returncode=exit_code)
    monkeypatch.setattr(runner.subprocess, "run", execute)
    return runner.run(root, pairs_ref, "777")


def test_joined_prepare_is_read_only_and_distinct(joined):
    root, ref, stage, _ = joined
    report = runner.run(root, ref, "777", check_only=True)
    assert report["status"] == "prepared_unrun" and report["accuracy_admitted"] is False
    assert report["work"].endswith("/nq07") and report["results"].endswith("/full_native_07")
    assert report["fas_protocol"]["newly_computed_pair_cap"] == 9000
    assert "not all submitted" in report["fas_protocol"]["population"]
    assert not Path(report["cwd"]).exists()


@pytest.mark.parametrize("problem", ["source", "input", "mapping", "request", "owner", "output", "native_failed", "comment"])
def test_joined_binding_rejects_changes(joined, monkeypatch, problem):
    root, ref, stage, _ = joined
    if problem in {"source", "input", "mapping", "request", "owner", "output"}:
        modified = deepcopy(stage)
        if problem == "source":
            modified["source"] = modified["native_input"]
        elif problem == "input":
            modified["input_fastas"] = []
        elif problem == "mapping":
            modified["mapping"] = modified["native_input"]
        elif problem == "request":
            modified["native_job_id"] = 1
        elif problem == "owner":
            modified["gene_ownership_sha256"] = "bad"
        else:
            modified["native_input"] = modified["pairs"]
        ref = put(Path(ref["path"]), modified)
    else:
        monkeypatch.setattr(runner, "verify_terminal", lambda job: dict(source="live_controller",
            verified=dict(fields=dict(JobState="FAILED" if problem == "native_failed" else "COMPLETED",
                ExitCode="0:0", Comment="wrong"))))
    with pytest.raises(ValueError):
        runner.prepare(root, ref, "777")


def test_live_conversion_rejected_before_reading(joined, monkeypatch):
    root, ref, _, _ = joined
    monkeypatch.setattr(runner, "accounting", lambda job: ("live", dict(State="RUNNING")))
    with pytest.raises(ValueError):
        runner.prepare(root, ref, "777")


def test_empty_submission_failure_retained_without_zero_score(joined, monkeypatch):
    root, ref, stage, _ = joined
    stage.update(total_pairs=0, retained_pairs=0, expected_pairs=0, empty_predictions=True,
        pair_coverage=dict(pair_rows=0, input_accessions=3, accessions_in_any_pair=0, fraction_inputs_in_any_pair=0.))
    stage["pairs"] = put(Path(stage["pairs"]["path"]), "")
    stage["filtered_pairs"] = put(Path(stage["filtered_pairs"]["path"]), "")
    review = json.loads(Path(stage["terminal_review"]["path"]).read_text())
    output_ref = review["reviews"]["outputs_or_failure"]
    output_review = json.loads(Path(output_ref["path"]).read_text())
    output_review["phylogeny"]["native_pair_rows"] = 0
    review["reviews"]["outputs_or_failure"] = put(Path(output_ref["path"]), output_review)
    stage["terminal_review"] = put(Path(stage["terminal_review"]["path"]), review)
    ref = put(Path(ref["path"]), stage)
    monkeypatch.setattr(runner.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=1))
    with pytest.raises(ValueError, match="assessment failed"):
        runner.run(root, ref, "777")
    output = root / "benchmarks/results/full_native_qfo_assessment_v1" / stage["cell"]
    report = json.loads((output / "results.json").read_text())
    assert report["stage"]["empty_predictions"] is True and report["status"] == "failed"
    assert "assessment" not in report and report["outputs"] == []


def test_invalid_fas_cap_not_accepted(joined):
    root, _, _, manifest = joined
    source = Path(manifest["pipeline_files"][1]["path"])
    manifest["pipeline_files"][1] = put(source, "MAX_PAIRS_COMPUTE = 9001\n")
    with pytest.raises(ValueError, match="sampling cap"):
        runner.fas_protocol(manifest)


def test_long_darwin_path_still_rejected(joined):
    root, ref, stage, manifest = joined
    with pytest.raises(ValueError, match="path exceeds"):
        runner.execution_spec(root / ("x" * 130), ref, stage, manifest, {})


@pytest.mark.parametrize("allocation", [None, "2", "32"])
def test_runner_requires_eight_cpu_own_allocation(joined, monkeypatch, allocation):
    root, ref, _, _ = joined
    if allocation is None:
        monkeypatch.delenv("SLURM_JOB_ID")
    else:
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", allocation)
    with pytest.raises(ValueError):
        runner.run(root, ref, "777")


@pytest.mark.parametrize("problem", ["exit", "tamper", "interrupt"])
def test_failed_execution_retained_without_admission(joined, monkeypatch, problem):
    root, ref, stage, _ = joined
    if problem == "interrupt":
        def interrupt(*a, **k):
            raise KeyboardInterrupt()
        monkeypatch.setattr(runner.subprocess, "run", interrupt)
        with pytest.raises(KeyboardInterrupt):
            runner.run(root, ref, "777")
    else:
        with pytest.raises(ValueError):
            execute_fixture(joined, monkeypatch, exit_code=1 if problem == "exit" else 0, tamper=problem == "tamper")
    output = root / "benchmarks/results/full_native_qfo_assessment_v1" / stage["cell"]
    report = json.loads((output / "results.json").read_text())
    assert report["status"] == "failed" and report["accuracy_admitted"] is False
    assert (output / "preflight.json").exists() and (output / "scoring.log").exists()
    assert report["finished_monotonic_ns"] >= report["started_monotonic_ns"]


def assessment_scheduler(monkeypatch, state="COMPLETED"):
    monkeypatch.setattr(admission, "accounting", lambda job: ("synthetic assessment accounting",
        dict(JobIDRaw="778", State=state, ExitCode="0:0", NodeList="bizon", AllocCPUS="8")))


def destination(root):
    return root / "benchmarks/results/full_native_qfo_admission_v1/test_admission"


def test_joined_execution_and_native_admission(joined, monkeypatch):
    root, ref, _, _ = joined
    execution = execute_fixture(joined, monkeypatch)
    assert execution["status"] == "process_succeeded_pending_independent_admission"
    assert execution["accuracy_admitted"] is False
    with pytest.raises(FileExistsError):
        runner.run(root, ref, "777")
    assessment_scheduler(monkeypatch)
    result = admission.admit(root, ref, "777", "778", destination(root))
    assert result["accuracy_admitted"] is True and result["publication_ready"] is False
    assert set(result["assessment"]["endpoints"]) == set(AXES)
    assert result["assessment"]["secondary_six_metric_mean"] == .5
    assert result["fas_sample"]["sample_membership_verified"] is True
    assert result["fas_sample"]["sample_pairs"] == 2 and result["fas_sample"]["sample_fraction"] == 1
    assert result["fas_sample"]["proteins_in_multiple_sample_pairs"] == 1
    assert result["next_identity_authorized"] is False
    with pytest.raises(FileExistsError):
        admission.admit(root, ref, "777", "778", destination(root))


@pytest.mark.parametrize("problem", ["scheduler", "failed", "command", "preflight", "outputs", "trace",
    "metrics", "fas_membership", "fas_mean", "inventory", "source"])
def test_independent_admission_rejects_incomplete_or_changed_execution(joined, monkeypatch, problem):
    root, ref, _, _ = joined
    execution = execute_fixture(joined, monkeypatch)
    assessment_scheduler(monkeypatch, state="RUNNING" if problem == "scheduler" else "COMPLETED")
    output = Path(execution["cwd"])
    if problem in {"failed", "command", "outputs", "source"}:
        if problem == "failed":
            execution["status"] = "failed"
        elif problem == "command":
            execution["command"].append("-resume")
        elif problem == "outputs":
            execution["outputs"] = []
        else:
            execution["source"] = ref
        put(output / "results.json", execution)
    elif problem == "preflight":
        preflight = json.loads((output / "preflight.json").read_text())
        preflight["native_index"] = 6
        put(output / "preflight.json", preflight)
    elif problem != "scheduler":
        results = Path(execution["results"])
        if problem == "inventory":
            put(results / "extra.json", "{}")
        elif problem == "trace":
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
        if problem != "inventory":
            # Test content validators independently of the unchanged-inventory gate.
            execution["outputs"] = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
            put(output / "results.json", execution)
    with pytest.raises(ValueError):
        admission.admit(root, ref, "777", "778", destination(root))
    report = json.loads((destination(root) / "results.json").read_text())
    assert report["status"] == "full_native_factorial_qfo_admission_failed_retained"
    assert report["accuracy_admitted"] is False


def test_admission_cannot_write_inside_native_execution(joined):
    root, ref, stage, _ = joined
    with pytest.raises(ValueError, match="separate direct"):
        admission.admit(root, ref, "777", "778", root / "native/assessment")
    assert not (root / "native/assessment").exists()
