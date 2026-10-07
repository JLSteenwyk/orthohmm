"""Reporting metadata fixtures; not native inference or repeated raw admission."""

from copy import deepcopy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import export_allocated_native_qfo_scientific_scores as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_export_native_qfo_factorial_scores import fixture, write
from tests.unit.test_export_native_qfo_scientific_scores import recovered


def allocated(tmp_path, monkeypatch, index=10, changes=None):
    changes = changes or {}
    def change(name, value):
        if name in changes:
            changes[name](value)
        return value
    fasta = tmp_path / "input.fasta"
    fasta.write_text("".join(f">gene{i}\nAAAA\n" for i in range(10)))
    output_root = tmp_path / "native"
    (output_root / "measurement").mkdir(parents=True)
    ready_ref = write(output_root / "measurement/ready.json", dict(metadata_fixture=True))
    suffix = "native/orthohmm_working_res/orthohmm_edges_clustered.txt" if index == 11 else \
        "native/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv"
    prediction_path = output_root / suffix
    prediction_path.parent.mkdir(parents=True)
    prediction_ref = write(prediction_path, dict(metadata_fixture=True))
    plan_ref, report_ref = fixture(tmp_path, index=index, changes={"plan": lambda value:
        value["runs"][index].update(repeat=0, output_root=str(output_root), inputs=[record(fasta)])})
    plan = json.loads(Path(plan_ref["path"]).read_text())
    report = json.loads(Path(report_ref["path"]).read_text())
    root = module.contract.ROOT
    execution_amendment = dict(root=str(root), historical_plan=plan_ref, historical_prefix=[{}]*10,
        new_sources=module.contract.sources())
    amendment_ref = write(tmp_path / "amendment.json", execution_amendment)
    # Contract tests/live amendment cover its transitive integrity; this fixture
    # deliberately tests the direct report handoff without claiming that closure.
    monkeypatch.setattr(module.contract, "validate_amendment", lambda value: plan)
    stage = report["conversion"]
    cpu_ids = list(range(52,84))
    request = dict(schema="allocated_native_factorial_request_v1", execution_authorized=True,
        job_id=10, index=index, plan=plan_ref, amendment=amendment_ref, history=[{}]*index,
        scheduler_command=str(module.contract.SCRIPT), allocation_cwd=str(root), automatic_retry=False)
    request_ref = write(tmp_path / "request.json", change("request", request))
    outputs = dict(schema="allocated_native_factorial_output_review_v1", native_outputs_validated=True,
        request=request_ref, plan=plan_ref, amendment=amendment_ref, index=index, cell=stage["cell"], job_id=10,
        execution_scope=module.contract.SCOPE, native_cpu_ids=cpu_ids, allocated_ready=ready_ref,
        source=record(root / "benchmark_tools/validate_allocated_native_factorial_outputs.py"),
        semantic_validator_source=record(root / "benchmark_tools/validate_native_factorial_outputs.py"),
        gene_ownership_sha256="ownership", checked_files=[prediction_ref], phylogeny=dict(native_pair_rows=3))
    outputs_ref = write(tmp_path / "outputs.json", change("outputs", outputs))
    metadata_ref = write(tmp_path / "metadata.json", {})
    review = dict(schema="allocated_native_factorial_terminal_review_v1", request=request_ref,
        plan=plan_ref, amendment=amendment_ref, job_id=10, index=index, dataset="qfo_corrected",
        cell=stage["cell"], repeat=0, status="native_success", scheduler_state="COMPLETED", scheduler_exit_code="0:0",
        terminal_reviewed=True, native_outputs_validated=True, primary_resources_replayed=True,
        shared_host_resources_reviewed=True, execution_scope=module.contract.SCOPE, resource_scopes=module.original.SCOPES,
        uncontended_timing=False, automatic_retry=False, accuracy_evaluated=False,
        scientific_timings_admitted=False, publication_ready=False, native_cpu_ids=cpu_ids,
        source=record(root / "benchmark_tools/review_allocated_native_factorial_attempt.py"),
        common_reviewer_source=record(root / "benchmark_tools/review_native_factorial_attempt.py"),
        reviews=dict(runtime=metadata_ref, resources=metadata_ref, environment=metadata_ref, outputs_or_failure=outputs_ref),
        resources=dict(wall_seconds=10., cpu_seconds=100., peak_memory_bytes=1024),
        whole_run_maximum_foreign_average_cores=50.)
    review_ref = write(tmp_path / "review.json", change("review", review))
    stage.update(schema="allocated_native_factorial_qfo_conversion_v1",
        status="allocated_native_factorial_qfo_pairs_prepared_unscored", request=request_ref, terminal_review=review_ref,
        amendment=amendment_ref, input_fastas=[record(fasta)], native_cpu_ids=cpu_ids, allocated_ready=ready_ref,
        source=record(root / "benchmark_tools/prepare_allocated_native_factorial_qfo_pairs.py"),
        conversion_kernel_source=record(root / "benchmark_tools/prepare_native_factorial_qfo_pairs.py"),
        pairs=record(fasta), filtered_pairs=record(fasta),
        conversion_started_monotonic_ns=10, conversion_finished_monotonic_ns=20,
        gene_ownership_sha256="ownership", native_input=prediction_ref)
    stage = change("stage", stage)
    env_ref = record(root / "benchmark_tools/results/qfo_assessment_environment_20260917.json")
    stage["environment_manifest"] = env_ref
    manifest = json.loads(Path(env_ref["path"]).read_text())
    pairs_ref = write(tmp_path / "pairs.json", stage)
    raw = "10|COMPLETED|0:0|64|128G|bizon|gpu|1560|2026-10-06T20:00:00|2026-10-06T20:01:00|2026-10-07T10:01:00|orthohmm_allocated_factorial\n"
    terminal = dict(source="fresh_accounting_after_controller_expiry", observation=dict(stdout=raw,
        command=["sacct","-X","-j","10","-n","-P",
            "--format=JobIDRaw,State,ExitCode,AllocCPUS,ReqMem,NodeList,Partition,TimelimitRaw,Submit,Start,End,JobName"]),
        verified=module.contract.terminal_accounting(raw,10),
        controller_observation=dict(command=["scontrol","show","job","10","--oneliner"],
            returncode=1,stdout="",stderr="slurm_load_jobs error: Invalid job id specified"))
    preflight = dict(module.execution_spec(root,pairs_ref,stage,manifest,{}),status="running",job_id="30",
        native_scheduler=deepcopy(terminal))
    preflight_ref = write(tmp_path / "preflight.json", change("preflight", preflight))
    execution = dict(preflight, status="process_succeeded_pending_independent_admission", exit_code=0)
    execution_ref = write(tmp_path / "execution.json", change("execution", execution))
    report.update(schema="allocated_native_factorial_qfo_admission_v1",
        status="allocated_native_factorial_qfo_assessment_admitted", amendment=amendment_ref, conversion=stage,
        pairs_manifest=pairs_ref, preflight=preflight_ref, execution_report=execution_ref,
        source=record(root / "benchmark_tools/admit_allocated_native_factorial_qfo_assessment.py"),
        environment_manifest=env_ref,native_scheduler=deepcopy(terminal))
    return plan_ref, write(tmp_path / "admission.json", change("admission", report))


def collect(refs):
    plan, report = refs
    return module.collect(plan["path"], plan["sha256"], [], [], [(report["path"], report["sha256"])])


@pytest.mark.parametrize("index", [10,11,12])
def test_new_scientific_route_retains_each_endpoint_and_measured_scope(tmp_path, monkeypatch, index):
    result = collect(allocated(tmp_path, monkeypatch, index))
    row = result["rows"][index-6]
    assert result["supplied_allocated_admissions"] == result["supplied_admissions"] == 1
    assert row["status"] == "supplied_allocated_native_admission" and row["accuracy_admitted"] is True
    assert row["scores"]["VGNC"] == .6 and row["endpoint_details"]["VGNC"]["precision"] == .75
    assert row["endpoint_details"]["VGNC"]["recall"] == .5
    assert row["scores"]["GO"] == .2 and row["endpoint_details"]["GO"]["statistic"] != "F1"
    assert row["resources"]["peak_memory_bytes"] == 1024 and row["scientific_timings_admitted"] is False
    assert row["whole_run_maximum_foreign_average_cores"] == 50
    assert row["relation_coverage"] == .4 and len(row["native_cpu_ids"]) == 32
    assert row["prediction_semantics"] == ("cross-species group-derived clique pairs" if index == 11 else
        "native phylogenetically inferred pairs")
    assert result["rows"][0]["scores"]["VGNC"] is None and result["publication_ready"] is False


@pytest.mark.parametrize("target,key,value", [
    ("admission", "schema", "full_native_factorial_qfo_admission_v1"), ("admission", "status", "running"),
    ("admission", "accuracy_admitted", False), ("admission", "publication_ready", True),
    ("admission", "automatic_retry", True), ("admission", "native_index", 9),
    ("admission", "native_index", True), ("admission", "native_job_id", True),
    ("admission", "source", {}), ("admission", "amendment", {}), ("admission", "conversion", {}),
    ("admission", "native_scheduler", dict(verified=dict(State="RUNNING", ExitCode="0:0"))),
    ("request", "schema", "native_factorial_cost_request_v1"), ("request", "index", 11),
    ("request", "amendment", {}), ("request", "automatic_retry", True),
    ("review", "schema", "native_factorial_terminal_review_v1"), ("review", "status", "native_failure_retained"),
    ("review", "source", {}), ("review", "common_reviewer_source", {}), ("review", "repeat", True),
    ("review", "native_outputs_validated", False), ("review", "shared_host_resources_reviewed", False),
    ("review", "scientific_timings_admitted", True), ("review", "resources", None),
    ("review", "whole_run_maximum_foreign_average_cores", float("nan")),
    ("review", "native_cpu_ids", [0]), ("review", "native_cpu_ids", list(range(31))+[True]),
    ("outputs", "schema", "native_factorial_output_review_v1"), ("outputs", "source", {}),
    ("outputs", "semantic_validator_source", {}), ("outputs", "native_cpu_ids", [0]),
    ("outputs", "gene_ownership_sha256", "different"), ("outputs", "allocated_ready", {}),
    ("stage", "schema", "full_native_factorial_qfo_conversion_v1"), ("stage", "conversion_kind", "group"),
    ("stage", "semantics", "preclustering edges"), ("stage", "removed_mapping_pairs", 1),
    ("stage", "total_pairs", True), ("stage", "input_fastas", []), ("stage", "native_cpu_ids", [0]),
    ("stage", "conversion_kernel_source", {}), ("stage", "source", {}),
    ("execution", "schema", "full_native_factorial_qfo_execution_v1"), ("execution", "exit_code", True),
    ("execution", "source", {}), ("execution", "accuracy_admitted", True), ("execution", "amendment", {}),
    ("preflight", "status", "prepared_unrun")])
def test_startup_historical_or_mismatched_report_not_scored(tmp_path, monkeypatch, target, key, value):
    refs = allocated(tmp_path, monkeypatch, changes={target: lambda item: item.update({key:value})})
    with pytest.raises((ValueError, KeyError)):
        collect(refs)


@pytest.mark.parametrize("change", [
    lambda value: value["assessment"]["endpoints"].pop("GO"),
    lambda value: value["assessment"]["endpoints"]["GO"].update(score=.3),
    lambda value: value["assessment"]["endpoints"]["GO"].update(score_semantics="F1"),
    lambda value: value["assessment"].update(secondary_six_metric_mean=0.),
    lambda value: value["assessment"]["endpoints"]["VGNC"]["native_participant"].update(metric_x=True),
    lambda value: value["assessment"]["endpoints"]["VGNC"]["native_participant"].update(metric_y=float("inf")),
    lambda value: value["assessment"]["endpoints"]["FAS"]["native_participant"].update(participant_id="cached"),
    lambda value: value["fas_sample"].update(sample_membership_verified=False),
    lambda value: value["conversion_scheduler"].update(AllocCPUS="8"),
    lambda value: value["scheduler"].update(State="RUNNING")])
def test_actual_statistic_source_and_native_validation_required(tmp_path, monkeypatch, change):
    with pytest.raises(ValueError):
        collect(allocated(tmp_path, monkeypatch, changes={"admission":change}))


@pytest.mark.parametrize("key,value", [("command", ["nextflow","-resume"]), ("cwd", "/old/assessment"),
    ("work", "/old/w"), ("results", "/old/scoring"), ("environment_overrides", {"NXF_OFFLINE":"false"})])
def test_coherently_resealed_preflight_does_not_override_frozen_command_or_namespace(tmp_path, monkeypatch, key, value):
    with pytest.raises(ValueError,match="command, namespace or environment"):
        collect(allocated(tmp_path,monkeypatch,changes={"preflight":lambda item:item.update({key:value})}))


@pytest.mark.parametrize("change", [
    lambda value:value["native_scheduler"].update(source="cached"),
    lambda value:value["native_scheduler"]["observation"].update(command=["sacct","-j","11"]),
    lambda value:value["native_scheduler"]["observation"].update(stdout="incomplete"),
    lambda value:value["native_scheduler"]["verified"].update(AllocCPUS="32"),
    lambda value:value["native_scheduler"]["verified"].update(JobIDRaw="11"),
    lambda value:value["native_scheduler"]["controller_observation"].update(stderr="Connection refused"),
    lambda value:value["native_scheduler"]["controller_observation"].update(returncode=0),
    lambda value:value["native_scheduler"]["controller_observation"].update(returncode=True)])
def test_retained_accounting_not_relabelled_as_new_native_success(tmp_path,monkeypatch,change):
    with pytest.raises(ValueError):
        collect(allocated(tmp_path,monkeypatch,changes={"admission":change}))


def live_terminal():
    from benchmark_tools.capture_array_scheduler import REQUIRED
    fields=dict.fromkeys(REQUIRED - {"ArrayJobId","ArrayTaskId"},"unused")
    fields.update(JobId="10",JobName="orthohmm_allocated_factorial",JobState="COMPLETED",Partition="gpu",
        NodeList="bizon",NumNodes="1",NumCPUs="64",NumTasks="1",OverSubscribe="OK",MinMemoryNode="128G",
        Requeue="0",Restarts="0",Command=str(module.contract.SCRIPT),WorkDir=str(module.contract.ROOT),
        TimeLimit="1-02:00:00",ExitCode="0:0",Comment="a"*64)
    fields["CPUs/Task"]="64"
    raw=" ".join(key+"="+value for key,value in fields.items())
    parsed=module.contract.validate_controller(raw,10,"terminal",command=str(module.contract.SCRIPT),
        cwd=str(module.contract.ROOT),time_limit="1-02:00:00",allocation_mode="shared")
    return dict(source="live_controller",observation=dict(command=["scontrol","show","job","10","--oneliner"],
        returncode=0,stdout=raw,stderr=""),verified=parsed)


def test_retained_successful_live_controller_is_reparsed():
    module.native_terminal(live_terminal(),10,dict(sha256="a"*64),module.contract.ROOT)


@pytest.mark.parametrize("change", ["comment","name","state","parser","command","returncode"])
def test_retained_live_controller_wrong_binding_refused(change):
    value=live_terminal()
    if change == "parser": value["verified"]["fields"]["NumCPUs"]="32"
    elif change == "command": value["observation"]["command"][3]="11"
    elif change == "returncode": value["observation"]["returncode"]=False
    else:
        old,new={"comment":("Comment="+"a"*64,"Comment="+"b"*64),
            "name":("JobName=orthohmm_allocated_factorial","JobName=old"),
            "state":("JobState=COMPLETED","JobState=RUNNING")}[change]
        value["observation"]["stdout"]=value["observation"]["stdout"].replace(old,new)
    with pytest.raises(ValueError):
        module.native_terminal(value,10,dict(sha256="a"*64),module.contract.ROOT)


def test_historical_successful_and_recovered_rows_are_not_rewritten(tmp_path):
    normal_path, recovery_path = tmp_path/"normal", tmp_path/"recovered"
    normal_path.mkdir(); recovery_path.mkdir()
    for refs, legacy, recovered_refs in ((fixture(normal_path), True, False),
            (recovered(recovery_path), False, True)):
        plan, ref = refs
        args = [(ref["path"],ref["sha256"])]
        normal, recovery = (args,[]) if legacy else ([],args)
        old = module.historical.collect(plan["path"],plan["sha256"],normal,recovery)
        new = module.collect(plan["path"],plan["sha256"],normal,recovery,[])
        assert old["rows"] == new["rows"]
        assert new["supplied_allocated_admissions"] == 0
        if recovered_refs:
            row = new["rows"][1]
            assert row["resources"] is None and row["timing_eligible"] is row["timing_admitted"] is False


def test_duplicate_cross_route_identity_refused(tmp_path, monkeypatch):
    old_path, new_path = tmp_path/"old", tmp_path/"new"
    old_path.mkdir(); new_path.mkdir()
    plan, old_ref = fixture(old_path,index=10)
    _, new_ref = allocated(new_path,monkeypatch)
    with pytest.raises(ValueError,match="duplicate"):
        module.collect(plan["path"],plan["sha256"],[(old_ref["path"],old_ref["sha256"])],[],
            [(new_ref["path"],new_ref["sha256"])])


def test_export_keeps_missing_values_disclosure_and_no_overwrite(tmp_path, monkeypatch):
    plan, ref = allocated(tmp_path,monkeypatch)
    output = tmp_path/"table"
    args = (plan["path"],plan["sha256"],[],[],[(ref["path"],ref["sha256"])],output)
    result = module.export(*args)
    with (output/"scores.tsv").open() as stream:
        rows = list(csv.DictReader(stream,delimiter="\t"))
    assert len(rows) == 7 and rows[0]["VGNC F1"] == "" and rows[4]["VGNC F1"] == "0.6"
    text = (output/"scores.md").read_text()
    assert "Unavailable" in text and "GO similarity" in text and "Secondary mean" in text
    assert result["timing_disclosure"] in text
    assert result["new_scoring_or_admission"] is False
    with pytest.raises(ValueError,match="already exists"):
        module.export(*args)


def test_modified_supplied_receipt_fails_before_writing_export(tmp_path, monkeypatch):
    plan, ref = allocated(tmp_path,monkeypatch)
    Path(ref["path"]).write_text("{}\n")
    with pytest.raises(ValueError,match="checksum"):
        module.export(plan["path"],plan["sha256"],[],[],[(ref["path"],ref["sha256"])],tmp_path/"table")
    assert not (tmp_path/"table").exists()
