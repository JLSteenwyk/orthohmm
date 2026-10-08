"""New review lineage gates, not production review success or scores."""

import json
import os
from pathlib import Path

import pytest

from benchmark_tools import native12_composed_review_binding as current
from benchmark_tools.prepare_allocated_native_factorial_qfo_pairs import admit_conversion
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture():
    request_ref = {"request": 12}
    request = dict(schema=current.executor.REQUEST_SCHEMA, job_id=24036,
        plan={"plan": 1}, amendment={"amendment": 1})
    run = dict(index=12, cell="p1_c1_r1", dataset="qfo_corrected", repeat=0)
    review = dict(schema=current.reviewer.SCHEMA, status="native_success", request=request_ref,
        plan=request["plan"], amendment=request["amendment"], job_id=24036,
        index=12, cell=run["cell"], dataset=run["dataset"], repeat=0,
        scheduler_state="COMPLETED", scheduler_exit_code="0:0", execution_scope=current.SCOPE,
        resource_scopes=current.SCOPES, terminal_reviewed=True, native_outputs_validated=True,
        primary_resources_replayed=True, shared_host_resources_reviewed=True,
        prospective_current_inventory_equality=True, current_original_os_inventory_equality=False,
        original_review_translated=False, continuous_runtime_integrity_established=False,
        next_identity_authorized=False, accuracy_evaluated=False, scientific_timings_admitted=False,
        uncontended_timing=False, automatic_retry=False, publication_ready=False,
        historical_review_failures_retained=[23986, 24033],
        reviews={key: {} for key in ("runtime", "resources", "environment", "outputs_or_failure")})
    producer = dict(JobIDRaw="25000", State="COMPLETED", ExitCode="0:0", NodeList="bizon",
        AllocCPUS="2", ReqMem="128G")
    return review, request_ref, request, run, producer, 25000


def test_new_successful_review_accepts_native_pairs_not_old_schema():
    values = fixture()
    assert current.admit_review(*values) == "native"
    with pytest.raises(ValueError):
        admit_conversion(*values[:4])


@pytest.mark.parametrize("field", ["schema", "status", "request", "plan", "amendment", "job_id",
    "index", "cell", "dataset", "repeat", "scheduler_state", "scheduler_exit_code", "execution_scope",
    "resource_scopes", "terminal_reviewed", "native_outputs_validated", "primary_resources_replayed",
    "shared_host_resources_reviewed", "prospective_current_inventory_equality",
    "current_original_os_inventory_equality", "original_review_translated", "continuous_runtime_integrity_established",
    "next_identity_authorized", "accuracy_evaluated", "scientific_timings_admitted", "uncontended_timing",
    "automatic_retry", "publication_ready", "historical_review_failures_retained", "reviews"])
def test_changed_review_identity_or_proof_is_rejected(field):
    values = fixture()
    review = values[0]
    review[field] = not review[field] if type(review[field]) is bool else "wrong"
    with pytest.raises(ValueError):
        current.admit_review(*values)


@pytest.mark.parametrize("field", ["JobIDRaw", "State", "ExitCode", "NodeList", "AllocCPUS", "ReqMem"])
def test_producer_accounting_must_be_actual_successful_envelope(field):
    values = fixture()
    values[4][field] = "wrong"
    with pytest.raises(ValueError, match="producer"):
        current.admit_review(*values)


@pytest.mark.parametrize("object_index,field", [(0, "repeat"), (0, "index"), (3, "index"),
    (2, "schema"), (2, "job_id"), (3, "cell")])
def test_boolean_indices_or_wrong_request_or_group_cell_refuse(object_index, field):
    values = fixture()
    values[object_index][field] = True if field in {"repeat", "index"} else "wrong"
    with pytest.raises(ValueError):
        current.admit_review(*values)


def raw_held(job, **changes):
    fields = dict(JobId=str(job), JobName="ohmm_native12_review", JobState="PENDING", Reason="JobHeldUser",
        Partition="gpu", ReqNodeList="bizon", NumCPUs="2", NumTasks="1", NumNodes="1-1",
        MinMemoryNode="128G", TimeLimit="06:00:00", Requeue="0", Restarts="0",
        Command=str(current.REVIEW_BATCH), WorkDir=str(current.ROOT), Comment=current.REVIEWER_SHA,
        UserId=f"fixture({os.getuid()})", Dependency="(null)")
    fields["CPUs/Task"] = "2"
    fields.update(changes)
    return " ".join(key + "=" + value for key, value in fields.items()) + "\n"


def test_exact_owned_held_producer_envelope():
    assert current.producer_envelope(raw_held(25000), 25000)["NumCPUs"] == "2"


@pytest.mark.parametrize("field", ["JobId", "JobName", "JobState", "Reason", "Partition", "ReqNodeList",
    "NumCPUs", "NumTasks", "NumNodes", "MinMemoryNode", "TimeLimit", "Requeue", "Restarts", "Command",
    "WorkDir", "Comment", "UserId", "Dependency", "CPUs/Task", "ArrayJobId", "HetJobId"])
def test_wrong_or_repeated_producer_envelope_cannot_admit(field):
    with pytest.raises(ValueError, match="envelope"):
        current.producer_envelope(raw_held(25000, **{field: "wrong"}), 25000)


def test_duplicate_scheduler_fields_refuse():
    with pytest.raises(ValueError):
        current.producer_envelope(raw_held(25000).strip() + " NumCPUs=2\n", 25000)


@pytest.fixture
def joined(tmp_path, monkeypatch):
    def save(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value, sort_keys=True))
        return record(path)

    review, _, request, run, producer, producer_job = fixture()
    root = tmp_path / "repo"
    tools = root / "benchmark_tools"
    tools.mkdir(parents=True)
    destination = root / "review"
    source_ref = record(current.reviewer.__file__)
    batch = tools / "review.sh"
    batch.write_text("fixture batch\n")
    batch_ref = record(batch)
    monkeypatch.setattr(current, "ROOT", root)
    monkeypatch.setattr(current, "REVIEW_BATCH", batch)
    monkeypatch.setattr(current, "REVIEW_BATCH_SHA", batch_ref["sha256"])
    monkeypatch.setattr(current.reviewer, "DESTINATION", destination)
    request_path = root / "request.json"
    monkeypatch.setattr(current.executor, "REQUEST", request_path)
    run["output_root"] = str(root / "native")
    ready = save(root / "native/measurement/ready.json", {"fixture": "placement"})
    request["new_sources"] = [source_ref, batch_ref]
    request["runtime_basis"] = save(root / "basis.json", {"fixture": "runtime"})
    request["held_scheduler"] = {"stdout": "stubbed native held receipt"}
    request_ref = save(request_path, request)
    monkeypatch.setattr(current, "REQUEST_SHA", request_ref["sha256"])
    semantic = tools / "validate_native_factorial_outputs.py"
    semantic.write_text("fixture original validator\n")
    replay_source = tools / "replay_allocated_threadripper_scaling.py"
    replay_source.write_text("fixture original accounting\n")
    output = dict(schema="native12_composed_output_review_v1", source=source_ref,
        semantic_validator_source=record(semantic), native_outputs_validated=True, accuracy_evaluated=False,
        index=12, job_id=24036, cell="p1_c1_r1", request=request_ref,
        plan=request["plan"], amendment=request["amendment"], native_cpu_ids=[1], allocated_ready=ready,
        phylogeny=dict(native_pair_rows=3), checked_files=[], evidence=[])
    runtime = dict(schema="native12_composed_runtime_review_v1", status="fresh_runtime_brackets_and_lookup_replayed",
        runtime_basis=request["runtime_basis"], current_original_os_inventory_equality=False,
        prospective_current_inventory_equality=True, continuous_runtime_integrity_established=False,
        phases={"before": {}, "after": {}})
    resources = dict(primary={"fixture": "measurement"}, primary_scopes=current.SCOPES)
    environment = dict(sampled_environment_evidence_valid=True, amendment=request["amendment"],
        uncontended_timing=False, background_cpu_used_for_eligibility=False, pressure_thresholds_used_for_eligibility=False)
    replay = dict(schema="native12_composed_resource_replay_summary_v1", full_replay_executed=True,
        measured_matches_retained_wrapper=True, source=record(replay_source), evidence=[])
    for name, data in (("runtime", runtime), ("resources", resources), ("environment", environment),
                       ("outputs_or_failure", output)):
        review["reviews"][name] = save(destination / (name + ".json"), data)
    review.update(request=request_ref, source=source_ref, resources=resources["primary"], native_cpu_ids=[1],
        resource_replay=save(destination / "resource_replay_summary.json", replay),
        scheduler=save(destination / "scheduler.json", {"fixture": "native scheduler"}), evidence=[])
    review_ref = save(destination / "review.json", review)
    stdout = tmp_path / "producer.out"
    stdout.write_text(repr(review_ref) + "\n")
    stderr = tmp_path / "producer.err"
    stderr.write_text("")
    held = dict(schema="native12_composed_review_submission_v1", job_id=producer_job, submission_count=1,
        held_comparison_passed=True, destination=str(destination), native_inference_reexecuted=False,
        automatic_retry=False, references=dict(request=request_ref, worker=source_ref, batch=batch_ref),
        controller=dict(stdout=raw_held(producer_job, StdOut=str(stdout), StdErr=str(stderr))))
    held_ref = save(tmp_path / "held.json", held)
    release = dict(schema="native12_composed_review_release_v1", job_id=producer_job, returncode=0,
        release_count=1, held=held_ref)
    release_ref = save(tmp_path / "release.json", release)
    plan = dict(runs=[{}] * 12 + [run], baseline=ready, helper_sources=[], evidence=[])
    context = ({}, dict(new_sources=[]), plan)
    monkeypatch.setattr(current.executor, "execution_binding", lambda *args: (request, context, dict(evidence=[])))
    monkeypatch.setattr(current.executor, "held_gate", lambda *args: None)
    monkeypatch.setattr(current, "accounting", lambda *args, **kwargs: ("fixture accounting", producer))
    terminal = dict(source="live_controller", verified=dict(fields=dict(
        JobState="COMPLETED", ExitCode="0:0", Comment=request_ref["sha256"])))
    monkeypatch.setattr(current.reviewer, "verify_terminal", lambda *args: terminal)
    monkeypatch.setattr(current, "check", lambda ref: None if "path" not in ref else check_real(ref))
    return dict(args=(request_ref, review_ref, producer_job, held_ref, release_ref), save=save,
        output=output, review=review, runtime=runtime, environment=environment, resources=resources,
        replay=replay, held=held, release=release, terminal=terminal, stdout=stdout)


def check_real(ref):
    from benchmark_tools.prepare_ob_candidate_neighborhood import check
    check(ref)


def test_joined_binding_keeps_new_type_and_returns_checked_evidence(joined):
    result = current.native_binding(*joined["args"])
    assert result[6] == "native"
    assert result[4]["schema"] == current.reviewer.SCHEMA
    assert result[-1]["review_producer_job_id"] == 25000
    assert result[-1]["composed_schema_preserved"] is True
    assert result[-1]["original_review_translated"] is False
    assert result[-1]["next_identity_authorized"] is False
    assert record(joined["stdout"]) in result[-2]


@pytest.mark.parametrize("component,field,new", [
    ("runtime", "prospective_current_inventory_equality", False),
    ("runtime", "current_original_os_inventory_equality", True),
    ("environment", "sampled_environment_evidence_valid", False),
    ("environment", "background_cpu_used_for_eligibility", True),
    ("resources", "primary", {}), ("replay", "full_replay_executed", False),
    ("output", "schema", "ordinary"), ("output", "cell", "p1_c1_r0"),
    ("output", "native_cpu_ids", [2]), ("output", "accuracy_evaluated", True),
    ("output", "phylogeny", {"native_pair_rows": True}),
    ("held", "submission_count", 2), ("release", "release_count", 2),
    ("release", "returncode", 1), ("held", "native_inference_reexecuted", True)])
def test_resealed_invalid_components_still_refuse(joined, component, field, new):
    joined[component][field] = new
    args = list(joined["args"])
    save = joined["save"]
    if component in {"held", "release"}:
        offset = 3 if component == "held" else 4
        args[offset] = save(Path(args[offset]["path"]), joined[component])
        if component == "held":
            joined["release"]["held"] = args[3]
            args[4] = save(Path(args[4]["path"]), joined["release"])
    else:
        key = "outputs_or_failure" if component == "output" else component
        ref = joined["review"]["resource_replay"] if key == "replay" else joined["review"]["reviews"][key]
        updated = save(Path(ref["path"]), joined[component])
        if key == "replay":
            joined["review"]["resource_replay"] = updated
        else:
            joined["review"]["reviews"][key] = updated
        args[1] = save(Path(args[1]["path"]), joined["review"])
        joined["stdout"].write_text(repr(args[1]) + "\n")
    with pytest.raises(ValueError):
        current.native_binding(*args)


def test_changed_output_bytes_without_resealing_are_rejected(joined):
    Path(joined["review"]["reviews"]["outputs_or_failure"]["path"]).write_text("{}\n")
    with pytest.raises(ValueError):
        current.native_binding(*joined["args"])


def test_producer_stdout_must_name_this_review(joined):
    joined["stdout"].write_text(repr({"wrong": "review"}) + "\n")
    with pytest.raises(ValueError, match="stdout"):
        current.native_binding(*joined["args"])
