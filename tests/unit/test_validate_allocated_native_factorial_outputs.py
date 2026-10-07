from copy import deepcopy
import json

import pytest

from benchmark_tools import validate_allocated_native_factorial_outputs as validator
from benchmark_tools import run_allocated_native_factorial_cost as controller
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_run_allocated_native_factorial_cost import native_case


@pytest.fixture
def output_case(native_case, monkeypatch):
    d = native_case
    controller.native(d["amendment"], 10)
    execution, plan = d["execution"], d["plan"]
    execution["new_sources"] = []
    plan["helper_sources"] = []
    d["run"]["inputs"] = []
    d["run"]["proteomes"] = 1
    request = dict(amendment=d["amendment"], plan=execution["historical_plan"], job_id=42, index=10)
    request_path = d["root"].parent / "request.json"
    request_path.write_text(json.dumps(request))
    request_ref = record(request_path)
    semantic_source = record(validator.ROOT / "benchmark_tools/validate_native_factorial_outputs.py")
    semantic = dict(schema="native_factorial_output_review_v1", status="native_outputs_validated", cell=d["run"]["cell"],
        factors=d["flags"], gene_ownership_sha256="ownership-digest", per_species_counts={"s0.fa":8},
        source=semantic_source, checked_files=[], native_outputs_validated=True, accuracy_evaluated=False,
        resource_measurements_admitted=False, next_identity_authorized=False)
    preparation = dict(status="fresh_factorial_inputs_prepared", gene_ownership_sha256="ownership-digest",
        per_species_counts={"s0.fa":8}, genes=8)
    (d["root"] / "preparation.json").write_text(json.dumps(preparation))
    terminal = dict(source="live_controller",verified=dict(fields=dict(JobState="COMPLETED",ExitCode="0:0",
        Comment=request_ref["sha256"])))
    monkeypatch.setattr(validator, "amendment", lambda ref:(execution,plan))
    monkeypatch.setattr(validator, "validate_request", lambda *a:None)
    monkeypatch.setattr(validator, "verify_terminal", lambda job:deepcopy(terminal))
    original_check = validator.check
    monkeypatch.setattr(validator, "check", lambda ref:None if ref["path"]=="/plan" else original_check(ref))
    contexts = []
    def validate_semantics(context):
        contexts.append(context)
        return deepcopy(semantic)
    monkeypatch.setattr(validator, "validate_semantics", validate_semantics)
    return dict(d, request=request_ref, semantic=semantic, contexts=contexts, terminal=terminal)


def test_exact_semantic_kernel_and_new_provenance(output_case):
    d = output_case
    result = validator.validate(d["request"])
    assert result["schema"] == "allocated_native_factorial_output_review_v1"
    assert result["semantic_validator_source"] == d["semantic"]["source"]
    assert result["source"] == record(validator.__file__)
    assert result["amendment"] == d["amendment"] and result["native_cpu_ids"] == list(range(32,64))
    assert result["native_outputs_validated"] is True and result["accuracy_evaluated"] is False
    context = d["contexts"][0]
    assert context["cpu"]==32 and context["threads_per_worker"]==4
    assert "--amendment" in context["command"] and "--plan" not in context["command"]
    assert "-B" not in context["command"]
    assert context["cell"]=="p1_c0_r1" and context["genes"]==8


@pytest.mark.parametrize("change", ["live", "failed", "comment", "old_native", "amendment", "cpu", "parent", "source",
                                    "preparation", "factors"])
def test_output_gate_rejects_route_and_provenance_tampering(output_case, change):
    d = output_case
    if change in {"live", "failed", "comment"}:
        fields=d["terminal"]["verified"]["fields"]
        if change=="live":fields["JobState"]="RUNNING"
        elif change=="failed":fields.update(JobState="FAILED",ExitCode="1:0")
        else:fields["Comment"]="wrong"
    elif change=="preparation":
        path=d["root"] / "preparation.json"
        value=json.loads(path.read_text());value["gene_ownership_sha256"]="changed"
        path.write_text(json.dumps(value))
    else:
        path=d["root"] / "native_execution.json"
        value=json.loads(path.read_text())
        if change=="old_native":value["schema"]="historical"
        elif change=="amendment":value["amendment"]={}
        elif change=="cpu":value["native_cpu_ids"]=list(range(32))
        elif change=="parent":value["parent_pid"]=999
        elif change=="source":value["source"]={}
        else:value["factors"]={}
        path.write_text(json.dumps(value))
    with pytest.raises(ValueError):
        validator.validate(d["request"])
    if change not in {"preparation","factors"}:
        assert d["contexts"] == []
