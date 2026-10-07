from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import native_factorial_allocated_execution as contract
from benchmark_tools import prepare_allocated_native_factorial_request as preparation


def ref(path, digest="a"*64):
    return dict(path=str(path), bytes=5, sha256=digest)


@pytest.fixture
def execution(monkeypatch):
    prefix = [ref(f"/synthetic/review_{i}.json") for i in range(10)]
    sources = [ref(contract.ROOT / "benchmark_tools" / name) for name in contract.SOURCE_NAMES]
    plan = dict(helper_sources=[], evidence=[])
    value = dict(schema="allocated_native_factorial_execution_amendment_v1", root=str(contract.ROOT),
        execution_scope=contract.SCOPE, allowed_indices=[10, 11, 12], placement_policy="actual_allocated_physical_cores_v1",
        resources=deepcopy(contract.RESOURCES), automatic_retry=False, scientific_method_modified=False,
        historical_evidence_modified=False, production_execution_authorized=True, qfo_conversion_scoring_ready=True,
        scientific_timings_admitted=False, publication_ready=False,
        historical_plan=ref(contract.PLAN, contract.PLAN_SHA), placement_fixture=ref(contract.FIXTURE, contract.FIXTURE_SHA),
        new_sources=sources, historical_prefix=deepcopy(prefix))
    def record(path):
        return ref(path, contract.PLAN_SHA if path == contract.PLAN else
                   contract.FIXTURE_SHA if path == contract.FIXTURE else "a"*64)
    def read(r):
        if r == value["historical_plan"]:
            return plan
        if r == value["placement_fixture"]:
            return dict(status="engineering_fixture_terminal_and_replay_verified", job_id=23897,
                independent_replay_matches=True, independent_resource_arithmetic_matches=True,
                current_sysfs_topology_matches=True)
        return dict(index=prefix.index(r), terminal_reviewed=True, next_identity_authorized=True)
    monkeypatch.setattr(contract, "record", record)
    monkeypatch.setattr(contract, "read", read)
    monkeypatch.setattr(contract, "check", lambda *a: None)
    monkeypatch.setattr(contract, "sources", lambda: deepcopy(sources))
    monkeypatch.setattr(contract, "historical_prefix", lambda: deepcopy(prefix))
    monkeypatch.setattr(contract, "validate_plan", lambda p: [])
    return value, plan


def test_only_prospective_placement_changes(execution):
    value, plan = execution
    assert contract.validate_amendment(value) is plan
    assert value["resources"]["native_physical_cores"] == 32
    assert value["historical_plan"]["sha256"] == contract.PLAN_SHA


@pytest.mark.parametrize("key,value", [("schema", "historical"), ("root", "/other"),
    ("execution_scope", "isolated"), ("allowed_indices", [9, 10, 11, 12]),
    ("placement_policy", "fixed_mask"), ("automatic_retry", True), ("scientific_method_modified", True),
    ("historical_evidence_modified", True), ("production_execution_authorized", False),
    ("qfo_conversion_scoring_ready", False), ("scientific_timings_admitted", True), ("publication_ready", True),
    ("historical_prefix", []), ("new_sources", [])])
def test_wrong_amendment_scope_refused(execution, key, value):
    execution[0][key] = value
    with pytest.raises(ValueError):
        contract.validate_amendment(execution[0])


@pytest.mark.parametrize("key,value", [("native_physical_cores", 64), ("slurm_slots", 32),
    ("memory_bytes", 64*1024**3), ("timeout_s", 900), ("sample_period_s", True),
    ("host_period_s", 60), ("minimum_available_memory_bytes", 0)])
def test_resource_mutation_refused(execution, key, value):
    execution[0]["resources"][key] = value
    with pytest.raises(ValueError):
        contract.validate_amendment(execution[0])


@pytest.mark.parametrize("key", ["historical_plan", "placement_fixture"])
def test_stale_historical_receipt_refused(execution, key):
    execution[0][key]["sha256"] = "0"*64
    with pytest.raises(ValueError):
        contract.validate_amendment(execution[0])


def test_reordered_passing_prefix_is_not_substitutable(execution):
    execution[0]["historical_prefix"].reverse()
    with pytest.raises(ValueError):
        contract.validate_amendment(execution[0])


def test_missing_downstream_file_refuses_binding(monkeypatch):
    def record(path):
        if path.name == "admit_allocated_native_factorial_qfo_assessment.py":
            raise FileNotFoundError(path)
        return ref(path)
    monkeypatch.setattr(contract, "record", record)
    with pytest.raises(FileNotFoundError):
        contract.sources()


def request(execution, index=10):
    value = execution[0]
    amendment = ref("/amendment.json")
    return dict(schema="allocated_native_factorial_request_v1", execution_authorized=True,
        job_id=30000, index=index, plan=value["historical_plan"], amendment=amendment,
        scheduler_command=str(contract.SCRIPT), allocation_cwd=str(contract.ROOT),
        history=deepcopy(value["historical_prefix"])+[ref(f"/new_review_{i}.json") for i in range(10,index)],
        automatic_retry=False), amendment


@pytest.mark.parametrize("index", [10, 11, 12])
def test_remaining_requests_preserve_prefix(execution, index):
    value, amendment = request(execution, index)
    contract.validate_request(value, amendment, execution[0], 30000)


@pytest.mark.parametrize("key,value", [("schema", "native_factorial_cost_request_v1"), ("index", 9),
    ("index", True), ("index", 13), ("job_id", True), ("job_id", 29999), ("execution_authorized", False),
    ("automatic_retry", True), ("plan", {}), ("amendment", {}), ("scheduler_command", "/old.sh"),
    ("allocation_cwd", "/other"), ("history", [])])
def test_request_rejects_retries_or_mismatched_route(execution, key, value):
    result, amendment = request(execution)
    result[key] = value
    with pytest.raises(ValueError):
        contract.validate_request(result, amendment, execution[0], 30000)


def test_native_command_and_metrics_command_distinction():
    amendment = ref("/amendment.json", "b"*64)
    baseline = dict(tool_entrypoints=dict(orthohmm_python=dict(absolute_path="/private/python")))
    args = contract.native_command(amendment, dict(index=10), baseline)
    assert args[:2] == ["/private/python", "-B"]
    assert args[2] == str(contract.ROOT / "benchmark_tools/run_allocated_native_factorial_cost.py")
    assert "--plan" not in args and args[-2:] == ["--index", "10"]
    assert contract.native_command(amendment, dict(index=10), baseline, metrics=True) == [args[0], *args[2:]]


def accounting():
    return "30000|COMPLETED|0:0|64|128G|bizon|gpu|1560|2026-10-06T20:00:00|2026-10-06T20:01:00|2026-10-07T10:01:00|orthohmm_allocated_factorial\n"


def test_terminal_accounting_new_route_only():
    assert contract.terminal_accounting(accounting(), 30000)["State"] == "COMPLETED"
    with pytest.raises(ValueError):
        contract.terminal_accounting(accounting().replace("orthohmm_allocated_factorial", "orthohmm_factorial_cost"), 30000)


@pytest.mark.parametrize("old,new", [("30000|", "30001|"), ("COMPLETED", "RUNNING"), ("|64|", "|32|"),
    ("128G", "64G"), ("bizon", "other"), ("1560", "1200"), ("2026-10-07T10:01:00", "Unknown")])
def test_accounting_changed_identity_refused(old, new):
    with pytest.raises(ValueError):
        contract.terminal_accounting(accounting().replace(old,new), 30000)


def held(nodes="1-1"):
    return (f"JobId=30000 JobName=orthohmm_allocated_factorial JobState=PENDING Reason=JobHeldUser "
        f"Partition=gpu ReqNodeList=bizon NumCPUs=64 NumTasks=1 MinMemoryNode=128G Requeue=0 Restarts=0 "
        f"TimeLimit=1-02:00:00 CPUs/Task=64 Command={contract.SCRIPT} WorkDir={contract.ROOT} "
        f"NumNodes={nodes} UserId=bizon({preparation.os.getuid()})")


@pytest.mark.parametrize("nodes", ["1", "1-1"])
def test_held_single_node_representations(nodes):
    assert preparation.held_job(held(nodes), 30000)["NumNodes"] == nodes


@pytest.mark.parametrize("change", ["node", "owner", "array", "duplicate", "script", "requeue"])
def test_held_owner_and_envelope_refusals(change):
    raw = held()
    if change == "node": raw = held("2-2")
    elif change == "owner": raw = raw.replace(f"({preparation.os.getuid()})", "(99999)")
    elif change == "array": raw += " ArrayJobId=30000"
    elif change == "duplicate": raw += " JobId=30000"
    elif change == "script": raw = raw.replace(str(contract.SCRIPT), "/old.sh")
    else: raw = raw.replace("Requeue=0", "Requeue=1")
    with pytest.raises(ValueError):
        preparation.held_job(raw, 30000)


def test_poll_failure_is_not_controller_expiry(monkeypatch):
    monkeypatch.setattr(contract.subprocess, "run", lambda *a, **k:
        SimpleNamespace(returncode=1, stdout="", stderr="Connection failure"))
    with pytest.raises(ValueError, match="not expired"):
        contract.verify_terminal(30000)


@pytest.fixture
def mixed_history(monkeypatch):
    prefix = [ref(f"/prior_{i}.json") for i in range(10)]
    new_ref = ref("/new_10.json")
    plan_ref, amendment_ref = ref("/plan.json"), ref("/amendment.json")
    plan = dict(runs=[{}]*10+[dict(dataset="qfo_corrected",cell="p1_c0_r1",repeat=0)])
    reviews = {r["path"]:dict(index=i,job_id=10000+i,plan=plan_ref,terminal_reviewed=True,
        next_identity_authorized=True,scheduler_state="COMPLETED",scheduler_exit_code="0:0")
        for i,r in enumerate(prefix)}
    prior = dict(schema="allocated_native_factorial_terminal_review_v1",status="native_success",index=10,
        job_id=30000,plan=plan_ref,amendment=amendment_ref,dataset="qfo_corrected",cell="p1_c0_r1",repeat=0,
        execution_scope=contract.SCOPE,automatic_retry=False,scientific_timings_admitted=False,
        terminal_reviewed=True,next_identity_authorized=True,scheduler_state="COMPLETED",scheduler_exit_code="0:0",
        source=contract.record(contract.ROOT/"benchmark_tools/review_allocated_native_factorial_attempt.py"))
    reviews[new_ref["path"]]=prior
    terminal=dict(verified=dict(fields=dict(JobState="COMPLETED",ExitCode="0:0")))
    calls=[]
    monkeypatch.setattr(contract,"read",lambda r:deepcopy(reviews[r["path"]]))
    monkeypatch.setattr(contract,"historical_terminal",lambda job:(calls.append(("old",job)) or deepcopy(terminal)))
    monkeypatch.setattr(contract,"verify_terminal",lambda job:(calls.append(("new",job)) or deepcopy(terminal)))
    return dict(request=dict(index=11,history=prefix+[new_ref]),amendment=amendment_ref,
        execution=dict(historical_plan=plan_ref,historical_prefix=prefix),plan=plan,prior=prior,calls=calls)


def test_mixed_prefix_uses_actual_matching_terminal_routes(mixed_history):
    d=mixed_history
    evidence=contract.reviewed_history(d["request"],d["amendment"],d["execution"],d["plan"])
    assert len(evidence)==11
    assert d["calls"][:10]==[("old",10000+i) for i in range(10)]
    assert d["calls"][-1]==("new",30000)


@pytest.mark.parametrize("key,value",[("schema","native_factorial_terminal_review_v1"),("amendment",{}),
    ("dataset","orthobench"),("cell","p1_c1_r0"),("repeat",True),("execution_scope","isolated"),
    ("automatic_retry",True),("scientific_timings_admitted",True),("next_identity_authorized",False),
    ("status","unclassified_infrastructure_failure"),("source",{}),("scheduler_state","FAILED")])
def test_mixed_history_rejects_forged_or_unresolved_success(mixed_history,key,value):
    d=mixed_history;d["prior"][key]=value
    with pytest.raises(ValueError):
        contract.reviewed_history(d["request"],d["amendment"],d["execution"],d["plan"])
