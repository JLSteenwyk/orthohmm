"""Prospective allocated-core execution contract; historical plan stays immutable."""

import csv
from pathlib import Path
import re
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_native_factorial_cost import (
    ROOT, SCOPE, MEMORY, read, validate_plan, verify_terminal as historical_terminal)
from benchmark_tools.verify_lineage_native_provenance import same
from benchmark_tools.verify_threadripper_controller import validate as validate_controller


SCRIPT = ROOT / "benchmark_tools/run_allocated_native_factorial_cost.sh"
PLAN = ROOT / "benchmark_tools/results/native_factorial_receipt_amendment_20261004/plan.json"
PLAN_SHA = "6c87babcbb5581830e0b9e7b9bf9aaba30a85bde1c4ab465e561017e67e9c89b"
FIXTURE = ROOT / "benchmark_tools/results/allocated_threadripper_fixture_terminal_23897.json"
FIXTURE_SHA = "1bf2389a1a9761a16e5ad423f5f167367a2addebd41cc2c708b123629597aa7c"
PREFIX_REQUEST = ROOT / "benchmarks/work/native_factorial_launch_20261004/request_09_receipt_amended.json"
PREFIX_REQUEST_SHA = "ece9eb34c8087cecb3ee0d563e1a8db8f1f12ce9a928127419ee04c262dbba6e"
PREFIX_FAILURE = ROOT / "benchmarks/work/native_factorial_launch_failure_review_23891/review.json"
PREFIX_FAILURE_SHA = "1f6167cf997f37dcff5f795c9479475bd065aa7359a68a73994f7f53abf4ed9a"
RESOURCES = dict(native_physical_cores=32, slurm_slots=64, memory_bytes=MEMORY,
    timeout_s=85800, sample_period_s=1., host_period_s=30., minimum_available_memory_bytes=MEMORY)
SOURCE_NAMES = (
    "native_factorial_allocated_execution.py", "run_allocated_native_factorial_cost.py",
    "run_allocated_native_factorial_cost.sh", "prepare_allocated_native_factorial_request.py",
    "prepare_allocated_native_factorial_amendment.py",
    "validate_allocated_native_factorial_outputs.py", "review_allocated_native_factorial_attempt.py",
    "prepare_allocated_native_factorial_qfo_pairs.py", "run_allocated_native_factorial_qfo_assessment.py",
    "admit_allocated_native_factorial_qfo_assessment.py", "native_factorial_allocated_placement.py",
    "native_factorial_cpu_selection.py", "measure_allocated_threadripper_scaling.py",
    "replay_allocated_threadripper_scaling.py")


def sources():
    # Missing downstream modules fail before any production amendment can freeze.
    return [record(ROOT / "benchmark_tools" / name) for name in SOURCE_NAMES]


def historical_prefix():
    request_ref, failure_ref = record(PREFIX_REQUEST), record(PREFIX_FAILURE)
    if request_ref["sha256"] != PREFIX_REQUEST_SHA or failure_ref["sha256"] != PREFIX_FAILURE_SHA:
        raise ValueError("Pinned historical request or explicit failure disposition changed")
    return read(request_ref)["history"] + [failure_ref]


def validate_amendment(value):
    if (value.get("schema") != "allocated_native_factorial_execution_amendment_v1"
            or value.get("root") != str(ROOT) or value.get("execution_scope") != SCOPE
            or value.get("allowed_indices") != [10, 11, 12]
            or value.get("placement_policy") != "actual_allocated_physical_cores_v1"
            or not same(value.get("resources"), RESOURCES)
            or value.get("automatic_retry") is not False
            or value.get("scientific_method_modified") is not False
            or value.get("historical_evidence_modified") is not False
            or value.get("production_execution_authorized") is not True
            or value.get("qfo_conversion_scoring_ready") is not True
            or value.get("scientific_timings_admitted") is not False
            or value.get("publication_ready") is not False):
        raise ValueError("Wrong prospective placement scope, resources or readiness")
    if (not same(value.get("historical_plan"), record(PLAN))
            or value["historical_plan"]["sha256"] != PLAN_SHA
            or not same(value.get("placement_fixture"), record(FIXTURE))
            or value["placement_fixture"]["sha256"] != FIXTURE_SHA
            or not same(value.get("new_sources"), sources())):
        raise ValueError("Historical plan, actual fixture or complete route bindings differ")
    prefix = value.get("historical_prefix")
    if not isinstance(prefix, list) or len(prefix) != 10 or prefix != historical_prefix():
        raise ValueError("Require the ten reviewed historical identities")
    for index, ref in enumerate(prefix):
        prior = read(ref)
        if (type(prior.get("index")) is not int or prior["index"] != index
                or prior.get("terminal_reviewed") is not True
                or prior.get("next_identity_authorized") is not True):
            raise ValueError("Historical prefix is unresolved or reordered")
    fixture = read(value["placement_fixture"])
    if (fixture.get("status") != "engineering_fixture_terminal_and_replay_verified"
            or fixture.get("job_id") != 23897
            or any(fixture.get(k) is not True for k in (
                "independent_replay_matches", "independent_resource_arithmetic_matches",
                "current_sysfs_topology_matches"))):
        raise ValueError("Missing actual allocation-aware fixture verification")
    plan = read(value["historical_plan"])
    validate_plan(plan)
    for ref in [value["historical_plan"], value["placement_fixture"], *value["new_sources"],
                *prefix, *plan["helper_sources"], *plan["evidence"]]:
        check(ref)
    return plan


def amendment(ref):
    value = read(ref)
    return value, validate_amendment(value)


def validate_request(value, amendment_ref, execution, job):
    index = value.get("index")
    if (value.get("schema") != "allocated_native_factorial_request_v1"
            or value.get("execution_authorized") is not True
            or type(job) is not int or job <= 0 or type(value.get("job_id")) is not int
            or value["job_id"] != job or type(index) is not int or index not in (10, 11, 12)
            or value.get("amendment") != amendment_ref
            or value.get("plan") != execution["historical_plan"]
            or value.get("scheduler_command") != str(SCRIPT)
            or value.get("allocation_cwd") != str(ROOT)
            or value.get("automatic_retry") is not False
            or not isinstance(value.get("history"), list) or len(value["history"]) != index
            or value["history"][:10] != execution["historical_prefix"]):
        raise ValueError("Wrong allocated request, frozen identity, job or historical prefix")


def native_command(amendment_ref, run, baseline, *, metrics=False):
    args = [baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"]]
    if not metrics:
        args.append("-B")
    return args + [str(ROOT / "benchmark_tools/run_allocated_native_factorial_cost.py"), "--native",
        "--amendment", amendment_ref["path"], "--amendment-sha256", amendment_ref["sha256"],
        "--index", str(run["index"])]


def terminal_accounting(raw, job):
    rows = list(csv.reader(raw.strip().splitlines(), delimiter="|"))
    names = ("JobIDRaw", "State", "ExitCode", "AllocCPUS", "ReqMem", "NodeList", "Partition",
             "TimelimitRaw", "Submit", "Start", "End", "JobName")
    if len(rows) != 1 or len(rows[0]) != len(names):
        raise ValueError("Require one full allocated-native accounting row")
    fields = dict(zip(names, rows[0]))
    from benchmark_tools.capture_array_scheduler import TERMINAL
    if (fields["JobIDRaw"] != str(job) or fields["State"] not in TERMINAL
            or not re.fullmatch(r"[0-9]+:[0-9]+", fields["ExitCode"])
            or fields["AllocCPUS"] != "64"
            or fields["ReqMem"] not in {"128G", "128Gn", "131072M", "131072Mn"}
            or fields["NodeList"] != "bizon" or fields["Partition"] != "gpu"
            or fields["TimelimitRaw"] != "1560"
            or fields["JobName"] != "orthohmm_allocated_factorial"
            or any(not re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}", fields[k])
                   for k in ("Submit", "Start", "End"))):
        raise ValueError("Allocated-native accounting identity/envelope differs")
    return fields


def verify_terminal(job):
    query = ["scontrol", "show", "job", str(job), "--oneliner"]
    result = subprocess.run(query, capture_output=True, text=True, timeout=5)
    observation = dict(command=query, returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
    if result.returncode == 0:
        verified = validate_controller(result.stdout, job, "terminal", command=str(SCRIPT), cwd=str(ROOT),
                                       time_limit="1-02:00:00", allocation_mode="shared")
        if verified["fields"].get("JobName") != "orthohmm_allocated_factorial":
            raise ValueError("Allocated-native job name differs")
        return dict(source="live_controller", observation=observation, verified=verified)
    if "Invalid job id specified" not in result.stderr:
        raise ValueError("Controller query failed, not expired evidence")
    query = ["sacct", "-X", "-j", str(job), "-n", "-P",
        "--format=JobIDRaw,State,ExitCode,AllocCPUS,ReqMem,NodeList,Partition,TimelimitRaw,Submit,Start,End,JobName"]
    result = subprocess.run(query, capture_output=True, text=True, check=True, timeout=5)
    return dict(source="fresh_accounting_after_controller_expiry", controller_observation=observation,
        observation=dict(command=query, stdout=result.stdout, stderr=result.stderr),
        verified=terminal_accounting(result.stdout, job),
        limitation="Terminal accounting lacks historical command/comment; bound request and held submission remain required.")


def reviewed_history(request, amendment_ref, execution, plan):
    from benchmark_tools.amend_native_factorial_receipt import adopt
    evidence = []
    plan_ref = execution["historical_plan"]
    for index, ref in enumerate(request["history"]):
        prior = read(ref)
        if (type(prior.get("index")) is not int or prior["index"] != index
                or prior.get("terminal_reviewed") is not True
                or prior.get("next_identity_authorized") is not True):
            raise ValueError("Previous native identity is unresolved")
        adoption = None
        if index < 10:
            if ref != execution["historical_prefix"][index]:
                raise ValueError("Historical prefix changed")
            if prior.get("plan") != plan_ref:
                adoption = adopt(ref, prior, plan)
            terminal = historical_terminal(prior["job_id"])
        else:
            run = plan["runs"][index]
            if (prior.get("schema") != "allocated_native_factorial_terminal_review_v1"
                    or prior.get("amendment") != amendment_ref or prior.get("plan") != plan_ref
                    or prior.get("dataset") != run["dataset"] or prior.get("cell") != run["cell"]
                    or type(prior.get("repeat")) is not int or prior["repeat"] != run["repeat"]
                    or prior.get("execution_scope") != SCOPE or prior.get("automatic_retry") is not False
                    or prior.get("scientific_timings_admitted") is not False
                    or prior.get("status") not in {"native_success", "native_failure_retained", "native_timeout_retained"}
                    or prior.get("source") != record(ROOT / "benchmark_tools/review_allocated_native_factorial_attempt.py")):
                raise ValueError("Previous allocated-native review route differs")
            terminal = verify_terminal(prior["job_id"])
        fields = terminal["verified"].get("fields", terminal["verified"])
        if (fields.get("JobState", fields.get("State")) != prior.get("scheduler_state")
                or fields.get("ExitCode") != prior.get("scheduler_exit_code")):
            raise ValueError("Historical scheduler outcome changed")
        evidence.append(dict(review=ref, scheduler=terminal, adoption=adoption))
    return evidence
