"""Per-run measurement composition for the frozen replacement scaling plan.

This library does not authorize or submit runs. Its caller must first verify
execution authorization, environmental policy, allocation and recipe identity.
There is deliberately no standalone execution CLI.
"""

import os
from pathlib import Path
import re

from benchmark_tools.launch_dgx_native_run import read_pinned, native_enumerator
from benchmark_tools.measure_scaling_root_context import measure as measure_native_run
from benchmark_tools.measure_native_scaling_run import measure_run
from benchmark_tools.prepare_root_context_scaling import OUTPUT_ROOT
from benchmark_tools.replay_scaling_root_context import native_outcome
from benchmark_tools.verify_lineage_native_provenance import same

PLAN_SHA = "f54790499a48f95e2a866750ebc78433195265591e12c3dc6aff1389e29ed084"


def classify_measurement(verification, replayed, expected_job):
    """Classify retained outcomes; never authorize continuation from wrapper exit alone."""
    if type(expected_job) is not int or expected_job <= 0:
        raise ValueError("Require positive expected scheduler job identity")
    if not isinstance(verification, dict):
        raise ValueError("Require a retained wrapper record")
    reported = verification.get("measurement")
    result = dict(status="measurement_evidence_incomplete", wrapper_status=verification.get("status"),
        reported_native=reported.get("native") if isinstance(reported, dict) else None,
        corroborated_native_outcome=None, next_submission_authorized=False, automatic_retry=False,
        native_outputs_validated=False, runtime_identity_verified=False,
        scientific_timings_admitted=False, environmental_validity_established=False,
        required_followup=["Terminal scheduler and bounded-session verification", "Frozen task/runtime/recipe provenance",
                           "Environmental policy and whole-run evidence", "Native outputs or retained failure audit"])
    if verification.get("status") in {"verified_wrapper_failed", "runtime_changed_or_unverifiable"}:
        return dict(result, status="infrastructure_or_provenance_failure",
                    reason="Wrapper or runtime verification failed; retain any observed native result and pause")
    if replayed is None:
        return dict(result, reason="Independent raw replay is missing or failed; do not infer a native-only failure")
    try:
        if verification["scientific_results_admitted"] is not False:
            raise ValueError("Unexpected wrapper admission")
        if not isinstance(verification["before"], dict) or not isinstance(verification["after"], dict):
            raise ValueError("Missing before/after verification records")
        measured = verification["measurement"]
        if (replayed["status"] != "scaling_root_context_measurement_replayed"
                or not same(replayed["measured"], measured)):
            raise ValueError("Replay does not bind the wrapper measurement")
        if type(measured["job_id"]) is not int or measured["job_id"] != expected_job:
            raise ValueError("Measured scheduler identity differs")
        if any(replayed[k] is not False for k in (
                "scientific_timings_admitted", "environmental_validity_established", "native_outputs_validated", "publication_ready")):
            raise ValueError("Unexpected replay admission")
        outcome, wall = native_outcome(measured, measured["native"])
        if (verification["status"] != measured["status"] or replayed["native_outcome"] != outcome
                or not same(replayed["native_exit_code"], measured["native"]["exit_code"])
                or not same(replayed["native_wall_s"], wall)):
            raise ValueError("Wrapper, native outcome and replay summary disagree")
    except (KeyError, TypeError, ValueError) as error:
        return dict(result, status="infrastructure_or_provenance_failure", reason=str(error))
    return dict(result, status="native_" + outcome, corroborated_native_outcome=outcome,
                reason="Native outcome corroborated by raw replay; remaining checks still gate continuation")


def load_task(plan_path, index):
    if type(index) is not int or not 0 <= index < 27:
        raise ValueError("Require original task index in [0, 26]")
    plan = read_pinned(plan_path, PLAN_SHA)
    if (plan["collector"]["module"] != measure_native_run.__module__
            or plan["collector"]["entry"] != measure_native_run.__name__):
        raise ValueError("Collector import differs from frozen specification")
    task = plan["runs"][index]
    if type(task["index"]) is not int or task["index"] != index:
        raise ValueError("Task identity differs")
    matching = [order for order in plan["orders"]
                if order["input_directory"] == task["run"]["dataset"]["input_directory"]]
    if len(matching) != 1:
        raise ValueError("Require one native input order for this dataset")
    return plan, task, matching[0]


def measure_task(plan_path, index, recipe_path, recipe_sha, job):
    if type(job) is not int or job <= 0:
        raise ValueError("Require positive scheduler job identity")
    if not isinstance(recipe_sha, str) or not re.fullmatch(r"[0-9a-f]{64}", recipe_sha):
        raise ValueError("Require recipe SHA-256")
    plan, task, order = load_task(plan_path, index)
    recipe_path = Path(recipe_path).resolve()
    run = task["run"]
    cache = OUTPUT_ROOT / f"cache_{index:02d}"
    if cache.exists() or cache.is_symlink():
        raise ValueError("Task cache prefix must be absent, including dangling links")
    overrides = dict(plan["environment_overrides"],
        PATH=os.pathsep.join(plan["environment_paths"][run["environment_role"]]),
        PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(cache))
    keys = set(overrides) | set(plan["unset_environment"])
    previous, cwd = {key: os.environ.get(key) for key in keys}, Path.cwd()
    try:
        for key in plan["unset_environment"]:
            os.environ.pop(key, None)
        os.environ.update(overrides)
        os.chdir(run["cwd"])
        specifications = [(row["path"], row["sha256"]) for row in plan["runtime_manifests"]]
        specifications.append((str(recipe_path), recipe_sha))
        return measure_run(run, order, specifications, native_enumerator(plan["enumerator"]),
            measure_native_run, job, cpus=plan["allocation"]["cpus"],
            memory_gib=plan["allocation"]["memory_gib"], timeout_s=plan["native_timeout_s"])
    finally:
        os.chdir(cwd)
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
