"""Preserve reviewed index zero while fixing only prospective receipt/control code."""

import argparse
import ast
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import ROOT, read, validate_plan


OLD_PLAN_SHA = "61756d74d1778fd268c9c4fe4f85b32cc130121f9bf02793f470b202bc4c2b1a"
OLD_REVIEW_SHA = "c3e55f042a230e54b531f429e207a9f2928709a66695a308869dd7feff90ef13"
RECOVERY_SHA = "0c919cfb97f497a98ca7ad6e614cb370731ddb58a63bdc4550fa0e640f4cbc15"
OLD_EXECUTOR_SHA = "07a926e38ba49144b56bf9a7d23478d3dd050fc3a14b98e1ce8889d1f98a8d20"
EXECUTOR = ROOT / "benchmark_tools/run_native_factorial_cost.py"
ARCHIVE = ROOT / "benchmark_tools/results/native_factorial_failed_source_22427/run_native_factorial_cost.py"
CHANGED_FIELDS = {"source_commit", "helper_sources", "evidence", "history_adoption"}
ADDED_HELPERS = {str(ROOT / "benchmark_tools" / name) for name in (
    "amend_native_factorial_receipt.py", "review_native_factorial_attempt.py",
    "recover_native_factorial_receipt_22427.py", "validate_native_factorial_outputs.py")}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def historical_check(ref):
    if ref["path"] == str(EXECUTOR):
        archived = record(ARCHIVE)
        require(ref["sha256"] == archived["sha256"] == OLD_EXECUTOR_SHA
            and ref["bytes"] == archived["bytes"], "Historical executor archive differs")
    else:
        check(ref)


def scientific_source_unchanged():
    old = {node.name: node for node in ast.parse(ARCHIVE.read_text()).body if isinstance(node, ast.FunctionDef)}
    current = {node.name: node for node in ast.parse(EXECUTOR.read_text()).body if isinstance(node, ast.FunctionDef)}
    for name in old:
        if name not in {"execute", "native"}:
            require(name in current and ast.dump(old[name]) == ast.dump(current[name]),
                    "Unexpected executor change outside receipt/history control: " + name)
    native = deepcopy(current["native"])
    initial = ast.parse("initial_receipt = dict(receipt)").body[0]
    require(sum(ast.dump(node) == ast.dump(initial) for node in native.body) == 1,
            "Unexpected running-receipt snapshot")
    native.body = [node for node in native.body if ast.dump(node) != ast.dump(initial)]
    final = native.body[-1]
    require(isinstance(final, ast.Try) and len(final.finalbody) == 1
        and ast.dump(final.finalbody[0]) == ast.dump(ast.parse(
            "update_native_receipt(report, initial_receipt, receipt)").body[0]), "Unexpected native terminal control")
    final.finalbody = ast.parse("save(report, receipt)").body
    require(ast.dump(native) == ast.dump(old["native"]), "Scientific native invocation changed")


def compare_plans(old, new):
    validate_plan(old)
    validate_plan(new)
    require({k: v for k, v in old.items() if k not in CHANGED_FIELDS}
        == {k: v for k, v in new.items() if k not in CHANGED_FIELDS},
        "Amendment changed frozen identities, paths, inputs, settings or resources")
    before = {r["path"]: r for r in old["helper_sources"]}
    after = {r["path"]: r for r in new["helper_sources"]}
    require(set(before) <= set(after) and set(after) - set(before) <= ADDED_HELPERS
        and {p for p in before if before[p] != after[p]} == {str(EXECUTOR)},
        "Unexpected historical helper removal/change")
    require(new["evidence"][:len(old["evidence"])] == old["evidence"], "Historical plan evidence removed")
    for ref in old["helper_sources"]:
        historical_check(ref)
    for ref in new["helper_sources"]:
        check(ref)
    scientific_source_unchanged()


def adopt(prior_ref, prior, plan):
    amendment = plan.get("history_adoption", {})
    require(amendment.get("schema") == "native_factorial_receipt_history_adoption_v1"
        and amendment.get("prior_review") == prior_ref and prior_ref["sha256"] == OLD_REVIEW_SHA
        and amendment.get("prior_plan") == prior.get("plan")
        and prior["plan"]["sha256"] == OLD_PLAN_SHA and amendment.get("recovery", {}).get("sha256") == RECOVERY_SHA
        and amendment.get("executed_source") == record(ARCHIVE), "Unapproved cross-plan history adoption")
    require(read(prior_ref) == prior and prior.get("index") == 0 and prior.get("job_id") == 22427
        and prior.get("dataset") == "orthobench" and prior.get("cell") == "p0_c0_r0"
        and prior.get("status") == "native_failure_retained" and prior.get("scheduler_state") == "FAILED"
        and prior.get("scheduler_exit_code") == "1:0" and prior.get("terminal_reviewed") is True
        and prior.get("next_identity_authorized") is True and prior.get("primary_resources_replayed") is True
        and prior.get("shared_host_resources_reviewed") is True, "Unreviewed or different historical attempt")
    old = read(prior["plan"])
    compare_plans(old, plan)
    recovery = read(amendment["recovery"])
    require(recovery.get("terminal_review") == prior_ref and recovery.get("plan") == prior["plan"]
        and recovery.get("request") == prior.get("request") and recovery.get("job_id") == 22427
        and recovery.get("index") == 0 and recovery.get("cell") == "p0_c0_r0"
        and recovery.get("status") == "scientific_outputs_recovered_from_failed_wrapper"
        and recovery.get("native_outputs_validated") is True and recovery.get("accuracy_evaluated") is True
        and all(recovery.get(k) is False for k in ("scheduler_success", "native_command_success",
            "timing_success_established", "original_receipts_rewritten", "inference_reexecuted", "automatic_retry")),
        "Recovery changes failed-attempt identity/scope")
    references = [prior_ref, prior["plan"], prior["request"], prior["source"], prior["scheduler"],
        prior["resource_replay"], *prior["reviews"].values(), *prior["evidence"],
        amendment["recovery"], recovery["source"], recovery["outputs"], recovery["score"], *recovery["evidence"]]
    unique = {}
    for ref in references:
        require(unique.get(ref["path"], ref) == ref, "Conflicting historical evidence pins")
        unique[ref["path"]] = ref
    for ref in unique.values():
        historical_check(ref)
    require(read(recovery["outputs"]).get("input_genes") == old["runs"][0]["genes"], "Recovered universe differs")
    return dict(status="reviewed_failed_index_zero_adopted_without_retry", prior_review=prior_ref,
        prior_plan=prior["plan"], recovery=amendment["recovery"], executed_source=amendment["executed_source"],
        retained_evidence_checked=len(unique), scientific_settings_unchanged=True,
        scheduler_success=False, timing_success_established=False, inference_reexecuted=False)


def prepare(output_directory):
    output_directory = Path(output_directory)
    require(output_directory.is_absolute() and output_directory.resolve() == output_directory
        and output_directory.is_relative_to(ROOT / "benchmark_tools/results"), "Require a direct results destination")
    if output_directory.exists():
        raise FileExistsError(output_directory)
    recovery_ref = record(ROOT / "benchmark_tools/results/native_factorial_recovery_22427.json")
    require(recovery_ref["sha256"] == RECOVERY_SHA, "Recovery evidence differs")
    recovery = read(recovery_ref)
    old_ref, review_ref = recovery["plan"], recovery["terminal_review"]
    old, prior = read(old_ref), read(review_ref)
    plan = deepcopy(old)
    plan["source_commit"] = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    paths = {r["path"] for r in old["helper_sources"]}
    paths.update(ADDED_HELPERS)
    plan["helper_sources"] = [record(p) for p in sorted(paths)]
    plan["history_adoption"] = dict(schema="native_factorial_receipt_history_adoption_v1",
        prior_plan=old_ref, prior_review=review_ref, recovery=recovery_ref, executed_source=record(ARCHIVE))
    plan["evidence"].extend([old_ref, review_ref, recovery_ref, record(ARCHIVE),
        record(ROOT / "benchmark_tools/PUBLICATION_GOAL_20261003.txt"),
        record(ROOT / "benchmark_tools/NATIVE_FACTORIAL_RECEIPT_AMENDMENT_20261004.md")])
    adoption = adopt(review_ref, prior, plan)
    for run in plan["runs"][1:]:
        require(not Path(run["output_root"]).exists() and not (Path(plan["panel_root"]) / "sessions" /
            f"run_{run['index']:02d}").exists(), "A remaining identity already has an attempt")
    policy_ref = record(ROOT / "benchmark_tools/results/native_factorial_cost_policy_repaired_20261004.json")
    policy = read(policy_ref)
    require(policy["plan_sha256"] == old_ref["sha256"], "Historical policy binds a different plan")
    output_directory.mkdir(exist_ok=False)
    plan_path = output_directory / "plan.json"
    save(plan_path, plan)
    new_ref = record(plan_path)
    policy["plan_sha256"] = new_ref["sha256"]
    policy["review_reference"] = "Receipt/control amendment only; failed index zero and recovered science retained; same shared-host resources"
    policy["evidence"].extend([policy_ref, new_ref, recovery_ref])
    save(output_directory / "policy.json", policy)
    save(output_directory / "preparation.json", dict(status="receipt_amendment_prepared_not_launched",
        plan=new_ref, policy=record(output_directory / "policy.json"), adoption=adoption,
        source=record(__file__), first_remaining_index=1, automatic_retry=False,
        scientific_execution_authorized=False, timing_scope="shared_host_matched_resources"))
    return record(output_directory / "preparation.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.output_directory.absolute()), sort_keys=True))
