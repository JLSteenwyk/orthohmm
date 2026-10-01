"""Integrate the recovered private CPM arm with unchanged parameter statistics."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import audit_qfo_parameter_swiss as legacy
from benchmark_tools import admit_private_helper_cpm_assessment as gate
from benchmark_tools import run_qfo_parameter_uncertainty as frozen
from benchmark_tools.bootstrap_qfo_parameter_neighborhood import ARMS, calculate, REPLICATES, SEED, MULTIPLICITY
from benchmark_tools.cpm_replay_context import PLAN_SHA, PROTOCOL_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

PROTOCOL = "benchmark_tools/results/QFO_PRIVATE_CPM_PARAMETER_INTEGRATION_PROTOCOL_20261001.md"
PRIOR = "benchmark_tools/results/qfo_parameter_uncertainty_partial_inventory_20260923.json"
PRIOR_SHA = "7d866297598f90bac2f58703746abcbb3726c5fe9e09cb41e11d1c1e885d1d8c"
GATE_EXECUTOR = "benchmarks/work/qfo_private_cpm_score_admission_executor_20261001"
GATE_COMMIT = "2e99680c590810e2bd853f8428a7301cbba1e2c1"
GATE_SHA = "a0b9a01adf7bf802cda37b2767d60ab1a41136b6786d96931950617ed5196de8"
GATE_PROTOCOL_SHA = "d44805159699728ff0cef6a4e4414912c36b2ffb194382ae6b91d951c6dfd071"
ADMISSION = "benchmarks/work/qfo_private_cpm_score_admission_20261001.json"
EXTRA_SOURCES = {
    "run_qfo_parameter_uncertainty.py": "2c551e6a961ce130f7973c064096056beef159ca2b1212bd0a414169b962c719",
    "admit_private_helper_cpm_assessment.py": GATE_SHA,
    "run_private_helper_cpm_assessment.py": gate.RUNNER_SHA,
    "prepare_private_helper_cpm_pairs.py": gate.runner.CONVERTER_SHA,
    "reproduce_qfo_parameter_uncertainty.py": "6b337c936bf1cf6896fea06de6e1924092c1dc6e5b3871ad621edbdc2138ddf8",
}


def validate_inventory(inventory, prior):
    if (inventory["status"] != "qfo_parameter_score_admission_inventory"
            or prior["status"] != inventory["status"]
            or [row["arm"] for row in inventory["arms"]] != list(ARMS)
            or [row["arm"] for row in prior["arms"]] != list(ARMS)):
        raise ValueError("Require unchanged ordered seven-arm inventory")
    for row, old in zip(inventory["arms"], prior["arms"]):
        if row["arm"] != "cpm_high" and row != old:
            raise ValueError("Previously admitted parameter arm changed")
        if row["arm"] != "cpm_high" and row["status"] != "admitted":
            raise ValueError("Previously admitted parameter arm omitted")
    high = inventory["arms"][2]
    if high["status"] == "not_admitted":
        if (set(high) != {"arm", "status", "reason"} or not isinstance(high["reason"], str)
                or not high["reason"].strip()):
            raise ValueError("Unavailable high-CPM arm requires reason without imputed evidence")
    elif high["status"] != "admitted" or set(high) != {"arm", "status", "admission", "admission_job", "admission_submission"}:
        raise ValueError("Recovered high-CPM arm requires terminal admission and actual submission")


def private_submission(root, row, report):
    executor = root / GATE_EXECUTOR
    job = row["admission_job"]
    receipt_path = root / f"benchmark_tools/results/qfo_private_cpm_score_admission_submission_{job}.json"
    if row["admission_submission"] != record(receipt_path):
        raise ValueError("Private score-admission submission path or bytes differ")
    receipt = read_frozen(receipt_path, row["admission_submission"]["sha256"])
    script = "benchmark_tools/results/qfo_private_cpm_score_admission_20261001.sh"
    names = ("benchmark_tools/admit_private_helper_cpm_assessment.py", gate.PROTOCOL,
             script, "tests/unit/test_admit_private_helper_cpm_assessment.py")
    sources = [record(executor / name) for name in names]
    if (receipt["status"] != "private_recovered_qfo_high_cpm_score_admission_submitted"
            or receipt["job_id"] != job or receipt["executor"] != str(executor)
            or receipt["executor_commit"] != GATE_COMMIT or receipt["executor_clean"] is not True
            or receipt["source_records"] != sources
            or any(receipt[key] is not False for key in
                   ("native_inference", "accuracy_evaluated", "controlled_timing", "publication_ready"))
            or receipt["assessment_report"] != report["execution_report"]
            or receipt["parent_submission"] != report["submission"]
            or receipt["submission_argv"] != ["sbatch", "--parsable", str(executor / script), str(executor),
                GATE_COMMIT, GATE_PROTOCOL_SHA, report["scheduler"]["JobIDRaw"],
                report["execution_report"]["sha256"], report["submission"]["sha256"]]):
        raise ValueError("Actual private score-admission submission differs")
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != GATE_COMMIT:
        raise ValueError("Private score-admission executor revision changed")
    if subprocess.check_output(["git", "-C", str(executor), "status", "--porcelain", "--untracked-files=all"], text=True).strip():
        raise ValueError("Private score-admission executor is not clean")
    for ref, name in zip(sources, names):
        blob = subprocess.check_output(["git", "-C", str(executor), "show", GATE_COMMIT + ":" + name])
        if len(blob) != ref["bytes"] or hashlib.sha256(blob).hexdigest() != ref["sha256"]:
            raise ValueError("Private score-admission submitted Git source differs")
        check(ref)
    if (report["source"] != sources[0] or sources[0]["sha256"] != GATE_SHA
            or report["protocol"] != record(root / gate.PROTOCOL)
            or report["protocol"]["sha256"] != GATE_PROTOCOL_SHA
            or sources[1]["sha256"] != GATE_PROTOCOL_SHA):
        raise ValueError("Wrong frozen private score-admission source or protocol")
    return [row["admission_submission"], *sources]


def validate_private(report, conversion, candidate):
    arm = candidate["candidate_arm"]
    context = dict(candidate_admission=conversion["candidate_admission"],
        private_native_admission=conversion["private_native_admission"],
        recovered_arm=dict(label="cpm_high", partition=arm["candidate_partition"],
            constraints=arm["membership_constraints"], seed_partition=arm["seed_partition"]))
    if (report["status"] != "private_recovered_cpm_assessment_admitted"
            or report["arm"] != "cpm_high" or type(report["index"]) is not int or report["index"] != 1
            or report["source"]["sha256"] != GATE_SHA
            or report["accuracy_admitted"] is not True or report["publication_ready"] is not False
            or report["controlled_timing"] is not False or report["conversion"] != conversion
            or candidate["status"] != "cpm_helper_recovered_candidates_admitted_unscored"
            or candidate["candidate_admitted"] is not True or candidate["arm"] != "cpm_high"
            or type(candidate["index"]) is not int or candidate["index"] != 1
            or report["context"] != context
            or report["assessment"]["participant"] != conversion["participant"]):
        raise ValueError("Wrong recovered private high-CPM score admission")
    gate.runner.validate_stage(conversion, report["conversion_scheduler"], conversion["source"])
    if conversion["source"]["sha256"] != gate.runner.CONVERTER_SHA:
        raise ValueError("Wrong recovered private conversion source")
    for key in ("execution_report", "preflight", "submission", "pairs_manifest", "environment_manifest", "source", "protocol"):
        if report[key] not in report["checked_records"]:
            raise ValueError("Missing independently checked private score evidence")


def audit(root, inventory, baseline, checked):
    def verify(items):
        for ref in items:
            old = checked.get(ref["path"])
            if old is not None and old != ref:
                raise ValueError("Conflicting parameter evidence file identities")
            if old is None:
                check(ref)
                checked[ref["path"]] = ref

    verify([baseline["reference"], baseline["stages"][0]["raw_file"]])
    entries, private_state = [], None
    for row in inventory["arms"]:
        arm = row["arm"]
        if row["status"] == "not_admitted":
            entries.append(dict(row))
            continue
        if arm == "cpm_high":
            job = row["admission_job"]
            gate.assessment_job(job)
            text = gate.runner.converter.accounting(job)
            scheduler = gate.runner.converter.completed_admission(text, job)
            if row["admission"]["path"] != str(root / ADMISSION):
                raise ValueError("Wrong recovered private score-admission output path")
            private_state = (job, scheduler)
        verify([row["admission"]])
        report = read_frozen(Path(row["admission"]["path"]), row["admission"]["sha256"])
        verify(legacy.file_records(report))
        pair, execution_pin = report["pairs_manifest"], report["execution_report"]
        if execution_pin not in report["checked_records"]:
            raise ValueError("Parameter execution was not independently checked")
        conversion = read_frozen(Path(pair["path"]), pair["sha256"])
        if arm == "cpm_high":
            verify(private_submission(root, row, report))
            candidate_pin = conversion["candidate_admission"]
            verify([candidate_pin])
            candidate = read_frozen(Path(candidate_pin["path"]), candidate_pin["sha256"])
            validate_private(report, conversion, candidate)
        else:
            legacy.validate_admission(arm, report, conversion)
        execution = read_frozen(Path(execution_pin["path"]), execution_pin["sha256"])
        raw = legacy.raw_from_execution(arm, report, execution)
        verify([raw])
        entries.append(dict(arm=arm, status="admitted", assessment=report["assessment"], raw_file=raw))
    counts = legacy.assemble(entries, baseline)
    if private_state is not None:
        job, scheduler = private_state
        if gate.runner.converter.completed_admission(gate.runner.converter.accounting(job), job) != scheduler:
            raise ValueError("Private score-admission completion changed during count audit")
    counts.update(source=record(__file__), checked_inputs=list(checked.values()))
    return counts


def run(root, inventory_path, inventory_sha, baseline_path, protocol_sha, output):
    if not output.is_absolute() or output.resolve() != output:
        raise ValueError("Require direct absolute private parameter-analysis output")
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    source, protocol = record(__file__), record(root / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Recovered parameter integration protocol changed")
    helpers = [record(Path(__file__).with_name(name)) for name in {**frozen.SOURCES, **EXTRA_SOURCES}]
    for ref, digest in zip(helpers, {**frozen.SOURCES, **EXTRA_SOURCES}.values()):
        if ref["sha256"] != digest:
            raise ValueError("Frozen parameter count/statistics/admission source changed")
    plan, scientific_protocol = root / "benchmark_tools/results/qfo_parameter_neighborhood_plan_20260919.json", root / "benchmark_tools/results/QFO_PARAMETER_NEIGHBORHOOD_PROTOCOL_20260919.md"
    read_frozen(plan, PLAN_SHA)
    if record(scientific_protocol)["sha256"] != PROTOCOL_SHA:
        raise ValueError("Frozen scientific parameter protocol changed")
    prior = read_frozen(root / PRIOR, PRIOR_SHA)
    inventory = read_frozen(inventory_path, inventory_sha)
    validate_inventory(inventory, prior)
    baseline = read_frozen(baseline_path, legacy.BASE_COUNTS_SHA)
    refs = [source, *helpers, protocol, record(plan), record(scientific_protocol), record(root / PRIOR),
            record(inventory_path), record(baseline_path)]
    for ref in refs:
        check(ref)
    checked = {ref["path"]: ref for ref in refs}
    result = dict(status="validating", source=source, helpers=helpers, protocol=record(scientific_protocol),
        execution_protocol=protocol, plan=record(plan), admission_inventory=record(inventory_path),
        baseline_audit=record(baseline_path), scientific_inputs_admitted=False, uncertainty_admitted=False,
        complete_panel=False, controlled_timing=False, publication_ready=False)
    try:
        counts = audit(root, inventory, baseline, checked)
        calculated = calculate(counts)
        if (calculated["protocol_controls_match"] is not True or calculated["replicates"] != REPLICATES
                or calculated["seed"] != SEED or calculated["multiplicity_endpoints"] != MULTIPLICITY):
            raise ValueError("Parameter bootstrap controls differ from frozen protocol")
        for ref in checked.values():
            check(ref)
        estimated = sum(row["status"] == "estimated" for row in calculated["comparisons"])
        result.update(calculated, status="corrected_qfo_parameter_uncertainty_audited", reconstructed_counts=counts,
            checked_inputs=list(checked.values()), scientific_inputs_admitted=True, uncertainty_admitted=estimated > 0,
            estimated_contrasts=estimated, complete_panel=estimated == 6)
        result["limitations"] = [text for text in result["limitations"] if not text.startswith("Numerical calculation only;")]
        result["limitations"].extend([
            "Recovered high-CPM score admission has separate private-runtime provenance; other six admitted arms are unchanged.",
            "Raw SwissTrees families are reconstructed against the same reference anchor; no inference/scoring is rerun.",
            "Shared-host recovery is not controlled timing; QfO is development-exposed.",
            "Independent numerical reproduction, figure/table export and other publication requirements remain separate."])
    except BaseException as error:
        result.update(status="validation_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "inventory", "baseline", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    for name in ("inventory-sha256", "protocol-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.inventory.resolve(), args.inventory_sha256, args.baseline.resolve(),
        args.protocol_sha256, args.output)
