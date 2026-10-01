import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import run_private_cpm_parameter_uncertainty as module
from benchmark_tools import bootstrap_qfo_parameter_neighborhood as kernel
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.reproduce_qfo_parameter_uncertainty import verify as reproduce
from tests.unit.test_audit_qfo_parameter_swiss import fixture as family_fixture
from tests.unit.test_run_private_helper_cpm_assessment import fixture as conversion_fixture


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True))
    return record(path)


def prior_inventory():
    return json.loads((Path(module.__file__).resolve().parents[1] / module.PRIOR).read_text())


def test_frozen_sources_and_prior_six_rows():
    sources = {**module.frozen.SOURCES, **module.EXTRA_SOURCES}
    assert all(record(Path(module.__file__).with_name(name))["sha256"] == digest for name, digest in sources.items())
    root = Path(module.__file__).resolve().parents[1]
    assert record(root / module.PRIOR)["sha256"] == module.PRIOR_SHA
    module.validate_inventory(prior_inventory(), prior_inventory())


@pytest.mark.parametrize("problem", ["order", "missing", "changed_prior", "omitted_prior", "status", "reason", "imputed", "missing_job"])
def test_exact_inventory_and_missing_arm_contract(problem):
    prior = prior_inventory()
    inventory = copy.deepcopy(prior)
    if problem == "order":
        inventory["arms"].reverse()
    elif problem == "missing":
        inventory["arms"].pop()
    elif problem == "changed_prior":
        inventory["arms"][0]["admission"]["sha256"] = "wrong"
    elif problem == "omitted_prior":
        inventory["arms"][1] = dict(arm="cpm_low", status="not_admitted", reason="omit")
    elif problem == "status":
        inventory["status"] = "unverified"
    elif problem == "reason":
        inventory["arms"][2]["reason"] = " "
    elif problem == "imputed":
        inventory["arms"][2]["metrics"] = {}
    else:
        inventory["arms"][2] = dict(arm="cpm_high", status="admitted", admission={})
    with pytest.raises(ValueError):
        module.validate_inventory(inventory, prior)


def private_report(stage, candidate):
    arm = candidate["candidate_arm"]
    report = dict(status="private_recovered_cpm_assessment_admitted", arm="cpm_high", index=1,
        source={"sha256": module.GATE_SHA}, conversion=stage, conversion_scheduler={"JobIDRaw": stage["job_id"]},
        context=dict(candidate_admission=stage["candidate_admission"], private_native_admission=stage["private_native_admission"],
            recovered_arm=dict(label="cpm_high", partition=arm["candidate_partition"],
                constraints=arm["membership_constraints"], seed_partition=arm["seed_partition"])),
        accuracy_admitted=True, publication_ready=False,
        controlled_timing=False, assessment={"participant": stage["participant"]}, checked_records=[])
    for key in ("execution_report", "preflight", "submission", "pairs_manifest", "environment_manifest", "protocol"):
        report[key] = {"path": "/fixture/" + key, "bytes": 0, "sha256": "fixture"}
    report["checked_records"] = [report[key] for key in
        ("execution_report", "preflight", "submission", "pairs_manifest", "environment_manifest", "source", "protocol")]
    return report


@pytest.mark.parametrize("problem", [None, "status", "arm", "index", "source", "accuracy", "ready", "timing",
    "conversion", "context", "participant", "checked", "counts", "semantics", "converter_source", "candidate", "partition"])
def test_private_arm_requires_original_converted_stage(tmp_path, monkeypatch, problem):
    data = conversion_fixture(tmp_path, monkeypatch)
    stage = data.stage
    candidate = dict(status="cpm_helper_recovered_candidates_admitted_unscored", candidate_admitted=True,
        arm="cpm_high", index=1, candidate_arm=dict(candidate_partition=stage["candidate_admission"],
            membership_constraints=stage["private_native_admission"], seed_partition=stage["native_input"]))
    report = private_report(stage, candidate)
    if problem in ("status", "arm", "index", "source", "accuracy", "ready", "timing"):
        key, value = {"status": ("status", "cpm_assessment_admitted"), "arm": ("arm", "cpm_low"),
            "index": ("index", True), "source": ("source", {"sha256": "wrong"}),
            "accuracy": ("accuracy_admitted", False), "ready": ("publication_ready", True),
            "timing": ("controlled_timing", 0)}[problem]
        report[key] = value
    elif problem == "conversion":
        report["conversion"] = {}
    elif problem == "context":
        report["context"]["candidate_admission"] = {}
    elif problem == "participant":
        report["assessment"]["participant"] = "wrong"
    elif problem == "checked":
        report["checked_records"].pop()
    elif problem == "counts":
        stage["retained_pairs"] -= 1
    elif problem == "semantics":
        stage["semantics"] = "RootHOG cliques"
    elif problem == "converter_source":
        stage["source"]["sha256"] = "wrong"
    elif problem == "candidate":
        candidate["candidate_admitted"] = False
    elif problem == "partition":
        report["context"]["recovered_arm"]["partition"] = {}
    if problem is None:
        module.validate_private(report, stage, candidate)
    else:
        with pytest.raises((ValueError, KeyError)):
            module.validate_private(report, stage, candidate)


def submission_fixture(tmp_path, monkeypatch):
    root = Path(module.__file__).resolve().parents[1]
    executor = tmp_path / module.GATE_EXECUTOR
    names = ("benchmark_tools/admit_private_helper_cpm_assessment.py", module.gate.PROTOCOL,
        "benchmark_tools/results/qfo_private_cpm_score_admission_20261001.sh",
        "tests/unit/test_admit_private_helper_cpm_assessment.py")
    blobs = {}
    for name in names:
        path = executor / name
        path.parent.mkdir(parents=True, exist_ok=True)
        blobs[name] = (root / name).read_bytes()
        path.write_bytes(blobs[name])
    protocol = tmp_path / module.gate.PROTOCOL
    protocol.parent.mkdir(parents=True, exist_ok=True)
    protocol.write_bytes((root / module.gate.PROTOCOL).read_bytes())
    report = dict(source=record(executor / names[0]), protocol=record(protocol),
        scheduler={"JobIDRaw": "27000"}, execution_report=dump(tmp_path / "execution.json", {}),
        submission=dump(tmp_path / "assessment_receipt.json", {}))
    receipt = dict(status="private_recovered_qfo_high_cpm_score_admission_submitted", job_id="27001",
        executor=str(executor), executor_commit=module.GATE_COMMIT, executor_clean=True,
        source_records=[record(executor / name) for name in names], native_inference=False,
        accuracy_evaluated=False, controlled_timing=False, publication_ready=False,
        assessment_report=report["execution_report"], parent_submission=report["submission"],
        submission_argv=["sbatch", "--parsable", str(executor / names[2]), str(executor), module.GATE_COMMIT,
            module.GATE_PROTOCOL_SHA, "27000", report["execution_report"]["sha256"], report["submission"]["sha256"]])
    path = tmp_path / "benchmark_tools/results/qfo_private_cpm_score_admission_submission_27001.json"
    row = dict(admission_job="27001", admission_submission=dump(path, receipt))

    def git(argv, **kwargs):
        if "rev-parse" in argv:
            return module.GATE_COMMIT + "\n"
        if "status" in argv:
            return ""
        if "show" in argv:
            return blobs[argv[-1].split(":", 1)[1]]
        pytest.fail("Unexpected command " + repr(argv))

    monkeypatch.setattr(module.subprocess, "check_output", git)
    return SimpleNamespace(report=report, receipt=receipt, path=path, row=row, executor=executor)


@pytest.mark.parametrize("problem", [None, "status", "job_id", "executor", "executor_commit", "executor_clean",
    "source_records", "assessment_report", "parent_submission", "submission_argv", "native_inference",
    "accuracy_evaluated", "controlled_timing", "publication_ready", "git_blob", "dirty"])
def test_actual_private_validator_submission_binding(tmp_path, monkeypatch, problem):
    data = submission_fixture(tmp_path, monkeypatch)
    if problem in data.receipt:
        data.receipt[problem] = True if problem in ("native_inference", "accuracy_evaluated", "controlled_timing", "publication_ready") else "wrong"
        data.row["admission_submission"] = dump(data.path, data.receipt)
    elif problem == "git_blob":
        path = Path(data.receipt["source_records"][0]["path"])
        path.write_text("# changed source\n")
        data.receipt["source_records"][0] = record(path)
        data.row["admission_submission"] = dump(data.path, data.receipt)
    elif problem == "dirty":
        original = module.subprocess.check_output
        monkeypatch.setattr(module.subprocess, "check_output", lambda argv, **kwargs:
            "?? file.py\n" if "status" in argv else original(argv, **kwargs))
    if problem is None:
        assert len(module.private_submission(tmp_path, data.row, data.report)) == 5
    else:
        with pytest.raises(ValueError):
            module.private_submission(tmp_path, data.row, data.report)


def setup(tmp_path, monkeypatch, missing_high=False):
    entries, baseline = family_fixture(tmp_path)
    reference = tmp_path / "reference.drw"
    reference.write_text("# frozen reference fixture\n")
    baseline["reference"] = record(reference)
    monkeypatch.setattr(module.legacy, "REFERENCE_SHA", baseline["reference"]["sha256"])
    monkeypatch.setattr(kernel, "REFERENCE_SHA", baseline["reference"]["sha256"])
    baseline_path = tmp_path / "baseline.json"
    monkeypatch.setattr(module.legacy, "BASE_COUNTS_SHA", dump(baseline_path, baseline)["sha256"])
    inventory = dict(status="qfo_parameter_score_admission_inventory", arms=[])
    for index, entry in enumerate(entries):
        arm = entry["arm"]
        native_index = 7 if arm == "control" else 0 if arm == "cpm_low" else 1 if arm == "cpm_high" else index - 3
        raw = tmp_path / arm / "SwissTrees/participant.raw.txt.gz"
        raw.parent.mkdir(parents=True)
        raw.write_bytes(Path(entry["raw_file"]["path"]).read_bytes())
        conversion = {"participant": entry["assessment"]["participant"], "fixture": "converted " + arm}
        if arm == "cpm_high":
            conversion["candidate_admission"] = dump(tmp_path / "candidate.json", {})
        pairs = dump(tmp_path / f"pairs_{arm}.json", conversion)
        execution = dict(status="process_succeeded_pending_independent_admission", exit_code=0,
            job_id=str(28000 + index), index=native_index, pairs_manifest=pairs, stage=conversion,
            outputs=[record(raw)])
        execution["cell" if arm == "control" else "arm" if arm.startswith("cpm") else "variant"] = "p1_c1_r1" if arm == "control" else arm
        context = {"arm": arm}
        if arm.startswith("cpm"):
            execution["context"] = context
        execution_pin = dump(tmp_path / f"execution_{arm}.json", execution)
        report = dict(index=native_index, pairs_manifest=pairs, conversion=conversion, context=context,
            scheduler=dict(State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="8", JobIDRaw=execution["job_id"]),
            assessment=entry["assessment"], execution_report=execution_pin, checked_records=[execution_pin, pairs, record(raw)])
        path = tmp_path / module.ADMISSION if arm == "cpm_high" else tmp_path / f"admission_{arm}.json"
        row = dict(arm=arm, status="admitted", admission=dump(path, report))
        if arm == "cpm_high":
            row.update(admission_job="28010", admission_submission=dump(tmp_path / "gate_submission.json", {}))
        inventory["arms"].append(row)
    prior = copy.deepcopy(inventory)
    prior["arms"][2] = dict(arm="cpm_high", status="not_admitted", reason="historical failure retained")
    monkeypatch.setattr(module, "PRIOR_SHA", dump(tmp_path / module.PRIOR, prior)["sha256"])
    if missing_high:
        inventory["arms"][2] = copy.deepcopy(prior["arms"][2])
    root = Path(module.__file__).resolve().parents[1]
    for name in (module.PROTOCOL, "benchmark_tools/results/qfo_parameter_neighborhood_plan_20260919.json",
                 "benchmark_tools/results/QFO_PARAMETER_NEIGHBORHOOD_PROTOCOL_20260919.md"):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes((root / name).read_bytes())
    inventory_path = tmp_path / "inventory.json"
    digest = dump(inventory_path, inventory)["sha256"]
    seen = []
    monkeypatch.setattr(module.legacy, "validate_admission", lambda arm, *args: seen.append(arm))
    monkeypatch.setattr(module, "validate_private", lambda *args: seen.append("private_cpm_high"))
    monkeypatch.setattr(module, "private_submission", lambda root, row, report: [row["admission_submission"]])
    monkeypatch.setattr(module.gate.runner.converter, "accounting", lambda job:
        "JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n"
        f"{job}|{job}|COMPLETED|0:0|00:01:00|bizon|2|64G\n"
        f"{job}.batch|{job}.batch|COMPLETED|0:0|00:01:00|bizon|2|\n")
    args = (tmp_path, inventory_path, digest, baseline_path, record(tmp_path / module.PROTOCOL)["sha256"], tmp_path / "result.json")
    return SimpleNamespace(args=args, seen=seen, inventory=inventory, prior=prior, baseline=baseline)


@pytest.mark.parametrize("missing_high", [False, True])
def test_actual_raw_assembly_kernel_and_independent_reproduction(tmp_path, monkeypatch, missing_high):
    data = setup(tmp_path, monkeypatch, missing_high)
    result = module.run(*data.args)
    assert result["status"] == "corrected_qfo_parameter_uncertainty_audited"
    assert result["estimated_contrasts"] == (5 if missing_high else 6)
    assert result["complete_panel"] is (not missing_high)
    assert result["uncertainty_admitted"] is True and result["publication_ready"] is False
    assert result["replicates"] == 100000 and result["seed"] == 20260925 and result["multiplicity_endpoints"] == 18
    assert result["reconstructed_counts"]["reference_relation_count"] == 270
    assert reproduce(result) == (15 if missing_high else 18)
    assert data.seen == ["control", "cpm_low", *([] if missing_high else ["private_cpm_high"]),
        "norm_low", "norm_high", "margin_low", "margin_high"]
    assert json.loads(data.args[-1].read_text()) == result


@pytest.mark.parametrize("problem", ["input_pin", "plan", "protocol", "sources", "inventory", "live", "raw", "controls", "changed_at_close"])
def test_fail_closed_without_uncertainty_promotion(tmp_path, monkeypatch, problem):
    data = setup(tmp_path, monkeypatch)
    args = list(data.args)
    if problem == "input_pin":
        args[2] = "0" * 64
    elif problem == "plan":
        dump(tmp_path / "benchmark_tools/results/qfo_parameter_neighborhood_plan_20260919.json", {})
    elif problem == "protocol":
        args[4] = "0" * 64
    elif problem == "sources":
        monkeypatch.setattr(module, "EXTRA_SOURCES", {"admit_private_helper_cpm_assessment.py": "wrong"})
    elif problem == "inventory":
        data.inventory["arms"][0]["admission"]["sha256"] = "wrong"
        args[2] = dump(args[1], data.inventory)["sha256"]
    elif problem == "live":
        monkeypatch.setattr(module.gate.runner.converter, "accounting", lambda job:
            "JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n"
            f"{job}|{job}|RUNNING|0:0|00:01:00|bizon|2|64G\n")
    elif problem == "raw":
        Path(data.baseline["stages"][0]["raw_file"]["path"]).write_text("changed raw\n")
    else:
        original = module.calculate
        def changed(counts):
            value = original(counts)
            if problem == "controls":
                value["multiplicity_endpoints"] = 15
            else:
                args[1].write_text("changed during calculation")
            return value
        monkeypatch.setattr(module, "calculate", changed)
    with pytest.raises((ValueError, json.JSONDecodeError)):
        module.run(*args)
    if args[-1].exists():
        failed = json.loads(args[-1].read_text())
        assert failed["status"] == "validation_failed" and failed["uncertainty_admitted"] is False
        assert failed["scientific_inputs_admitted"] is False and failed["complete_panel"] is False


@pytest.mark.parametrize("symlink", [False, True])
def test_no_output_overwrite(tmp_path, monkeypatch, symlink):
    output = tmp_path / "result.json"
    if symlink:
        output.symlink_to(tmp_path / "missing")
    else:
        output.write_text("preserved")
    monkeypatch.setattr(module, "record", lambda *args: pytest.fail("Read inputs before overwrite gate"))
    with pytest.raises((ValueError, FileExistsError)):
        module.run(tmp_path, None, None, None, None, output)


@pytest.mark.parametrize("job", [True, 28010, "", "028010", "28010_1", "\u0661"])
def test_invalid_private_job_never_queries_scheduler(tmp_path, monkeypatch, job):
    data = setup(tmp_path, monkeypatch)
    data.inventory["arms"][2]["admission_job"] = job
    args = list(data.args)
    args[2] = dump(args[1], data.inventory)["sha256"]
    monkeypatch.setattr(module.gate.runner.converter, "accounting", lambda *args: pytest.fail("Queried invalid job"))
    with pytest.raises(ValueError, match="standalone"):
        module.run(*args)
    assert json.loads(args[-1].read_text())["scientific_inputs_admitted"] is False


def test_live_private_admission_is_not_read(tmp_path, monkeypatch):
    data = setup(tmp_path, monkeypatch)
    original = module.read_frozen
    private_path = tmp_path / module.ADMISSION
    def read(path, sha):
        if path == private_path:
            pytest.fail("Read active private score-admission output")
        return original(path, sha)
    monkeypatch.setattr(module, "read_frozen", read)
    monkeypatch.setattr(module.gate.runner.converter, "accounting", lambda job:
        "JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n"
        f"{job}|{job}|RUNNING|0:0|00:01:00|bizon|2|64G\n")
    with pytest.raises(ValueError, match="completed"):
        module.run(*data.args)


def test_private_completion_stable_at_close(tmp_path, monkeypatch):
    data = setup(tmp_path, monkeypatch)
    original = module.legacy.assemble
    def changed(entries, baseline):
        value = original(entries, baseline)
        monkeypatch.setattr(module.gate.runner.converter, "accounting", lambda job:
            "JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n"
            f"{job}|{job}|FAILED|1:0|00:01:00|bizon|2|64G\n")
        return value
    monkeypatch.setattr(module.legacy, "assemble", changed)
    monkeypatch.setattr(module, "calculate", lambda *args: pytest.fail("Calculated changed admission"))
    with pytest.raises(ValueError):
        module.run(*data.args)
    assert json.loads(data.args[-1].read_text())["uncertainty_admitted"] is False
