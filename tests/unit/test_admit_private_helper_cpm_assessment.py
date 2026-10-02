import copy
from contextlib import nullcontext
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile

import pytest

from benchmark_tools import admit_private_helper_cpm_assessment as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.validate_qfo_native_assessment import AXES
from tests.unit.test_admit_qfo_recovered_assessment import trace
from tests.unit.test_validate_qfo_native_assessment import records as native_records

JOB = "25001"
CONVERSION_JOB = "25000"


@pytest.fixture(name="tmp_path")
def short_tmp_path():
    # Exercise the unchanged Darwin command builder within its real path limit.
    with tempfile.TemporaryDirectory(prefix="oh-qg-", dir="/tmp") as directory:
        yield Path(directory).resolve()


def test_short_fixture_resolves_symlink_root(tmp_path, monkeypatch):
    alias = tmp_path / "temporary_alias"
    alias.symlink_to(tmp_path, target_is_directory=True)
    monkeypatch.setattr(tempfile, "TemporaryDirectory", lambda **kwargs: nullcontext(str(alias)))
    fixture = short_tmp_path.__wrapped__()
    assert next(fixture) == tmp_path
    with pytest.raises(StopIteration):
        next(fixture)


def accounting(state="COMPLETED", code="0:0", cpus="8", memory="64G", node="bizon"):
    return ("JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n"
        f"{JOB}|{JOB}|{state}|{code}|00:10:00|{node}|{cpus}|{memory}\n"
        f"{JOB}.batch|{JOB}.batch|{state}|{code}|00:10:00|{node}|{cpus}|\n")


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True))


def fixture(tmp_path, monkeypatch):
    root = Path(module.__file__).resolve().parents[1]
    executor = tmp_path / module.EXECUTOR
    names = ("benchmark_tools/run_private_helper_cpm_assessment.py", module.runner.PROTOCOL,
        "benchmark_tools/results/qfo_private_cpm_assessment_20261001.sh",
        "tests/unit/test_run_private_helper_cpm_assessment.py")
    blobs = {}
    for name in [*names, *("benchmark_tools/" + name for name in sorted(module.REQUIRED_HELPERS))]:
        path = executor / name
        path.parent.mkdir(parents=True, exist_ok=True)
        blobs[name] = (root / name).read_bytes()
        path.write_bytes(blobs[name])
    scoring_protocol = tmp_path / module.runner.PROTOCOL
    scoring_protocol.parent.mkdir(parents=True, exist_ok=True)
    scoring_protocol.write_bytes((root / module.runner.PROTOCOL).read_bytes())
    protocol = tmp_path / module.PROTOCOL
    protocol.write_text("# Prospective score gate fixture\n")
    mapping = tmp_path / "mapping.json.gz"
    mapping.write_text("fixture mapping\n")
    pairs = tmp_path / "native_pairs.qfo.tsv"
    pairs.write_text("A\tB\n")
    stage = dict(mapping=record(mapping), filtered_pairs=record(pairs),
        participant="ohmm_qfo_parameter_cpm_high", arm="cpm_high", index=1,
        semantics="native phylogenetically inferred pairs")
    conversion = tmp_path / "conversion.json"
    conversion_submission = tmp_path / "conversion_submission.json"
    dump(conversion, stage)
    dump(conversion_submission, {"fixture": "conversion submitted"})
    verified = dict(stage=stage, pairs_manifest=record(conversion), submission=record(conversion_submission),
        context=dict(candidate_admission=record(conversion), private_native_admission=record(conversion_submission),
            recovered_arm={"index": 1}), scheduler=dict(JobIDRaw=CONVERSION_JOB),
        accounting="completed conversion fixture", checked_records=[record(pairs), record(mapping)])
    monkeypatch.setattr(module.runner, "verify_conversion", lambda *args: copy.deepcopy(verified))
    monkeypatch.setattr(module.runner.converter, "accounting", lambda job: accounting())
    pipeline = tmp_path / "qfo_benchmark/benchmark-webservice"
    reference = pipeline / "reference_data"
    family = reference / "2020/ReconciledTrees_SwissTrees.drw"
    family.parent.mkdir(parents=True, exist_ok=True)
    family.write_text("ReconciledTrees['X'] := RecTreeCase('X',...)\n")
    results = tmp_path / module.runner.RESULTS
    rows = native_records()
    for row in rows:
        row["participant_id"] = stage["participant"]
    dump(results / "assessment_out/Assessment_datasets.json", rows)
    for challenge, axes in AXES.items():
        data = {"type": "aggregation", "challenge_ids": [challenge], "datalink": {"inline_data": {
            "visualization": {"x_axis": axes[0], "y_axis": axes[1], "type": "2D-plot"},
            "challenge_participants": [{"participant_id": stage["participant"],
                "metric_x": 4 if axes[0] == "NR_ORTHOLOGS" else .5, "metric_y": .5,
                "stderr_x": .01, "stderr_y": .01}]}}}
        dump(results / "results" / challenge / (challenge + ".json"), data)
        dump(reference / "data" / (challenge + ".json"), data)
    task_trace = results / "stats/trace_fixture.txt"
    task_trace.parent.mkdir(parents=True)
    task_trace.write_text("task_id\tname\tstatus\texit\n" + "\n".join("\t".join(row) for row in trace()))
    infrastructure = tmp_path / "runtime.cfg"
    infrastructure.write_text("# frozen scoring config fixture\n")
    environment = dict(status="local_qfo_assessment_environment_frozen", accuracy_evaluated=False,
        source=record(infrastructure), execution_config=record(infrastructure), singularity_config=record(infrastructure),
        pipeline=str(pipeline), pipeline_files=[], java_files=[], images=[], executables=[], singularity_support=[],
        reference_files=[record(mapping), *[record(path) for path in sorted(reference.rglob("*")) if path.is_file()]],
        environment_overrides={"JAVA_HOME": "/frozen/java", "NXF_OFFLINE": "true"})
    env_path = tmp_path / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    dump(env_path, environment)
    monkeypatch.setattr(module.runner.converter, "ENV_SHA", record(env_path)["sha256"])
    submitted = dict(status="private_recovered_qfo_high_cpm_assessment_submitted", job_id=JOB,
        executor=str(executor), executor_commit=module.COMMIT, executor_clean=True,
        native_inference=False, accuracy_evaluated=False, controlled_timing=False, publication_ready=False,
        pair_conversion=verified["pairs_manifest"], parent_submission=verified["submission"],
        source_records=[record(executor / name) for name in names],
        submission_argv=["sbatch", "--parsable", str(executor / names[2]), str(executor), module.COMMIT,
            module.RUNNER_PROTOCOL_SHA, CONVERSION_JOB, verified["pairs_manifest"]["sha256"], verified["submission"]["sha256"]])
    submission = tmp_path / f"benchmark_tools/results/qfo_private_cpm_assessment_submission_{JOB}.json"
    dump(submission, submitted)
    helpers = [record(executor / "benchmark_tools" / name) for name in sorted(module.REQUIRED_HELPERS)]
    base = [record(scoring_protocol), verified["pairs_manifest"], verified["submission"], record(env_path),
        *module.environment_records(environment), *verified["checked_records"]]
    directory, work = tmp_path / module.runner.OUTPUT, tmp_path / module.runner.WORK
    directory.mkdir(parents=True)
    log = directory / "scoring.log"
    log.write_text("fixture scoring succeeded\n")
    report = dict(status="process_succeeded_pending_independent_admission", exit_code=0, job_id=JOB,
        arm="cpm_high", index=1, context=verified["context"], stage=stage, source=record(executor / names[0]),
        protocol=record(scoring_protocol), pairs_manifest=verified["pairs_manifest"], conversion_submission=verified["submission"],
        environment_manifest=record(env_path), conversion_scheduler=verified["scheduler"], conversion_accounting=verified["accounting"],
        conversion_job=CONVERSION_JOB, conversion_sha256=verified["pairs_manifest"]["sha256"],
        conversion_submission_sha256=verified["submission"]["sha256"], converter_commit=module.runner.CONVERTER_COMMIT,
        command=module.command_for(tmp_path, stage, environment, work, results), cwd=str(directory), work=str(work),
        results=str(results), verified_records=[*base, *helpers], helpers=helpers,
        environment_overrides=environment["environment_overrides"], accuracy_admitted=False,
        controlled_timing=False, publication_ready=False, log=record(log),
        outputs=[record(path) for path in sorted(results.rglob("*")) if path.is_file()])
    preflight = {key: value for key, value in report.items() if key not in ("exit_code", "log", "outputs")}
    preflight["status"] = "running"
    report_path, preflight_path = directory / "results.json", directory / "preflight.json"
    dump(report_path, report)
    dump(preflight_path, preflight)

    def git(argv, **kwargs):
        if "rev-parse" in argv:
            return module.COMMIT + "\n"
        if "status" in argv:
            return ""
        if "show" in argv:
            return blobs[argv[-1].split(":", 1)[1]]
        pytest.fail("Unexpected fixture command: " + repr(argv))

    monkeypatch.setattr(module.subprocess, "check_output", git)
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=0))
    return SimpleNamespace(verified=verified, protocol=protocol, submission=submission, submitted=submitted,
        report=report, report_path=report_path, preflight=preflight, preflight_path=preflight_path,
        results=results, trace=task_trace, executor=executor, environment=environment, env_path=env_path,
        base=base, output=tmp_path / "admitted.json")


def invoke(root, data):
    return module.admit(root, JOB, record(data.report_path)["sha256"], record(data.submission)["sha256"],
        record(data.protocol)["sha256"], data.output)


def save_execution(data):
    data.preflight = {key: value for key, value in data.report.items() if key not in ("exit_code", "log", "outputs")}
    data.preflight["status"] = "running"
    dump(data.report_path, data.report)
    dump(data.preflight_path, data.preflight)


def test_full_gate_uses_real_native_metrics_and_trace(tmp_path, monkeypatch):
    data = fixture(tmp_path, monkeypatch)
    result = invoke(tmp_path, data)
    assert result["status"] == "private_recovered_cpm_assessment_admitted"
    assert result["accuracy_admitted"] is True
    assert result["controlled_timing"] is False and result["publication_ready"] is False
    assert len(result["native_tasks"]) == 15 and len(result["metric_files"]) == 14
    assert result["assessment"]["secondary_six_metric_mean"] == .5
    assert set(result["assessment"]["endpoints"]) == set(AXES)
    assert result["assessment"]["endpoints"]["GO"]["score_semantics"] == "avg Schlicker"
    assert result["conversion"]["semantics"] == "native phylogenetically inferred pairs"
    assert json.loads(data.output.read_text()) == result


@pytest.mark.parametrize("change", [{"state": "RUNNING"}, {"state": "FAILED"}, {"code": "1:0"},
    {"cpus": "2"}, {"memory": "32G"}, {"node": "dgx"}])
def test_wrong_allocation_before_output_reads(tmp_path, monkeypatch, change):
    monkeypatch.setattr(module.runner.converter, "accounting", lambda job: accounting(**change))
    monkeypatch.setattr(module, "record", lambda *args: pytest.fail("Read active/invalid scoring outputs"))
    with pytest.raises(ValueError, match="completed"):
        module.admit(tmp_path, JOB, "", "", "", tmp_path / "admitted.json")


@pytest.mark.parametrize("problem", ["missing_batch", "duplicate_batch", "bad_step", "foreign", "duplicate_parent"])
def test_all_reported_steps_required(problem):
    rows = accounting().splitlines()
    if problem == "missing_batch":
        rows.pop()
    elif problem == "duplicate_batch":
        rows.append(rows[-1])
    elif problem == "bad_step":
        rows.append(f"{JOB}.0|{JOB}.0|FAILED|1:0|00:00:01|bizon|8|")
    elif problem == "foreign":
        rows.append("123|123|COMPLETED|0:0|00:00:01|bizon|8|64G")
    else:
        rows.append(rows[1])
    with pytest.raises(ValueError):
        module.completed_assessment("\n".join(rows), JOB)


@pytest.mark.parametrize("job", [True, 25001, "025001", "25001_1", "", "-1", "\u0661"])
def test_job_identity(tmp_path, monkeypatch, job):
    monkeypatch.setattr(module.runner.converter, "accounting", lambda *args: pytest.fail("Queried invalid job identity"))
    with pytest.raises(ValueError):
        module.completed_assessment(accounting(), job)
    with pytest.raises(ValueError):
        module.admit(tmp_path, job, "", "", "", tmp_path / "admitted.json")


@pytest.mark.parametrize("kind", ["file", "symlink", "indirect", "relative"])
def test_fresh_direct_destination_before_scheduler(tmp_path, monkeypatch, kind):
    output = tmp_path / "admitted.json"
    if kind == "file":
        output.write_text("preserve")
    elif kind == "symlink":
        output.symlink_to(tmp_path / "missing")
    elif kind == "indirect":
        alias = tmp_path / "alias"
        alias.symlink_to(tmp_path, target_is_directory=True)
        output = alias / "admitted.json"
    else:
        output = Path("relative.json")
    monkeypatch.setattr(module.runner.converter, "accounting", lambda *args: pytest.fail("Read scheduler"))
    with pytest.raises((ValueError, FileExistsError)):
        module.admit(tmp_path, JOB, "", "", "", output)


@pytest.mark.parametrize("key,value", [("status", "failed"), ("job_id", "999"), ("exit_code", 1),
    ("exit_code", False),
    ("accuracy_admitted", True), ("controlled_timing", 0), ("publication_ready", True),
    ("index", True), ("arm", "cpm_low"), ("context", {}), ("stage", {}), ("source", {}),
    ("pairs_manifest", {}), ("conversion_submission", {}), ("command", ["nextflow", "-resume"]),
    ("environment_overrides", {}), ("work", "old"), ("results", "old"), ("cwd", "old"),
    ("converter_commit", "wrong"), ("conversion_scheduler", {}), ("conversion_accounting", "wrong")])
def test_changed_execution_retained_without_admission(tmp_path, monkeypatch, key, value):
    data = fixture(tmp_path, monkeypatch)
    data.report[key] = value
    save_execution(data)
    with pytest.raises((ValueError, KeyError)):
        invoke(tmp_path, data)
    failed = json.loads(data.output.read_text())
    assert failed["status"] == "validation_failed" and failed["accuracy_admitted"] is False
    assert "assessment" not in failed and failed["publication_ready"] is False


@pytest.mark.parametrize("key,value", [("status", "wrong"), ("job_id", "999"), ("executor", "other"),
    ("executor_commit", "wrong"), ("executor_clean", False), ("source_records", []), ("pair_conversion", {}),
    ("parent_submission", {}), ("submission_argv", []), ("native_inference", True),
    ("accuracy_evaluated", 0), ("controlled_timing", True), ("publication_ready", True)])
def test_actual_submission_binding(tmp_path, monkeypatch, key, value):
    data = fixture(tmp_path, monkeypatch)
    data.submitted[key] = value
    dump(data.submission, data.submitted)
    with pytest.raises(ValueError, match="submission"):
        invoke(tmp_path, data)


@pytest.mark.parametrize("problem", ["base", "missing", "duplicate", "foreign", "suffix", "changed_blob", "dirty"])
def test_helper_and_git_source_rejection(tmp_path, monkeypatch, problem):
    data = fixture(tmp_path, monkeypatch)
    if problem == "base":
        data.report["verified_records"][0] = {}
    elif problem == "missing":
        data.report["helpers"].pop()
        data.report["verified_records"] = data.base + data.report["helpers"]
    elif problem == "duplicate":
        data.report["helpers"].append(data.report["helpers"][0])
        data.report["verified_records"] = data.base + data.report["helpers"]
    elif problem == "foreign":
        data.report["helpers"][0] = {"path": str(tmp_path / "foreign.py")}
        data.report["verified_records"] = data.base + data.report["helpers"]
    elif problem == "suffix":
        data.report["verified_records"].append(data.report["helpers"][0])
    elif problem == "changed_blob":
        path = Path(data.report["helpers"][0]["path"])
        path.write_text("# changed Git source\n")
        data.report["helpers"][0] = record(path)
        data.report["verified_records"] = data.base + data.report["helpers"]
    else:
        original = module.subprocess.check_output
        monkeypatch.setattr(module.subprocess, "check_output", lambda argv, **kw:
            "?? foreign.py\n" if "status" in argv else original(argv, **kw))
    save_execution(data)
    with pytest.raises(ValueError):
        invoke(tmp_path, data)


@pytest.mark.parametrize("problem", ["missing", "cached", "failed", "metric", "inventory", "reference", "preflight", "report_pin", "submission_pin"])
def test_trace_arithmetic_inventory_and_pins(tmp_path, monkeypatch, problem):
    data = fixture(tmp_path, monkeypatch)
    if problem in ("missing", "cached", "failed"):
        text = data.trace.read_text()
        if problem == "missing":
            text = "\n".join(text.splitlines()[:-1])
        else:
            text = text.replace("COMPLETED", problem.upper(), 1)
        data.trace.write_text(text)
        data.report["outputs"] = [record(path) for path in sorted(data.results.rglob("*")) if path.is_file()]
        save_execution(data)
    elif problem == "metric":
        path = data.results / "results/GO/GO.json"
        content = json.loads(path.read_text())
        content["datalink"]["inline_data"]["challenge_participants"][0]["metric_y"] = .9
        dump(path, content)
        data.report["outputs"] = [record(path) for path in sorted(data.results.rglob("*")) if path.is_file()]
        save_execution(data)
    elif problem == "inventory":
        (data.results / "unrecorded.txt").write_text("new\n")
    elif problem == "reference":
        data.environment["reference_files"].pop()
        dump(data.env_path, data.environment)
        monkeypatch.setattr(module.runner.converter, "ENV_SHA", record(data.env_path)["sha256"])
    elif problem == "preflight":
        data.preflight["context"] = {}
        dump(data.preflight_path, data.preflight)
    if problem in ("report_pin", "submission_pin"):
        with pytest.raises(ValueError):
            module.admit(tmp_path, JOB, "0" * 64 if problem == "report_pin" else record(data.report_path)["sha256"],
                "0" * 64 if problem == "submission_pin" else record(data.submission)["sha256"],
                record(data.protocol)["sha256"], data.output)
    else:
        with pytest.raises(ValueError):
            invoke(tmp_path, data)
    assert json.loads(data.output.read_text())["accuracy_admitted"] is False


@pytest.mark.parametrize("problem", ["context", "accounting", "source", "output"])
def test_evidence_stability_at_close(tmp_path, monkeypatch, problem):
    data = fixture(tmp_path, monkeypatch)
    original = module.validate_directory

    def change_after_metrics(*args):
        value = original(*args)
        if problem == "context":
            data.verified["context"] = {}
        elif problem == "accounting":
            monkeypatch.setattr(module.runner.converter, "accounting", lambda job: accounting(state="FAILED"))
        elif problem == "source":
            Path(data.report["source"]["path"]).write_text("# changed source\n")
        else:
            (data.results / "new.txt").write_text("changed\n")
        return value

    monkeypatch.setattr(module, "validate_directory", change_after_metrics)
    with pytest.raises(ValueError):
        invoke(tmp_path, data)
    result = json.loads(data.output.read_text())
    assert result["status"] == "validation_failed" and result["accuracy_admitted"] is False


@pytest.mark.parametrize("problem", ["revision", "source_pin", "protocol_pin", "log_path", "helper_reuse", "uninventoried_reference"])
def test_additional_source_and_reference_contracts(tmp_path, monkeypatch, problem):
    data = fixture(tmp_path, monkeypatch)
    if problem == "revision":
        original = module.subprocess.check_output
        monkeypatch.setattr(module.subprocess, "check_output", lambda argv, **kw:
            "wrong\n" if "rev-parse" in argv else original(argv, **kw))
    elif problem == "source_pin":
        monkeypatch.setattr(module, "RUNNER_SHA", "0" * 64)
    elif problem == "protocol_pin":
        with pytest.raises(ValueError, match="protocol"):
            module.admit(tmp_path, JOB, record(data.report_path)["sha256"], record(data.submission)["sha256"],
                "0" * 64, data.output)
        return
    elif problem == "log_path":
        data.report["log"] = record(data.protocol)
        save_execution(data)
    elif problem == "helper_reuse":
        original = module.record
        def changed_local(path):
            value = original(path)
            if str(path) == str(Path(module.__file__).with_name("prepare_private_helper_cpm_pairs.py")):
                value["sha256"] = "0" * 64
            return value
        monkeypatch.setattr(module, "record", changed_local)
    else:
        data.environment["reference_files"].pop()
        dump(data.env_path, data.environment)
        monkeypatch.setattr(module.runner.converter, "ENV_SHA", record(data.env_path)["sha256"])
        data.base = [record(tmp_path / module.runner.PROTOCOL), data.verified["pairs_manifest"], data.verified["submission"],
            record(data.env_path), *module.environment_records(data.environment), *data.verified["checked_records"]]
        data.report["verified_records"] = data.base + data.report["helpers"]
        data.report["environment_manifest"] = record(data.env_path)
        save_execution(data)
    with pytest.raises(ValueError):
        invoke(tmp_path, data)
