import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import run_publication_completion_audits as module


def stage(tmp_path, name="restored_orthobench", value=None, exit_code=0):
    output=tmp_path / "audit-result.json"
    value = value or dict(status="restored_archive_orthobench_independently_audited", job_id=22377,
                         reproduction_equal=True)
    code = ("import json;from pathlib import Path;Path(" + repr(str(output)) + ").write_text(" +
            repr(json.dumps(value)) + ");raise SystemExit(" + str(exit_code) + ")")
    return dict(name=name,command=[sys.executable,"-B","-c",code],result_path=str(output),
                parent_job=22377 if name=="restored_orthobench" else 22378)


def test_success_requires_matching_result_and_scheduler(tmp_path):
    result=module.execute(stage(tmp_path),tmp_path,parent_success=True)
    assert result["passed"] is True
    assert result["reproduction_equal"] is True
    assert result==json.loads((tmp_path/"restored_orthobench.json").read_text())


def test_zero_exit_mismatch_retained_as_failure(tmp_path):
    result=module.execute(stage(tmp_path,value=dict(status="restored_archive_orthobench_independently_audited",
        job_id=22377,reproduction_equal=False)),tmp_path,parent_success=True)
    assert result["returncode"]==0 and result["passed"] is False
    assert result["reproduction_equal"] is False


@pytest.mark.parametrize("defect",["scheduler","exit","wrong_job","wrong_status"])
def test_false_success_refused(tmp_path,defect):
    value=dict(status="restored_archive_orthobench_independently_audited",job_id=22377,reproduction_equal=True)
    if defect=="wrong_job":value["job_id"]=22376
    if defect=="wrong_status":value["status"]="unrecognized"
    result=module.execute(stage(tmp_path,value=value,exit_code=3 if defect=="exit" else 0),
                          tmp_path,parent_success=defect!="scheduler")
    assert result["passed"] is False


def test_failed_calibration_retained(tmp_path):
    result=module.execute(stage(tmp_path,"observer_calibration",dict(status="calibration_checks_failed",job_id=22378),1),
                          tmp_path,parent_success=True)
    assert result["reported_status"]=="calibration_checks_failed"
    assert result["passed"] is False


def test_successful_calibration(tmp_path):
    result=module.execute(stage(tmp_path,"observer_calibration",dict(status="calibration_checks_passed",job_id=22378)),
                          tmp_path,parent_success=True)
    assert result["passed"] is True


def test_existing_result_not_reexecuted(tmp_path):
    value=stage(tmp_path)
    path=Path(value["result_path"])
    path.write_text("original")
    result=module.execute(value,tmp_path,parent_success=True)
    assert result["error_type"]=="FileExistsError"
    assert result["returncode"] is None and result["passed"] is False
    assert path.read_text()=="original"


def test_timeout_owned_child_reaped(tmp_path):
    value=stage(tmp_path)
    value["command"]=[sys.executable,"-B","-c","import os,time;print(os.getpid(),flush=True);time.sleep(60)"]
    result=module.execute(value,tmp_path,timeout=1,parent_success=True)
    assert result["timed_out"] is True and result["returncode"]!=0
    assert result["passed"] is False
    pid=int((tmp_path/"restored_orthobench.log").read_text().strip())
    assert not Path(f"/proc/{pid}").exists()


@pytest.mark.parametrize("defect",["missing_batch","failed","signal","duplicate","unrelated_batch"])
def test_parent_scheduler_checks(defect):
    lines=["JobIDRaw|State|ExitCode", "22377|COMPLETED|0:0", "22377.batch|COMPLETED|0:0"]
    assert module.scheduler_success("\n".join(lines),22377) is True
    if defect=="missing_batch":lines.pop()
    elif defect=="failed":lines[-1]="22377.batch|FAILED|1:0"
    elif defect=="signal":lines[-1]="22377.batch|COMPLETED|0:9"
    elif defect=="duplicate":lines.append(lines[-1])
    else:lines[-1]="22378.batch|COMPLETED|0:0"
    assert module.scheduler_success("\n".join(lines),22377) is False


def test_protocol_digest_checked_first(tmp_path):
    path=tmp_path/"protocol.json"
    path.write_text("not JSON")
    with pytest.raises(ValueError,match="differs"):
        module.protocol(path,"wrong")


@pytest.mark.parametrize("raw",["", "wrong-header\nvalue", "JobIDRaw|State|ExitCode\n22377|COMPLETED",
                                    "JobIDRaw|State|ExitCode|State\n22377|COMPLETED|0:0|COMPLETED"])
def test_incomplete_accounting_is_not_success(raw):
    assert module.scheduler_success(raw,22377) is False


def test_failed_stage_does_not_suppress_other_audit(tmp_path,monkeypatch):
    first_dir,second_dir=tmp_path/"first",tmp_path/"second"
    first_dir.mkdir()
    second_dir.mkdir()
    stages=[stage(first_dir,exit_code=2),stage(second_dir,"observer_calibration",
             dict(status="calibration_checks_passed",job_id=22378))]
    value=dict(output=str(tmp_path/"execution"),interpreter=module.record(sys.executable),
               stages=stages,stage_timeout_seconds=10)
    monkeypatch.setattr(module,"protocol",lambda *_:value)
    monkeypatch.setattr(module.os,"uname",lambda:SimpleNamespace(nodename="bizon"))
    monkeypatch.setenv("SLURM_CPUS_PER_TASK","4")
    monkeypatch.setenv("SLURM_MEM_PER_NODE","32768")
    monkeypatch.setenv("SLURM_JOB_ID","99999")
    accounting="JobIDRaw|State|ExitCode\n22377|COMPLETED|0:0\n22377.batch|COMPLETED|0:0\n22378|COMPLETED|0:0\n22378.batch|COMPLETED|0:0\n"
    monkeypatch.setattr(module.subprocess,"run",lambda *_args,**_kwargs:SimpleNamespace(stdout=accounting,stderr=""))
    proto=tmp_path/"protocol.json"
    proto.write_text("synthetic unit-test protocol")
    result=module.run(proto,"synthetic")
    assert result["status"]=="completion_audits_failed"
    assert [s["passed"] for s in result["stages"]]==[False,True]
    assert result["publication_ready"] is False


def test_valid_protocol_checks_runner_and_sources(tmp_path):
    path=tmp_path/"protocol.json"
    pin=module.record(module.__file__)
    value=dict(schema="publication_completion_audits_v1",jobs=[22377,22378],production_timing=False,
        attempts_per_stage=1,stage_timeout_seconds=7200,interpreter=module.record(sys.executable),
        submission_script=pin,sources=[pin],output=str(tmp_path/"output"),
        stages=[dict(name="restored_orthobench",parent_job=22377),dict(name="observer_calibration",parent_job=22378)])
    path.write_text(json.dumps(value))
    assert module.protocol(path,module.record(path)["sha256"])==value
    value["sources"]=[]
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError,match="runner"):
        module.protocol(path,module.record(path)["sha256"])
