from copy import deepcopy
import subprocess
from types import SimpleNamespace

import pytest

from benchmark_tools import observe_threadripper_services as service


def test_commands_are_local_read_only_and_exclude_sensitive_properties():
    commands = service.commands()
    assert len(commands) == 5
    for argv in commands.values():
        assert argv[0] in {"squeue", "systemctl"}
        assert not set(argv) & {"stop", "start", "restart", "mask", "ssh", "--system"}
        assert not any("ExecStart" in arg or "Environment" in arg for arg in argv)


@pytest.mark.parametrize("failure", [None, "nonzero", "timeout", "missing_command", "parse", "unreadable"])
def test_observation_retains_errors_without_approval(monkeypatch, failure):
    monkeypatch.setattr(service.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    calls=[]
    def run(argv, **kwargs):
        calls.append(argv)
        assert kwargs["timeout"] == 10 and kwargs["env"]["LC_ALL"] == "C"
        if failure == "timeout": raise subprocess.TimeoutExpired(argv, 10)
        if failure == "missing_command": raise FileNotFoundError("synthetic")
        return SimpleNamespace(returncode=1 if failure == "nonzero" else 0, stdout="synthetic", stderr="")
    def fingerprint(raw):
        if failure == "parse": raise ValueError("synthetic parse failure")
        return {"fixture.service": dict(properties={}, files=[dict(path="/synthetic", error_type="PermissionError")]
                                        if failure == "unreadable" else [])}
    result=service.collect(run=run, fingerprint=fingerprint)
    assert len(calls)==5 and not result["policy_approved"] and not result["scientific_timings_admitted"]
    assert bool(result["errors"]) == (failure is not None)


def test_wrong_host_refused(monkeypatch):
    monkeypatch.setattr(service.os, "uname", lambda: SimpleNamespace(nodename="other"))
    with pytest.raises(ValueError): service.collect()


@pytest.mark.parametrize("change", ["added", "removed", "changed", "missing_inventory", None])
def test_comparison_keeps_changes_and_missing_evidence(change):
    a=dict(configurations={scope:{"fixture.service":dict(properties={},files=[])}
                          for scope in ("system_configuration","user_configuration")},errors=[])
    b=deepcopy(a); scope="user_configuration"
    if change=="added": b["configurations"][scope]["other.service"]={}
    elif change=="removed": del b["configurations"][scope]["fixture.service"]
    elif change=="changed": b["configurations"][scope]["fixture.service"]["files"]=[{"sha256":"changed"}]
    elif change=="missing_inventory": del b["configurations"][scope]
    result=service.compare(a,b)
    assert not result["policy_approved"]
    assert [x["kind"] for x in result["changes"]] == ([] if change is None else [change])
