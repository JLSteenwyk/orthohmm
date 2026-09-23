import hashlib
import json
from types import SimpleNamespace
import subprocess

import pytest

from benchmark_tools import dgx_service_guard as module


@pytest.fixture
def guard(tmp_path, monkeypatch):
    source = tmp_path / "user/service"
    source.parent.mkdir()
    source.write_text("frozen service\n")
    enabled = tmp_path / "enabled"
    enabled.symlink_to(source)
    monkeypatch.setattr(module, "SOURCE", source)
    monkeypatch.setattr(module, "SOURCE_SHA", hashlib.sha256(source.read_bytes()).hexdigest())
    monkeypatch.setattr(module, "ENABLED", enabled)
    monkeypatch.setattr(module, "MASK", tmp_path / "runtime/user.control/service")
    monkeypatch.setattr(module, "LOWER_MASK", tmp_path / "runtime/user/service")
    monkeypatch.setenv("ORTHOHMM_APPROVED_SAMWISE_STOP", "20260920")
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="spark-7ff0"))
    monkeypatch.setattr(module.os, "getuid", lambda: 1000)
    state = {"prior": "activating", "accounting": "42|COMPLETED|2026-09-23T15:00:00\n",
             "commands": [], "timeout": False}

    def run(argv, **kwargs):
        assert kwargs["timeout"] == 15
        state["commands"].append(argv)
        if state["timeout"]:
            raise subprocess.TimeoutExpired(argv, 15)
        if argv[0] == "sacct":
            return SimpleNamespace(returncode=0, stdout=state["accounting"], stderr="")
        stdout = ""
        if argv[2] == "show":
            masked = module.MASK.is_symlink()
            values = dict(LoadState="masked" if masked else "loaded",
                ActiveState="inactive" if masked else state["prior"], SubState="dead" if masked else "auto-restart",
                MainPID="0", FragmentPath=str(module.MASK if masked else source))
            stdout = "\n".join(f"{k}={v}" for k, v in values.items()) + "\n"
        return SimpleNamespace(returncode=0, stdout=stdout, stderr="")

    monkeypatch.setattr(module.subprocess, "run", run)
    return module.ServiceGuard(tmp_path / "receipts"), state


@pytest.mark.parametrize("prior", ["active", "activating", "failed", "inactive"])
def test_normal_lifecycle_preserves_configuration(guard, prior):
    g, state = guard
    state["prior"] = prior
    original = module.SOURCE.read_bytes()
    g.begin()
    g.check()
    g.before_submission()
    g.bind_job(42)
    g.restore()
    assert not module.MASK.is_symlink()
    assert module.SOURCE.read_bytes() == original
    assert (["systemctl", "--user", "start", module.UNIT] in state["commands"]) == (prior != "inactive")
    assert json.loads((g.evidence / "restored.json").read_text())["application_health_verified"] is False


@pytest.mark.parametrize("accounting", ["", "42|RUNNING|Unknown\n", "42|COMPLETING|Unknown\n",
    "42|COMPLETED|Unknown\n", "42|COMPLETED|garbage\n",
    "43|COMPLETED|2026-09-23\n", "42.batch|COMPLETED|2026-09-23\n",
    "42|COMPLETED|2026-09-23\n42|FAILED|2026-09-23\n"])
def test_unverified_terminal_state_keeps_mask(guard, accounting):
    g, state = guard
    g.begin()
    g.bind_job(42)
    state["accounting"] = accounting
    with pytest.raises(ValueError):
        g.restore()
    assert module.MASK.is_symlink()
    assert not g.restored


def test_poll_timeout_keeps_mask_and_records_error(guard):
    g, state = guard
    g.begin()
    g.bind_job(42)
    state["timeout"] = True
    with pytest.raises(subprocess.TimeoutExpired):
        g.restore()
    assert module.MASK.is_symlink()
    receipt = json.loads((g.evidence / f"command_{g.sequence:06d}.json").read_text())
    assert receipt["error_type"] == "TimeoutExpired"


def test_unresolved_submission_does_not_restore(guard):
    g, _ = guard
    g.begin()
    g.before_submission()
    with pytest.raises(ValueError, match="identity unresolved"):
        g.restore()
    with pytest.raises(ValueError, match="repeat"):
        g.before_submission()
    assert module.MASK.is_symlink()


@pytest.mark.parametrize("change", ["source", "mask", "missing_mask", "enabled"])
def test_drift_prevents_check_and_restoration(guard, change):
    g, _ = guard
    g.begin()
    if change == "source":
        module.SOURCE.write_text("changed")
    elif change == "enabled":
        module.ENABLED.unlink()
    else:
        module.MASK.unlink()
        if change == "mask":
            module.MASK.symlink_to("/foreign")
    with pytest.raises(ValueError):
        g.check()
    with pytest.raises(ValueError):
        g.restore()
    assert not g.restored


def test_binding_survives_failed_mask_observation(guard):
    g, _ = guard
    g.begin()
    module.MASK.unlink()
    with pytest.raises(ValueError):
        g.bind_job(42)
    assert g.job == 42


def test_no_job_prelaunch_restoration(guard):
    g, state = guard
    g.begin()
    g.restore()
    assert not any(command[0] == "sacct" for command in state["commands"])


def test_requires_explicit_approval(guard, monkeypatch):
    g, state = guard
    monkeypatch.delenv("ORTHOHMM_APPROVED_SAMWISE_STOP")
    with pytest.raises(ValueError):
        g.begin()
    assert not state["commands"]


def test_refuses_preexisting_mask(guard):
    g, state = guard
    module.MASK.parent.mkdir(parents=True)
    module.MASK.symlink_to("/dev/null")
    with pytest.raises(ValueError):
        g.begin()
    assert not state["commands"]


@pytest.mark.parametrize("status", sorted(module.TERMINAL) + ["CANCELLED by 1000"])
def test_terminal_failures_also_allow_restoration(guard, status):
    g, state = guard
    g.begin()
    g.bind_job(42)
    state["accounting"] = f"42|{status}|2026-09-23T15:00:00\n"
    g.restore()
    assert g.restored


def test_stop_observation_failure_allows_prelaunch_restoration(guard, monkeypatch):
    g, state = guard
    original = g.ctl

    def ctl(*args):
        if args[0] == "stop":
            raise subprocess.TimeoutExpired(args, 15)
        return original(*args)

    monkeypatch.setattr(g, "ctl", ctl)
    with pytest.raises(subprocess.TimeoutExpired):
        g.begin()
    assert g.stopped and not g.owned
    g.restore()
    assert g.restored
