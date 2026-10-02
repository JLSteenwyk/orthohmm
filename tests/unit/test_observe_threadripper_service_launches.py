import copy
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from benchmark_tools import observe_threadripper_service_launches as module

SECRET = "never-retain-this-secret-value"

pytestmark = pytest.mark.usefixtures("synthetic_linux_boot_id")


@pytest.fixture(autouse=True)
def local_host(monkeypatch):
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="bizon"))


def reply(values):
    return dict(type="a{sv}", data=[values])


def property_values(fields):
    defaults = {"s": "", "as": [], "b": False, "u": 0, "a(sb)": []}
    return {name: dict(type=kind, data=copy.deepcopy(defaults[kind])) for name, kind in fields.items()}


def launch_values():
    return {key: dict(type="a(sasbttttuii)", data=[]) for key in module.EXEC_FIELDS}


def unit_list(name="sample.service"):
    return dict(type="a(ssssssouso)", data=[[[name, "Example", "loaded", "active", "running", "",
        "/org/freedesktop/systemd1/unit/sample_2eservice", 0, "", "/"]]])


def runner_factory(tmp_path, defect=None):
    executable = tmp_path / "declared_executable"
    executable.write_bytes(b"synthetic declared executable; never launched\n")
    unit = property_values(module.UNIT_FIELDS)
    unit["Id"]["data"] = "sample.service"
    path = tmp_path / "sample.service"
    path.write_text(f"[Service]\nExecStart={executable}\n")
    unit["FragmentPath"]["data"] = str(path)
    service = property_values(module.SERVICE_FIELDS)
    service.update(launch_values())
    service["ExecStart"]["data"] = [[str(executable), [str(executable), SECRET], False, 0, 0, 0, 0, 0, 0, 0]]
    service["Environment"] = dict(type="as", data=["CREDENTIAL=" + SECRET])
    if defect == "missing_file":
        unit["FragmentPath"]["data"] = str(tmp_path / "missing")
    if defect == "pending_reload":
        unit["NeedDaemonReload"]["data"] = True
    seen, calls = {"system": 0, "user": 0}, []
    def runner(command, **kwargs):
        calls.append(command)
        assert command[:2] == ["busctl", "--json=short"]
        assert kwargs["capture_output"] is True and 0 < kwargs["timeout"] <= 5
        assert command[6:8] == ["org.freedesktop.systemd1.Manager", "ListUnits"] or command[6:8] == ["org.freedesktop.DBus.Properties", "GetAll"]
        scope = command[2][2:]
        if defect == "nonzero":
            return SimpleNamespace(returncode=1, stdout=SECRET, stderr=SECRET)
        if defect == "timeout":
            raise subprocess.TimeoutExpired(command, 5, output=SECRET)
        if command[-1] == "ListUnits":
            seen[scope] += 1
            value = unit_list("changed.service" if defect == "changed" and seen[scope] > 1 else "sample.service")
        else:
            value = reply(unit if command[-1].endswith(".Unit") else service)
        return SimpleNamespace(returncode=0, stdout=json.dumps(value), stderr=SECRET)
    return runner, calls


def test_real_file_hash_without_content(tmp_path):
    path = tmp_path / "file"
    path.write_text(SECRET)
    result = module.file_identity(str(path))
    assert result["observed_stable"] is True
    assert result["bytes"] == len(SECRET)
    assert SECRET not in json.dumps(result)


def test_missing_or_nonregular_file_retained(tmp_path):
    for path in (tmp_path / "missing", tmp_path):
        result = module.file_identity(str(path))
        assert result["observed_stable"] is False
        assert "error" in result


def test_relative_file_rejected():
    with pytest.raises(ValueError, match="absolute"):
        module.file_identity("relative")


def test_large_declared_file_not_read(tmp_path):
    path = tmp_path / "large"
    with path.open("wb") as handle:
        handle.truncate(128 * 1024**2 + 1)
    assert module.file_identity(str(path))["observed_stable"] is False


def test_capture_is_readonly_and_never_approves_or_exposes_values(tmp_path, synthetic_linux_boot_id):
    runner, calls = runner_factory(tmp_path)
    result = module.collect(run=runner)
    assert result["boot_unchanged"] is True
    assert result["boot_id"] == synthetic_linux_boot_id
    assert not result["errors"]
    assert result["policy_approved"] is False and result["scientific_timings_admitted"] is False
    assert SECRET not in json.dumps(result)
    for scope in result["scopes"].values():
        row = scope["units"]["sample.service"]
        assert row["service_approved"] is False
        assert row["launches"]["ExecStart"][0]["argc"] == 2
        assert row["files"]
        assert not scope["inventory_changes"]
    assert len(calls) == 8


@pytest.mark.parametrize("defect", ["nonzero", "timeout", "missing_file", "pending_reload", "changed"])
def test_failures_and_changes_remain_explicit_without_secrets(tmp_path, defect):
    runner, calls = runner_factory(tmp_path, defect)
    result = module.collect(run=runner)
    assert result["policy_approved"] is False
    assert SECRET not in json.dumps(result)
    assert result["errors"] or any(v["inventory_changes"] for v in result["scopes"].values())


def test_deadline_prevents_bus_calls(tmp_path):
    runner, calls = runner_factory(tmp_path)
    values = iter([0, 181, 182])
    result = module.collect(run=runner, clock=lambda: next(values))
    assert not calls
    assert len(result["errors"]) == 2


@pytest.mark.parametrize("defect", ["outer", "missing", "kind", "data", "environment_file"])
def test_malformed_properties_rejected(defect):
    value = reply(property_values(module.SERVICE_FIELDS))
    if defect == "outer":
        value["type"] = "s"
    elif defect == "missing":
        del value["data"][0]["MainPID"]
    elif defect == "kind":
        value["data"][0]["MainPID"]["type"] = "s"
    elif defect == "data":
        value["data"][0]["MainPID"]["data"] = True
    else:
        value["data"][0]["EnvironmentFiles"]["data"] = [["/path", "not a bool"]]
    with pytest.raises(ValueError):
        module.properties(value, module.SERVICE_FIELDS)


@pytest.mark.parametrize("defect", ["missing", "kind", "relative", "argv", "length"])
def test_malformed_launches_rejected(defect):
    value = reply(launch_values())
    value["data"][0]["ExecStart"]["data"] = [["/bin/true", ["true"], False, 0, 0, 0, 0, 0, 0, 0]]
    command = value["data"][0]["ExecStart"]["data"][0]
    if defect == "missing":
        del value["data"][0]["ExecStop"]
    elif defect == "kind":
        value["data"][0]["ExecStart"]["type"] = "as"
    elif defect == "relative":
        command[0] = "relative/path"
    elif defect == "argv":
        command[1] = "not a list"
    else:
        command.pop()
    with pytest.raises(ValueError):
        module.launches(value)


def test_bare_command_retained_without_guessed_resolution():
    value = reply(launch_values())
    value["data"][0]["ExecStart"]["data"] = [["systemctl", ["systemctl", SECRET], False, 0, 0, 0, 0, 0, 0, 0]]
    result = module.launches(value)
    assert result["ExecStart"][0]["executable"] == "systemctl"
    assert SECRET not in json.dumps(result)


def test_bare_command_collection_has_explicit_resolution_gap(tmp_path):
    base, calls = runner_factory(tmp_path)
    def runner(command, **kwargs):
        result = base(command, **kwargs)
        value = json.loads(result.stdout)
        if command[-1].endswith(".Service"):
            value["data"][0]["ExecStart"]["data"][0][0] = "true"
            result.stdout = json.dumps(value)
        return result
    result = module.collect(run=runner)
    for scope in result["scopes"].values():
        row = scope["units"]["sample.service"]
        assert row["errors"] == [dict(kind="bare_executable_resolution_unverified", executable="true")]
        assert all(file["role"] != "declared_executable" for file in row["files"])


def test_remote_host_refused(monkeypatch):
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="other"))
    with pytest.raises(ValueError, match="Threadripper"):
        module.collect()


@pytest.mark.parametrize("defect", ["duplicate", "empty", "path", "job_type"])
def test_malformed_unit_lists_rejected(defect):
    value = unit_list()
    if defect == "duplicate":
        value["data"][0].append(value["data"][0][0])
    elif defect == "empty":
        value["data"][0].clear()
    elif defect == "path":
        value["data"][0][0][6] = "/other"
    else:
        value["data"][0][0][7] = True
    with pytest.raises(ValueError):
        module.units(value)


def test_reboot_between_service_reads_remains_unapproved(tmp_path, monkeypatch, synthetic_linux_boot_id):
    runner, calls = runner_factory(tmp_path)
    values = iter([synthetic_linux_boot_id, "00000000-0000-4000-8000-000000000043"])
    read_text = Path.read_text
    def changed(path, *args, **kwargs):
        if path == Path("/proc/sys/kernel/random/boot_id"):
            return next(values) + "\n"
        return read_text(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", changed)
    result = module.collect(run=runner)
    assert result["boot_id"] == synthetic_linux_boot_id
    assert result["boot_unchanged"] is False
    assert result["policy_approved"] is result["scientific_timings_admitted"] is False
    assert len(calls) == 8 and SECRET not in json.dumps(result)


def test_missing_boot_fails_before_any_service_query(tmp_path, monkeypatch):
    runner, calls = runner_factory(tmp_path)
    read_text = Path.read_text
    def missing(path, *args, **kwargs):
        if path == Path("/proc/sys/kernel/random/boot_id"):
            raise FileNotFoundError("synthetic missing boot counter")
        return read_text(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", missing)
    with pytest.raises(FileNotFoundError):
        module.collect(run=runner)
    assert not calls
