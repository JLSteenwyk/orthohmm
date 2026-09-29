"""Read effective local systemd launch/file evidence without approving services."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import time

from benchmark_tools.run_integrated_full_job import record

DESTINATION = "org.freedesktop.systemd1"
MANAGER = "/org/freedesktop/systemd1"
UNIT_FIELDS = dict(Id="s", LoadState="s", ActiveState="s", SubState="s", FragmentPath="s",
                   DropInPaths="as", NeedDaemonReload="b")
SERVICE_FIELDS = dict(WorkingDirectory="s", RootDirectory="s", User="s", Group="s",
                      MainPID="u", ControlGroup="s", Restart="s", EnvironmentFiles="a(sb)")
EXEC_FIELDS = ("ExecStart", "ExecStartPre", "ExecStartPost", "ExecReload", "ExecStop", "ExecStopPost")


def error_record(error):
    return dict(type=type(error).__name__, errno=getattr(error, "errno", None))


def file_identity(name):
    path = Path(name)
    if not path.is_absolute():
        raise ValueError("Require absolute declared file paths")
    try:
        before = path.stat()
        if not stat.S_ISREG(before.st_mode):
            raise ValueError("Declared path is not a regular file")
        if before.st_size > 128 * 1024**2:
            raise ValueError("Declared file exceeds 128-MiB read envelope")
        digest, size = hashlib.sha256(), 0
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024**2), b""):
                size += len(chunk)
                if size > 128 * 1024**2:
                    raise ValueError("Declared file grew beyond the read envelope")
                digest.update(chunk)
        pin = dict(path=str(path), bytes=size, sha256=digest.hexdigest())
        after = path.stat()
        if ((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns)
                != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns)
                or pin["bytes"] != before.st_size):
            raise ValueError("Declared file changed while reading")
        return dict(pin, resolved_path=str(path.resolve()), observed_stable=True)
    except (OSError, ValueError) as error:
        return dict(path=str(path), error=error_record(error), observed_stable=False)


def properties(reply, fields):
    if reply.get("type") != "a{sv}" or not isinstance(reply.get("data"), list) or len(reply["data"]) != 1:
        raise ValueError("Invalid property reply")
    values = reply["data"][0]
    if not isinstance(values, dict):
        raise ValueError("Invalid property dictionary")
    result = {}
    for name, kind in fields.items():
        value = values.get(name)
        if not isinstance(value, dict) or value.get("type") != kind or "data" not in value:
            raise ValueError("Missing or mistyped required property: " + name)
        item = value["data"]
        valid = ((kind == "s" and isinstance(item, str)) or (kind == "b" and type(item) is bool)
                 or (kind == "u" and type(item) is int and item >= 0)
                 or (kind in {"as", "a(sb)"} and isinstance(item, list)))
        if not valid:
            raise ValueError("Invalid property value: " + name)
        if kind == "as" and not all(isinstance(v, str) for v in item):
            raise ValueError("Invalid string list")
        if kind == "a(sb)" and not all(isinstance(v, list) and len(v) == 2
                                      and isinstance(v[0], str) and type(v[1]) is bool for v in item):
            raise ValueError("Invalid declared environment-file list")
        result[name] = item
    return result


def launches(reply):
    if reply.get("type") != "a{sv}" or not isinstance(reply.get("data"), list) or len(reply["data"]) != 1:
        raise ValueError("Invalid launch-property reply")
    values, result = reply["data"][0], {}
    for field in EXEC_FIELDS:
        value = values.get(field)
        if (not isinstance(value, dict) or value.get("type") != "a(sasbttttuii)"
                or not isinstance(value.get("data"), list)):
            raise ValueError("Missing or mistyped launch property: " + field)
        commands = []
        for command in value["data"]:
            if (not isinstance(command, list) or len(command) != 10
                    or not isinstance(command[0], str) or not command[0]
                    or (not Path(command[0]).is_absolute() and "/" in command[0])
                    or not isinstance(command[1], list) or not all(isinstance(v, str) for v in command[1])
                    or type(command[2]) is not bool):
                raise ValueError("Invalid structured launch command")
            encoded = json.dumps(command[1], ensure_ascii=True, separators=(",", ":")).encode()
            commands.append(dict(executable=command[0], argc=len(command[1]), ignore_failure=command[2],
                                 argv_bytes=len(encoded), argv_sha256=hashlib.sha256(encoded).hexdigest()))
        result[field] = commands
    return result


def units(reply):
    if reply.get("type") != "a(ssssssouso)" or not isinstance(reply.get("data"), list) or len(reply["data"]) != 1:
        raise ValueError("Invalid ListUnits reply")
    result = {}
    for row in reply["data"][0]:
        if (not isinstance(row, list) or len(row) != 10
                or not all(isinstance(v, str) for v in row[:7] + row[8:])
                or type(row[7]) is not int or not 0 <= row[7] < 2**32):
            raise ValueError("Invalid structured unit entry")
        name = row[0]
        if name.endswith((".service", ".timer")):
            if name in result or not row[6].startswith(MANAGER + "/unit/"):
                raise ValueError("Duplicate unit or unexpected object path")
            result[name] = dict(object_path=row[6], load_state=row[2], active_state=row[3], sub_state=row[4])
    if not result or len(result) > 512:
        raise ValueError("Unit inventory outside bounded capture scope")
    return result


def collect(*, run=subprocess.run, fingerprint=file_identity, clock=time.monotonic):
    if os.uname().nodename != "bizon":
        raise ValueError("Require the approved local Threadripper")
    boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    deadline = clock() + 180
    result = dict(schema="threadripper_effective_service_inventory_v1", host="bizon", boot_id=boot,
                  started_unix_ns=time.time_ns(), scopes={}, errors=[], policy_approved=False,
                  scientific_timings_admitted=False, source=record(__file__))
    def call(scope, object_path, interface, method, *args):
        remaining = deadline - clock()
        if remaining <= 0:
            raise TimeoutError("Global service capture deadline")
        command = ["busctl", "--json=short", "--" + scope, "call", DESTINATION, object_path,
                   interface, method, *args]
        completed = run(command, capture_output=True, text=True, timeout=min(5., remaining))
        if completed.returncode:
            raise RuntimeError("Read-only bus query failed")
        # GetAll can include secrets. Never retain its stdout/stderr or unused fields.
        return json.loads(completed.stdout)
    for scope in ("system", "user"):
        observations = dict(units={}, errors=[], inventory_changes=[])
        result["scopes"][scope] = observations
        try:
            before = units(call(scope, MANAGER, DESTINATION + ".Manager", "ListUnits"))
            observations["before"] = before
            for name, listed in sorted(before.items()):
                row = dict(listed=listed, files=[], errors=[], service_approved=False)
                observations["units"][name] = row
                try:
                    object_path = listed["object_path"]
                    unit_reply = call(scope, object_path, "org.freedesktop.DBus.Properties", "GetAll", "s", DESTINATION + ".Unit")
                    row["unit"] = properties(unit_reply, UNIT_FIELDS)
                    files = [(v, "unit_fragment") for v in [row["unit"]["FragmentPath"]] if v]
                    files.extend((v, "unit_dropin") for v in row["unit"]["DropInPaths"])
                    if name.endswith(".service"):
                        service_reply = call(scope, object_path, "org.freedesktop.DBus.Properties", "GetAll", "s", DESTINATION + ".Service")
                        row["service"] = properties(service_reply, SERVICE_FIELDS)
                        row["launches"] = launches(service_reply)
                        files.extend((v[0], "declared_environment_file") for v in row["service"]["EnvironmentFiles"])
                        files.extend((command["executable"], "declared_executable")
                                     for commands in row["launches"].values() for command in commands
                                     if Path(command["executable"]).is_absolute())
                        for commands in row["launches"].values():
                            for command in commands:
                                if not Path(command["executable"]).is_absolute():
                                    row["errors"].append(dict(kind="bare_executable_resolution_unverified",
                                                              executable=command["executable"]))
                    for path, role in sorted(set(files)):
                        if clock() >= deadline:
                            raise TimeoutError("Global service capture deadline")
                        pin = fingerprint(path)
                        row["files"].append(dict(role=role, observation=pin))
                        if not pin.get("observed_stable"):
                            row["errors"].append(dict(kind="declared_file_unverified", role=role, path=path))
                    if row["unit"]["NeedDaemonReload"]:
                        row["errors"].append(dict(kind="pending_daemon_reload"))
                except (OSError, ValueError, KeyError, TypeError, RuntimeError, subprocess.TimeoutExpired) as error:
                    row["errors"].append(error_record(error))
                if row["errors"]:
                    observations["errors"].append(dict(unit=name, errors=row["errors"]))
            after = units(call(scope, MANAGER, DESTINATION + ".Manager", "ListUnits"))
            observations["after"] = after
            observations["inventory_changes"] = [name for name in sorted(before.keys() | after.keys())
                                                  if before.get(name) != after.get(name)]
        except (OSError, ValueError, KeyError, TypeError, RuntimeError, subprocess.TimeoutExpired) as error:
            observations["errors"].append(error_record(error))
        result["errors"].extend(dict(scope=scope, detail=v) for v in observations["errors"])
    result["finished_unix_ns"] = time.time_ns()
    result["boot_unchanged"] = boot == Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    result["limitations"] = ["Effective loaded-unit properties are sequential observations, not policy approval.",
        "Only unit fragments/drop-ins, declared executables and EnvironmentFiles are fingerprinted.",
        "Command arguments are represented by digest/length only; environment values and query errors are not retained.",
        "Interpreted scripts, dependencies, runtime overrides and unlisted/transient scopes require separate review.",
        "Bare executable names are retained, not resolved using the observer's PATH or fingerprinted as guessed files.",
        "Bus calls use a 180-second cooperative budget; file reads are capped at 128 MiB, not a hard wall-time guarantee.",
        "Missing or unreadable declared files, pending reloads and inventory changes remain explicit.",
        "No service action, signal, remote connection, quiet-host certification or timing admission."]
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as handle:
        json.dump(collect(), handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
