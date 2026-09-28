"""Read-only local service inventory for policy review, never policy approval."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from benchmark_tools.observe_dgx_environment import configuration_fingerprints
from benchmark_tools.prepare_ob_candidate_neighborhood import record


PROPERTIES = "Id,LoadState,FragmentPath,DropInPaths,NeedDaemonReload"


def commands():
    result = {"scheduler": ["squeue", "--noheader", "--format=%i,%T,%j,%C,%N"]}
    for scope, prefix in (("system", ["systemctl"]), ("user", ["systemctl", "--user"])):
        result[scope + "_running"] = prefix + ["list-units", "--type=service", "--state=running",
                                               "--no-pager", "--plain", "--no-legend"]
        result[scope + "_configuration"] = prefix + ["show", "*.service", "*.timer",
            "--property=" + PROPERTIES, "--no-pager"]
    return result


def collect(*, run=subprocess.run, fingerprint=configuration_fingerprints):
    if os.uname().nodename != "bizon":
        raise ValueError("Require the authorized local Threadripper host")
    started = time.time_ns()
    result = dict(schema="threadripper_service_inventory_v1", host="bizon",
        started_unix_ns=started, commands={}, configurations={}, errors=[],
        policy_approved=False, scientific_timings_admitted=False,
        sources=[record(Path(__file__)), record(Path(__file__).with_name("observe_dgx_environment.py"))])
    for name, argv in commands().items():
        row = dict(argv=argv, started_monotonic_ns=time.monotonic_ns())
        try:
            completed = run(argv, capture_output=True, text=True, timeout=10,
                env=dict(os.environ, LC_ALL="C", SYSTEMD_COLORS="0", SYSTEMD_PAGER=""))
            row.update(returncode=completed.returncode, stdout=completed.stdout, stderr=completed.stderr)
            if completed.returncode != 0:
                result["errors"].append(dict(command=name, type="nonzero_exit"))
        except (OSError, subprocess.TimeoutExpired) as error:
            row.update(error_type=type(error).__name__, error=str(error))
            result["errors"].append(dict(command=name, type=type(error).__name__))
        row["finished_monotonic_ns"] = time.monotonic_ns()
        result["commands"][name] = row
        if name.endswith("_configuration") and row.get("returncode") == 0:
            try:
                units = fingerprint(row["stdout"])
                result["configurations"][name] = units
                for unit, entry in units.items():
                    for item in entry["files"]:
                        if "error_type" in item:
                            result["errors"].append(dict(command=name, unit=unit,
                                path=item["path"], type=item["error_type"]))
            except (ValueError, OSError) as error:
                result["errors"].append(dict(command=name, type=type(error).__name__, error=str(error)))
    result.update(finished_unix_ns=time.time_ns(), limitations=[
        "Loaded service/timer units only; scopes, transient processes and unloaded units are not a complete configuration inventory.",
        "Fingerprints cover fragments/drop-ins, not external scripts, EnvironmentFiles, libraries or runtime manager overrides.",
        "No unit contents, command lines or environment values are retained; only selected properties and file digests.",
        "Commands and file reads are sequential, not an atomic snapshot or an approved ordinary-background policy.",
        "The imported fingerprint parser is reused locally; no DGX collector or remote connection is invoked."])
    return result


def compare(before, after):
    changes = []
    for scope in ("system_configuration", "user_configuration"):
        a, b = before["configurations"].get(scope), after["configurations"].get(scope)
        if a is None or b is None:
            changes.append(dict(scope=scope, kind="missing_inventory"))
            continue
        for unit in sorted(set(a) | set(b)):
            if a.get(unit) != b.get(unit):
                changes.append(dict(scope=scope, unit=unit,
                    kind="added" if unit not in a else "removed" if unit not in b else "changed"))
    return dict(changes=changes, observation_errors=len(before["errors"]) + len(after["errors"]),
                policy_approved=False, scientific_timings_admitted=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as handle:
        first, second = collect(), collect()
        json.dump(dict(schema="threadripper_service_pair_v1", snapshots=[first, second],
                       comparison=compare(first, second)), handle, indent=2, sort_keys=True)
        handle.write("\n")
