"""Check observed native Perl modules and memory mappings against a runtime snapshot."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify
from benchmark_tools.probe_orthomcl_bpo_parity import TOOL
from benchmark_tools.run_qfo_corrected_blast import environment

VERSIONS = {"perl": "v5.26.2", "bioperl_searchio": "1.007002", "storable": "3.15", "orthomcl": "1.4"}
EXTERNAL = ["/bin/sh", "/lib/x86_64-linux-gnu/libcrypt.so.1", "/lib64/ld-linux-x86-64.so.2"]


def validate(runtime, observed, cwd):
    if observed["versions"] != VERSIONS:
        raise ValueError("Unexpected native Perl/module versions")
    if runtime["external_symlinks"] != EXTERNAL:
        raise ValueError("Unreviewed external runtime symlinks")
    files = {}
    for item in runtime["records"]:
        if item["kind"] == "file":
            files[str(Path(item["path"]).resolve())] = (item["bytes"], item["sha256"])
        elif item["kind"] == "symlink" and "target_sha256" in item:
            files[item["resolved"]] = (item["target_bytes"], item["target_sha256"])
    if not observed["search_path"] or any(not isinstance(p, str) or (p != "." and not Path(p).is_absolute())
                                           for p in observed["search_path"]):
        raise ValueError("Relative or executable-hook Perl search path")
    if "." in observed["search_path"] and any(
            p.name not in {"probe.json", "probe.log"} or p.is_symlink() or not p.is_file()
            for p in cwd.iterdir()):
        raise ValueError("Unreviewed file in current-directory Perl search path")
    if not observed["loaded_modules"] or not observed["mapped_files"]:
        raise ValueError("Empty loaded runtime inventory")
    paths = [*observed["loaded_modules"].values(), *observed["mapped_files"]]
    if any(not isinstance(path, str) or not Path(path).is_absolute() for path in paths):
        raise ValueError("Invalid loaded runtime path")
    checked = []
    for path in sorted(set(paths)):
        item = record(Path(path).resolve())
        if files.get(item["path"]) != (item["bytes"], item["sha256"]):
            raise ValueError("Loaded module/library outside snapshot or changed: " + path)
        checked.append(item)
    return checked


def inspect(runtime_path, runtime_sha, output):
    if output.exists():
        raise FileExistsError(output)
    runtime = read_frozen(runtime_path, runtime_sha)
    verify(runtime)
    driver = Path(__file__).with_suffix(".pl")
    checked = [record(runtime_path), record(__file__), record(driver),
               record(Path(__file__).with_name("snapshot_runtime_trees.py")),
               record(Path(__file__).with_name("run_qfo_corrected_blast.py"))]
    output.mkdir(parents=True, exist_ok=False)
    argv = [str(TOOL / "venv_orthomcl/bin/perl"), "-I" + str(TOOL), str(driver)]
    done = subprocess.run(argv, cwd=output, env=environment(), capture_output=True, timeout=60)
    (output / "probe.json").write_bytes(done.stdout)
    (output / "probe.log").write_bytes(done.stderr)
    if done.returncode or done.stderr:
        raise ValueError("Native runtime probe failed or emitted diagnostics")
    observed = json.loads(done.stdout)
    loaded = validate(runtime, observed, output)
    verify(runtime)
    for item in [*checked, *loaded]:
        check(item)
    report = {"status": "native_orthomcl_perl_runtime_observed_and_bound", "runtime": record(runtime_path),
              "checked_records": checked, "loaded_records": loaded, "observed": observed,
              "command": argv, "cwd": str(output), "environment": environment(),
              "outputs": [record(output / "probe.json"), record(output / "probe.log")],
              "accuracy_admitted": False, "execution_authorized": False, "publication_ready": False,
              "limitations": [
                  "Observed imports/maps for a module-load probe, not every future process or syscall.",
                  "Snapshot includes native sources, complete environment, MCL, selected system libraries, shell and date.",
                  "Shell/date/MCL are inventoried but not executed by this Perl probe.",
                  "Legacy Perl includes dot in its search path; this probe uses a fresh directory restricted to its two output files.",
                  "Production working-directory module lookup requires an explicit policy; this probe does not authorize it.",
                  "Before/after hashes cannot detect all temporary changes; this is not a hermetic operating-system image.",
                  "Configured native module and pair-parallel script still require separate source freezing and validation."]}
    with (output / "report.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime", type=Path, required=True)
    parser.add_argument("--runtime-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = inspect(args.runtime.resolve(), args.runtime_sha256, args.output.resolve())
    print(json.dumps({"status": report["status"], "loaded_records": len(report["loaded_records"])}))
