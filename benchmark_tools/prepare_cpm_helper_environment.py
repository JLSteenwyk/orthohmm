"""Prepare a separate offline runtime for the frozen high-CPM benchmark helper."""

import argparse
import json
import os
from pathlib import Path
import shutil

from benchmark_tools.audit_recovery_install import installed_payload
from benchmark_tools.build_private_timing_environment import canonical, hash_lock, wheel_metadata
from benchmark_tools.probe_cpm_private_runtime import (
    NATIVE, PINS, RESTORED, bound_base_records, execute, package_records, unique_records,
)
from benchmark_tools.probe_cpm_partition_parser import check, record

PROTOCOL = "benchmark_tools/results/QFO_CPM_HELPER_ENVIRONMENT_PROTOCOL_20260930.md"
BASE = Path("/tmp/orthohmm-base-reconstruction-20260927/python-runtime")
IMPORTS = ("benchmark_tools.replay_high_sensitivity",
           "benchmark_tools.audit_historical_profile_ablation",
           "benchmark_tools.orthobench_stage_diagnostics",
           "orthohmm.accuracy", "orthohmm.refinement", "numpy")


def select_wheels(plan, audit):
    command = plan["command"]
    assets = Path(command[command.index("--assets") + 1]) / "wheels"
    readers = Path(command[command.index("--reader-wheels") + 1])
    selected = {canonical(row["name"]): row["version"]
                for row in audit["inference"]["inventory"]["package"]}
    if len(selected) != 11 or "biopython" in selected:
        raise ValueError("Unexpected restored inference inventory")
    selected["biopython"] = "1.87"
    wheels, seen = [], set()
    for row in plan["checked_records"]:
        path = Path(row["path"])
        if path.suffix != ".whl" or path.parent not in (assets, readers):
            continue
        check([row])
        name, version = wheel_metadata(path)
        key = canonical(name)
        if path.parent == readers and key != "biopython":
            continue
        if key not in selected or selected[key] != version or key in seen:
            raise ValueError("Unexpected, duplicate or wrong-version helper wheel")
        seen.add(key)
        wheels.append(dict(row, name=name, version=version))
    if seen != set(selected):
        raise ValueError("Missing helper wheel")
    return selected, wheels


def validate_imports(observed, selected, prefix, launcher, base):
    inventory = {canonical(k): v for k, v in observed["packages"].items()}
    if inventory != selected:
        raise ValueError("Helper package inventory differs")
    if Path(observed["base_prefix"]).resolve() != base.resolve():
        raise ValueError("Helper interpreter base differs")
    if not observed["gc_enabled"] or observed["gc_thresholds"] != [700, 10, 10]:
        raise ValueError("Helper GC defaults differ")
    for name in IMPORTS:
        if name not in observed["modules"]:
            raise ValueError("Missing worker import")
    for name, value in observed["modules"].items():
        path = Path(value).resolve()
        if not any(path.is_relative_to(p.resolve()) for p in (prefix, launcher, base)):
            raise ValueError("Helper module outside declared roots: " + name)
        if name.startswith(("benchmark_tools.", "orthohmm.")) and not path.is_relative_to(launcher):
            raise ValueError("Nonfrozen scientific/helper import")


def prepare(root, output, protocol_sha):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    prerequisites = [record(root / name) for name in PINS]
    if [row["sha256"] for row in prerequisites] != list(PINS.values()):
        raise ValueError("Changed restored-runtime prerequisite")
    protocol = record(root / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Unreviewed helper-environment protocol")
    plan = json.loads((root / RESTORED / "plan.json").read_bytes())
    audit = json.loads((root / RESTORED / "independent_admission/package_audits.json").read_bytes())
    selected, wheels = select_wheels(plan, audit)
    base_records = bound_base_records(plan, BASE)
    base_python = record(BASE / "bin/python")
    if base_python not in base_records:
        raise ValueError("Unbound helper base interpreter")
    modules = json.loads((root / NATIVE / "refinement.json").read_bytes())["modules"]
    sources = [record(root / ("benchmark_tools/" + name)) for name in (
        "prepare_cpm_helper_environment.py", "probe_cpm_private_runtime.py",
        "probe_cpm_partition_parser.py", "audit_recovery_install.py",
        "build_private_timing_environment.py", "audit_frozen_overlay_install.py",
        "audit_leiden_recovery_wheel.py", "prepare_ob_candidate_neighborhood.py",
        "verify_frozen_source_archive.py", "prepare_frozen_build_overlay.py")]
    inputs = unique_records([*prerequisites, protocol, *sources, *base_records, *modules,
                             *[{k: row[k] for k in ("path", "bytes", "sha256")} for row in wheels]])
    check(inputs)
    output.mkdir(parents=True)
    result = dict(status="helper_environment_started", checked_records=inputs, selected=selected,
                  source=record(__file__), stages=[], refinement_attempts=0,
                  seed_admitted=False, accuracy_evaluated=False, publication_ready=False)
    try:
        wheelhouse = output / "wheels"
        wheelhouse.mkdir()
        for row in wheels:
            target = wheelhouse / Path(row["path"]).name
            shutil.copyfile(row["path"], target)
            if any(record(target)[key] != row[key] for key in ("bytes", "sha256")):
                raise ValueError("Copied helper wheel differs")
        lock_text, copied = hash_lock(wheelhouse, selected)
        lock = output / "requirements.txt"
        lock.write_text(lock_text)
        result.update(lock=record(lock), wheels=copied)
        env = dict(PATH="/usr/bin:/bin", HOME=str(output / "home"), LANG="C.UTF-8",
                   PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0", PYTHONMALLOC="debug",
                   PYTHONFAULTHANDLER="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                   MKL_NUM_THREADS="1")
        result["environment"] = env
        result["stage_limits"] = dict(cpu_seconds=300, wall_seconds=180,
                                      address_space_bytes=64 * 1024**3, affinity_cpus=1)
        Path(env["HOME"]).mkdir()

        def stage(name, command):
            code, timed_out = execute(command, output, env, output / (name + ".stdout"),
                                      output / (name + ".stderr"), 180)
            result["stages"].append(dict(name=name, command=command, returncode=code,
                                          timed_out=timed_out,
                                          stdout=record(output / (name + ".stdout")),
                                          stderr=record(output / (name + ".stderr"))))
            if code or timed_out:
                raise RuntimeError("Helper preparation stage failed: " + name)

        prefix = output / "venv"
        python = prefix / "bin/python"
        stage("venv", [str(BASE / "bin/python"), "-I", "-B", "-m", "venv",
                       "--without-pip", "--copies", str(prefix)])
        stage("install", [str(BASE / "bin/python"), "-I", "-B", "-m", "pip", "--isolated",
            "--python", str(python), "install", "--no-index", "--no-deps", "--no-compile",
            "--only-binary=:all:", "--require-hashes", "--find-links", str(wheelhouse), "-r", str(lock)])
        stage("dependencies", [str(python), "-I", "-B", "-m", "pip", "check"])
        launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
        probe = ("import sys,json,gc,importlib,importlib.metadata as m; sys.path.insert(0,sys.argv[1]); "
                 "[importlib.import_module(n) for n in json.loads(sys.argv[2])]; "
                 "print(json.dumps(dict(base_prefix=sys.base_prefix,gc_enabled=gc.isenabled(),"
                 "gc_thresholds=list(gc.get_threshold()),packages={d.metadata['Name']:d.version for d in m.distributions()},"
                 "modules={n:v.__file__ for n,v in sys.modules.items() if getattr(v,'__file__',None)}),sort_keys=True))")
        stage("imports", [str(python), "-I", "-B", "-c", probe, str(launcher), json.dumps(IMPORTS)])
        observed = json.loads((output / "imports.stdout").read_bytes())
        validate_imports(observed, selected, prefix, launcher, BASE)
        site = prefix / "lib/python3.10/site-packages"
        installed = [dict(name=row["name"], version=row["version"],
                          **installed_payload(Path(row["path"]), site)) for row in copied]
        import_records = [record(path) for path in observed["modules"].values()]
        payload_records = package_records(dict(packages=installed), site)
        cfg = record(prefix / "pyvenv.cfg")
        if "include-system-site-packages = false" not in (prefix / "pyvenv.cfg").read_text():
            raise ValueError("Helper environment permits system site packages")
        interpreter = record(python)
        if any(interpreter[key] != base_python[key] for key in ("bytes", "sha256")):
            raise ValueError("Copied helper interpreter differs")
        records = unique_records([*inputs, record(lock), cfg, interpreter, *import_records,
            *payload_records, *[{k: row[k] for k in ("path", "bytes", "sha256")} for row in copied]])
        check(records)
        result.update(status="helper_complete_private_environment_prepared", checked_records=records,
            prefix=str(prefix), base_prefix=str(BASE), base_python=base_python, interpreter=interpreter,
            import_report=record(output / "imports.stdout"), observed_imports=observed,
            package_audits=installed, limitations=[
                "Restored inference wheels plus retained reader Biopython 1.87; not the historical shared runtime.",
                "Observed worker-import closure, not branch-complete execution or ELF/OS dependency closure.",
                "Wheel comparison excludes generated RECORD/bytecode and relocated non-site payloads.",
                "No refinement, optimizer, score, controlled timing or scientific admission executed."])
    except BaseException as error:
        result.update(status="helper_environment_preparation_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "result.json").open("x") as handle:
            json.dump(result, handle, indent=2, sort_keys=True)
            handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output.absolute(), args.protocol_sha256)
