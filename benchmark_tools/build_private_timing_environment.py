"""Stage a private timing candidate; never changes shared installed packages."""

import argparse
from email.parser import BytesParser
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import zipfile

PACKAGES = ("Cython", "pillow", "colorama", "coverage", "cycler", "python-dateutil",
    "defusedxml", "DendroPy", "igraph", "python-igraph", "leidenalg", "llvmlite",
    "matplotlib", "numba", "numpy", "packaging", "plotly", "pyparsing", "scipy",
    "six", "texttable", "wcwidth", "PyYAML", "kiwisolver", "contourpy", "fonttools",
    "narwhals", "psutil", "setuptools")


def record(path):
    return dict(path=str(path.resolve()), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def write(path, data):
    with path.open("x") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def selected_versions(baseline, patched_deployment=False):
    packages = baseline["environments"]["orthohmm"]["packages"]
    selected = {name: packages[name] for name in PACKAGES}
    if any(not isinstance(v, str) or not v or any(c.isspace() for c in v) for v in selected.values()):
        raise ValueError("Malformed frozen version")
    # Installer only, outside the frozen scientific package selection.
    selected["pip"] = "26.2.1"
    if patched_deployment:
        # Explicit deployment amendment; historical selections remain available.
        selected.update(packaging="26.1", pillow="12.3.0", setuptools="83.0.0")
    return selected


def wheel_metadata(path):
    with zipfile.ZipFile(path) as archive:
        names = [n for n in archive.namelist() if len(n.split("/")) == 2 and n.endswith(".dist-info/METADATA")]
        if len(names) != 1:
            raise ValueError("Ambiguous wheel metadata")
        meta = BytesParser().parsebytes(archive.read(names[0]))
    return meta["Name"], meta["Version"]


def canonical(name):
    return name.lower().replace("_", "-").replace(".", "-")


def hash_lock(wheels, selected):
    wanted = {canonical(k): (k, v) for k, v in selected.items()}
    found, records, lines = set(), [], []
    for path in sorted(wheels.glob("*.whl")):
        name, version = wheel_metadata(path)
        key = canonical(name)
        if key in found or key not in wanted or wanted[key][1] != version:
            raise ValueError("Unexpected or duplicate wheel")
        found.add(key)
        item = record(path)
        records.append(dict(name=name, version=version, **item))
        lines.append(f"{name}=={version} --hash=sha256:{item['sha256']}")
    if found != set(wanted):
        raise ValueError("Missing pinned wheels")
    return "\n".join(lines) + "\n", records


def build(baseline_path, expected_sha, output, patched_deployment=False):
    if output.exists():
        raise FileExistsError(output)
    raw = baseline_path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected_sha:
        raise ValueError("Baseline checksum differs")
    baseline = json.loads(raw)
    selected = selected_versions(baseline, patched_deployment)
    original = selected_versions(baseline)
    changes = {k: dict(previous=original[k], selected=v)
               for k, v in selected.items() if original[k] != v}
    output.mkdir(parents=True)
    wheels = output / "wheels"
    wheels.mkdir()
    requirements = output / "versions.txt"
    requirements.write_text("".join(f"{k}=={v}\n" for k, v in sorted(selected.items())))
    env = {"PATH": "/usr/bin:/bin", "HOME": str(output / "home"), "LANG": "C.UTF-8",
           "PYTHONDONTWRITEBYTECODE": "1", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
           "MKL_NUM_THREADS": "1", "PYTHONHASHSEED": "0"}
    Path(env["HOME"]).mkdir()
    stages = []
    def execute(label, command):
        write(output / f"{label}_started.json", dict(command=command, environment=env))
        started = time.monotonic()
        with (output / f"{label}.log").open("x") as log:
            result = subprocess.run(command, env=env, cwd=output, stdout=log,
                                    stderr=subprocess.STDOUT, timeout=1200)
        stage = dict(command=command, exit_code=result.returncode, wall_s=time.monotonic()-started,
                     log=record(output / f"{label}.log"))
        write(output / f"{label}_finished.json", stage)
        stages.append(stage)
        if result.returncode:
            raise RuntimeError(f"{label} failed; retain the attempt")
    pip = [sys.executable, "-I", "-B", "-m", "pip", "--isolated"]
    write(output / "started.json", dict(baseline=record(baseline_path), source=record(Path(__file__)),
          base_python=record(Path(sys.executable)), selected=selected, deployment_changes=changes,
          scientific_execution_authorized=False))
    execute("download", [*pip, "download", "--index-url", "https://pypi.org/simple", "--no-deps",
        "--only-binary=:all:", "--dest", str(wheels), "--requirement", str(requirements)])
    lock_text, wheel_records = hash_lock(wheels, selected)
    lock = output / "requirements.txt"
    lock.write_text(lock_text)
    prefix = output / "venv"
    execute("venv", [sys.executable, "-I", "-B", "-m", "venv", "--without-pip", "--copies", str(prefix)])
    python = prefix / "bin/python"
    execute("install", [*pip, "--python", str(python), "install", "--no-index", "--no-deps",
        "--only-binary=:all:", "--require-hashes", "--find-links", str(wheels), "-r", str(lock)])
    execute("dependencies", [str(python), "-I", "-B", "-m", "pip", "check"])
    probe = ("import sys,json,importlib.metadata as m; sys.path.insert(0,sys.argv[1]); "
             "import orthohmm.orthohmm,orthohmm.search.engine,igraph,leidenalg,dendropy; "
             "print(json.dumps(dict(executable=sys.executable,base_prefix=sys.base_prefix,path=sys.path,"
             "packages={d.metadata['Name']:d.version for d in m.distributions()},"
             "modules={n:getattr(v,'__file__',None) for n,v in sorted(sys.modules.items()) "
             "if getattr(v,'__file__',None)}),sort_keys=True))")
    execute("imports", [str(python), "-I", "-B", "-c", probe, baseline["core_root"]])
    observed = json.loads((output / "imports.log").read_text())
    if {canonical(k): v for k, v in observed["packages"].items()} != {canonical(k): v for k, v in selected.items()}:
        raise ValueError("Installed package inventory differs")
    leaks = [p for p in [*observed["path"], *observed["modules"].values()] if "/home/bizon/anaconda3" in p]
    if leaks:
        raise ValueError("Shared prefix in private Python lookup")
    for row in wheel_records:
        if record(Path(row["path"]))["sha256"] != row["sha256"]:
            raise ValueError("Wheel changed during installation")
    result = dict(status="private_timing_environment_candidate_installed", selected=selected,
        deployment_changes=changes,
        baseline=record(baseline_path), source=record(Path(__file__)), lock=record(lock), wheels=wheel_records,
        stages=stages, import_report=record(output / "imports.log"), base_prefix=observed["base_prefix"],
        shared_python_prefix_observed=False, scientific_execution_authorized=False,
        limitations=["Explicit deployment_changes supersede historical metadata; scientific package versions otherwise preserved.",
                     "Unrelated editable startup hooks omitted; pip installer is 26.2.1. No scientific inference or timing run.",
                     "Import-path check is not a file-access sandbox, native-library audit or complete branch coverage."])
    write(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--patched-deployment", action="store_true",
                        help="Apply the documented packaging alignment and Pillow/setuptools security revision")
    args = parser.parse_args()
    result = build(args.baseline.resolve(), args.sha256, args.output.resolve(), args.patched_deployment)
    print(json.dumps(dict(status=result["status"], packages=len(result["selected"]))))
