"""Install an offline private controller without changing scientific environments."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.build_private_timing_environment import canonical, hash_lock, record, write

SELECTED = dict(biopython="1.87", DendroPy="5.1.0", numpy="2.2.6", pip="26.2.1",
                setuptools="83.0.0", psutil="7.2.2")
MODULES = ["benchmark_tools.run_threadripper_scaling", "benchmark_tools.probe_verified_threadripper_fixture",
           "benchmark_tools.validate_threadripper_outputs", "benchmark_tools.replay_threadripper_scaling"]


def require_inventory(observed):
    if {canonical(k): v for k, v in observed.items()} != {canonical(k): v for k, v in SELECTED.items()}:
        raise ValueError("Controller distribution inventory differs")


def build(repo, output):
    if output.exists():
        raise FileExistsError(output)
    reader_lock = repo / "benchmark_tools/results/publication_reader_requirements_20260927_v2.txt"
    candidate_path = repo / "benchmark_tools/results/threadripper_patched_runtime_20260928.json"
    candidate = json.loads(candidate_path.read_text())
    psutil = [r for r in candidate["wheels"] if canonical(r["name"]) == "psutil"]
    if len(psutil) != 1 or psutil[0]["version"] != SELECTED["psutil"]:
        raise ValueError("Missing pinned psutil wheel")
    wheel = Path(psutil[0]["path"])
    if record(wheel)["sha256"] != psutil[0]["sha256"]:
        raise ValueError("psutil wheel changed")
    output.mkdir(parents=True)
    wheels = output / "wheels"
    wheels.mkdir()
    for source in [*Path("/tmp/orthohmm-reader-patched-20260927/wheels").glob("*.whl"), wheel]:
        shutil.copyfile(source, wheels / source.name)
    lock_text, records = hash_lock(wheels, SELECTED)
    expected = {line.strip() for line in reader_lock.read_text().splitlines() if line and not line.startswith("#")}
    expected.add(f"psutil==7.2.2 --hash=sha256:{psutil[0]['sha256']}")
    if {line.lower() for line in lock_text.splitlines()} != {line.lower() for line in expected}:
        raise ValueError("Controller wheel hashes differ from retained locks")
    lock = output / "requirements.txt"
    lock.write_text(lock_text)
    env = dict(PATH="/usr/bin:/bin", HOME=str(output / "home"), LANG="C.UTF-8",
               PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    Path(env["HOME"]).mkdir()
    result = dict(status="started", source=record(Path(__file__)), base_python=record(Path(sys.executable)),
                  reader_lock=record(reader_lock), candidate=record(candidate_path),
                  selected=SELECTED, wheels=records, lock=record(lock), stages=[],
                  scientific_execution_authorized=False)
    write(output / "started.json", result)
    def execute(name, command):
        write(output / f"{name}_started.json", dict(command=command, environment=env))
        with (output / f"{name}.log").open("x") as log:
            process = subprocess.run(command, env=env, cwd=output, stdout=log,
                                     stderr=subprocess.STDOUT, timeout=180)
        stage = dict(name=name, command=command, exit_code=process.returncode, log=record(output / f"{name}.log"))
        write(output / f"{name}_finished.json", stage)
        result["stages"].append(stage)
        if process.returncode:
            raise RuntimeError(name + " failed")
    try:
        python = output / "venv/bin/python"
        execute("venv", [sys.executable, "-I", "-B", "-m", "venv", "--without-pip", "--copies", str(python.parent.parent)])
        execute("install", [sys.executable, "-I", "-B", "-m", "pip", "--isolated", "--python", str(python),
                            "install", "--no-index", "--no-deps", "--only-binary=:all:", "--require-hashes",
                            "--find-links", str(wheels), "-r", str(lock)])
        execute("dependencies", [str(python), "-I", "-B", "-m", "pip", "check"])
        probe = ("import sys,json,importlib,importlib.metadata as m; sys.path.insert(0,sys.argv[1]); "
                 "[importlib.import_module(n) for n in json.loads(sys.argv[2])]; "
                 "print(json.dumps(dict(base_prefix=sys.base_prefix,packages={d.metadata['Name']:d.version for d in m.distributions()},"
                 "modules={n:v.__file__ for n,v in sys.modules.items() if getattr(v,'__file__',None)})))")
        execute("imports", [str(python), "-I", "-B", "-c", probe, str(repo), json.dumps(MODULES)])
        observed = json.loads((output / "imports.log").read_text())
        # The checkout adds source distribution metadata, not an installed wheel.
        packages = dict(observed["packages"])
        result["checkout_distribution"] = packages.pop("orthohmm", None)
        require_inventory(packages)
        if any("/home/bizon/anaconda3" in p for p in observed["modules"].values()) or "/home/bizon/anaconda3" in observed["base_prefix"]:
            raise ValueError("Shared Python prefix observed")
        result.update(status="private_controller_installed_imports_passed", import_report=record(output / "imports.log"),
                      base_prefix=observed["base_prefix"], interpreter=record(python),
                      limitations=["Declared controller imports only, not branch-complete execution or ELF closure.",
                                   "Controller reader dependencies are separate from native inference packages.",
                                   "Runtime binding and collector fixture validation remain required."])
    except Exception as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        write(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(build(args.repo.resolve(), args.output.resolve())["status"])
