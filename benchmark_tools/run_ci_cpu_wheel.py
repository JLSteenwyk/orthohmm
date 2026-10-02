"""Build and verify an isolated Linux CPU wheel, not benchmark equivalence."""

import argparse
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.verify_cpu_wheel_install import verify

SCOPES = ("setup.py", "README.md", "LICENSE.md", "orthohmm")
REQUIRED = {"setup.py", "README.md", "LICENSE.md", "orthohmm/__init__.py", "orthohmm/version.py"}


def stage(root, destination):
    commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    names = subprocess.check_output(
        ["git", "-C", str(root), "ls-files", "-z", "--", *SCOPES], text=True).split("\0")
    names = [name for name in names if name]
    if not REQUIRED <= set(names) or len(set(names)) != len(names):
        raise ValueError("Require complete, unique tracked package inputs")
    prepared = []
    for name in names:
        relative = PurePosixPath(name)
        if (relative.is_absolute() or ".." in relative.parts
                or not (name in SCOPES or name.startswith("orthohmm/"))
                or relative.suffix in {".so", ".dylib", ".dll", ".pyc"}):
            raise ValueError("Unsafe or inherited binary package input")
        source = root / name
        if source.resolve() != source or not source.is_file():
            raise ValueError("Require direct regular source files")
        content = source.read_bytes()
        committed = subprocess.check_output(["git", "-C", str(root), "show", commit + ":" + name])
        if content != committed:
            raise ValueError("Package source differs from the committed revision")
        prepared.append((name, content, record(source)))
    destination.mkdir(parents=True, exist_ok=False)
    files = []
    for name, content, original in prepared:
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
        staged = record(target)
        if any(staged[key] != original[key] for key in ("bytes", "sha256")):
            raise ValueError("Staged package bytes differ")
        files.append(dict(relative_path=name, original=original, staged=staged))
    return dict(commit=commit, files=files)


def execute(argv, log, *, cwd, env):
    with log.open("x") as stream:
        subprocess.run(argv, cwd=cwd, env=env, stdout=stream,
                       stderr=subprocess.STDOUT, check=True, timeout=900)


def run(root, output):
    if sys.platform != "linux":
        raise NotImplementedError("CPU wheel linkage verification requires Linux")
    root = root.resolve()
    if not output.is_absolute() or output.resolve() != output:
        raise ValueError("Require a direct absolute output directory")
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    output.mkdir(parents=True, exist_ok=False)
    result = dict(status="cpu_wheel_attempt_started", controlled_timing=False,
                  accuracy_evaluated=False, publication_ready=False, commands=[])
    try:
        result["driver"] = record(__file__)
        result["package_source"] = stage(root, output / "source")
        result["requirements"] = [record(root / name) for name in ("requirements.txt", "tests/requirements.txt")]
        env = {k: v for k, v in os.environ.items() if k not in
               {"PYTHONPATH", "PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"}}
        env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        build_env = dict(env, PATH="/usr/bin:/bin", ORTHOHMM_CPU_TARGET="baseline")
        result["build_policy"] = dict(PATH=build_env["PATH"], ORTHOHMM_CPU_TARGET="baseline",
                                     inherited_shared_libraries=False)
        python = output / "venv/bin/python"

        def command(argv, name, environment=env):
            log = output / (name + ".log")
            result["commands"].append(dict(argv=[str(arg) for arg in argv], log=str(log)))
            execute([str(arg) for arg in argv], log, cwd=output, env=environment)
            result["commands"][-1]["log_record"] = record(log)

        command([sys.executable, "-m", "venv", output / "venv"], "create_environment")
        command([python, "-m", "pip", "install", "-r", root / "requirements.txt",
                 "-r", root / "tests/requirements.txt"], "install_dependencies")
        wheels = output / "wheels"
        command([python, "-m", "pip", "wheel", "--no-deps", "--no-build-isolation",
                 "--wheel-dir", wheels, output / "source"], "build_wheel", build_env)
        candidates = list(wheels.glob("*.whl"))
        if len(candidates) != 1:
            raise ValueError("Require exactly one freshly built project wheel")
        result["wheel"] = record(candidates[0])
        command([python, "-m", "pip", "install", "--no-deps", candidates[0]], "install_wheel")
        command([python, "-m", "pip", "check"], "check_dependencies")
        result["verification"] = verify(root, python, candidates[0], output / "verification")
        for row in result["package_source"]["files"]:
            check(row["original"])
            check(row["staged"])
        for pin in [result["driver"], result["wheel"], *result["requirements"]]:
            check(pin)
        result.update(status="linux_cpu_wheel_installation_verified", installed_runtime_verified=True,
            limitations=["Current committed development package, not a frozen benchmark executor.",
                         "Fixture runs do not establish biological accuracy or controlled efficiency.",
                         "No phylogenetic pipeline, cross-architecture, hermetic runtime or rights clearance.",
                         "Dependency versions and logs are retained, not a complete transitive hash lock."])
    except BaseException as error:
        result.update(status="cpu_wheel_attempt_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root, args.output)
