"""Install the frozen recovery lock offline into a fresh, recorded environment."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


LOCK_SHA = "813751262c9405155cc2f1974a2c32187b352f88307e2332e0cec9df35708ae3"


def record(path):
    path = Path(path).absolute()
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def commands(base_python, installer_python, wheels, lock, output):
    return [
        [str(base_python), "-I", "-m", "venv", "--without-pip", str(output / "venv")],
        [str(installer_python), "-I", "-m", "pip", "--isolated", "--disable-pip-version-check",
         "--python", str(output / "venv/bin/python"), "install", "--no-index", "--require-hashes",
         "--only-binary=:all:", "--no-cache-dir", "--find-links", str(wheels),
         "--report", str(output / "install_report.json"), "-r", str(lock)],
        [str(output / "venv/bin/python"), "-I", "-m", "pip", "check"],
    ]


def install(base_python, installer_python, wheels, lock, output):
    base_python, installer_python, wheels, lock, output = map(
        lambda p: Path(p).absolute(), (base_python, installer_python, wheels, lock, output))
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    records = [record(lock), *[record(p) for p in sorted(wheels.glob("*.whl"))]]
    if records[0]["sha256"] != LOCK_SHA or len(records) != 12:
        raise ValueError("Require the frozen recovery lock and eleven wheels")
    plan = commands(base_python, installer_python, wheels, lock, output)
    output.mkdir(parents=True)
    preparation = dict(commands=plan, inputs=records, source=record(__file__),
        interpreters=[record(base_python), record(installer_python)], publication_ready=False)
    (output / "preparation.json").write_text(json.dumps(preparation, indent=2) + "\n")
    outcomes = []
    try:
        for i, command in enumerate(plan):
            path = output / f"command_{i}.log"
            with path.open("x") as stream:
                result = subprocess.run(command, cwd=output, stdout=stream,
                    stderr=subprocess.STDOUT, timeout=600)
            outcomes.append(dict(command=command, returncode=result.returncode, log=record(path)))
            if result.returncode:
                raise RuntimeError(f"Installation command {i} failed; no retry")
        for item in records:
            if record(item["path"]) != item:
                raise ValueError("Installation input changed")
    except BaseException as error:
        (output / "failure.json").write_text(json.dumps(dict(error=str(error),
            type=type(error).__name__, outcomes=outcomes, retry=False), indent=2) + "\n")
        raise
    result = dict(status="offline_recovery_install_commands_complete", outcomes=outcomes,
        preparation=record(output / "preparation.json"),
        install_report=record(output / "install_report.json"),
        publication_ready=False, package_bytes_audited=False, native_inference_executed=False)
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("base-python", "installer-python", "wheels", "lock", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(install(args.base_python, args.installer_python,
                            args.wheels, args.lock, args.output), indent=2))
