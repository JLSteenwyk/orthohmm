"""Record OrthoFinder's own subprocess PATH resolution, not historical exec traces."""

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

FLAGS = {"diamond": ["version"], "FastTree": ["-help"], "mafft": ["--version"],
         "mcl": ["--version"], "famsa": ["-h"], "fastme": ["-V"]}


def file_record(path):
    path = Path(path).resolve()
    return {"path": str(path), "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def inspect_tools(path):
    result = {}
    for name, flags in FLAGS.items():
        executable = shutil.which(name, path=path)
        if executable is None:
            result[name] = {"status": "missing"}
            continue
        record = file_record(executable)
        env = os.environ.copy()
        env["PATH"] = path
        completed = subprocess.run([executable, *flags], env=env, text=True,
                                   capture_output=True, timeout=30)
        result[name] = {"status": "resolved", "file": record, "argv": [executable, *flags],
                        "exit_code": completed.returncode, "stdout": completed.stdout,
                        "stderr": completed.stderr}
        if file_record(executable) != record:
            raise ValueError("Executable changed during inspection")
    return result


def inspect(baseline=None):
    source = None
    if baseline:
        source = file_record(baseline)
        manifest = json.loads(Path(baseline).read_text())
        os.environ.update(manifest["environment_overrides"])
        os.environ["PATH"] = os.pathsep.join([*manifest["prepend_path"], os.environ.get("PATH", "")])
    outer = os.environ["PATH"]
    from orthofinder.utils import parallel_task_manager
    package_root = Path(parallel_task_manager.__file__).resolve().parents[1]
    return {"status": "current_resolution_observed", "historical_execution_proven": False,
            "source": file_record(__file__), "baseline": source, "machine": platform.machine(),
            "python": sys.version, "python_executable": file_record(sys.executable),
            "packages": {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()},
            "conda_prefix": os.environ.get("CONDA_PREFIX"),
            "outer_tools": inspect_tools(outer),
            "child_tools": inspect_tools(parallel_task_manager.my_env["PATH"]),
            "package_sources": [file_record(p) for p in sorted(package_root.rglob("*"))
                                if p.is_file() and p.suffix in {".py", ".json"}],
            "limitations": ["Current resolution/version probes, not historical execution tracing.",
                            "No inference, accuracy or timing admission."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = inspect(args.baseline)
    with args.output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
