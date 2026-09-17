"""Read-only SonicParanoid dependency probe; never call installation helpers."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.snapshot_runtime_trees import inventory
from benchmark_tools.snapshot_orthohmm_input_order import record

FLAGS = {"blastp": ["-version"], "makeblastdb": ["-version"],
         "diamond": ["version"], "mmseqs": ["version"], "mcl": ["--version"]}


def resolve_mcl(root, system, path):
    if system["is_conda"] or system["is_mamba"]:
        result = shutil.which("mcl", path=path)
        if result is None:
            raise ValueError("Missing PATH-selected MCL")
        return Path(result)
    return root / "bin/mcl"


def probe(path, flags):
    before = record(path)
    done = subprocess.run([str(path), *flags], capture_output=True, text=True, timeout=30)
    if done.returncode or record(path) != before:
        raise ValueError("Dependency probe failed or binary changed: " + str(path))
    return {"file": before, "argv": [str(path), *flags], "exit_code": done.returncode,
            "stdout": done.stdout, "stderr": done.stderr}


def inspect():
    from sonicparanoid import sonic_paranoid, sys_tools, workers
    root = Path(sonic_paranoid.__file__).resolve().parent
    before = inventory([root])
    system = sys_tools.get_sys_info()
    paths = {"blastp": Path(workers.get_blastp_path()),
             "makeblastdb": Path(workers.get_makeblastdb_path()),
             "diamond": Path(workers.get_dmnd_path()), "mmseqs": Path(workers.get_mmseqs_path()),
             "mcl": resolve_mcl(root, system, os.environ["PATH"])}
    mode = sonic_paranoid.map_mode2sensitivity(
        argparse.Namespace(mode="default", diamond="", mmseqs=0.0, blast=False))
    probes = {name: probe(path, FLAGS[name]) for name, path in paths.items()}
    if inventory([root]) != before:
        raise ValueError("SonicParanoid package changed during inspection")
    return {"status": "read_only_current_resolution_observed", "execution_authorized": False,
            "source": record(__file__), "version": importlib.metadata.version("sonicparanoid"),
            "python": record(sys.executable), "system": system, "path": os.environ["PATH"],
            "default_mode": list(mode), "tools": probes,
            "declared_dependency_versions": sys_tools.get_binaries_info(),
            "package_inventory": before,
            "limitations": ["Current resolver inspection, not historical child execution proof.",
                            "No installer, download or inference entrypoint is called.",
                            "Interpreter/transitive/system libraries require a separate inventory.",
                            "Default DIAMOND mode does not imply absence of later MMseqs profile searches."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = inspect()
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
