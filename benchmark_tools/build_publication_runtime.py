"""Build and verify a CPU-only native runtime without changing frozen sources."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.validate_profile_runtime import require_profile_runtime

COMMIT = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
KERNELS = ("hmm_viterbi", "kmer_prefilter", "pair_align")


def record(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size}


def verify_runtime(manifest, root):
    root = Path(root).resolve()
    data = json.loads(Path(manifest).read_text())
    if data.get("status") != "complete" or data.get("commit") != COMMIT or data.get("root") != str(root):
        raise ValueError("Wrong or incomplete native runtime manifest")
    if data.get("profile_probe", {}).get("status") != "passed":
        raise ValueError("Missing successful native profile probe")
    csrc = root / "orthohmm/search/csrc"
    expected = {csrc / (name + ".so") for name in KERNELS}
    if set(csrc.glob("*.so")) != expected:
        raise ValueError("Unexpected native runtime library set")
    records = data["sources"] + data["binaries"]
    if {Path(item["path"]) for item in data["binaries"]} != expected:
        raise ValueError("Incomplete native binary manifest")
    expected_sources = {csrc / (name + ".c") for name in KERNELS} | {root / "setup.py"}
    if {Path(item["path"]) for item in data["sources"]} != expected_sources:
        raise ValueError("Incomplete native source manifest")
    for item in records:
        if record(item["path"]) != item:
            raise ValueError("Changed native runtime file: " + item["path"])
    return data


def build(root, output, python=sys.executable):
    root = Path(root).resolve()
    output = Path(output).resolve()
    if output.exists():
        raise ValueError("Refusing to overwrite build provenance")
    actual = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    if actual != COMMIT:
        raise ValueError("Unexpected frozen source revision")
    subprocess.run(["git", "-C", str(root), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "setup.py"], check=True, capture_output=True)
    csrc = root / "orthohmm/search/csrc"
    if list(csrc.glob("*.so")):
        raise ValueError("Refusing to modify an existing native runtime")
    compiler = shutil.which("gcc")
    if compiler is None:
        raise ValueError("GCC required for the frozen Linux CPU build")
    compiler = str(Path(compiler).resolve())
    data = {"schema_version": 1, "status": "building", "root": str(root), "commit": actual,
            "started_at": datetime.now(timezone.utc).isoformat(), "builder": record(__file__),
            "compiler": record(compiler),
            "compiler_version": subprocess.check_output([compiler, "--version"], text=True),
            "platform": platform.platform(), "machine": platform.machine(),
            "sources": [record(csrc / (name + ".c")) for name in KERNELS] + [record(root / "setup.py")],
            "commands": [], "binaries": [], "cuda_enabled": False,
            "limitations": ["CPU-only build; -march=native binaries are machine-specific.",
                            "Synthetic profile smoke is not historical partition equivalence."]}
    with output.open("x") as handle:
        handle.write(json.dumps(data, indent=2, sort_keys=True) + "\n")
    try:
        for name in KERNELS:
            command = [compiler, "-O3", "-fopenmp", "-shared", "-fPIC", "-march=native"]
            if name == "hmm_viterbi":
                command.append("-mavx2")
            command.extend(["-o", str(csrc / (name + ".so")), str(csrc / (name + ".c"))])
            completed = subprocess.run(command, capture_output=True, text=True)
            data["commands"].append({"argv": command, "exit_code": completed.returncode,
                                     "stdout": completed.stdout, "stderr": completed.stderr})
            completed.check_returncode()
            binary = record(csrc / (name + ".so"))
            data["binaries"].append(binary)
        data["profile_probe"] = require_profile_runtime(root, python)
        data["status"] = "complete"
    except Exception as error:
        data.update(status="failed", error=str(error))
        raise
    finally:
        data["finished_at"] = datetime.now(timezone.utc).isoformat()
        output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    verify_runtime(output, root)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    if args.verify_only:
        verify_runtime(args.output, args.root)
    else:
        build(args.root, args.output, args.python)


if __name__ == "__main__":
    main()
