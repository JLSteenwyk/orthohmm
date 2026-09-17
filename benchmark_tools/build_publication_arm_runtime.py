"""Build unchanged frozen kernels on AArch64; not a numerical equivalence gate."""

import argparse
import ctypes
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.build_publication_runtime import COMMIT, KERNELS, record, verify_runtime
from benchmark_tools.validate_profile_runtime import require_profile_runtime

SYMBOLS = {
    "hmm_viterbi": ("batch_hmm_viterbi_c", "batch_hmm_viterbi_xdrop_c",
                    "hmm_set_num_threads", "hmm_have_avx2"),
    "kmer_prefilter": ("batch_prefilter_c", "prefilter_set_num_threads"),
    "pair_align": ("batch_pair_align_c", "pair_align_set_num_threads"),
}


def compile_command(compiler, csrc, name, machine):
    if machine != "aarch64" or name not in KERNELS:
        raise ValueError("Expected AArch64 host and frozen kernel")
    # Baseline ISA runs on both core types; do not assume one native core type.
    return [compiler, "-O3", "-fopenmp", "-shared", "-fPIC", "-march=armv8-a",
            "-o", str(csrc / (name + ".so")), str(csrc / (name + ".c"))]


def inspect_library(path, name):
    lib = ctypes.CDLL(str(path))
    for symbol in SYMBOLS[name]:
        getattr(lib, symbol)
    result = {"symbols": list(SYMBOLS[name])}
    if name == "hmm_viterbi":
        result["multipair_avx2_symbol"] = hasattr(lib, "batch_hmm_viterbi_multipair_avx2_c")
        lib.hmm_have_avx2.restype = ctypes.c_int32
        lib.hmm_have_avx2.argtypes = []
        result["hmm_have_avx2"] = lib.hmm_have_avx2()
        if result["hmm_have_avx2"] != 0:
            raise ValueError("Unexpected AVX2 runtime on ARM")
    return result


def build(root, output, python=sys.executable):
    root, output = Path(root).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    machine = platform.machine()
    if machine != "aarch64":
        raise ValueError("ARM build must execute on AArch64")
    revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    if revision != COMMIT:
        raise ValueError("Unexpected frozen source revision")
    subprocess.run(["git", "-C", str(root), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "setup.py"], check=True, capture_output=True)
    csrc = root / "orthohmm/search/csrc"
    if list(csrc.glob("*.so")):
        raise ValueError("Refusing to modify existing native runtime")
    compiler = shutil.which("gcc")
    if not compiler:
        raise ValueError("GCC required")
    compiler = str(Path(compiler).resolve())
    data = {"status": "building", "commit": revision, "root": str(root),
            "machine": machine, "platform": platform.platform(), "builder": record(__file__),
            "compiler": record(compiler),
            "compiler_version": subprocess.check_output([compiler, "--version"], text=True),
            "sources": [record(csrc / (name + ".c")) for name in KERNELS] + [record(root / "setup.py")],
            "binaries": [], "commands": [], "library_probes": {}, "cuda_enabled": False,
            "hmm_backend": "scalar_fallback_no_AVX2", "numerical_equivalence_verified": False,
            "started_at": datetime.now(timezone.utc).isoformat(),
            "limitations": ["Build and profile smoke only; numerical and end-to-end output checks remain required.",
                            "ARM timings cannot be pooled with x86 AVX2 results."]}
    with output.open("x") as handle:
        handle.write(json.dumps(data, indent=2, sort_keys=True) + "\n")
    try:
        for name in KERNELS:
            command = compile_command(compiler, csrc, name, machine)
            completed = subprocess.run(command, capture_output=True, text=True)
            data["commands"].append({"argv": command, "exit_code": completed.returncode,
                                     "stdout": completed.stdout, "stderr": completed.stderr})
            completed.check_returncode()
            binary = csrc / (name + ".so")
            data["binaries"].append(record(binary))
            data["library_probes"][name] = inspect_library(binary, name)
        data["profile_probe"] = require_profile_runtime(root, python)
        for source in data["sources"]:
            if record(source["path"]) != source:
                raise ValueError("Source changed during build")
        data["status"] = "complete"
    except Exception as error:
        data.update(status="failed", error=str(error))
        raise
    finally:
        data["finished_at"] = datetime.now(timezone.utc).isoformat()
        output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    verify_runtime(output, root)
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()
    build(args.root, args.output, args.python)
