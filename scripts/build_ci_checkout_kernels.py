"""Prepare native kernels for checkout tests, not frozen benchmark execution."""

import argparse
import ctypes
import json
import os
from pathlib import Path
import subprocess
import sys

KERNELS = ("hmm_viterbi.so", "kmer_prefilter.so", "pair_align.so")


def build(root, output):
    root = root.resolve()
    source = root / "orthohmm/search/csrc"
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if list(source.glob("*.so")):
        raise ValueError("Refuse inherited checkout shared libraries")
    inputs = {name: (source / name.replace(".so", ".c")).read_bytes() for name in KERNELS}
    output.mkdir(parents=True, exist_ok=False)
    subprocess.run(
        [sys.executable, "setup.py", "build_py", "--build-lib", str(output)],
        cwd=root, env=dict(os.environ, ORTHOHMM_CPU_TARGET="baseline"),
        check=True, timeout=300,
    )
    compiled = output / "orthohmm/search/csrc"
    if {path.name for path in compiled.glob("*.so")} != set(KERNELS):
        raise ValueError("Require all three freshly built CPU kernels")
    # Validate the entire build before exposing any library to checkout imports.
    libraries = {name: ctypes.CDLL(str(compiled / name)) for name in KERNELS}
    libraries["hmm_viterbi.so"].hmm_have_avx2.restype = ctypes.c_int32
    if libraries["hmm_viterbi.so"].hmm_have_avx2() != 0:
        raise ValueError("Require scalar baseline Viterbi build")
    if any((source / name.replace(".so", ".c")).read_bytes() != content
           for name, content in inputs.items()):
        raise ValueError("Kernel source changed during build")
    for name in KERNELS:
        with (source / name).open("xb") as handle:
            handle.write((compiled / name).read_bytes())
    return dict(status="checkout_cpu_kernels_load_verified", cpu_target="baseline",
                kernels=list(KERNELS), controlled_timing=False,
                frozen_benchmark_runtime=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.root, args.output), sort_keys=True))
