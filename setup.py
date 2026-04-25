"""Setup for OrthoHMM.

The package includes several performance-critical C/OpenMP kernels under
``orthohmm/search/csrc/`` and an optional CUDA kernel for NVIDIA GPUs.
``pip install`` compiles the CPU kernels automatically with ``gcc``; the
CUDA kernel is compiled with ``nvcc`` when available and silently skipped
otherwise (the Python runtime falls back to the CPU multipair / scalar
kernels seamlessly).

All kernels are ctypes-loaded shared libraries (``.so``), NOT Python
extension modules. A custom ``build_py`` subclass handles the builds.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path
from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop

from orthohmm.version import __version__

HERE = Path(__file__).resolve().parent
CSRC = HERE / "orthohmm" / "search" / "csrc"

# ─── Kernel build configuration ────────────────────────────────────────

# (source filename, needs_avx2, needs_math_lib)
# Production kernels only. Experimental kernels under
# orthohmm/search/_experimental/csrc/ are intentionally not built —
# they are retained as reference implementations for the paper's
# methods section and are not loaded at runtime.
CPU_KERNELS = [
    ("hmm_viterbi.c",     True,  False),  # scalar + multipair AVX2 + xdrop
    ("kmer_prefilter.c",  False, False),  # CSR k-mer prefilter
    ("pair_align.c",      False, False),  # NW for center-star MSA
]
CUDA_KERNELS = [
    ("hmm_viterbi_cuda.cu",),
]


def _run(cmd, cwd=None):
    try:
        subprocess.check_call(cmd, cwd=cwd)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"[orthohmm] kernel build step failed: {cmd[0]} -> {e}",
              file=sys.stderr)
        return False


def _have_gcc():
    return shutil.which("gcc") is not None or shutil.which("cc") is not None


def _have_nvcc():
    return shutil.which("nvcc") is not None


def _gcc_supports(flag, cc="gcc"):
    """Probe whether the compiler accepts a flag."""
    import tempfile
    try:
        with tempfile.TemporaryDirectory() as td:
            src = Path(td) / "probe.c"
            src.write_text("int main(void) { return 0; }\n")
            out = Path(td) / "probe"
            return subprocess.call(
                [cc, flag, "-o", str(out), str(src)],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            ) == 0
    except (OSError, FileNotFoundError):
        return False


def build_cpu_kernels(csrc_dir: Path) -> list[str]:
    """Compile CPU kernels to .so. Returns list of successfully-built names."""
    if not _have_gcc():
        print("[orthohmm] no C compiler found — skipping CPU kernel build. "
              "The package will fall back to Numba (slower).", file=sys.stderr)
        return []

    cc = "gcc" if shutil.which("gcc") else "cc"
    base = [cc, "-O3", "-fopenmp", "-shared", "-fPIC"]
    # -march=native and -mavx2 aren't universally available; probe.
    march_native = _gcc_supports("-march=native", cc=cc)
    mavx2 = _gcc_supports("-mavx2", cc=cc)
    if march_native:
        base.append("-march=native")

    built = []
    for fname, needs_avx2, needs_math in CPU_KERNELS:
        if needs_avx2 and not mavx2:
            print(f"[orthohmm] skipping {fname} (no AVX2)", file=sys.stderr)
            continue
        src = csrc_dir / fname
        dst = csrc_dir / fname.replace(".c", ".so")
        cmd = list(base)
        if needs_avx2:
            cmd.append("-mavx2")
        cmd += ["-o", str(dst), str(src)]
        if needs_math:
            cmd.append("-lm")
        if _run(cmd):
            built.append(fname.replace(".c", ".so"))
            print(f"[orthohmm] built {dst.name}")
    return built


def build_cuda_kernels(csrc_dir: Path) -> list[str]:
    """Compile CUDA kernels via nvcc. Silently skip if nvcc missing."""
    if not _have_nvcc():
        print("[orthohmm] no nvcc — skipping CUDA kernel build. "
              "Runtime will use CPU kernels.", file=sys.stderr)
        return []
    built = []
    for (fname,) in CUDA_KERNELS:
        src = csrc_dir / fname
        dst = csrc_dir / fname.replace(".cu", ".so")
        # Target a broad arch list so the build works on most modern GPUs.
        # Ada (sm_89), Ampere (sm_80/86), Turing (sm_75), Volta (sm_70).
        cmd = [
            "nvcc", "-O3", "-Xcompiler", "-fPIC", "-shared",
            "-gencode", "arch=compute_70,code=sm_70",
            "-gencode", "arch=compute_80,code=sm_80",
            "-gencode", "arch=compute_86,code=sm_86",
            "-gencode", "arch=compute_89,code=sm_89",
            "-o", str(dst), str(src),
        ]
        if _run(cmd):
            built.append(fname.replace(".cu", ".so"))
            print(f"[orthohmm] built {dst.name}")
    return built


class BuildKernelsMixin:
    def _build_kernels(self):
        build_cpu_kernels(CSRC)
        build_cuda_kernels(CSRC)


class CustomBuildPy(build_py, BuildKernelsMixin):
    def run(self):
        self._build_kernels()
        super().run()
        # Ensure .so files land in the install tree too
        pkg_csrc_out = Path(self.build_lib) / "orthohmm" / "search" / "csrc"
        pkg_csrc_out.mkdir(parents=True, exist_ok=True)
        for so in CSRC.glob("*.so"):
            shutil.copy2(so, pkg_csrc_out / so.name)


class CustomDevelop(develop, BuildKernelsMixin):
    def run(self):
        self._build_kernels()
        super().run()


# ─── Package metadata ──────────────────────────────────────────────────

with open(HERE / "README.md", encoding="utf-8") as f:
    long_description = f.read()

CLASSIFIERS = [
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

REQUIRES = [
    "numpy>=1.23",
    "numba>=0.58",
    "python-igraph",
    "leidenalg",
]

setup(
    name="orthohmm",
    description="HMM-based orthogroups",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Jacob L. Steenwyk",
    author_email="jlsteenwyk@gmail.com",
    url="https://github.com/jlsteenwyk/orthohmm",
    packages=["orthohmm", "orthohmm.search"],
    package_data={
        "orthohmm.search": ["csrc/*.c", "csrc/*.cu", "csrc/*.so"],
    },
    classifiers=CLASSIFIERS,
    entry_points={"console_scripts": ["orthohmm = orthohmm.orthohmm:main"]},
    version=__version__,
    include_package_data=True,
    install_requires=REQUIRES,
    cmdclass={
        "build_py": CustomBuildPy,
        "develop": CustomDevelop,
    },
)
