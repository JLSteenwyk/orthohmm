"""Build a distinct frozen-source CPU wheel offline, without historical admission."""

import argparse
from email.parser import BytesParser
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import zipfile

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools import bundle_publication_source as source
from benchmark_tools.audit_frozen_overlay_install import scientific_members, local_install_wheels
from benchmark_tools.audit_recovery_install import installed_payload
from benchmark_tools.run_integrated_publication_workflow import record, check, save, stage

BUILD_WHEELS = {
    "pip-26.2.1-py3-none-any.whl": dict(bytes=1816632,
        sha256="71138adf1f4ca900cdb7d289c21b7494329f2332b6d85f0e1c42108c0384ed3e"),
    "setuptools-83.0.0-py3-none-any.whl": dict(bytes=1008090,
        sha256="29b23c360f22f414dc7336bb39178cc7bcbf6021ed2733cde173f09dba19abb3"),
}
WHEEL_NAME = "orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl"
KERNELS = {
    "hmm_viterbi.so": ["hmm_set_num_threads", "hmm_have_avx2", "batch_hmm_viterbi_c"],
    "kmer_prefilter.so": ["prefilter_set_num_threads", "batch_prefilter_c"],
    "pair_align.so": ["pair_align_set_num_threads", "batch_pair_align_c"],
}
RUNTIME_PROBE = """import hashlib,json,platform,sys,sysconfig
from pathlib import Path
from importlib import metadata
site=Path(sysconfig.get_paths()['purelib'])
files=[]
for path in sorted(site.rglob('*')):
    if path.is_symlink(): raise ValueError('Symlinked base-site payload')
    if path.is_file():
        data=path.read_bytes()
        files.append(dict(path=path.relative_to(site).as_posix(),bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
print(json.dumps(dict(version=platform.python_version(),implementation=platform.python_implementation(),
    machine=platform.machine(),prefix=sys.prefix,site=str(site),files=files,
    distributions=sorted((d.metadata['Name'].lower(),d.version) for d in metadata.distributions(path=[str(site)])))))
"""
BUILD_PROBE = ("import json,sysconfig; from importlib import metadata; "
    "site=sysconfig.get_paths()['purelib']; print(json.dumps(dict(site=site,distributions="
    "sorted((d.metadata['Name'].lower(),d.version) for d in metadata.distributions(path=[site])))))")
KERNEL_PROBE = """import ctypes,json,sys
from pathlib import Path
names=json.loads(sys.argv[2])
results=[]
for name,symbols in sorted(names.items()):
    lib=ctypes.CDLL(str(Path(sys.argv[1])/name))
    for symbol in symbols:
        getattr(lib,symbol)
    setter=getattr(lib,symbols[0]); setter.argtypes=[ctypes.c_int32]; setter.restype=None; setter(1)
    if name=='hmm_viterbi.so':
        lib.hmm_have_avx2.argtypes=[]; lib.hmm_have_avx2.restype=ctypes.c_int32
        if lib.hmm_have_avx2()!=0: raise ValueError('Expected baseline, not AVX2 build')
    results.append(dict(kernel=name,symbols=symbols))
print(json.dumps(dict(status='baseline_cpu_kernels_load_verified',kernels=results)))
"""


def preflight(args):
    output = args.output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if (output.resolve() != output or type(args.timeout) is not int or args.timeout < 1
            or args.acknowledge_historical_runtime is not True):
        raise ValueError("Require canonical fresh output, positive timeout and historical acknowledgement")
    if platform.system() != "Linux" or platform.machine() != "x86_64":
        raise ValueError("Frozen CPU candidate requires Linux x86-64")
    component = args.component.resolve(strict=True)
    if output.is_relative_to(component):
        raise ValueError("Build output must be outside the immutable source component")
    verified = source.verify(component, args.manifest_sha256)
    index = json.loads((component / "SOURCE_INDEX.json").read_bytes())
    if index.get("profile") != "native-build" or verified.get("build_revision") != source.BUILD_REVISION:
        raise ValueError("Require the explicit native-build source profile")
    python = args.base_python.absolute()
    if python.parent.name == "bin" and output.is_relative_to(python.parent.parent):
        raise ValueError("Build output must be outside the supplied base prefix")
    if not python.is_file() or not os.access(python, os.X_OK) or record(python)["sha256"] != args.base_python_sha256:
        raise ValueError("Supplied base interpreter identity differs")
    wheels = [args.pip_wheel.absolute(), args.setuptools_wheel.absolute()]
    for wheel, name in zip(wheels, BUILD_WHEELS):
        ref = record(wheel) if wheel.is_file() else None
        if (wheel.is_symlink() or not wheel.is_file() or wheel.name != name
                or wheel.stat().st_size != BUILD_WHEELS[name]["bytes"]
                or {k: ref[k] for k in ("bytes", "sha256")} != BUILD_WHEELS[name]):
            raise ValueError("Require exact regular supplied build wheels")
    compiler_name = shutil.which("gcc", path="/usr/bin:/bin")
    if compiler_name is None:
        raise ValueError("Require the separately supplied system GCC toolchain")
    if shutil.which("nvcc", path="/usr/bin:/bin") is not None:
        raise ValueError("The CPU-only build PATH must not expose optional nvcc")
    compiler = Path(compiler_name).resolve(strict=True)
    watched = [record(p) for p in [python, compiler, *wheels, Path(__file__), Path(source.__file__),
        Path(sys.modules[stage.__module__].__file__), Path(sys.modules[scientific_members.__module__].__file__),
        Path(sys.modules[installed_payload.__module__].__file__), component / "SOURCE_INDEX.json"]]
    return output, component, verified, index, python, wheels, compiler, watched


def inspect_wheel(wheel, files, output):
    if wheel.name != WHEEL_NAME or wheel.stat().st_size > 10 * 1024 ** 2:
        raise ValueError("Unexpected rebuilt wheel identity/platform/size")
    members = scientific_members(wheel, files)
    if len(members) != 33:
        raise ValueError("Incomplete frozen scientific source membership")
    with zipfile.ZipFile(wheel) as archive:
        metadata_name = "orthohmm-0.5.0.dist-info/METADATA"
        message = BytesParser().parsebytes(archive.read(metadata_name))
        if (message.get_all("Name") != ["orthohmm"] or message.get_all("Version") != ["0.5.0"]
                or message.get_all("Requires-Python") != [">=3.10"]):
            raise ValueError("Unexpected rebuilt distribution metadata")
        wheel_meta = BytesParser().parsebytes(archive.read("orthohmm-0.5.0.dist-info/WHEEL"))
        if wheel_meta.get_all("Tag") != ["cp310-cp310-linux_x86_64"] or wheel_meta.get_all("Root-Is-Purelib") != ["false"]:
            raise ValueError("Unexpected rebuilt wheel tag/purity")
        payloads = []
        output.mkdir()
        for name in sorted(KERNELS):
            member = "orthohmm/search/csrc/" + name
            if archive.getinfo(member).file_size > 10 * 1024 ** 2:
                raise ValueError("Unexpected native member size")
            data = archive.read(member)
            if data[:4] != b"\x7fELF":
                raise ValueError("Required CPU kernel is not ELF")
            path = output / name
            with path.open("xb") as stream:
                stream.write(data)
            path.chmod(0o644)
            payloads.append(record(path))
    return dict(scientific_members=members, kernels=payloads)


def run(args):
    output, component, verified, index, python, wheels, compiler, watched = preflight(args)
    output.mkdir(parents=True)
    for name in ("home", "tmp", "build_wheels", "wheels"):
        (output / name).mkdir()
    environment = dict(HOME=str(output / "home"), PATH="/usr/bin:/bin", TMPDIR=str(output / "tmp"),
        LANG="C.UTF-8", PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1", PYTHONHASHSEED="0",
        PIP_CONFIG_FILE="/dev/null", ORTHOHMM_CPU_TARGET="baseline", SOURCE_DATE_EPOCH="1789568686",
        OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    save(output / "started.json", dict(inputs=watched, source=verified, environment=environment,
        attempts=1, historical_admission=False, publication_ready=False))
    outcomes = []
    try:
        staged, files = [], {}
        for row in index["files"]:
            if row["path"].startswith("scientific/"):
                original = component / row["path"]
                name = row["git_path"]
                data = original.read_bytes()
                files[name] = dict(content=data, git_blob=row["git_blob"])
                target = output / "source" / name
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open("xb") as stream:
                    stream.write((component / "build/setup.py").read_bytes() if name == "setup.py" else data)
                target.chmod(row["mode"])
                staged.append(record(target))
        for wheel in wheels:
            target = output / "build_wheels" / wheel.name
            shutil.copyfile(wheel, target)
            check(dict(record(wheel), path=str(target)))
            staged.append(record(target))
        requirements = output / "build_requirements.txt"
        with requirements.open("x") as stream:
            for name, pin in BUILD_WHEELS.items():
                project, version = name.split("-")[:2]
                stream.write(project + "==" + version + " --hash=sha256:" + pin["sha256"] + "\n")
        staged.append(record(requirements))
        venv = output / "venv"
        installed_python = str(venv / "bin/python")
        bootstrap = ("import sys,runpy;sys.path.insert(0," + repr(str(output / "build_wheels" / wheels[0].name))
                     + ");runpy.run_module('pip',run_name='__main__')")
        commands = [
            ("base_runtime", [str(python), "-I", "-B", "-c", RUNTIME_PROBE]),
            ("compiler_version", [str(compiler), "--version"]),
            ("create_build_environment", [str(python), "-I", "-B", "-m", "venv", "--without-pip", str(venv)]),
            ("install_build_dependencies", [installed_python, "-I", "-S", "-B", "-c", bootstrap,
                "--isolated", "--disable-pip-version-check", "install", "--prefix", str(venv),
                "--ignore-installed", "--no-index", "--no-deps",
                "--require-hashes", "--only-binary=:all:", "--no-cache-dir", "--find-links",
                str(output / "build_wheels"), "--report", str(output / "build_install.json"), "-r", str(requirements)]),
            ("build_dependency_check", [installed_python, "-I", "-B", "-m", "pip", "check"]),
            ("build_environment", [installed_python, "-I", "-B", "-c", BUILD_PROBE]),
            ("build_wheel", [installed_python, "-I", "-B", "-m", "pip", "--isolated",
                "--disable-pip-version-check", "wheel", "--no-index", "--no-deps", "--no-build-isolation",
                "--no-cache-dir", "--wheel-dir", str(output / "wheels"), str(output / "source")]),
        ]
        for name, command in commands:
            for ref in watched + staged:
                check(ref)
            outcomes.append(stage(output, name, command, environment, args.timeout))
            if name == "base_runtime":
                runtime = json.loads((output / "base_runtime.log").read_bytes())
                if (any(runtime[k] != v for k, v in dict(version="3.10.13", implementation="CPython", machine="x86_64").items())
                        or runtime["distributions"] != [["pip", "26.2.1"]]):
                    raise ValueError("Supplied base runtime is not the frozen Python version/platform")
        build_runtime = json.loads((output / "build_environment.log").read_bytes())
        site = Path(build_runtime["site"])
        if (build_runtime["distributions"] != [["pip", "26.2.1"], ["setuptools", "83.0.0"]]
                or site.resolve() != site or not site.is_relative_to(venv)):
            raise ValueError("Installed build dependency inventory/site differs")
        installed = local_install_wheels(json.loads((output / "build_install.json").read_bytes()), output / "build_wheels")
        if {(r["name"].lower(), r["version"]) for r in installed} != {("pip", "26.2.1"), ("setuptools", "83.0.0")}:
            raise ValueError("Installed build dependency report differs")
        payloads = [dict(wheel=record(w), matched_files=installed_payload(w, site)["matched_files"]) for w in wheels]
        candidates = list((output / "wheels").iterdir())
        if len(candidates) != 1 or candidates[0].is_symlink() or not candidates[0].is_file():
            raise ValueError("Require exactly one fresh project wheel")
        wheel = candidates[0]
        built_wheel = record(wheel)
        inspected = inspect_wheel(wheel, files, output / "kernels")
        outcomes.append(stage(output, "native_load", [installed_python, "-I", "-B", "-c", KERNEL_PROBE,
            str(output / "kernels"), json.dumps(KERNELS, sort_keys=True)], environment, args.timeout))
        loaded = json.loads((output / "native_load.log").read_bytes())
        if (loaded["status"] != "baseline_cpu_kernels_load_verified" or loaded["kernels"] !=
                [dict(kernel=name, symbols=symbols) for name, symbols in sorted(KERNELS.items())]):
            raise ValueError("Incomplete native load check")
        outcomes.append(stage(output, "base_unchanged", [str(python), "-I", "-B", "-c", RUNTIME_PROBE],
            environment, args.timeout))
        if json.loads((output / "base_unchanged.log").read_bytes()) != runtime:
            raise ValueError("Supplied base runtime/site payload changed during build")
        for ref in watched + staged + [built_wheel] + inspected["kernels"]:
            check(ref)
        if source.verify(component, args.manifest_sha256) != verified or shutil.which("gcc", path="/usr/bin:/bin") is None:
            raise ValueError("Source component or system compiler changed")
        if Path(shutil.which("gcc", path="/usr/bin:/bin")).resolve() != compiler:
            raise ValueError("System compiler alias changed")
        result = dict(status="frozen_source_cpu_wheel_candidate_built", source=verified, inputs=watched,
            staged_inputs=staged, outcomes=outcomes, wheel=record(wheel), environment=environment,
            base_unchanged=True, base_site_payload_files=len(runtime["files"]),
            build_runtime=build_runtime, build_payloads=payloads, inspection=inspected, native_load=loaded,
            attempts=1, retry=False, historical_wheel_reproduced=(built_wheel["bytes"] == 144444 and
                built_wheel["sha256"] == "cfdfde5ed1be29e4080dd3571f5c0fc5fc5f45b57ebe5ee9b3c7559096c3b93d"),
            historical_admission=False,
            scientific_inference_executed=False, controlled_timing=False, publication_ready=False,
            security_clearance=False, redistribution_clearance=False,
            limitations=["A new baseline CPU wheel, not a replacement in historical locks or scientific admissions.",
                "Frozen source parity and three-library load/ABI checks, not full numerical or inference equivalence.",
                "Supplied system GCC entrypoint is pinned, not its compiler/OS/library/source-rights closure.",
                "Same-host private offline build, not cross-host restoration or syscall sandbox."])
        save(output / "complete.json", result)
        return result
    except BaseException as error:
        save(output / "failed.json", dict(status="project_wheel_build_failed", type=type(error).__name__,
            error=str(error), outcomes=outcomes, attempts=1, retry=False, publication_ready=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("component", "base-python", "pip-wheel", "setuptools-wheel", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--base-python-sha256", required=True)
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--acknowledge-historical-runtime", action="store_true")
    args = parser.parse_args()
    result = run(args)
    print(json.dumps(dict(status=result["status"], complete=record(args.output / "complete.json")), indent=2))
