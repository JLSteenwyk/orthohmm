import json
from types import SimpleNamespace
import zipfile

import pytest

from benchmark_tools import build_publication_project_wheel as module


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    component = tmp_path / "component"
    (component / "scientific/orthohmm").mkdir(parents=True)
    (component / "build").mkdir()
    (component / "build/setup.py").write_text("BUILD = 'synthetic'\n")
    names = ["orthohmm/source" + str(i) + ".py" for i in range(33)] + ["setup.py"]
    rows = []
    for name in names:
        path = component / "scientific" / name
        path.write_text("pass\n")
        rows.append(dict(path="scientific/" + name, git_path=name, git_blob="a" * 40, mode=0o644))
    (component / "SOURCE_INDEX.json").write_text(json.dumps(dict(profile="native-build", files=rows)))
    verified = dict(status="publication_source_components_verified", build_revision=module.source.BUILD_REVISION)
    monkeypatch.setattr(module.source, "verify", lambda *args: verified)
    python, compiler = tmp_path / "python", tmp_path / "gcc"
    for path in (python, compiler):
        path.write_text("Synthetic executable identity, never executed")
        path.chmod(0o755)
    monkeypatch.setattr(module.shutil, "which", lambda name, path: str(compiler) if name == "gcc" else None)
    monkeypatch.setattr(module.platform, "system", lambda: "Linux")
    monkeypatch.setattr(module.platform, "machine", lambda: "x86_64")
    wheels, pins = [], {}
    for name in module.BUILD_WHEELS:
        wheel = tmp_path / name
        wheel.write_text("Synthetic build wheel, never installed")
        ref = module.record(wheel)
        pins[name] = {k: ref[k] for k in ("bytes", "sha256")}
        wheels.append(wheel)
    monkeypatch.setattr(module, "BUILD_WHEELS", pins)
    args = SimpleNamespace(component=component, manifest_sha256="a" * 64, base_python=python,
        base_python_sha256=module.record(python)["sha256"], pip_wheel=wheels[0], setuptools_wheel=wheels[1],
        output=tmp_path / "fresh with spaces", timeout=60, acknowledge_historical_runtime=True)
    return args


def write_wheel(args, fault=None):
    path = args.output / "wheels" / module.WHEEL_NAME
    rows = json.loads((args.component / "SOURCE_INDEX.json").read_bytes())["files"]
    with zipfile.ZipFile(path, "w") as archive:
        for row in rows:
            if row["git_path"].startswith("orthohmm/"):
                archive.writestr(row["git_path"], b"changed" if fault == "source" else b"pass\n")
        for name in module.KERNELS:
            if fault == "fallback" and name == "hmm_viterbi.so": continue
            archive.writestr("orthohmm/search/csrc/" + name, b"not ELF" if fault == "elf" else b"\x7fELFsynthetic")
        archive.writestr("orthohmm-0.5.0.dist-info/METADATA",
            "Name: orthohmm\nVersion: 0.5.0\nRequires-Python: >=3.10\n")
        archive.writestr("orthohmm-0.5.0.dist-info/WHEEL",
            "Tag: cp310-cp310-linux_x86_64\nRoot-Is-Purelib: false\n")


def execution(args, monkeypatch, fail=None, wheel_fault=None):
    calls = []
    monkeypatch.setattr(module, "installed_payload", lambda *args: dict(matched_files=1))
    def stage(output, name, command, env, timeout):
        calls.append((name, command, env, timeout))
        if name == fail: raise RuntimeError("Injected build stage failure")
        value = {}
        if name in {"base_runtime", "base_unchanged"}:
            value = dict(version="3.10.13", implementation="CPython", machine="x86_64", files=[],
                         distributions=[["pip", "26.2.1"]])
        elif name == "build_environment":
            site = output / "venv/lib/python3.10/site-packages"
            site.mkdir(parents=True)
            value = dict(site=str(site), distributions=[["pip", "26.2.1"], ["setuptools", "83.0.0"]])
        elif name == "install_build_dependencies":
            rows = []
            for project, version, wheel in (("pip", "26.2.1", args.pip_wheel), ("setuptools", "83.0.0", args.setuptools_wheel)):
                local = output / "build_wheels" / wheel.name
                rows.append(dict(metadata=dict(name=project, version=version), download_info=dict(url=local.as_uri(),
                    archive_info=dict(hashes=dict(sha256=module.record(local)["sha256"])))))
            (output / "build_install.json").write_text(json.dumps(dict(install=rows)))
        elif name == "build_wheel": write_wheel(args, wheel_fault)
        elif name == "native_load":
            value = dict(status="baseline_cpu_kernels_load_verified",
                kernels=[dict(kernel=k, symbols=v) for k, v in sorted(module.KERNELS.items())])
        (output / (name + ".log")).write_text(json.dumps(value))
        return dict(returncode=0, command=command)
    monkeypatch.setattr(module, "stage", stage)
    return calls


def test_mocked_offline_build_and_source_parity(inputs, monkeypatch):
    calls = execution(inputs, monkeypatch)
    result = module.run(inputs)
    assert result["status"] == "frozen_source_cpu_wheel_candidate_built"
    assert len(result["inspection"]["scientific_members"]) == 33 and len(result["inspection"]["kernels"]) == 3
    assert len(calls) == 9 and result["base_unchanged"] is True
    env = calls[0][2]
    assert env["PATH"] == "/usr/bin:/bin" and env["ORTHOHMM_CPU_TARGET"] == "baseline"
    assert env["HOME"] == str(inputs.output / "home") and env["TMPDIR"] == str(inputs.output / "tmp")
    assert "PYTHONPATH" not in env and "LD_LIBRARY_PATH" not in env
    assert (inputs.output / "source/setup.py").read_text() == "BUILD = 'synthetic'\n"
    assert (inputs.component / "scientific/setup.py").read_text() == "pass\n"
    install = next(c[1] for c in calls if c[0] == "install_build_dependencies")
    assert all(f in install for f in ("--no-index", "--no-deps", "--require-hashes", "--only-binary=:all:", "--ignore-installed"))
    assert install[install.index("--prefix") + 1] == str(inputs.output / "venv")
    build = next(c[1] for c in calls if c[0] == "build_wheel")
    assert all(f in build for f in ("--no-index", "--no-deps", "--no-build-isolation", "--no-cache-dir"))
    assert all(result[k] is False for k in ("historical_wheel_reproduced", "historical_admission", "retry",
        "scientific_inference_executed", "controlled_timing", "publication_ready", "security_clearance", "redistribution_clearance"))


@pytest.mark.parametrize("fault", ["ack", "timeout", "existing", "inside_component", "inside_base", "host", "profile", "python_hash", "python_execute", "wheel_hash", "wheel_symlink", "cuda", "compiler"])
def test_preflight_precedes_output_and_execution(inputs, monkeypatch, tmp_path, fault):
    args = inputs
    calls = execution(args, monkeypatch)
    if fault == "ack": args.acknowledge_historical_runtime = False
    elif fault == "timeout": args.timeout = True
    elif fault == "existing": args.output.mkdir()
    elif fault == "inside_component": args.output = args.component / "build-output"
    elif fault == "inside_base":
        prefix = tmp_path / "base"
        (prefix / "bin").mkdir(parents=True)
        python = prefix / "bin/python"
        python.write_bytes(args.base_python.read_bytes()); python.chmod(0o755)
        args.base_python = python; args.output = prefix / "build-output"
    elif fault == "host": monkeypatch.setattr(module.platform, "machine", lambda: "aarch64")
    elif fault == "profile":
        p = args.component / "SOURCE_INDEX.json"
        data = json.loads(p.read_bytes()); data["profile"] = "native-wheels"; p.write_text(json.dumps(data))
    elif fault == "python_hash": args.base_python_sha256 = "0" * 64
    elif fault == "python_execute": args.base_python.chmod(0o644)
    elif fault == "wheel_hash": args.pip_wheel.write_text("changed")
    elif fault == "wheel_symlink":
        link = tmp_path / "link"; link.symlink_to(args.pip_wheel); args.pip_wheel = link
    elif fault == "cuda": monkeypatch.setattr(module.shutil, "which", lambda name, path: str(args.base_python))
    else: monkeypatch.setattr(module.shutil, "which", lambda *args, **kwargs: None)
    with pytest.raises((ValueError, FileExistsError)): module.run(args)
    assert not calls
    if fault != "existing": assert not args.output.exists()


@pytest.mark.parametrize("name", ["base_runtime", "compiler_version", "create_build_environment", "install_build_dependencies",
    "build_dependency_check", "build_environment", "build_wheel", "native_load", "base_unchanged"])
def test_stage_failure_retained_without_retry(inputs, monkeypatch, name):
    calls = execution(inputs, monkeypatch, name)
    with pytest.raises(RuntimeError): module.run(inputs)
    assert calls[-1][0] == name and len(calls) == len({c[0] for c in calls})
    assert not (inputs.output / "complete.json").exists()


def test_all_stdlib_probes_compile():
    for name in ("RUNTIME_PROBE", "BUILD_PROBE", "KERNEL_PROBE"):
        compile(getattr(module, name), name, "exec")
    assert json.loads((inputs.output / "failed.json").read_bytes())["retry"] is False


@pytest.mark.parametrize("fault", ["fallback", "source", "elf"])
def test_rebuilt_wheel_rejects_fallback_or_changed_sources(inputs, monkeypatch, fault):
    calls = execution(inputs, monkeypatch, wheel_fault=fault)
    with pytest.raises(ValueError): module.run(inputs)
    assert calls[-1][0] == "build_wheel" and not (inputs.output / "complete.json").exists()


@pytest.mark.parametrize("fault", ["runtime", "base_dependencies", "base_changed", "dependencies", "site", "native_status", "native_symbols", "changed_kernel"])
def test_installed_and_native_validation_guards(inputs, monkeypatch, fault):
    execution(inputs, monkeypatch)
    original = module.stage
    def stage(output, name, command, env, timeout):
        result = original(output, name, command, env, timeout)
        path = output / (name + ".log")
        value = json.loads(path.read_bytes())
        if name == "base_runtime" and fault == "runtime": value["version"] = "0.0.0"
        if name == "base_runtime" and fault == "base_dependencies": value["distributions"].append(["setuptools", "83.0.0"])
        if name == "base_unchanged" and fault == "base_changed": value["files"].append(dict(path="new", bytes=1, sha256="a"*64))
        if name == "build_environment":
            if fault == "dependencies": value["distributions"].append(["unknown", "1"])
            elif fault == "site": value["site"] = str(output.parent)
        if name == "native_load":
            if fault == "native_status": value["status"] = "unknown"
            elif fault == "native_symbols": value["kernels"][0]["symbols"] = []
            elif fault == "changed_kernel": (output / "kernels/hmm_viterbi.so").write_text("changed")
        path.write_text(json.dumps(value))
        return result
    monkeypatch.setattr(module, "stage", stage)
    with pytest.raises(ValueError): module.run(inputs)
    assert not (inputs.output / "complete.json").exists()
