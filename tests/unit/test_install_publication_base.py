import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import install_publication_base as module


@pytest.fixture
def arguments(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "supported_host", lambda: True)
    cache = tmp_path / "cache"
    cache.mkdir()
    packages = []
    for i in range(19):
        filename = f"package{i}-1.0-main.conda"
        path = cache / filename
        path.write_bytes(str(i).encode())
        hashes = module.archives.identity(path)
        packages.append(dict(package=dict(name=f"package{i}", fn=filename, version="1.0", build="main",
            subdir="linux-64", sha256=hashes["sha256"], md5=hashes["md5"],
            url="https://repo.anaconda.com/pkgs/main/linux-64/" + filename),
            archive=dict(bytes=hashes["bytes"], sha256=hashes["sha256"])))
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps(dict(acquisition=dict(packages=packages))))
    monkeypatch.setattr(module, "RECEIPT_SHA", module.record(receipt)["sha256"])
    wheel = tmp_path / module.PIP_NAME
    wheel.write_bytes(b"Synthetic wheel, never executed")
    monkeypatch.setattr(module, "PIP_SHA", module.record(wheel)["sha256"])
    monkeypatch.setattr(module, "PIP_BYTES", wheel.stat().st_size)
    conda = tmp_path / "conda"
    conda.write_bytes(b"Synthetic Conda identity, never executed")
    conda.chmod(0o755)
    args = SimpleNamespace(receipt=receipt, cache=cache, conda=conda, pip_wheel=wheel,
        conda_sha256=module.record(conda)["sha256"], output=tmp_path / "fresh with spaces",
        timeout=60, acknowledge_historical_runtime=True)
    return args


def install_metadata(prefix, packages):
    directory = prefix / "conda-meta"
    directory.mkdir(parents=True)
    for item in packages:
        package = item["package"]
        (directory / (package["name"] + ".json")).write_text(json.dumps(package))


def mock_execution(args, monkeypatch, failure=None):
    calls = []
    packages = module.archives.packages_from_receipt(json.loads(args.receipt.read_bytes()))
    prefix = args.output / "python-runtime"
    def execute(output, name, command, environment, timeout):
        calls.append((name, command, environment, timeout))
        if name == failure:
            raise RuntimeError("Injected stage failure")
        if name == "conda_install":
            install_metadata(prefix, packages)
            (prefix / "bin").mkdir()
            (prefix / "bin/python").write_text("Synthetic Python identity, never executed")
        if name == "pip_bootstrap":
            module.save(output / "bootstrap_install.json", dict(install=[dict(
                metadata=dict(name="pip", version="26.2.1"),
                download_info=dict(url=(output / "bootstrap_wheels" / module.PIP_NAME).as_uri(),
                                   archive_info=dict(hashes=dict(sha256=module.PIP_SHA))))]))
        if name == "runtime_snapshot":
            site = prefix / "lib/python3.10/site-packages"
            site.mkdir(parents=True)
            module.save(output / "runtime.json", dict(prefix=str(prefix), version="3.10.13",
                implementation="CPython", machine="x86_64", site=str(site),
                python=module.record(prefix / "bin/python"), distributions=[dict(name="pip", version="26.2.1")]))
        return dict(command=command, returncode=0)
    monkeypatch.setattr(module, "stage", execute)
    monkeypatch.setattr(module, "installed_payload", lambda *args: dict(matched_files=1, excluded=[]))
    return calls


def test_mocked_offline_composition(arguments, monkeypatch):
    calls = mock_execution(arguments, monkeypatch)
    result = module.run(arguments)
    assert [r[0] for r in calls] == ["conda_install", "pip_bootstrap", "pip_check", "runtime_snapshot"]
    assert len(result["installed_conda_packages"]) == 19
    assert result["historical_runtime"] is True
    assert result["native_inference_executed"] is False
    assert result["security_clearance"] is False
    assert result["publication_ready"] is False
    command = calls[0][1]
    assert "--offline" in command and "--copy" in command and "--no-default-packages" in command
    assert command[command.index("--prefix") + 1] == str(arguments.output / "python-runtime")
    bootstrap = calls[1][1]
    assert all(flag in bootstrap for flag in ("--no-index", "--require-hashes", "--no-deps", "--only-binary=:all:"))
    assert calls[0][2]["CONDA_PKGS_DIRS"] == str(arguments.output / "cache")
    assert calls[0][2]["HOME"] == str(arguments.output / "home")
    assert "PYTHONPATH" not in calls[0][2] and "LD_PRELOAD" not in calls[0][2]
    assert not (arguments.output / "failed.json").exists()


@pytest.mark.parametrize("stage", ["conda_install", "pip_bootstrap", "pip_check", "runtime_snapshot"])
def test_stage_failure_retained_without_retry(arguments, monkeypatch, stage):
    calls = mock_execution(arguments, monkeypatch, stage)
    with pytest.raises(RuntimeError, match="Injected"):
        module.run(arguments)
    assert [r[0] for r in calls][-1] == stage
    assert len({r[0] for r in calls}) == len(calls)
    failed = json.loads((arguments.output / "failed.json").read_text())
    assert failed["retry"] is False and failed["attempts"] == 1
    assert not (arguments.output / "complete.json").exists()


def test_unsupported_host_rejected(arguments, monkeypatch):
    monkeypatch.setattr(module, "supported_host", lambda: False)
    with pytest.raises(ValueError, match="Linux x86-64"):
        module.run(arguments)
    assert not arguments.output.exists()


@pytest.mark.parametrize("defect", ["package_inventory", "report_hash", "report_source", "pip_payload"])
def test_post_install_failure_not_admitted(arguments, monkeypatch, defect):
    calls = mock_execution(arguments, monkeypatch)
    original = module.stage
    def mutate(output, name, *args):
        result = original(output, name, *args)
        if defect == "package_inventory" and name == "conda_install":
            path = output / "python-runtime/conda-meta/extra.json"
            path.write_text(json.dumps(dict(name="extra", version="1", build="extra")))
        if defect.startswith("report") and name == "pip_bootstrap":
            path = output / "bootstrap_install.json"
            report = json.loads(path.read_text())
            download = report["install"][0]["download_info"]
            if defect == "report_hash":
                download["archive_info"]["hashes"]["sha256"] = "0" * 64
            else:
                download["url"] = "https://example.invalid/pip.whl"
            path.write_text(json.dumps(report))
        return result
    monkeypatch.setattr(module, "stage", mutate)
    if defect == "pip_payload":
        def reject(*args):
            raise ValueError("Installed payload differs")
        monkeypatch.setattr(module, "installed_payload", reject)
    with pytest.raises(ValueError):
        module.run(arguments)
    assert (arguments.output / "failed.json").exists()
    assert not (arguments.output / "complete.json").exists()
    assert len(calls) == len({r[0] for r in calls})


@pytest.mark.parametrize("defect", ["acknowledgement", "timeout", "boolean_timeout", "receipt", "conda_pin",
    "conda_mode", "wheel", "wheel_name", "cache_changed", "cache_missing", "cache_symlink", "output"])
def test_preflight_rejects_before_execution(arguments, monkeypatch, defect):
    args = arguments
    if defect == "acknowledgement":
        args.acknowledge_historical_runtime = False
    elif defect == "timeout":
        args.timeout = 0
    elif defect == "boolean_timeout":
        args.timeout = True
    elif defect == "receipt":
        args.receipt.write_text("{}")
    elif defect == "conda_pin":
        args.conda_sha256 = "0" * 64
    elif defect == "conda_mode":
        args.conda.chmod(0o644)
    elif defect == "wheel":
        args.pip_wheel.write_text("Changed")
    elif defect == "wheel_name":
        args.pip_wheel = args.pip_wheel.rename(args.pip_wheel.with_name("wrong.whl"))
    elif defect.startswith("cache"):
        path = next(args.cache.iterdir())
        if defect == "cache_changed":
            path.write_text("changed")
        else:
            path.unlink()
            if defect == "cache_symlink":
                path.symlink_to(args.receipt)
    else:
        args.output.mkdir()
    monkeypatch.setattr(module, "stage", lambda *args: pytest.fail("No native stage allowed"))
    with pytest.raises((ValueError, FileExistsError)):
        module.run(args)
    if defect != "output":
        assert not args.output.exists()


def test_inventory_rejects_extra_and_duplicate(arguments):
    packages = module.archives.packages_from_receipt(json.loads(arguments.receipt.read_bytes()))
    prefix = arguments.output / "python-runtime"
    install_metadata(prefix, packages)
    assert len(module.conda_inventory(prefix, packages)) == 19
    path = prefix / "conda-meta/extra.json"
    path.write_text(json.dumps(dict(name="extra", version="1", build="extra")))
    with pytest.raises(ValueError, match="inventory differs"):
        module.conda_inventory(prefix, packages)
    path.write_text((prefix / "conda-meta/package0.json").read_text())
    with pytest.raises(ValueError, match="Duplicate"):
        module.conda_inventory(prefix, packages)


@pytest.mark.parametrize("field,value", [("version", "3.12.3"), ("implementation", "PyPy"),
    ("machine", "aarch64"), ("prefix", "/wrong"), ("site", "/outside"),
    ("distributions", []), ("distributions", [dict(name="pip", version="24.0")])])
def test_wrong_runtime_rejected(arguments, monkeypatch, field, value):
    mock_execution(arguments, monkeypatch)
    original = module.validate_snapshot
    def corrupt(snapshot, prefix):
        snapshot[field] = value
        return original(snapshot, prefix)
    monkeypatch.setattr(module, "validate_snapshot", corrupt)
    with pytest.raises(ValueError, match="runtime differs"):
        module.run(arguments)
    assert (arguments.output / "failed.json").exists()
    assert not (arguments.output / "complete.json").exists()


def test_input_mutation_stops_before_next_stage(arguments, monkeypatch):
    calls = mock_execution(arguments, monkeypatch)
    original = module.stage
    def mutate(output, name, *args):
        result = original(output, name, *args)
        if name == "conda_install":
            (output / "staging/explicit.txt").write_text("changed")
        return result
    monkeypatch.setattr(module, "stage", mutate)
    with pytest.raises(ValueError, match="Changed pinned file"):
        module.run(arguments)
    assert len(calls) == 1
    assert not (arguments.output / "complete.json").exists()


@pytest.mark.parametrize("escape", [False, True])
def test_installed_python_alias_must_stay_inside_prefix(arguments, monkeypatch, escape):
    mock_execution(arguments, monkeypatch)
    original = module.stage
    def replace_alias(output, name, *args):
        result = original(output, name, *args)
        if name == "runtime_snapshot":
            python = output / "python-runtime/bin/python"
            target = output / "outside-python" if escape else python.with_name("python3.10")
            python.rename(target)
            python.symlink_to(target if escape else target.name)
        return result
    monkeypatch.setattr(module, "stage", replace_alias)
    if escape:
        with pytest.raises(ValueError, match="escapes the new prefix"):
            module.run(arguments)
        assert not (arguments.output / "complete.json").exists()
    else:
        assert module.run(arguments)["status"] == "historical_publication_base_installed"
