import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools.run_integrated_publication_workflow import assembly_base, install_commands, record, relocate_assets, stage, validate_assets, validate_base_probe, validate_data


def test_private_base_anchor_and_output_scope(tmp_path):
    from types import SimpleNamespace
    base = tmp_path / "base/bin/python"
    base.parent.mkdir(parents=True)
    base.write_bytes(b"pinned interpreter fixture")
    base.chmod(0o755)
    args = SimpleNamespace(base_python=base, installer_python=base, base_python_sha256=record(base)["sha256"],
                           assets=tmp_path / "bundle/assets")
    assert assembly_base(args, tmp_path / "output") == record(base)
    with pytest.raises(ValueError):
        assembly_base(args, tmp_path / "base/output")
    with pytest.raises(ValueError):
        assembly_base(args, tmp_path / "bundle/output")
    args.base_python_sha256 = "0" * 64
    with pytest.raises(ValueError):
        assembly_base(args, tmp_path / "output")


@pytest.mark.parametrize("defect", [None, "version", "distributions", "prefix", "site", "machine"])
def test_assembly_base_probe_pins_private_runtime(tmp_path, defect):
    base = tmp_path / "base/bin/python"
    value = dict(version="3.10.13", implementation="CPython", machine="x86_64", prefix=str(base.parent.parent),
                 site=str(base.parent.parent / "lib/python3.10/site-packages"), distributions=[["pip", "26.2.1"]])
    if defect == "version":
        value["version"] = "3.10.14"
    elif defect == "distributions":
        value["distributions"].append(["setuptools", "83.0.0"])
    elif defect in {"prefix", "site"}:
        value[defect] = str(tmp_path / "unrelated")
    elif defect == "machine":
        value["machine"] = "aarch64"
    if defect is None:
        validate_base_probe(value, base)
    else:
        with pytest.raises(ValueError):
            validate_base_probe(value, base)


def fixture_manifest(tmp_path):
    def item(name):
        path = tmp_path / name
        path.write_text(name)
        return record(path)
    return dict(dataset="installation_fixture", genes=16,
                fasta=[item(f"S{i}.fa") for i in range(4)],
                references=[item(f"R{i}.txt") for i in range(3)], uncertain=[])


@pytest.mark.parametrize("mode,change_base,change_host", [
    ("run", False, False), ("run", True, False), ("run", False, True),
    ("preflight", False, False), ("reject_host", False, False)])
def test_assembled_execution_checks_base_without_changing_legacy_commands(tmp_path, monkeypatch, mode, change_base, change_host):
    from types import SimpleNamespace
    from benchmark_tools import run_integrated_publication_workflow as module
    bundle = tmp_path / "bundle"
    assets = bundle / "assets"
    locks = {"inference": assets / "benchmark_tools/results/publication_recovery_requirements_20260926.txt",
             "reader": bundle / "reader_requirements.txt"}
    for name, path in locks.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
        monkeypatch.setitem(module.LOCKS, name, record(path)["sha256"])
    for path in (assets / "wheels", assets / "mafft", bundle / "reader_wheels", bundle / "readers"):
        path.mkdir(parents=True)
    (assets / "FastTree").write_text("fake tool")
    base = tmp_path / "base/bin/python"
    base.parent.mkdir(parents=True)
    base.write_text("fake interpreter")
    base.chmod(0o755)
    data = tmp_path / "manifest.json"
    data.write_text(json.dumps(fixture_manifest(tmp_path)))
    abi = tmp_path / "abi.json"
    abi.write_text("synthetic ABI fixture")
    args = SimpleNamespace(assets=assets, readers=bundle / "readers", reader_wheels=bundle / "reader_wheels",
        reader_lock=locks["reader"], base_python=base, installer_python=base,
        base_python_sha256=record(base)["sha256"], assembly_manifest_sha256="assembly-anchor",
        abi_inventory=abi, abi_inventory_sha256=record(abi)["sha256"],
        preflight_only=mode == "preflight",
        data=data, data_sha256=record(data)["sha256"], output=tmp_path / "run", cpu=2, timeout=10)
    validations, commands = [], []
    monkeypatch.setattr(module, "validate_assets", lambda *values: validations.append(values))
    host_checks = []
    def fake_host(values):
        host_checks.append(values)
        if mode == "reject_host":
            raise ValueError("Host glibc is below declared floor")
        return dict(inventory=record(abi), source=record(module.__file__), loaders=[],
                    glibc_floor_check_passed=True, full_host_compatibility_verified=False,
                    controller_glibc="glibc 2.39" if not change_host or len(host_checks) == 1 else "glibc 2.40")
    monkeypatch.setattr(module, "validate_assembly_host", fake_host)
    def fake_stage(directory, name, command, environment, timeout):
        commands.append((name, command))
        if name.startswith("base_"):
            value = dict(version="3.10.13", implementation="CPython", machine="x86_64",
                prefix=str(base.parent.parent), site=str(base.parent.parent / "lib/python3.10/site-packages"),
                distributions=[["pip", "26.2.1"]], files=[dict(path="pip.py", sha256="original")])
            if change_base and name == "base_after":
                value["files"][0]["sha256"] = "changed"
            (directory / (name + ".log")).write_text(json.dumps(value))
        elif name == "readback_score":
            (directory / "score.json").write_text("{}")
        return dict(returncode=0, name=name)
    monkeypatch.setattr(module, "stage", fake_stage)
    if mode == "reject_host":
        with pytest.raises(ValueError, match="below declared floor"):
            module.run(args)
        assert commands == [] and not args.output.exists()
        return
    if mode == "preflight":
        result = module.run(args)
        assert result["status"] == "integrated_assembled_preflight_complete"
        assert len(host_checks) == 2 and commands == []
        assert not (args.output / "complete.json").exists()
        assert not (args.output / "inference").exists()
        assert (args.output / "preflight.json").exists()
        assert all(result[key] is False for key in ("installation_executed", "scientific_inference_executed",
            "private_base_runtime_probe_executed", "reviewed_native_code_executed", "inference_execution_permitted",
            "full_host_compatibility_verified", "controlled_timing", "publication_ready"))
        return
    if change_base or change_host:
        with pytest.raises(ValueError, match="site changed" if change_base else "preflight changed"):
            module.run(args)
        assert not (args.output / "complete.json").exists()
        assert json.loads((args.output / "failure.json").read_text())["retry"] is False
    else:
        result = module.run(args)
        assert result["base_unchanged"] is True and result["base_site_payload_files"] == 1
        assert result["glibc_preflight"] == record(args.output / "glibc_preflight.json")
        assert len(host_checks) == 2
        assert len(validations) == 2
    assert len(commands) == 10
    assert commands[0][0] == "base_before" and commands[-1][0] == "base_after"
    assert all(command[1] == "-B" for name, command in commands if "install" in name)
    assert "-B" not in install_commands(base, base, assets / "wheels", locks["inference"], tmp_path / "legacy")[0]
    with pytest.raises(FileExistsError):
        module.run(args)


@pytest.mark.parametrize("assembly,inventory,digest", [("anchor", None, None), ("anchor", "file", None),
    ("anchor", None, "digest"), (None, "file", "digest")])
def test_invalid_abi_arguments_fail_before_output_and_install(tmp_path, monkeypatch, assembly, inventory, digest):
    from types import SimpleNamespace
    from benchmark_tools import run_integrated_publication_workflow as module
    data = tmp_path / "data.json"
    data.write_text(json.dumps(fixture_manifest(tmp_path)))
    args = SimpleNamespace(output=tmp_path / "run", cpu=2, timeout=10, data=data,
        data_sha256=record(data)["sha256"], assembly_manifest_sha256=assembly,
        abi_inventory=inventory, abi_inventory_sha256=digest)
    monkeypatch.setattr(module, "validate_assets", lambda *a: pytest.fail("must reject before asset work"))
    monkeypatch.setattr(module, "stage", lambda *a: pytest.fail("must reject before installation"))
    with pytest.raises(ValueError, match="ABI inventory"):
        module.run(args)
    assert not args.output.exists()


def test_fixture_scope(tmp_path):
    validate_data(fixture_manifest(tmp_path))


def test_fixture_cannot_be_claimed_as_full_benchmark(tmp_path):
    data = fixture_manifest(tmp_path)
    data["dataset"] = "orthobench"
    with pytest.raises(ValueError, match="scope or dimensions"):
        validate_data(data)


def test_changed_input_rejected(tmp_path):
    data = fixture_manifest(tmp_path)
    Path(data["fasta"][0]["path"]).write_text("changed")
    with pytest.raises(ValueError, match="Changed pinned"):
        validate_data(data)


def test_duplicate_input_rejected(tmp_path):
    data = fixture_manifest(tmp_path)
    data["fasta"][1] = data["fasta"][0]
    with pytest.raises(ValueError, match="Duplicate input"):
        validate_data(data)


def test_installs_are_offline_and_hash_required(tmp_path):
    commands = install_commands(Path("/base"), Path("/installer"), Path("/wheels"), Path("/lock"), tmp_path / "env")
    assert "--without-pip" in commands[0]
    assert all(x in commands[1] for x in ("--no-index", "--require-hashes", "--only-binary=:all:"))
    assert commands[2][-2:] == ["pip", "check"]


def test_completed_stage_records_success(tmp_path):
    result = stage(tmp_path, "success", [sys.executable, "-I", "-c", "print('complete')"], dict(os.environ), 10)
    assert result["returncode"] == 0
    assert (tmp_path / "success.log").read_text().strip() == "complete"


def test_failed_stage_does_not_retry(tmp_path):
    with pytest.raises(RuntimeError, match="without retry"):
        stage(tmp_path, "failure", [sys.executable, "-I", "-c", "raise SystemExit(9)"], dict(os.environ), 10)
    assert json.loads((tmp_path / "failure_finished.json").read_text())["returncode"] == 9
    with pytest.raises(FileExistsError):
        stage(tmp_path, "failure", [sys.executable, "-I", "-c", "pass"], dict(os.environ), 10)


def test_timeout_is_preserved(tmp_path):
    with pytest.raises(subprocess.TimeoutExpired):
        stage(tmp_path, "timeout", [sys.executable, "-I", "-c", "import time; time.sleep(10)"], dict(os.environ), 0.05)
    assert json.loads((tmp_path / "timeout_failed.json").read_text()) == dict(timeout=True, retry=False)


def test_changed_assets_manifest_rejected(tmp_path):
    (tmp_path / "copied_assets.json").write_text("[]")
    with pytest.raises(ValueError, match="native-asset manifest"):
        validate_assets(tmp_path, tmp_path)


def test_asset_relocation_never_overwrites(tmp_path):
    with pytest.raises(FileExistsError):
        relocate_assets(tmp_path, tmp_path)


def test_absolute_mafft_links_become_relative_only_in_new_copy(tmp_path, monkeypatch):
    from benchmark_tools import run_integrated_publication_workflow as module
    source, output = tmp_path / "original", tmp_path / "relocated"
    for name in ("wheels", "benchmark_tools", "mafft/bin", "mafft/libexec/mafft"):
        (source / name).mkdir(parents=True)
    for name in ("FastTree", "copied_assets.json", "manifest.json"):
        (source / name).write_text("fixture")
    for name in ("mafft-distance", "mafft-profile"):
        target = source / "mafft/libexec/mafft" / name
        target.write_text(name)
        (source / "mafft/bin" / name).symlink_to(target)
    original_record = module.record
    def fake_record(path):
        value = original_record(path)
        if Path(path) == source / "copied_assets.json":
            value["sha256"] = "8f8be7f1609d549da79a4ec4231e937a08833a7b8acc415d03e66b16274d5067"
        return value
    monkeypatch.setattr(module, "record", fake_record)
    changes = relocate_assets(source, output)
    assert len(changes) == 2
    for name in ("mafft-distance", "mafft-profile"):
        assert os.path.isabs(os.readlink(source / "mafft/bin" / name))
        assert os.readlink(output / "mafft/bin" / name) == "../libexec/mafft/" + name
        assert (output / "mafft/bin" / name).read_text() == name
