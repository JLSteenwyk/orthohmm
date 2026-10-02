import hashlib
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest

from benchmark_tools import assemble_publication_runtime_assets as module
from benchmark_tools import run_integrated_publication_workflow as executor


def write(path, payload, mode=0o644):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    path.chmod(mode)
    return module.record(path)


@pytest.fixture
def assembly(tmp_path, monkeypatch):
    component, prepared, wheelhouse = [tmp_path / name for name in ("component", "prepared", "wheelhouse")]
    workflow = component / "workflow"
    names = set.union(*module.wheels.ROLES.values())
    inventory = dict(status="selected_wheel_elf_inventory", wheels=[])
    project = tmp_path / "project/orthohmm-1.0-py3-none-any.whl"
    for name in sorted(names):
        filename = name.replace("-", "_") + "-1.0-py3-none-any.whl"
        target = project if name == "orthohmm" else wheelhouse / filename
        row = write(target, ("Test wheel " + name).encode())
        inventory["wheels"].append(dict(name=name, version="1.0", wheel=row))
    inventory_path = workflow / "benchmark_tools/results/integrated_wheel_elf_20260927.json"
    row = write(inventory_path, json.dumps(inventory).encode())
    monkeypatch.setattr(module.wheels, "INVENTORY_SHA", row["sha256"])
    for role, filename in {"inference": "publication_recovery_requirements_20260926.txt",
                           "reader": "publication_reader_requirements_20260927_v2.txt"}.items():
        row = write(workflow / "benchmark_tools/results" / filename, ("Test lock " + role).encode())
        monkeypatch.setitem(module.LOCKS, role, row["sha256"])
    for name in module.HARNESS:
        row = write(workflow / "benchmark_tools" / name, ("# Frozen test harness " + name).encode())
        monkeypatch.setitem(module.HARNESS, name, row["sha256"])
    entrypoint = module.readers.ENTRYPOINT.replace(".", "/") + ".py"
    reader_data = b"pass\n"
    write(workflow / entrypoint, reader_data)
    reader_pins = {entrypoint: module.readers.identity(reader_data)}
    monkeypatch.setattr(module, "READER_PINS", reader_pins)
    reader_export = tmp_path / "reference-reader"
    module.readers.write_export({entrypoint: reader_data}, module.READER_REVISION, reader_export)
    monkeypatch.setattr(module, "READER_MANIFEST_SHA", module.record(reader_export / "manifest.json")["sha256"])
    science = "orthohmm/source.py"
    write(component / "scientific" / science, b"frozen source")
    write(component / "SOURCE_INDEX.json", json.dumps(dict(profile="native-build",
        files=[dict(path="scientific/" + science, git_path=science, git_blob="test-blob")])).encode())
    def verify(root, digest):
        if root != component or digest != "source-anchor":
            raise ValueError("Changed source anchor")
        return dict(scientific_revision=module.source.SCIENTIFIC_REVISION, workflow_revision="fixture-workflow")
    monkeypatch.setattr(module.source, "verify", verify)
    monkeypatch.setattr(module, "scientific_members", lambda *args: [dict(path=science)] * 33)
    helper_pins = {}
    for name in ("version", "mafft-distance", "mafft-profile"):
        row = write(prepared / "tools/mafft/libexec/mafft" / name, b"\x7fELF" + name.encode(), 0o755)
        helper_pins[name] = {key: row[key] for key in ("bytes", "sha256")}
    monkeypatch.setattr(module.tools, "HELPERS", helper_pins)
    write(prepared / "tools/mafft/bin/mafft", b"#!/bin/sh\necho test launcher\n", 0o755)
    for name in ("mafft-distance", "mafft-profile"):
        (prepared / "tools/mafft/bin" / name).symlink_to("../libexec/mafft/" + name)
    for name in ("license", "license.extensions", "README.md"):
        write(prepared / "tools/notices/mafft" / name, ("MAFFT notice " + name).encode())
    artifacts = module.tools.acquisition.artifacts()
    for item in artifacts[1:]:
        name = Path(item["relative"]).name
        row = write(prepared / "tools/fasttree" / name, ("FastTree artifact " + name).encode(),
                    0o755 if name == "FastTree" else 0o644)
        item.update({key: row[key] for key in ("bytes", "sha256")})
    monkeypatch.setattr(module.tools.acquisition, "artifacts", lambda: artifacts)
    stages = [dict(returncode=0, log=write(prepared / (str(i) + ".log"), b"completed")) for i in range(5)]
    report = dict(status="private_frozen_phylogeny_tools_prepared", old_installation_required=False,
                  source_unchanged=True, scientific_inference_executed=False, stages=stages,
                  inventory=module.tools.inventory(prepared / "tools"),
                  helper_files=module.tools.inspect_helpers(prepared / "tools/mafft"))
    row = write(prepared / "complete.json", json.dumps(report).encode())
    args = SimpleNamespace(component=component, manifest_sha256="source-anchor", prepared_tools=prepared,
                           prepared_tools_sha256=row["sha256"], wheels=wheelhouse,
                           project_wheel=project, output=tmp_path / "assembled")
    return args


def test_complete_assembly_and_frozen_reader_manifest(assembly):
    result = module.build(assembly)
    root = assembly.output / "bundle"
    digest = module.record(root / module.INDEX)["sha256"]
    assert result == json.loads((assembly.output / "complete.json").read_bytes())
    assert module.validate(root, digest) == result["assembly"]
    assert len(result["scientific_members"]) == 33
    assert len(list((root / "assets/wheels").iterdir())) == 11
    assert len(list((root / "reader_wheels").iterdir())) == 5
    assert (root / "assets/notices/mafft/license.extensions").is_file()
    assert (root / "assets/notices/fasttree/LICENSE").is_file()
    assert module.record(root / "readers/manifest.json")["sha256"] == module.READER_MANIFEST_SHA
    assert result["source_unchanged"] is True
    assert all(result[key] is False for key in ("installation_performed", "native_code_executed",
        "scientific_inference_executed", "controlled_timing", "publication_ready", "security_clearance",
        "redistribution_clearance", "retry"))


def test_relocated_verification_needs_no_original_inputs(assembly, tmp_path):
    module.build(assembly)
    root = assembly.output / "bundle"
    digest = module.record(root / module.INDEX)["sha256"]
    relocated = tmp_path / "relocated"
    shutil.copytree(root, relocated, symlinks=True)
    for name in ("component", "prepared", "wheelhouse", "project"):
        shutil.rmtree(tmp_path / name)
    assert module.validate(relocated, digest)["status"] == "publication_runtime_assets_verified"


@pytest.mark.parametrize("defect", ["existing", "canonical", "inside", "source", "tool_anchor", "tool_byte",
    "tool_mode", "tool_log", "project_name", "project_byte", "wheel", "harness", "reader", "missing_science"])
def test_preflight_rejects_without_creating_output(assembly, tmp_path, monkeypatch, defect):
    if defect == "existing":
        assembly.output.mkdir()
    elif defect == "canonical":
        alias = tmp_path / "alias"
        alias.symlink_to(tmp_path, target_is_directory=True)
        assembly.output = alias / "assembled"
    elif defect == "inside":
        assembly.output = assembly.prepared_tools / "assembled"
    elif defect == "source":
        assembly.manifest_sha256 = "wrong"
    elif defect == "tool_anchor":
        assembly.prepared_tools_sha256 = "0" * 64
    elif defect == "tool_log":
        (assembly.prepared_tools / "0.log").write_bytes(b"changed stage evidence")
    elif defect.startswith("tool_"):
        path = assembly.prepared_tools / "tools/mafft/libexec/mafft/version"
        if defect == "tool_byte":
            path.write_bytes(b"changed native bytes")
        else:
            path.chmod(0o644)
    elif defect == "project_name":
        renamed = assembly.project_wheel.with_name("wrong.whl")
        assembly.project_wheel.rename(renamed)
        assembly.project_wheel = renamed
    elif defect == "project_byte":
        assembly.project_wheel.write_bytes(b"changed wheel")
    elif defect == "wheel":
        next(assembly.wheels.iterdir()).write_bytes(b"changed wheel")
    elif defect == "harness":
        (assembly.component / "workflow/benchmark_tools" / next(iter(module.HARNESS))).write_bytes(b"changed")
    elif defect == "reader":
        path = assembly.component / "workflow" / next(iter(module.READER_PINS))
        path.write_bytes(b"raise RuntimeError('changed')")
    else:
        monkeypatch.setattr(module, "scientific_members", lambda *args: [])
    with pytest.raises((ValueError, FileExistsError)):
        module.build(assembly)
    if defect != "existing":
        assert not assembly.output.exists()


@pytest.mark.parametrize("defect", ["anchor", "payload", "extra", "mode", "link", "path", "schema",
                                    "reanchored_wheel", "reanchored_helper", "reanchored_reader"])
def test_verifier_fails_closed(assembly, defect):
    module.build(assembly)
    root = assembly.output / "bundle"
    index = root / module.INDEX
    digest = module.record(index)["sha256"]
    value = json.loads(index.read_bytes())
    if defect == "anchor":
        digest = "0" * 64
    elif defect == "payload":
        (root / "assets/FastTree").write_bytes(b"changed")
    elif defect == "extra":
        (root / "extra.py").write_bytes(b"pass")
    elif defect == "mode":
        (root / "assets/FastTree").chmod(0o644)
    elif defect == "link":
        path = root / "assets/mafft/bin/mafft-profile"
        path.unlink()
        path.symlink_to("/outside-prefix")
    elif defect == "path":
        value["entries"][0]["path"] = "../outside"
    elif defect == "schema":
        value["schema"] = "different"
    else:
        if defect == "reanchored_wheel":
            next((root / "assets/wheels").iterdir()).write_bytes(b"changed wheel")
        elif defect == "reanchored_helper":
            (root / "assets/mafft/libexec/mafft/version").write_bytes(b"changed helper")
        else:
            (root / "readers" / next(iter(module.READER_PINS))).write_bytes(b"changed reader")
        value["entries"] = module.entries(root)
    if defect in {"path", "schema"} or defect.startswith("reanchored_"):
        index.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
        digest = module.record(index)["sha256"]
    with pytest.raises((ValueError, FileNotFoundError)):
        module.validate(root, digest)


def test_failed_copy_retained_without_retry(assembly, monkeypatch):
    def fail(*args, **kwargs):
        raise OSError("Injected copy failure")
    monkeypatch.setattr(module, "copy_file", fail)
    with pytest.raises(OSError):
        module.build(assembly)
    failure = json.loads((assembly.output / "failed.json").read_bytes())
    assert failure["retry"] is False and failure["publication_ready"] is False
    assert not (assembly.output / "complete.json").exists()
    with pytest.raises(FileExistsError):
        module.build(assembly)


def test_executor_route_is_explicit_and_checks_new_paths(assembly, monkeypatch):
    module.build(assembly)
    root = assembly.output / "bundle"
    digest = module.record(root / module.INDEX)["sha256"]
    monkeypatch.setattr(module.tools, "cpu_compatible", lambda: dict(observed=True))
    result = executor.validate_assets(root / "assets", root / "readers", digest)
    assert result["status"] == "publication_runtime_assets_verified"
    with pytest.raises(FileNotFoundError):
        executor.validate_assets(root / "assets", root / "readers")
    with pytest.raises(ValueError):
        executor.validate_assets(root / "assets", root / "reader_wheels", digest)


def test_reader_writer_preserves_bound_content_and_canonical_modes(tmp_path):
    payload = {"benchmark_tools/audit_publication_pipeline.py": b"pass\n"}
    module.readers.write_export(payload, module.READER_REVISION, tmp_path / "reader")
    assert (tmp_path / "reader/manifest.json").stat().st_mode & 0o777 == 0o644
    assert (tmp_path / "reader/benchmark_tools/audit_publication_pipeline.py").stat().st_mode & 0o777 == 0o644
    with pytest.raises(ValueError):
        module.readers.write_export({"../outside.py": b"pass"}, module.READER_REVISION, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()
