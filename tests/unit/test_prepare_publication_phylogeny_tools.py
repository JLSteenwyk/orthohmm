import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import prepare_publication_phylogeny_tools as module


@pytest.fixture
def preparation(tmp_path, monkeypatch):
    component = tmp_path / "component"
    component.mkdir()
    (component / "SOURCE_INDEX.json").write_text("{}")
    def verify(path, digest):
        if path != component or digest != "anchored-index":
            raise ValueError("Changed source anchor")
        return dict(status="publication_source_components_verified", publication_ready=False)
    monkeypatch.setattr(module.source, "verify", verify)
    monkeypatch.setattr(module, "cpu_compatible", lambda: dict(system="Linux", machine="x86_64",
                                                               required_feature="avx2", observed=True))
    artifacts = tmp_path / "artifacts"
    rows = module.acquisition.artifacts()
    for row in rows:
        payload = ("Pinned test artifact: " + row["relative"]).encode()
        path = artifacts / row["relative"]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        row.update(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
    monkeypatch.setattr(module.acquisition, "artifacts", lambda: rows)
    def unpack(archive, output):
        assert archive == artifacts / rows[0]["relative"]
        source = output / module.acquisition.mafft.ROOT
        (source / "core").mkdir(parents=True)
        originals = []
        for name in ("core/Makefile", "core/input.c", "license", "license.extensions", "README.md"):
            path = source / name
            path.write_bytes(("Source " + name).encode())
            originals.append(module.record(path))
        return originals
    monkeypatch.setattr(module.acquisition.mafft, "unpack", unpack)
    monkeypatch.setattr(module, "SOURCE_FILES", 5)
    payloads = {name: (b"#!/usr/bin/perl\n" if name.endswith(".pl") else b"\x7fELF") + name.encode()
                for name in ("version", "mafft-distance", "mafft-profile", "script.pl")}
    monkeypatch.setattr(module, "HELPERS", {name: dict(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
                                           for name, payload in payloads.items()})
    calls = []
    def stage(output, name, command, environment, timeout):
        calls.append(dict(name=name, command=command, environment=environment, timeout=timeout))
        if name == "mafft_core_build":
            prefix = Path(next(argument.split("=", 1)[1] for argument in command if argument.startswith("PREFIX=")))
            helpers = prefix / "libexec/mafft"
            helpers.mkdir(parents=True)
            for helper, payload in payloads.items():
                path = helpers / helper
                path.write_bytes(payload)
                path.chmod(0o755)
            (prefix / "bin").mkdir()
            (prefix / "bin/mafft").write_bytes(b"#!/bin/sh\necho generated launcher\n")
            (prefix / "bin/mafft").chmod(0o755)
            for helper in ("mafft-distance", "mafft-profile"):
                (prefix / "bin" / helper).symlink_to(helpers / helper)
            (prefix / "bin/linsi").symlink_to("mafft")
        text = {"compiler_version": "gcc 13.3.0\n", "mafft_core_build": "build log\n",
                "mafft_version": "v7.525 (2024/Mar/13)\n", "mafft_helper_version": "7.525\n",
                "fasttree_help": "FastTree 2.2.0 Double precision:\n"}[name]
        log = output / (name + ".log")
        log.write_text(text)
        return dict(command=command, returncode=0, log=module.record(log))
    monkeypatch.setattr(module, "stage", stage)
    compiler = Path("/usr/bin/gcc").resolve()
    args = SimpleNamespace(component=component, manifest_sha256="anchored-index", artifacts=artifacts,
                           compiler_sha256=module.record(compiler)["sha256"], output=tmp_path / "prepared",
                           timeout=900, acknowledge_historical_runtime=True)
    return args, calls, rows, payloads, stage


def test_private_build_and_probes_without_old_installation(preparation, monkeypatch):
    args, calls, rows, _, _ = preparation
    monkeypatch.setenv("LD_PRELOAD", "/unrelated.so")
    monkeypatch.setenv("CC", "/unrelated-compiler")
    monkeypatch.setenv("MAFFT_BINARIES", "/unrelated-installation")
    result = module.run(args)
    assert result["status"] == "private_frozen_phylogeny_tools_prepared"
    assert result == json.loads((args.output / "complete.json").read_bytes())
    assert [row["name"] for row in calls] == ["compiler_version", "mafft_core_build", "mafft_version",
                                             "mafft_helper_version", "fasttree_help"]
    build = calls[1]["command"]
    assert "-j2" in build and "CFLAGS=-O3" in build and build[-1] == "install"
    assert build[build.index("-C") + 1].startswith(str(args.output / "source"))
    assert len(result["helper_files"]) == 4 and len(result["source_files"]) == 5
    assert result["helpers_historical_byte_equal"] is True and result["source_unchanged"] is True
    assert len(result["link_changes"]) == 2 and len(result["notices"]) == 3
    assert all(row["target"].startswith("../libexec/mafft/") for row in result["link_changes"])
    for call in calls:
        assert call["environment"]["PATH"] == "/usr/bin:/bin"
        assert not {"LD_PRELOAD", "CC", "PYTHONPATH"} & call["environment"].keys()
    assert calls[2]["environment"]["MAFFT_BINARIES"] == str(args.output / "tools/mafft/libexec/mafft")
    assert result["native_code_executed"] is True and result["private_tools_only"] is True
    assert all(result[key] is False for key in ("old_installation_required", "shared_environment_modified",
        "scientific_inference_executed", "historical_admission", "controlled_timing", "security_clearance",
        "redistribution_clearance", "publication_ready", "retry"))
    for row in rows:
        original = module.record(args.artifacts / row["relative"])
        assert all(original[key] == row[key] for key in ("bytes", "sha256"))
    assert (args.output / "tools/fasttree/FastTree").stat().st_mode & 0o777 == 0o755
    assert (args.artifacts / "fasttree/FastTree").stat().st_mode & 0o111 == 0


@pytest.mark.parametrize("defect", ["ack", "bool_timeout", "zero_timeout", "existing", "canonical",
    "spaces", "shell", "inside_component", "inside_artifacts", "artifacts_alias", "file_alias",
    "artifact_hash", "compiler", "anchor", "architecture"])
def test_preflight_refuses_before_creating_output(preparation, tmp_path, monkeypatch, defect):
    args, calls, rows, _, _ = preparation
    if defect == "ack":
        args.acknowledge_historical_runtime = False
    elif defect in {"bool_timeout", "zero_timeout"}:
        args.timeout = True if defect == "bool_timeout" else 0
    elif defect == "existing":
        args.output.mkdir()
    elif defect == "canonical":
        alias = tmp_path / "alias"
        alias.symlink_to(tmp_path, target_is_directory=True)
        args.output = alias / "prepared"
    elif defect in {"spaces", "shell"}:
        args.output = tmp_path / ("not quoted" if defect == "spaces" else "injected;command")
    elif defect == "inside_component":
        args.output = args.component / "prepared"
    elif defect == "inside_artifacts":
        args.output = args.artifacts / "prepared"
    elif defect == "artifacts_alias":
        alias = tmp_path / "artifacts-alias"
        alias.symlink_to(args.artifacts, target_is_directory=True)
        args.artifacts = alias
    elif defect == "file_alias":
        path = args.artifacts / rows[0]["relative"]
        saved = tmp_path / "saved"
        path.rename(saved)
        path.symlink_to(saved)
    elif defect == "artifact_hash":
        (args.artifacts / rows[0]["relative"]).write_bytes(b"changed artifact")
    elif defect == "compiler":
        args.compiler_sha256 = "0" * 64
    elif defect == "anchor":
        args.manifest_sha256 = "changed-index"
    else:
        def incompatible():
            raise ValueError("Unsupported host")
        monkeypatch.setattr(module, "cpu_compatible", incompatible)
    with pytest.raises((ValueError, FileExistsError)):
        module.run(args)
    assert not calls
    if defect != "existing":
        assert not args.output.exists()


@pytest.mark.parametrize("failure", ["compiler_version", "mafft_core_build", "mafft_version",
                                    "mafft_helper_version", "fasttree_help"])
def test_each_stage_failure_is_retained_without_retry(preparation, monkeypatch, failure):
    args, _, _, _, original = preparation
    def fail(output, name, *values):
        if name == failure:
            (output / (name + ".log")).write_text("Injected command failure")
            raise RuntimeError("Injected failure: " + name)
        return original(output, name, *values)
    monkeypatch.setattr(module, "stage", fail)
    with pytest.raises(RuntimeError, match="Injected failure"):
        module.run(args)
    report = json.loads((args.output / "failed.json").read_bytes())
    assert report["retry"] is False and report["attempts"] == 1
    assert report["publication_ready"] is False and report["scientific_inference_executed"] is False
    assert not (args.output / "complete.json").exists()
    assert (args.output / (failure + ".log")).is_file()
    with pytest.raises(FileExistsError):
        module.run(args)


@pytest.mark.parametrize("defect", ["helper_hash", "helper_mode", "helper_symlink", "extra_helper", "link",
                                     "mafft_version", "helper_version", "fasttree_version", "changed_source",
                                     "changed_input"])
def test_built_output_and_post_execution_guards(preparation, monkeypatch, defect):
    args, _, rows, _, original = preparation
    def change(output, name, *values):
        result = original(output, name, *values)
        if name == "mafft_core_build" and defect in {"helper_hash", "helper_mode", "helper_symlink", "extra_helper", "link"}:
            helpers = output / "tools/mafft/libexec/mafft"
            helper = helpers / "version"
            if defect == "helper_hash":
                helper.write_bytes(b"different native bytes")
            elif defect == "helper_mode":
                helper.chmod(0o644)
            elif defect == "helper_symlink":
                helper.unlink()
                helper.symlink_to(helpers / "mafft-profile")
            elif defect == "extra_helper":
                (helpers / "unexpected").write_bytes(b"x")
            else:
                link = output / "tools/mafft/bin/mafft-profile"
                link.unlink()
                link.symlink_to("/unrelated-prefix")
        bad_version = {"mafft_version": "mafft_version", "helper_version": "mafft_helper_version",
                       "fasttree_version": "fasttree_help"}
        if name == bad_version.get(defect):
            (output / (name + ".log")).write_text("Unexpected version\n")
        if name == "fasttree_help" and defect == "changed_source":
            (output / "source" / module.acquisition.mafft.ROOT / "core/input.c").write_bytes(b"changed source")
        if name == "fasttree_help" and defect == "changed_input":
            (args.artifacts / rows[0]["relative"]).write_bytes(b"changed input")
        return result
    monkeypatch.setattr(module, "stage", change)
    with pytest.raises(ValueError):
        module.run(args)
    assert (args.output / "failed.json").is_file() and not (args.output / "complete.json").exists()


def test_frozen_catalog_is_metadata_only():
    assert len(module.HELPERS) == 34
    assert len([name for name in module.HELPERS if name.endswith(".pl")]) == 2
    assert module.SOURCE_FILES == 173
    assert all(set(pin) == {"bytes", "sha256"} and len(pin["sha256"]) == 64 for pin in module.HELPERS.values())


@pytest.mark.parametrize("system,machine,flags,accepted", [
    ("Linux", "x86_64", "flags : sse avx2\n", True),
    ("Darwin", "x86_64", "flags : avx2\n", False),
    ("Linux", "aarch64", "flags : avx2\n", False),
    ("Linux", "x86_64", "flags : sse\n", False),
    ("Linux", "x86_64", "model : 1\n", False),
    ("Linux", "x86_64", "flags : avx2\nflags : sse\n", False),
])
def test_historical_fasttree_cpu_gate(monkeypatch, system, machine, flags, accepted):
    monkeypatch.setattr(module.platform, "system", lambda: system)
    monkeypatch.setattr(module.platform, "machine", lambda: machine)
    read = Path.read_text
    def cpu_text(path, *args, **kwargs):
        return flags if path == Path("/proc/cpuinfo") else read(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", cpu_text)
    if accepted:
        assert module.cpu_compatible()["observed"] is True
    else:
        with pytest.raises(ValueError):
            module.cpu_compatible()
