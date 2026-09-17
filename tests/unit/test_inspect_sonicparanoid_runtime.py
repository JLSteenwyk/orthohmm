from pathlib import Path

import pytest

from benchmark_tools.inspect_sonicparanoid_runtime import probe, resolve_mcl


def test_nonconda_uses_package_mcl_even_when_path_differs():
    assert resolve_mcl(Path("/package"), {"is_conda": False, "is_mamba": False}, "/other") == Path("/package/bin/mcl")


@pytest.mark.parametrize("field", ["is_conda", "is_mamba"])
def test_environment_uses_path_mcl(tmp_path, field):
    binary = tmp_path / "mcl"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(0o755)
    system = {"is_conda": False, "is_mamba": False, field: True}
    assert resolve_mcl(Path("/package"), system, str(tmp_path)) == binary
    with pytest.raises(ValueError, match="Missing"):
        resolve_mcl(Path("/package"), system, str(tmp_path / "absent"))


def test_probe_records_version_and_rejects_failure(tmp_path):
    binary = tmp_path / "tool"
    binary.write_text("#!/bin/sh\nprintf 'test-version\\n'\n")
    binary.chmod(0o755)
    assert probe(binary, ["version"])["stdout"] == "test-version\n"
    binary.write_text("#!/bin/sh\nexit 2\n")
    with pytest.raises(ValueError, match="probe failed"):
        probe(binary, ["version"])
