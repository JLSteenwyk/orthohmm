from pathlib import Path

import pytest

from benchmark_tools.prepare_qfo_corrected_primary import primary_commands


def baseline():
    return {"core_root": "/frozen", "tool_entrypoints": {
        "orthohmm_python": {"absolute_path": "/python"},
        "orthofinder": {"absolute_path": "/orthofinder"}}}


def test_native_commands_preserve_science_and_isolate_inputs():
    result = primary_commands(Path("/corrected"), Path("/new"), baseline())
    assert set(result) == {"orthohmm_high_sensitivity", "orthofinder_full"}
    hmm = result["orthohmm_high_sensitivity"]["native_argv"]
    assert hmm[:4] == ["/python", "-m", "orthohmm", "/corrected"]
    for flag, value in (("-c", "32"), ("--threads_per_worker", "4"), ("-x", "BLOSUM62"),
                        ("-e", "0.0001"), ("--clustering", "leiden"), ("--cpm_resolution", "0.1"),
                        ("--accuracy_profile", "high_sensitivity"), ("--stop", "infer")):
        assert hmm[hmm.index(flag) + 1] == value
    assert "--phylogeny" not in hmm and "--start" not in hmm
    of = result["orthofinder_full"]
    assert of["native_argv"] == ["/orthofinder", "-f", "/new/orthofinder_full/input", "-t", "32", "-a", "32", "-S", "diamond"]
    assert of["copy_inputs_from"] == "/corrected"
    assert all(r["search_reuse"] is False and r["accuracy_admitted"] is False for r in result.values())


@pytest.mark.parametrize("cpus", [0, -1, True, 1.5])
def test_invalid_resource_request(cpus):
    with pytest.raises(ValueError):
        primary_commands(Path("/corrected"), Path("/new"), baseline(), cpus)
