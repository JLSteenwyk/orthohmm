from pathlib import Path

import pytest

from benchmark_tools.run_ob_cpm_neighborhood import ARMS, replay_command


def baseline():
    return ["python", "replay.py", "--hits-pickle", "hits", "--fasta-directory", "fastas",
        "--cpm-resolution", "0.1", "--leiden-seed", "4", "--cpu", "32", "--matrix", "BLOSUM62",
        "--profile-iterations", "1", "--profile-min-species", "1", "--output-directory", "old",
        "--json", "old.json"]


@pytest.mark.parametrize("label,resolution", ARMS)
def test_only_resolution_and_destinations_change(label, resolution):
    original = baseline()
    expected = baseline()
    for flag, value in (("--cpm-resolution", resolution), ("--output-directory", "/out/replay"),
                        ("--json", "/out/replay.json")):
        expected[expected.index(flag) + 1] = value
    assert replay_command(original, label, Path("/out")) == expected
    assert original == baseline()


@pytest.mark.parametrize("flag", ["--cpm-resolution", "--leiden-seed", "--cpu", "--matrix",
    "--profile-iterations", "--profile-min-species"])
def test_fixed_parameters_cannot_change(flag):
    command = baseline()
    command[command.index(flag) + 1] = "changed"
    with pytest.raises(ValueError, match="parameters"):
        replay_command(command, "cpm_low", Path("/out"))


@pytest.mark.parametrize("change", ["arm", "scoring", "missing", "duplicate"])
def test_invalid_design_rejected(change):
    command = baseline()
    if change == "scoring":
        command += ["--official-benchmark", "benchmark.py"]
    elif change == "missing":
        command.remove("--json")
    elif change == "duplicate":
        command += ["--output-directory", "second"]
    with pytest.raises(ValueError):
        replay_command(command, "unknown" if change == "arm" else "cpm_low", Path("/out"))
