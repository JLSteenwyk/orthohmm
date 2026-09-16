from pathlib import Path

import pytest

from benchmark_tools.run_sequence_graph_control import check_profile_off, graph_command


def test_frozen_graph_command_has_no_profile_or_scoring_inputs():
    args = graph_command(Path("/launcher"), Path("/checkpoint"), "hash", Path("/output"))
    assert "--fasta-directory" not in args
    assert "--official-benchmark" not in args
    assert args[args.index("--leiden-seed") + 1] == "4"
    assert args[args.index("--cpm-resolution") + 1] == "0.1"
    assert args[args.index("--cpu") + 1] == "32"


@pytest.mark.parametrize("change", [None, "profile", "stages", "genes", "species"])
def test_profile_off_stage_admission(change):
    data = {"parameters": {"profile_expansion": False}, "counts": {"genes": 251378, "species": 12},
            "stages": [{"label": "multipass"}, {"label": "multipass_refined"}]}
    if change == "profile":
        data["parameters"]["profile_expansion"] = True
    elif change == "stages":
        data["stages"].append({"label": "strict_profiles"})
    elif change in {"genes", "species"}:
        data["counts"][change] -= 1
    if change is None:
        check_profile_off(data)
    else:
        with pytest.raises(ValueError):
            check_profile_off(data)
