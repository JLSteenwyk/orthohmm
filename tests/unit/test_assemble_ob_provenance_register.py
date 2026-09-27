import pytest

from benchmark_tools.assemble_ob_provenance_register import workflow_duration,details


def test_repeated_consistent_duration():
    assert workflow_duration(["Duration : 50m 23s","Duration    : 50m 23s"]) == 3023


@pytest.mark.parametrize("lines", [[],["Duration : 50m 23s","Duration : 51m 23s"],
                                    ["Duration : 5m 90s"],["Duration : unknown"]])
def test_ambiguous_duration_rejected(lines):
    with pytest.raises(ValueError):
        workflow_duration(lines)


def test_sequence_only_never_inherits_full_or_conversion_duration():
    r = details("orthofinder_3_1_5_sequence_only",{"of":{
        "checks":{"processed_sequences_match":True},"command":["full"],
        "checkpoint_conversion_command":["convert"],"checkpoint_conversion_resources":{"elapsed_seconds":.42},
        "full_inference_resources":{"elapsed_seconds":3421.07}}})
    assert r["inference_wall_seconds"] is None
    assert r["conversion_wall_seconds"] == .42
    assert r["command"] == ["convert"]


def test_orthomcl_uses_july_not_april():
    def run(name,seconds):
        return dict(run=name,sequences_exact=False,assigned_genes=10,logged_duration={"seconds":seconds},commands={})
    r = details("orthomcl_1_4",{"orthomcl":{"runs":[run("Apr_20",1),run("Jul_25",2)]}})
    assert r["inference_wall_seconds"] == 2
    assert r["run"] == "Jul_25"


def test_unknown_is_not_zero_or_proxy():
    r = details("orthohmm_high_sensitivity",{})
    assert r["inference_wall_seconds"] is None and r["input_sequences_equal"] is None


def test_unknown_method_rejected():
    with pytest.raises(ValueError):
        details("another_tool",{})
