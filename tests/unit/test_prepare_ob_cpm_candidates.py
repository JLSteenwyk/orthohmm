import pytest

from benchmark_tools.prepare_ob_cpm_candidates import LABELS, seed_arms


def report():
    return {"status": "cpm_replay_panel_verified_unscored", "accuracy_evaluated": False,
        "arms": [{"label": label, "stages": {"strict_profiles_refined": {"output": {"path": label + "_own_seed"}}}}
                 for label in LABELS]}


def test_each_variant_uses_own_refined_hmm_seed():
    assert seed_arms(report()) == [(label, {"path": label + "_own_seed"}) for label in LABELS]


@pytest.mark.parametrize("problem", ["status", "scored", "missing", "extra", "order"])
def test_unadmitted_or_incomplete_panel_rejected(problem):
    value = report()
    if problem == "status":
        value["status"] = "running"
    elif problem == "scored":
        value["accuracy_evaluated"] = True
    elif problem == "missing":
        value["arms"].pop()
    elif problem == "extra":
        value["arms"].append(value["arms"][0])
    else:
        value["arms"].reverse()
    with pytest.raises(ValueError):
        seed_arms(value)
