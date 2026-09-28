from copy import deepcopy

import pytest

from benchmark_tools.prepare_threadripper_baseline import amend


def test_only_known_cwd_difference_is_accepted():
    baseline = dict(core_root="/frozen", environments={
        "orthohmm": {"python": "3.10", "packages": {"orthohmm": "0.5.0"}},
        "orthofinder": {"python": "3.12", "packages": {"orthohmm": "0.5.0", "orthofinder": "3.1.5"}}})
    observed = deepcopy(baseline["environments"])
    del observed["orthofinder"]["packages"]["orthohmm"]
    changed = amend(baseline, observed)
    assert changed["environments"] == observed
    assert baseline["environments"]["orthofinder"]["packages"]["orthohmm"] == "0.5.0"
    for name in observed:
        bad = deepcopy(observed)
        bad[name]["packages"]["unexpected"] = "1.0"
        with pytest.raises(ValueError, match="beyond"):
            amend(baseline, bad)
    with pytest.raises(ValueError, match="beyond"):
        amend(baseline, baseline["environments"])
