from copy import deepcopy

import pytest

from benchmark_tools.admit_ob_cpm_candidates import check_arm


@pytest.mark.parametrize("problem", [None, "seed", "label", "disabled", "parameter", "profile"])
def test_fixed_expansion_uses_own_cpm_seed(problem):
    seed, parameters = {"path": "cpm_low_seed", "sha256": "a"}, {"min_norm": .03, "min_margin": 1.5}
    row = {"label": "cpm_low", "seed_partition": deepcopy(seed), "candidate_expansion": True,
           "expansion": {"parameters": deepcopy(parameters), "profile": "satellite_v2"}}
    if problem == "seed":
        row["seed_partition"]["path"] = "baseline_seed"
    elif problem == "label":
        row["label"] = "cpm_high"
    elif problem == "disabled":
        row["candidate_expansion"] = False
    elif problem == "parameter":
        row["expansion"]["parameters"]["min_norm"] = .024
    elif problem == "profile":
        row["expansion"]["profile"] = "other"
    if problem:
        with pytest.raises(ValueError):
            check_arm(row, "cpm_low", seed, parameters)
    else:
        check_arm(row, "cpm_low", seed, parameters)
