from pathlib import Path

import pytest

from benchmark_tools.admit_ob_cpm_replay import check_metrics, OUTPUTS


def fixture():
    row = {"resolution": .08, "command": ["replay", "--cpm-resolution", "0.08"]}
    cache, source = {"path": "cache"}, {"path": "source"}
    counts = {"genes": 251378, "species": 12, "significant_hits": 18235373, "profiles_built": 3,
        "profile_candidates": 10, "significant_profile_hits": 4, "strict_profile_edges": 2, "calibrated_profiles": 0}
    metrics = {"parameters": {"accuracy_profile": "high_sensitivity", "cpm_resolution": .08,
        "jackknife_profile_thresholds": False, "jackknife_single_copy_profiles": False,
        "leiden_seed": 4, "matrix": "BLOSUM62", "profile_expansion": True,
        "profile_iterations": 1, "profile_min_species": 1}, "command": row["command"],
        "input": cache, "source": source, "cwd": "/frozen",
        "git": {"commit": "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806", "dirty": False},
        "counts": counts, "profile_iterations": [{"iteration": 1, **{key: counts[key] for key in
            ("profiles_built", "profile_candidates", "significant_profile_hits", "strict_profile_edges", "calibrated_profiles")}}],
        "stages": [{"label": label, "clusters": 10, "output": {}} for label in OUTPUTS]}
    return metrics, row, cache, source, Path("/frozen")


def test_exact_metadata_accepted():
    check_metrics(*fixture())


@pytest.mark.parametrize("problem", ["resolution", "command", "cache", "git", "counts", "no_profiles", "iterations", "accounting", "stage", "scored"])
def test_changed_or_incomplete_profile_execution_rejected(problem):
    args = fixture()
    metrics = args[0]
    if problem == "resolution":
        metrics["parameters"]["cpm_resolution"] = .1
    elif problem == "command":
        metrics["command"] = ["other"]
    elif problem == "cache":
        metrics["input"] = {}
    elif problem == "git":
        metrics["git"]["dirty"] = True
    elif problem == "counts":
        metrics["counts"]["genes"] -= 1
    elif problem == "no_profiles":
        metrics["profile_iterations"][0]["profiles_built"] = 0
    elif problem == "iterations":
        metrics["profile_iterations"].append(dict(metrics["profile_iterations"][0]))
    elif problem == "accounting":
        metrics["counts"]["strict_profile_edges"] += 1
    elif problem == "stage":
        metrics["stages"].reverse()
    else:
        metrics["stages"][0]["f_score"] = 90
    with pytest.raises(ValueError):
        check_metrics(*args)
