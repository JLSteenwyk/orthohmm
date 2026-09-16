from types import SimpleNamespace

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import ARMS, controlled_expansion


def engine(parameters, fail=False, calls=1):
    seen = []
    def merge(*args, **kwargs):
        seen.append(kwargs)
        if fail:
            raise RuntimeError("native failure")
        return [], 0, 0, 0
    module = SimpleNamespace(merge_supported_satellite_candidate_clusters=merge)
    def expand(*args, profile):
        assert profile == "satellite_v2"
        for _ in range(calls):
            module.merge_supported_satellite_candidate_clusters(*args, **parameters, merge_trace=[])
        return {"parameters": dict(parameters), "_membership_constraints": []}
    module._expand_phylogeny_candidates = expand
    return module, merge, seen


@pytest.mark.parametrize("label,delta", ARMS)
def test_exact_one_parameter_override_and_restoration(label, delta):
    baseline = {"min_norm": .03, "min_margin": 1.5, "max_iterations": 2}
    module, original, seen = engine(baseline)
    result = controlled_expansion(module, baseline, label, ("fixture",))
    assert module.merge_supported_satellite_candidate_clusters is original
    assert len(seen) == 1
    assert {k: v for k, v in seen[0].items() if k != "merge_trace"} == {**baseline, **delta}
    assert result["applied_parameters"] == {**baseline, **delta}
    assert result["engine_fixed_profile_report"]["parameters"] == baseline
    assert baseline == {"min_norm": .03, "min_margin": 1.5, "max_iterations": 2}


def test_engine_restored_after_failure():
    baseline = {"min_norm": .03, "min_margin": 1.5}
    module, original, _ = engine(baseline, fail=True)
    with pytest.raises(RuntimeError, match="native failure"):
        controlled_expansion(module, baseline, "norm_low", ())
    assert module.merge_supported_satellite_candidate_clusters is original


@pytest.mark.parametrize("count", [0, 2])
def test_unexpected_call_count_rejected(count):
    baseline = {"min_norm": .03, "min_margin": 1.5}
    module, original, _ = engine(baseline, calls=count)
    with pytest.raises(ValueError, match="invocation"):
        controlled_expansion(module, baseline, "control", ())
    assert module.merge_supported_satellite_candidate_clusters is original


def test_changed_fixed_wrapper_parameters_rejected():
    module, original, _ = engine({"min_norm": .05, "min_margin": 1.5})
    with pytest.raises(ValueError, match="wrapper differs"):
        controlled_expansion(module, {"min_norm": .03, "min_margin": 1.5}, "control", ())
    assert module.merge_supported_satellite_candidate_clusters is original


def test_unplanned_arm_rejected():
    module, _, _ = engine({"min_norm": .03, "min_margin": 1.5})
    with pytest.raises(ValueError, match="Unexpected candidate arm"):
        controlled_expansion(module, {"min_norm": .03, "min_margin": 1.5}, "best_score", ())
