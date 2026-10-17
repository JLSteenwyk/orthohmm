import pytest

from benchmark_tools.results.hillclimb_campaign_completion_audit_20261016.supporting_readback import selected_module


@pytest.mark.parametrize("classname", [
    "tests.unit.test_profile", "tests.unit.test_profile.TestProfile", "tests.unit.test_profile.Outer.Inner",
])
def test_module_and_multiple_class_levels_map_to_one_selected_module(classname):
    assert selected_module(classname, {"tests.unit.test_profile", "tests.unit.test_other"}) == "tests.unit.test_profile"


@pytest.mark.parametrize("classname", ["tests.unit.test_profiles", "tests.unit.test_other", ""])
def test_unselected_or_similar_prefix_is_not_counted_as_a_selected_module(classname):
    with pytest.raises(ValueError, match="exactly one"):
        selected_module(classname, {"tests.unit.test_profile"})


def test_ambiguous_module_prefixes_cannot_double_count_a_case():
    with pytest.raises(ValueError, match="exactly one"):
        selected_module("tests.unit.test_profile.TestProfile", {"tests.unit.test_profile", "tests.unit.test_profile.TestProfile"})
