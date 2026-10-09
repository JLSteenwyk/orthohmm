from copy import deepcopy

import pytest

from benchmark_tools.controlled_fragment_runtime import CRITICAL, SCHEMA, adopt_inventory, freeze, package_version


@pytest.fixture
def inventories():
    parent = {"core_commit": "frozen", "argv": ["unchanged"], "environments": {
        "orthohmm": {"python": "same Python", "packages": {**{name: "1.0" for name in CRITICAL}, "notebook": "old"}},
        "orthofinder": {"python": "OF Python", "packages": {"of": "3.1.5"}}}}
    current = deepcopy(parent["environments"])
    current["orthohmm"]["packages"].pop("notebook")
    current["orthohmm"]["packages"]["new-unrelated-package"] = "2.0"
    return parent, current


def test_explicit_adoption_preserves_all_noninventory_scientific_bindings(inventories):
    parent, current = inventories
    original = deepcopy(parent)
    result = adopt_inventory(parent, current, {"sha256": "parent"})
    assert parent == original and result["schema"] == SCHEMA
    assert result["environments"] == current and result["parent_manifest"] == {"sha256": "parent"}
    for key in ("core_commit", "argv"):
        assert result[key] == parent[key]
    amendment = result["inventory_amendment"]
    assert amendment["historical_full_inventory_equal"] is False
    assert amendment["historical_output_equivalence_established"] is False
    assert amendment["differences"] == {"notebook": {"historical": "old", "prospective": None},
                                       "new-unrelated-package": {"historical": None, "prospective": "2.0"}}
    assert set(amendment["scientific_dependencies_unchanged"]) == set(CRITICAL)


@pytest.mark.parametrize("name", CRITICAL)
def test_any_changed_scientific_dependency_refuses_adoption(inventories, name):
    parent, current = inventories
    current["orthohmm"]["packages"][name] = "other"
    with pytest.raises(ValueError, match="Scientific dependency changed"):
        adopt_inventory(parent, current, {})


@pytest.mark.parametrize("name", CRITICAL)
def test_any_missing_scientific_dependency_refuses_adoption(inventories, name):
    parent, current = inventories
    current["orthohmm"]["packages"].pop(name)
    with pytest.raises(ValueError, match="required distribution"):
        adopt_inventory(parent, current, {})


@pytest.mark.parametrize("fault", ["python", "orthofinder", "missing_inventory"])
def test_other_runtime_changes_not_authorized(inventories, fault):
    parent, current = inventories
    if fault == "python":
        current["orthohmm"]["python"] = "other"
    elif fault == "orthofinder":
        current["orthofinder"]["packages"]["of"] = "other"
    else:
        current.pop("orthofinder")
    with pytest.raises(ValueError):
        adopt_inventory(parent, current, {})


def test_normalized_package_name_requires_unique_metadata():
    assert package_version({"Biopython": "1"}, "biopython") == "1"
    with pytest.raises(ValueError, match="ambiguous"):
        package_version({"Biopython": "1", "biopython": "1"}, "biopython")


def test_existing_provenance_never_replaced(tmp_path):
    path = tmp_path / "runtime.json"
    path.write_text("old")
    with pytest.raises(FileExistsError):
        freeze(tmp_path, path)
    assert path.read_text() == "old"
