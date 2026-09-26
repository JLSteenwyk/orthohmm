import hashlib

import pytest

from benchmark_tools.trace_native_residue_groups import scan


def test_membership_and_cross_species_pairs(tmp_path):
    path = tmp_path / "groups"
    data = b"ORTHOMCL0(3 genes,2 taxa): a(s1) b(s1) c(s2)\n"
    path.write_bytes(data)
    result = scan(path, hashlib.sha256(data).hexdigest(), {"a", "missing"})
    assert result["absent_targets"] == ["missing"]
    assert result["incident_cross_species_pairs"] == 1
    assert result["groups"][0]["targets"] == ["a"]


@pytest.mark.parametrize("data", [
    b"bad\n", b"ORTHOMCL0(2 genes,1 taxa): a(s1)\n",
    b"ORTHOMCL0(1 genes,1 taxa): a(s1)\nORTHOMCL1(1 genes,1 taxa): a(s1)\n",
])
def test_reject_malformed_or_duplicate(tmp_path, data):
    path = tmp_path / "groups"
    path.write_bytes(data)
    with pytest.raises(ValueError):
        scan(path, hashlib.sha256(data).hexdigest(), {"a"})


def test_reject_wrong_digest(tmp_path):
    path = tmp_path / "groups"
    path.write_text("ORTHOMCL0(1 genes,1 taxa): a(s1)\n")
    with pytest.raises(ValueError, match="digest"):
        scan(path, "0" * 64, {"a"})


def test_single_species_has_no_submitted_pairs(tmp_path):
    path = tmp_path / "groups"
    data = b"ORTHOMCL0(2 genes,1 taxa): a(s1) b(s1)\n"
    path.write_bytes(data)
    result = scan(path, hashlib.sha256(data).hexdigest(), {"a", "b"})
    assert result["incident_cross_species_pairs"] == 0
    assert result["absent_targets"] == []
