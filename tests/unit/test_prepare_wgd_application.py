import pytest

from benchmark_tools.prepare_wgd_application import SPECIES, collect_sequences, map_cohort


def mapping():
    assignments = {"a": {"species": "Scerevisiae", "pillar_lines": [1]},
                   "b": {"species": "Scerevisiae", "pillar_lines": [1]},
                   "c": {"species": "Smikatae", "pillar_lines": [1]}}
    owners = {gene: row["species"] for gene, row in assignments.items()}
    return [{"orf_pair": ["a", "b"]}], assignments, owners, {gene: "included" for gene in owners}


def test_complete_homology_mapping():
    mapped, ref, ambiguous = map_cohort(*mapping())
    assert mapped[0]["split_eligible"] and mapped[0]["reference_eligible"]
    assert mapped[0]["reference_pillar"] == "Pillar00001"
    assert ref == {"Pillar00001": ["a", "b", "c"]} and ambiguous == []


def test_ambiguous_reference_does_not_discard_input_pair():
    cohort, assignments, owners, status = mapping()
    assignments["c"]["pillar_lines"] = [1, 2]
    mapped, ref, _ = map_cohort(cohort, assignments, owners, status)
    assert mapped[0]["split_eligible"] and not mapped[0]["reference_eligible"]
    assert len(mapped) == 1 and not ref


def test_conflicting_pillars_not_silently_merged():
    cohort, assignments, owners, status = mapping()
    assignments["b"]["pillar_lines"] = [2]
    mapped, _, _ = map_cohort(cohort, assignments, owners, status)
    assert mapped[0]["reference_pillar"] is None
    assert mapped[0]["reference_reasons"][0]["reason"] == "anchors_in_different_pillars"


def test_missing_anchor_retained_with_reason():
    cohort, assignments, owners, status = mapping()
    del owners["b"]
    status["b"] = "internal_stop"
    mapped, _, _ = map_cohort(cohort, assignments, owners, status)
    assert not mapped[0]["split_eligible"]
    assert mapped[0]["input_reasons"] == [{"gene": "b", "reason": "internal_stop"}]


def test_collects_entire_species_sets_not_only_cohort(tmp_path):
    assignments = {str(i): {"species": species, "pillar_lines": [1]} for i, species in enumerate(SPECIES)}
    path = tmp_path / "input.fasta"
    path.write_text("".join(f">{gene} {{ON}}\nMAA*\n" for gene in assignments)
                    + ">off {OFF}\nMAA\n")
    rows, status, excluded = collect_sequences(path, assignments)
    assert all(len(values) == 1 for values in rows.values())
    assert str(rows["Scerevisiae"][0].seq) == "MAA"
    assert excluded["OFF"] == ["off"] and status["off"] == "OFF"


def test_duplicate_protein_rejected(tmp_path):
    path = tmp_path / "input.fasta"
    path.write_text(">a {OFF}\nMAA\n>a {OFF}\nMAA\n")
    with pytest.raises(ValueError, match="Duplicate"):
        collect_sequences(path, {})


@pytest.mark.parametrize("header", [">a", ">a {ON} {OFF}"])
def test_source_state_must_be_unambiguous(tmp_path, header):
    path = tmp_path / "input.fasta"
    path.write_text(header + "\nMAA\n")
    with pytest.raises(ValueError, match="ON/OFF"):
        collect_sequences(path, {})


def test_missing_complete_species_rejected(tmp_path):
    path = tmp_path / "input.fasta"
    path.write_text(">a {ON}\nMAA\n")
    with pytest.raises(ValueError, match="complete-proteome"):
        collect_sequences(path, {"a": {"species": "Scerevisiae"}})


def test_internal_stop_excluded_without_losing_other_proteins(tmp_path):
    assignments = {str(i): {"species": species} for i, species in enumerate(SPECIES)}
    assignments["bad"] = {"species": "Scerevisiae"}
    path = tmp_path / "input.fasta"
    path.write_text("".join(f">{gene} {{ON}}\nMAA*\n" for gene in assignments if gene != "bad")
                    + ">bad {ON}\nMA*AA*\n")
    rows, status, excluded = collect_sequences(path, assignments)
    assert sum(map(len, rows.values())) == 4
    assert status["bad"] == "internal_stop"
    assert excluded["internal_stop"] == ["bad"]
