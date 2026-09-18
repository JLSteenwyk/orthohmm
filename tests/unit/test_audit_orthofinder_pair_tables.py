from itertools import permutations

import pytest

from benchmark_tools.audit_orthofinder_pair_tables import audit_tables
from benchmark_tools.orthofinder_to_pairwise import iter_pairs


def fixture(tmp_path):
    owners = {"sp|a|A": "A", "sp|b|B": "B", "sp|c|C": "C", "sp|a2|A": "A"}
    for a, b in permutations("ABC", 2):
        directory = tmp_path / "Orthologues" / f"Orthologues_{a}"
        directory.mkdir(parents=True, exist_ok=True)
        rows = f"Orthogroup\t{a}\t{b}\n"
        if {a, b} == {"A", "B"}:
            genes = {"A": "sp|a|A, sp|a2|A", "B": "sp|b|B"}
            rows += f"OG0\t{genes[a]}\t{genes[b]}\n"
        (directory / f"{a}__v__{b}.tsv").write_text(rows)
    return owners, {g: "OG0" for g in owners}


def test_complete_native_tables_and_converter_agree(tmp_path):
    owners, groups = fixture(tmp_path)
    (tmp_path / "Orthologues/A.tsv").write_text("native per-species summary\n")
    result = audit_tables(tmp_path, owners, list("ABC"), groups)
    assert result["directed_tables"] == 6
    assert result["distinct_pairs"] == 2
    assert result["converter_emitted_pairs"] == len(list(iter_pairs(tmp_path))) == 2
    assert result["mcl_membership_checked"] is True
    assert result["accuracy_admitted"] is False


def test_duplicate_relations_counted_not_discarded(tmp_path):
    owners, groups = fixture(tmp_path)
    path = tmp_path / "Orthologues/Orthologues_A/A__v__B.tsv"
    with path.open("a") as stream:
        stream.write("OG0\tsp|a|A\tsp|b|B\n")
    result = audit_tables(tmp_path, owners, list("ABC"), groups)
    assert result["distinct_pairs"] == 2
    assert result["converter_emitted_pairs"] == 3
    assert result["converter_duplicate_pairs"] == 1


@pytest.mark.parametrize("mutation", ["missing", "extra", "header", "unknown", "species", "empty", "reverse", "cross_group", "incomplete_group", "alias", "duplicate_species"])
def test_invalid_evidence_rejected(tmp_path, mutation):
    owners, groups = fixture(tmp_path)
    species = list("ABC")
    path = tmp_path / "Orthologues/Orthologues_A/A__v__B.tsv"
    if mutation == "missing":
        path.unlink()
    elif mutation == "extra":
        path.with_name("extra.tsv").write_text(path.read_text())
    elif mutation == "header":
        path.write_text(path.read_text().replace("Orthogroup\tA\tB", "Orthogroup\tB\tA"))
    elif mutation in ("unknown", "species", "empty"):
        replacement = {"unknown": "absent", "species": "sp|c|C", "empty": ""}[mutation]
        path.write_text(path.read_text().replace("sp|a|A", replacement))
    elif mutation == "reverse":
        path.write_text("Orthogroup\tA\tB\n")
    elif mutation == "cross_group":
        groups["sp|a|A"] = "OTHER"
    elif mutation == "incomplete_group":
        groups.pop("sp|a|A")
    elif mutation == "alias":
        owners["tr|a|other"] = "A"
        groups["tr|a|other"] = "OG0"
    else:
        species.append("A")
    with pytest.raises(ValueError):
        audit_tables(tmp_path, owners, species, groups)
