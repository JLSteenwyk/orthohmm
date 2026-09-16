import pytest

from benchmark_tools.audit_ygob_overlap import audit, fingerprint, read_pillars


def setup_data(tmp_path):
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    fields = ["---"] * 33
    fields[0], fields[12], fields[32] = "a", "Anc_1", "b"
    (candidate / "Pillars.tab").write_text("\t".join(fields) + "\n")
    (candidate / "AA.fsa").write_text(">a {ON}\nAAAA*\n>b {ON}\nCCCC\n>off {OFF}\nAAAA\n")
    development = tmp_path / "development"
    development.mkdir()
    (development / "input.fa").write_text(">renamed\naaaa\n>another\nAAAA\n")
    return candidate, development


def test_overlap_uses_sequences_not_identifiers_and_counts_multiplicity(tmp_path):
    candidate, development = setup_data(tmp_path)
    result = audit(candidate, {"exposed": development})
    assert result["on_proteins"] == 2
    assert result["off_proteins"] == 1
    assert result["overlaps"]["exposed"]["matching_development_records"] == 2
    assert result["overlaps"]["exposed"]["matching_candidate_proteins"] == 1
    assert result["on_proteins_by_species"] == {"Vpolyspora": 2}
    assert not result["independence_established"]


def test_fingerprint_preserves_internal_stops():
    assert fingerprint("aa*") == fingerprint("AA")
    assert fingerprint("A*A") != fingerprint("AA")


def test_bad_pillar_width_rejected(tmp_path):
    path = tmp_path / "Pillars.tab"
    path.write_text("a\tb\n")
    with pytest.raises(ValueError, match="33 columns"):
        read_pillars(path)


def test_duplicate_pillar_gene_is_reported_not_silently_overwritten(tmp_path):
    candidate, _ = setup_data(tmp_path)
    path = candidate / "Pillars.tab"
    path.write_text(path.read_text() * 2)
    _, genes = read_pillars(path)
    assert genes["a"]["pillar_lines"] == [1, 2]


def test_unknown_on_protein_is_retained_as_unassigned(tmp_path):
    candidate, development = setup_data(tmp_path)
    (candidate / "AA.fsa").write_text(">missing {ON}\nAAAA\n")
    result = audit(candidate, {"exposed": development})
    assert result["on_proteins_without_pillar"] == ["missing"]
    assert result["on_proteins_by_species"] == {"unassigned": 1}
