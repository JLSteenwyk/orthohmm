import pytest

from benchmark_tools.audit_orthomcl_blast import audit, parse_diagnostics


def fixture_paths(tmp_path):
    log = tmp_path / "blast.log"
    log.write_text(
        "[blastall] WARNING:  [000.000]  sp|A|A_ONE: SetUpBlastSearch failed.\n"
        "[blastall] ERROR:  [000.000]  sp|A|A_ONE: BLASTSetUpSearch: Unable to calculate Karlin-Altschul params, check query sequence\n"
        "[blastall] WARNING:  [000.000]  sp|B|B_TWO: Blast: Selenocysteine (U) at position 2 replaced by X\n"
    )
    (tmp_path / "one.fasta").write_text(">sp|A|A_ONE description\nXXXX\n")
    (tmp_path / "two.fasta").write_text(">sp|B|B_TWO\nMUAA\n")
    groups = tmp_path / "groups.txt"
    groups.write_text("sp|A|A_ONE sp|B|B_TWO\n")
    return log, groups


def test_counts_queries_not_warning_lines_and_checks_final_membership(tmp_path):
    log, groups = fixture_paths(tmp_path)
    result = audit(log, tmp_path, groups, {"A": 1})
    assert result["failed_queries"] == 1
    assert result["failed_query_fraction"] == 0.5
    assert result["failed_queries_in_final_groups"] == 1
    assert result["failed_queries_qfo_mapping_valid"] == 1
    assert result["selenocysteine_replaced_queries"] == 1
    assert result["records"][0]["residue_counts"] == {"X": 4}
    assert result["records"][0]["final_group_size"] == 2


@pytest.mark.parametrize("text", ["unexpected output", "[blastall] ERROR: [0] gene: new failure"])
def test_unrecognized_diagnostics_fail_closed(tmp_path, text):
    path = tmp_path / "log"
    path.write_text(text)
    with pytest.raises(ValueError):
        parse_diagnostics(path)


def test_missing_diagnostic_gene_is_not_silently_ignored(tmp_path):
    log, groups = fixture_paths(tmp_path)
    (tmp_path / "one.fasta").unlink()
    with pytest.raises(ValueError, match="Diagnostic IDs absent"):
        audit(log, tmp_path, groups, {})


@pytest.mark.parametrize("groups_text", ["sp|A|A_ONE sp|A|A_ONE\n", "sp|A|A_ONE\nsp|A|A_ONE\n", "unknown\n"])
def test_invalid_groups_rejected(tmp_path, groups_text):
    log, groups = fixture_paths(tmp_path)
    groups.write_text(groups_text)
    with pytest.raises(ValueError):
        audit(log, tmp_path, groups, {})


def test_missing_failure_is_reported_without_inventing_subject_hit_evidence(tmp_path):
    log, groups = fixture_paths(tmp_path)
    groups.write_text("sp|B|B_TWO\n")
    result = audit(log, tmp_path, groups, {})
    assert result["failed_queries_in_final_groups"] == 0
    assert result["records"][0]["final_group_line"] is None
