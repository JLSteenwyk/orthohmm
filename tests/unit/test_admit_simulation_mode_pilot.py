import pytest

from benchmark_tools.admit_simulation_mode_pilot import retained_equivalence, admit
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("problem", [None, "matrix", "extra_comment", "missing_comment"])
def test_mcl_comment_exception_does_not_hide_membership_changes(tmp_path, problem):
    a, b = tmp_path / "a", tmp_path / "b"
    a.write_text("# cline: old\n(matrix\n0 1 2 $\n)\n")
    text = "# cline: new\n(matrix\n0 1 2 $\n)\n"
    if problem == "matrix":
        text = text.replace("0 1 2", "0 1 3")
    elif problem == "extra_comment":
        text += "# cline: extra\n"
    elif problem == "missing_comment":
        text = text.replace("# cline: new\n", "")
    b.write_text(text)
    key = "clusters_OrthoFinder_I1.2.txt"
    assert retained_equivalence("orthofinder_full", {key: record(a)}, {key: record(b)}) == (problem is None)


@pytest.mark.parametrize("key,allowed", [("Alignments_ids/SpeciesTreeAlignment.fa", True),
    ("Trees_ids/OG000001.txt", False), ("SequenceIDs.txt", False)])
def test_only_species_tree_alignment_may_be_absent(tmp_path, key, allowed):
    path = tmp_path / "x"
    path.write_text("data")
    assert retained_equivalence("orthofinder_full", {key: record(path)}, {}) == allowed


def test_orthohmm_requires_byte_identity(tmp_path):
    path = tmp_path / "x"
    path.write_text("data")
    before = {"x": record(path)}
    assert retained_equivalence("orthohmm_satellite_v2", before, before)
    assert not retained_equivalence("orthohmm_satellite_v2", before, {})


def test_existing_admission_refused(tmp_path):
    with pytest.raises(FileExistsError):
        admit(tmp_path, tmp_path)
