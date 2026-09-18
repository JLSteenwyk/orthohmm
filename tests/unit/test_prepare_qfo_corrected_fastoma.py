from copy import deepcopy
from pathlib import Path

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_qfo_corrected_fastoma import copy_inputs, prepare, validate_tree_binding


def bindings():
    primary = {"input_directory": "/corrected", "methods": {"orthofinder_full": {"output": "/of/corrected"}}}
    inputs = [{"path": f"/corrected/sp{i}.fasta", "bytes": 3, "sha256": "placeholder"} for i in range(78)]
    tree = {"path": "/of/corrected/Results/Species_Tree/SpeciesTree_rooted_node_labels.txt",
            "bytes": 3, "sha256": "tree"}
    content = {"species_tree": tree, "results_directory": "/of/corrected/Results"}
    admission = {"checked_records": [tree, *inputs]}
    return content, admission, primary, inputs


def test_corrected_tree_binding():
    content, admission, primary, inputs = bindings()
    assert validate_tree_binding(content, admission, primary, inputs) == content["species_tree"]


@pytest.mark.parametrize("change", ["old_tree", "wrong_tree_name", "unbound_tree", "unbound_input",
                                   "missing_species", "duplicate_species", "old_input", "wrong_suffix"])
def test_tree_binding_rejects_substitutions(change):
    content, admission, primary, inputs = deepcopy(bindings())
    if change == "old_tree":
        content["results_directory"] = "/of/original/Results"
    elif change == "wrong_tree_name":
        content["species_tree"]["path"] = "/of/corrected/Results/Species_Tree/other.txt"
    elif change == "unbound_tree":
        admission["checked_records"] = inputs[:]
    elif change == "unbound_input":
        admission["checked_records"].pop()
    elif change == "missing_species":
        inputs.pop()
    elif change == "duplicate_species":
        inputs[-1] = inputs[0]
    elif change == "old_input":
        inputs[0]["path"] = "/original/sp0.fasta"
    else:
        inputs[0]["path"] = "/corrected/sp0.fa"
    with pytest.raises(ValueError):
        validate_tree_binding(content, admission, primary, inputs)


def sources(tmp_path):
    first, second, tree = [tmp_path / name for name in ("sp1.fasta", "sp2.fasta", "tree.txt")]
    first.write_text(">a\nACDE\n")
    second.write_text(">b\nACDF\n")
    tree.write_text("(sp1:1,sp2:1)N0;\n")
    return [record(first), record(second)], record(tree)


def test_stage_copies_exact_bytes_without_links(tmp_path):
    inputs, tree = sources(tmp_path)
    directory = tmp_path / "new"
    copied = copy_inputs(inputs, tree, directory)
    assert len(copied) == 3
    assert {p.name for p in (directory / "proteome").iterdir()} == {"sp1.fa", "sp2.fa"}
    for row in copied:
        source, target = Path(row["source"]["path"]), Path(row["staged"]["path"])
        assert source.read_bytes() == target.read_bytes()
        assert source.stat().st_ino != target.stat().st_ino
        assert not target.is_symlink()


def test_source_tampering_is_rejected(tmp_path):
    inputs, tree = sources(tmp_path)
    Path(inputs[0]["path"]).write_text("changed")
    with pytest.raises((ValueError, RuntimeError)):
        copy_inputs(inputs, tree, tmp_path / "new")


def test_duplicate_output_names_rejected(tmp_path):
    inputs, tree = sources(tmp_path)
    with pytest.raises(ValueError, match="Duplicate"):
        copy_inputs([inputs[0], inputs[0]], tree, tmp_path / "new")


@pytest.mark.parametrize("existing", ["stage", "report"])
def test_no_reuse_or_overwrite(tmp_path, existing):
    (tmp_path / existing).touch()
    with pytest.raises(FileExistsError):
        prepare(tmp_path, tmp_path / "absent", "unused", 1,
                tmp_path / "stage", tmp_path / "report")


def test_stage_directory_not_overwritten(tmp_path):
    inputs, tree = sources(tmp_path)
    with pytest.raises(FileExistsError):
        copy_inputs(inputs, tree, tmp_path)
