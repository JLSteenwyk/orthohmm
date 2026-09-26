import json

import pytest

from benchmark_tools.audit_phylogeny_sequences import alignment_content, audit, sequences
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("text", ["", ">a\nAA\n>a\nAA\n", ">a\n"])
def test_bad_fasta(tmp_path, text):
    path = tmp_path / "input.fa"
    path.write_text(text)
    with pytest.raises(ValueError):
        sequences(path)


@pytest.mark.parametrize("text,valid", [
    (">g00000000\nAC-D\n>g00000001\na-ce\n", True),
    (">g00000000\nAC-D\n>g00000001\nA-CF\n", False),
    (">g00000000\nAC-D\n>g00000001\nACE\n", False),
    (">g00000000\nAC-D\n>foreign\nA-CE\n", False),
    (">g00000000\nAC-D\n", False),
])
def test_alignment_residues(tmp_path, text, valid):
    candidate, alignment = tmp_path / "candidate.fa", tmp_path / "alignment.fa"
    candidate.write_text(">g00000000\nACD\n>g00000001\nACE\n")
    alignment.write_text(text)
    if valid:
        assert len(alignment_content(candidate, alignment, ["b", "a"], {"a": "ACD", "b": "ACE"})) == 2
    else:
        with pytest.raises(ValueError):
            alignment_content(candidate, alignment, ["a", "b"], {"a": "ACD", "b": "ACE"})


@pytest.fixture
def evidence(tmp_path):
    directory = tmp_path / "phylo"
    directory.mkdir()
    species_dir = directory / "species_tree_inference"
    for base in (directory, species_dir):
        for name in ("candidate_fastas", "alignments"):
            (base / name).mkdir(parents=True)
            (base / name / "F.faa").write_text(">g00000000\nACD\n>g00000001\nACE\n")
    (directory / "checkpoints").mkdir()
    (directory / "checkpoints/F.json").write_text("{}")
    records, inputs = [], []
    for name, gene, seq in (("S1", "a", "ACD"), ("S2", "b", "ACE")):
        path = tmp_path / f"{name}.fa"
        path.write_text(f">{gene}\n{seq}\n")
        records.append(record(path))
        inputs.append(dict(filename=path.name, taxon=name, sha256=record(path)["sha256"]))
    roots = directory / "orthohmm_root_hogs.tsv"
    roots.write_text("root_hog\tsource_family\tgenes\nH\tF\ta,b\n")
    matrix = species_dir / "species_tree_alignment.faa"
    matrix.write_text(">s00000000\nACD\n>s00000001\nACE\n")
    (species_dir / "species_tree.raw.nwk").write_text("(s00000000,s00000001);\n")
    checkpoint = dict(status="complete", selected_family_ids=["F"], species_tree_sha256="test",
                      supermatrix_sha256=record(matrix)["sha256"])
    (species_dir / "checkpoint.json").write_text(json.dumps(checkpoint))
    manifest = dict(input_proteomes=inputs, species_tree_mode="infer", species_tree_sha256="test",
                    species_tree_taxa=["S1", "S2"], species_tree_inference=dict(selected_family_ids=["F"]))
    path = directory / "provenance_manifest.json"
    path.write_text(json.dumps(manifest))
    records.extend([record(path), record(roots), record(directory / "checkpoints/F.json")])
    structure = tmp_path / "structure.json"
    structure.write_text(json.dumps(dict(status="phylogeny_structure_verified", checked_records=records)))
    return directory, structure


def test_sequence_readback(evidence):
    result = audit(*evidence)
    assert result["gene_alignments"] == result["species_alignments"] == 1
    assert result["supermatrix_columns"] == 3
    assert not result["scientific_scores_admitted"]


@pytest.mark.parametrize("change", ["candidate", "alignment", "extra", "matrix", "marker"])
def test_sequence_readback_mutations(evidence, change):
    directory, structure = evidence
    if change in {"candidate", "alignment"}:
        part = "candidate_fastas" if change == "candidate" else "alignments"
        (directory / part / "F.faa").write_text(">g00000000\nACE\n>g00000001\nACD\n")
    elif change == "extra":
        (directory / "alignments/extra.faa").write_text(">extra\nAAA\n")
    elif change == "matrix":
        path = directory / "species_tree_inference/species_tree_alignment.faa"
        path.write_text(">s00000000\nACE\n>s00000001\nACD\n")
        cp = directory / "species_tree_inference/checkpoint.json"
        data = json.loads(cp.read_text())
        data["supermatrix_sha256"] = record(path)["sha256"]
        cp.write_text(json.dumps(data))
    else:
        cp = directory / "species_tree_inference/checkpoint.json"
        data = json.loads(cp.read_text())
        data["selected_family_ids"] = ["F", "F"]
        cp.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        audit(directory, structure)
