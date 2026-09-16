from pathlib import Path
import json

import pytest

from benchmark_tools import replay_phylogeny


def test_parser_accepts_checkpoint_source():
    args = replay_phylogeny.build_parser().parse_args([
        "--fasta-directory", "proteomes",
        "--candidate-clusters", "clusters.txt",
        "--output-directory", "output",
        "--json", "result.json",
        "--checkpoint-source", "previous",
        "--root-rule", "species_overlap",
        "--pair-rule", "positive_paralogy",
        "--species-tree-rooting", "min_variance",
    ])

    assert args.checkpoint_source == Path("previous")
    assert args.root_rule == "species_overlap"
    assert args.pair_rule == "positive_paralogy"
    assert args.species_tree_rooting == "min_variance"


def test_seed_checkpoint_output_copies_required_tree(tmp_path):
    source = tmp_path / "source" / "orthohmm_phylogeny"
    (source / "checkpoints").mkdir(parents=True)
    (source / "gene_trees").mkdir()
    (source / "checkpoints" / "Family0.json").write_text("{}\n")
    (source / "gene_trees" / "Family0.raw.nwk").write_text("(a,b);\n")
    output = tmp_path / "output"

    selected = replay_phylogeny.seed_checkpoint_output(source.parent, output)

    assert selected == source.resolve()
    assert (
        output / "orthohmm_phylogeny" / "checkpoints" / "Family0.json"
    ).read_text() == "{}\n"
    assert (
        output / "orthohmm_phylogeny" / "gene_trees" / "Family0.raw.nwk"
    ).read_text() == "(a,b);\n"


def test_seed_checkpoint_output_rejects_incomplete_source(tmp_path):
    source = tmp_path / "source"
    source.mkdir()

    with pytest.raises(SystemExit, match="missing required directories"):
        replay_phylogeny.seed_checkpoint_output(source, tmp_path / "output")


@pytest.mark.parametrize("constraints", [
    {}, [None], [{"source_genes": "a", "target_genes": ["b"]}],
    [{"source_genes": [], "target_genes": ["b"]}],
    [{"source_genes": ["a", "a"], "target_genes": ["b"]}],
    [{"source_genes": ["a"], "target_genes": ["a"]}],
    [{"source_genes": ["a"], "target_genes": ["unknown"]}],
    [{"source_genes": ["a"], "target_genes": ["c"]}],
])
def test_membership_constraints_reject_invalid_trace(tmp_path, constraints):
    candidates = tmp_path / "candidates.txt"
    candidates.write_text("a b\nc\n")
    trace = tmp_path / "trace.json"
    trace.write_text(json.dumps(constraints))
    with pytest.raises(ValueError):
        replay_phylogeny.load_membership_constraints(trace, candidates)


def test_membership_constraints_accept_empty_and_valid_trace(tmp_path):
    candidates = tmp_path / "candidates.txt"
    candidates.write_text("a b\nc\n")
    trace = tmp_path / "trace.json"
    for constraints in ([], [{"source_genes": ["a"], "target_genes": ["b"]}]):
        trace.write_text(json.dumps(constraints))
        assert replay_phylogeny.load_membership_constraints(trace, candidates) == constraints
    candidates.write_text("a b\na c\n")
    with pytest.raises(ValueError, match="Duplicate"):
        replay_phylogeny.load_membership_constraints(trace, candidates)


def test_main_passes_membership_to_production_stage(tmp_path, monkeypatch):
    candidates = tmp_path / "candidates.txt"
    candidates.write_text("a b\n")
    constraints = [{"source_genes": ["a"], "target_genes": ["b"]}]
    trace = tmp_path / "trace.json"
    trace.write_text(json.dumps(constraints))
    monkeypatch.setattr(replay_phylogeny, "fetch_fasta_files", lambda _: ["species.fa"])
    captured = {}
    class StageReached(Exception):
        pass
    def fake_stage(*args, **kwargs):
        captured.update(kwargs)
        raise StageReached()
    monkeypatch.setattr(replay_phylogeny, "run_phylogeny_stage", fake_stage)
    with pytest.raises(StageReached):
        replay_phylogeny.main([
            "--fasta-directory", str(tmp_path),
            "--candidate-clusters", str(candidates),
            "--output-directory", str(tmp_path / "out"),
            "--json", str(tmp_path / "metrics.json"),
            "--membership-constraints", str(trace),
        ])
    assert captured["membership_constraints"] == constraints
    metrics = json.loads((tmp_path / "metrics.json").read_text())
    assert metrics["metadata"]["membership_constraints"]["sha256"] == (
        replay_phylogeny.file_provenance(trace)["sha256"]
    )


def test_main_rejects_silently_omitted_satellite_trace(tmp_path):
    candidates = tmp_path / "candidates.txt"
    candidates.write_text("a b\n")
    (tmp_path / "phylogeny_candidate_merges.json").write_text("[]\n")
    with pytest.raises(SystemExit, match="contains a satellite merge trace"):
        replay_phylogeny.main([
            "--fasta-directory", str(tmp_path),
            "--candidate-clusters", str(candidates),
            "--output-directory", str(tmp_path / "out"),
            "--json", str(tmp_path / "metrics.json"),
        ])
    assert not (tmp_path / "out").exists()
