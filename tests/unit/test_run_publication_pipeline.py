from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from benchmark_tools.run_publication_pipeline import arguments, candidate_policy
from benchmark_tools.audit_publication_pipeline import constraint_argument


def test_candidate_policy_is_local_value_preserving_and_restored():
    seen = []
    def original(output, names, species, hits, profile):
        seen.append(hits)
        return "result"
    pipeline = SimpleNamespace(_expand_phylogeny_candidates=original)
    hits = (np.array([1, 0]), np.array([0, 1]), np.array([0.7, 0.3]))
    receipt = []
    with candidate_policy(pipeline, receipt):
        assert pipeline._expand_phylogeny_candidates("out", ["a", "b"], [0, 1], hits,
                                                     profile="satellite_v2") == "result"
    assert pipeline._expand_phylogeny_candidates is original
    assert list(seen[0][0]) == [0, 1]
    assert list(seen[0][2]) == [0.3, 0.7]
    assert list(hits[0]) == [1, 0]
    assert receipt[0]["scores_modified"] is False


@pytest.mark.parametrize("fault", ["exception", "wrong_profile", "twice"])
def test_restore_on_failure(fault):
    def original(*args, **kwargs):
        if fault == "exception":
            raise RuntimeError("fixture failure")
    pipeline = SimpleNamespace(_expand_phylogeny_candidates=original)
    with pytest.raises((ValueError, RuntimeError)):
        with candidate_policy(pipeline, []):
            call = lambda profile: pipeline._expand_phylogeny_candidates("out", ["a"], [0],
                (np.array([0]), np.array([0]), np.array([1.])), profile=profile)
            call("satellite_v1" if fault == "wrong_profile" else "satellite_v2")
            if fault == "twice":
                call("satellite_v2")
    assert pipeline._expand_phylogeny_candidates is original


def test_frozen_arguments():
    args = SimpleNamespace(input=Path("/input"), output=Path("/out"), cpu=2,
                           aligner=Path("/mafft"), tree_builder=Path("/fasttree"))
    argv = arguments(args)
    assert argv[argv.index("--cpm_resolution") + 1] == "0.1"
    assert argv[argv.index("--phylogeny_candidates") + 1] == "satellite_v2"
    assert argv[argv.index("--species_tree_mode") + 1] == "infer"
    assert argv[argv.index("-o") + 1] == "/out/inference"
    args.cpu = 0
    with pytest.raises(ValueError):
        arguments(args)


@pytest.mark.parametrize("text,expected", [("[]", False), ('[{"source_genes": ["a"]}]', True)])
def test_constraint_presence_matches_applied_semantics(tmp_path, text, expected):
    path = tmp_path / "constraints.json"
    path.write_text(text)
    assert constraint_argument(path) == (path if expected else None)


def test_reject_nonlist_constraints(tmp_path):
    path = tmp_path / "constraints.json"
    path.write_text("{}")
    with pytest.raises(ValueError):
        constraint_argument(path)
