import numpy as np

from orthohmm.orthohmm import _expand_phylogeny_candidates


def test_expand_phylogeny_candidates_writes_partition_and_seed_sidecar(tmp_path):
    working = tmp_path / "orthohmm_working_res"
    working.mkdir()
    clusters = working / "orthohmm_edges_clustered.txt"
    clusters.write_text("a b c\nd\n")

    result = _expand_phylogeny_candidates(
        str(tmp_path),
        ["a", "b", "c", "d"],
        np.asarray([0, 1, 2, 3], dtype=np.int32),
        (
            np.asarray([3, 0], dtype=np.int32),
            np.asarray([0, 3], dtype=np.int32),
            np.asarray([3.0, 3.0]),
        ),
    )

    assert clusters.read_text() == "a b c d\n"
    assert (working / "phylogeny_candidate_seeds.tsv").read_text() == (
        "candidate_family\tseed_families\n"
        "Family0000000\tSeed0000000,Seed0000001\n"
    )
    assert result["profile"] == "satellite_v1"
    assert result["merges"] == 1
