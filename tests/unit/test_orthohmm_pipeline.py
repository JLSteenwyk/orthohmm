import numpy as np

from orthohmm import orthohmm as pipeline


def test_broad_refinement_skips_unused_hit_and_edge_collection(tmp_path, monkeypatch):
    work_dir = tmp_path / "orthohmm_working_res"
    work_dir.mkdir()

    genes = [f"g{i}" for i in range(150)]
    species = [0] * 9
    for species_idx in range(1, 50):
        species.extend([species_idx] * 2)
    for species_idx in range(1, 44):
        species.append(species_idx)

    cluster_path = work_dir / "orthohmm_edges_clustered.txt"
    cluster_path.write_text(" ".join(genes) + "\n")

    gene_lengths = np.asarray(
        [(f"s{species_id}", gene, 100) for gene, species_id in zip(genes, species)],
        dtype=[("spp", object), ("name", object), ("length", int)],
    )

    def fail_hit_collection(*_args, **_kwargs):
        raise AssertionError("broad refinement should not collect directed hits")

    def collect_empty_edges(*_args, **_kwargs):
        return (
            np.asarray([], dtype=np.int32),
            np.asarray([], dtype=np.int32),
            np.asarray([], dtype=np.float32),
        )

    monkeypatch.setattr(pipeline, "_collect_search_hit_arrays", fail_hit_collection)
    monkeypatch.setattr(pipeline, "_collect_edge_arrays", collect_empty_edges)

    pipeline._refine_cluster_file(
        str(tmp_path),
        gene_lengths,
        search_results={("a", "b"): object()},
        edges={frozenset(("g0", "g1")): 1.0},
        evalue_threshold=0.0001,
    )

    refined = cluster_path.read_text().splitlines()
    assert len(refined) == 1
    assert set(refined) == {" ".join(sorted(genes))}


def test_refinement_profile_allows_explicit_overrides(tmp_path, monkeypatch):
    work_dir = tmp_path / "orthohmm_working_res"
    work_dir.mkdir()
    (work_dir / "orthohmm_edges_clustered.txt").write_text("g0 g1\n")

    gene_lengths = np.asarray(
        [("s0", "g0", 100), ("s1", "g1", 100)],
        dtype=[("spp", object), ("name", object), ("length", int)],
    )
    captured = {}

    def capture_refinement(*args, **kwargs):
        captured.update(kwargs)
        return args[0]

    monkeypatch.setattr(pipeline, "refine_cluster_indices", capture_refinement)

    pipeline._refine_cluster_file(
        str(tmp_path),
        gene_lengths,
        search_results={},
        edges={},
        evalue_threshold=0.0001,
        refinement_profile="qfo",
        copy_split_min_size=200,
    )

    assert captured["broad_copy_split_min_species_count"] == 24
    assert captured["copy_split_min_size"] == 200
