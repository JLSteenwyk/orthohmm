import numpy as np

from orthohmm import orthohmm as pipeline


def test_broad_refinement_skips_unused_hit_and_edge_collection(tmp_path, monkeypatch):
    work_dir = tmp_path / "orthohmm_working_res"
    work_dir.mkdir()

    genes = [f"g{i}" for i in range(150)]
    species = [0] * 10
    for species_idx in range(1, 50):
        species.extend([species_idx] * 2)
    for species_idx in range(1, 43):
        species.append(species_idx)

    cluster_path = work_dir / "orthohmm_edges_clustered.txt"
    cluster_path.write_text(" ".join(genes) + "\n")

    gene_lengths = np.asarray(
        [(f"s{species_id}", gene, 100) for gene, species_id in zip(genes, species)],
        dtype=[("spp", object), ("name", object), ("length", int)],
    )

    def fail_hit_collection(*_args, **_kwargs):
        raise AssertionError("broad refinement should not collect directed hits")

    def fail_edge_collection(*_args, **_kwargs):
        raise AssertionError("broad refinement should not collect RBNH edges")

    monkeypatch.setattr(pipeline, "_collect_search_hit_arrays", fail_hit_collection)
    monkeypatch.setattr(pipeline, "_collect_edge_arrays", fail_edge_collection)

    pipeline._refine_cluster_file(
        str(tmp_path),
        gene_lengths,
        search_results={("a", "b"): object()},
        edges={frozenset(("g0", "g1")): 1.0},
        evalue_threshold=0.0001,
    )

    refined = cluster_path.read_text().splitlines()
    assert len(refined) == len(genes)
    assert {line for line in refined} == set(genes)
