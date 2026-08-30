from io import StringIO

import pytest

from benchmark_tools.orthofinder_mcl_to_orthogroups import (
    iter_mcl_clusters,
    load_sequence_ids,
    write_orthogroups,
)


def test_sequence_only_conversion_preserves_wrapped_mcl_clusters(tmp_path):
    sequence_ids_path = tmp_path / "SequenceIDs.txt"
    sequence_ids_path.write_text(
        "0_0: sp|A1|A_ONE description\n"
        "0_1: sp|A2|A_TWO\n"
        "1_0: sp|B1|B_ONE description\n"
        "1_1: sp|B2|B_TWO\n"
    )
    clusters_path = tmp_path / "clusters.txt_id_pairs.txt"
    clusters_path.write_text(
        "(mclheader\n"
        "mcltype matrix\n"
        "dimensions 4x2\n"
        ")\n"
        "(mclmatrix\n"
        "begin\n"
        "0  0_0 1_0\n"
        "   1_1 $\n"
        "1  0_1 $\n"
        ")\n"
    )

    sequence_ids = load_sequence_ids(sequence_ids_path)
    clusters = list(iter_mcl_clusters(clusters_path))
    output = StringIO()

    assert clusters == [("0_0", "1_0", "1_1"), ("0_1",)]
    assert write_orthogroups(clusters, sequence_ids, output) == (2, 4)
    assert output.getvalue() == (
        "sp|A1|A_ONE sp|B1|B_ONE sp|B2|B_TWO\n"
        "sp|A2|A_TWO\n"
    )


def test_sequence_only_conversion_rejects_unknown_sequence_id(tmp_path):
    clusters_path = tmp_path / "clusters.txt"
    clusters_path.write_text("(mclmatrix\nbegin\n0  0_0 9_9 $\n)\n")

    with pytest.raises(ValueError, match="unknown sequence ID: 9_9"):
        write_orthogroups(
            iter_mcl_clusters(clusters_path), {"0_0": "gene_a"}, StringIO()
        )


def test_sequence_only_conversion_rejects_duplicate_cluster_member(tmp_path):
    clusters_path = tmp_path / "clusters.txt"
    clusters_path.write_text(
        "(mclmatrix\nbegin\n0  0_0 $\n1  0_0 $\n)\n"
    )

    with pytest.raises(ValueError, match="multiple clusters: 0_0"):
        write_orthogroups(
            iter_mcl_clusters(clusters_path), {"0_0": "gene_a"}, StringIO()
        )
