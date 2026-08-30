import pytest

from benchmark_tools.build_permissive_candidates import materialize_clusters


def test_materialize_clusters_maps_ids_and_restores_isolates(tmp_path):
    raw = tmp_path / "mcl.txt"
    raw.write_text("2\t0\n# comment\n")
    output = tmp_path / "candidates.txt"

    count = materialize_clusters(raw, output, ["a", "b", "c", "d"])

    assert count == 3
    assert output.read_text().splitlines() == ["a c", "b", "d"]


def test_materialize_clusters_rejects_repeated_members(tmp_path):
    raw = tmp_path / "mcl.txt"
    raw.write_text("0 1\n1 2\n")

    with pytest.raises(ValueError, match="repeated gene IDs"):
        materialize_clusters(raw, tmp_path / "out.txt", ["a", "b", "c"])
