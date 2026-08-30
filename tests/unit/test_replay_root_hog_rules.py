from benchmark_tools.replay_root_hog_rules import _read_clusters, write_groups


def test_root_hog_rule_replay_io(tmp_path):
    clusters = tmp_path / "clusters.txt"
    clusters.write_text("b a\n\nc\n")
    output = tmp_path / "output.txt"

    parsed = _read_clusters(clusters)
    write_groups(output, parsed)

    assert parsed == [("a", "b"), ("c",)]
    assert output.read_text() == "a b\nc\n"
