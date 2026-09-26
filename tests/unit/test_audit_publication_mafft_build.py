from benchmark_tools.audit_publication_mafft_build import compare_files


def test_content_comparison_not_filename(tmp_path):
    left, right = tmp_path / "left", tmp_path / "right"
    left.write_bytes(b"abc")
    right.write_bytes(b"abc")
    result = compare_files(left, right)
    assert result["byte_equal"]
    assert result["original"]["path"] != result["rebuilt"]["path"]
    right.write_bytes(b"abd")
    assert not compare_files(left, right)["byte_equal"]
