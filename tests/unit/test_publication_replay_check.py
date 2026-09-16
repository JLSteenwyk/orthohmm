import pytest

from benchmark_tools.run_publication_replay_check import OUTPUTS, compare_stages
from benchmark_tools.orthobench_stage_diagnostics import file_provenance


def test_stage_comparison_distinguishes_bytes_and_partition(tmp_path):
    expected = {}
    output = tmp_path / "output"
    output.mkdir()
    for name, filename in OUTPUTS.items():
        path = tmp_path / filename
        path.write_text("a b\nc\n")
        expected[name] = file_provenance(path)
        (output / filename).write_text("c\nb a\n")
    result = compare_stages(expected, output, {"a", "b", "c"})
    assert all(r["partition_equal"] and not r["byte_equal"] for r in result.values())
    (output / OUTPUTS["multipass"]).write_text("a c\nb\n")
    result = compare_stages(expected, output, {"a", "b", "c"})
    assert not result["multipass"]["partition_equal"]
    assert result["multipass"]["observed_only_groups"] == 2
    (output / OUTPUTS["multipass"]).write_text("a b\n")
    with pytest.raises(ValueError, match="cover"):
        compare_stages(expected, output, {"a", "b", "c"})
