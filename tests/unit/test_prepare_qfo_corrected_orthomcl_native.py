from pathlib import Path

import pytest

from benchmark_tools import prepare_qfo_corrected_orthomcl_native as module
from benchmark_tools.probe_orthomcl_native_inference import partition, weighted_graph


def source_fixture(tmp_path, monkeypatch):
    tool = tmp_path / "original"
    tool.mkdir()
    (tool / "orthomcl.pl").write_text("# original\n")
    (tool / "orthomcl_module.pm").write_text('our $PATH_TO_ORTHOMCL = "original/";\n'
        'our $ORTHOMCL_DATA_DIR = "data/";\nour $BLAST_NOCPU = 8;\nour $MCL_INFLATION_DEFAULT = 1.5;\n')
    monkeypatch.setattr(module, "TOOL", tool)
    return tool


def test_isolated_sources_preserve_originals_and_scientific_default(tmp_path, monkeypatch):
    tool = source_fixture(tmp_path, monkeypatch)
    original = (tool / "orthomcl_module.pm").read_bytes()
    output, data = tmp_path / "configured", tmp_path / "data"
    result = module.write_sources(output, data, 180, False)
    assert result["pair_parallel_patch"] is False
    assert (tool / "orthomcl_module.pm").read_bytes() == original
    assert (output / "orthomcl.pl").read_bytes() == (tool / "orthomcl.pl").read_bytes()
    text = (output / "orthomcl_module.pm").read_text()
    assert 'our $MCL_INFLATION_DEFAULT = 1.5;' in text
    assert 'our $BLAST_NOCPU = 180;' in text
    for item in result["originals"] + result["configured_sources"]:
        module.check(item)
    with pytest.raises(FileExistsError):
        module.write_sources(output, data, 180, False)


@pytest.mark.parametrize("threads", [0, -1, True, 1.5])
def test_invalid_resources(tmp_path, threads):
    with pytest.raises(ValueError):
        module.write_sources(tmp_path / "out", tmp_path / "data", threads)


@pytest.mark.parametrize("name", ["has space", 'has"quote', "has;command"])
def test_unsafe_native_shell_paths_rejected(tmp_path, name):
    with pytest.raises(ValueError, match="shell-safe"):
        module.write_sources(tmp_path / name, tmp_path / "data", 2)


def test_relative_native_path_rejected(tmp_path):
    with pytest.raises(ValueError, match="absolute"):
        module.write_sources(Path("relative"), tmp_path / "data", 2)


def test_parallel_patch_is_used_only_when_requested(tmp_path, monkeypatch):
    source_fixture(tmp_path, monkeypatch)
    calls = []
    monkeypatch.setattr(module, "parallelize_source", lambda text: calls.append(text) or text + "# patched\n")
    result = module.write_sources(tmp_path / "out", tmp_path / "data", 2, True)
    assert calls == ["# original\n"] and result["pair_parallel_patch"] is True


@pytest.mark.parametrize("text", ["", "OG0: A(sp) A(sp)\n", "OG0: A(sp)\nOG1: A(sp)\n", "OG0: X(sp)\n"])
def test_partition_rejects_invalid_groups(tmp_path, text):
    path = tmp_path / "groups"
    path.write_text(text)
    with pytest.raises(ValueError):
        partition(path, {"A", "B"})


def test_partition_ignores_group_label_and_member_order(tmp_path):
    path = tmp_path / "groups"
    path.write_text("OG8: B(sp2) A(sp1)\n")
    assert partition(path, {"A", "B"}) == [("A", "B")]


def test_weighted_graph_checks_scores_and_ignores_row_order(tmp_path):
    matrix, index, gg = (tmp_path / name for name in ("matrix", "index", "gg"))
    index.write_text("0\tA\n1\tB\n")
    gg.write_text("sp1: A\nsp2: B\n")
    matrix.write_text("(mclheader\nmcltype matrix\ndimensions 2x2\n)\n(mclmatrix\nbegin\n0 1:1.0 $\n1 0:1.0 $\n)\n")
    expected = weighted_graph(matrix, index, gg)
    assert expected == [("A", "B", "1.0"), ("B", "A", "1.0")]
    matrix.write_text(matrix.read_text().replace("0 1:1.0 $\n1 0:1.0 $", "1 0:1.0 $\n0 1:1.0 $"))
    assert weighted_graph(matrix, index, gg) == expected
    matrix.write_text(matrix.read_text().replace("1:1.0", "1:2.0"))
    with pytest.raises(ValueError):
        weighted_graph(matrix, index, gg)
