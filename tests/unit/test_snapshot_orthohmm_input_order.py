import pytest

from benchmark_tools.snapshot_orthohmm_input_order import record, snapshot


@pytest.fixture
def fixture(tmp_path):
    package = tmp_path / "orthohmm"
    package.mkdir()
    source = package / "files.py"
    source.write_text("def fetch_fasta_files(directory):\n    return ['z.fa', 'a.fa']\n")
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for name in ("a.fa", "z.fa"):
        (inputs / name).write_text(">gene\nMALW\n")
    data = {"datasets": [{"input_directory": str(inputs), "proteomes": 2,
                         "inputs": [record(inputs / n) for n in ("a.fa", "z.fa")]}]}
    runtime = {"records": [record(source)]}
    return tmp_path, source, data, runtime


def test_native_order_not_manifest_sort(fixture):
    core, source, data, runtime = fixture
    result = snapshot(core, data, runtime)
    assert result["datasets"][0]["native_order"] == ["z.fa", "a.fa"]
    assert snapshot(core, data, runtime) == result


def test_source_change_rejected(fixture):
    core, source, data, runtime = fixture
    source.write_text("def fetch_fasta_files(directory): return []")
    with pytest.raises(ValueError, match="frozen runtime"):
        snapshot(core, data, runtime)


def test_input_change_rejected(fixture):
    core, source, data, runtime = fixture
    (core / "inputs/a.fa").write_text(">gene\nCHANGED\n")
    with pytest.raises(ValueError, match="Input bytes"):
        snapshot(core, data, runtime)


@pytest.mark.parametrize("names", ["['a.fa']", "['a.fa', 'a.fa']", "['a.fa', 'other.fa']"])
def test_changed_membership_rejected(fixture, names):
    core, source, data, runtime = fixture
    source.write_text("def fetch_fasta_files(directory): return " + names)
    runtime["records"] = [record(source)]
    with pytest.raises(ValueError, match="membership"):
        snapshot(core, data, runtime)
