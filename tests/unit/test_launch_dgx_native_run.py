import json

import pytest

from benchmark_tools.launch_dgx_native_run import native_enumerator, read_pinned, select
from benchmark_tools.snapshot_runtime_trees import digest


@pytest.fixture
def spec(tmp_path):
    overhead = tmp_path / "overhead.json"
    overhead.write_text(json.dumps({"status": "verified_overhead_panel_evaluated", "all_protocol_gates_met": True}))
    return {"purpose": "launcher_smoke", "execution_authorized": True,
            "runs": [{"index": i} for i in range(3)],
            "overhead": {"path": str(overhead), "sha256": digest(overhead)}}


def test_exact_smoke_selection(spec):
    assert select(spec, 2) == {"index": 2}


@pytest.mark.parametrize("index", [-1, 3, 27])
def test_wrong_index(spec, index):
    with pytest.raises(ValueError):
        select(spec, index)


def test_unauthorized_and_reordered(spec):
    spec["execution_authorized"] = False
    with pytest.raises(ValueError, match="authorized"):
        select(spec, 0)
    spec["execution_authorized"] = True
    spec["runs"].reverse()
    with pytest.raises(ValueError, match="index/order"):
        select(spec, 0)


def test_scientific_runs_require_separate_smoke_proof(spec):
    spec.update(purpose="scientific_scaling", runs=[{"index": i} for i in range(27)])
    with pytest.raises(KeyError, match="validated_launcher_smoke"):
        select(spec, 0)


def test_changed_json_and_enumerator_rejected(tmp_path):
    source = tmp_path / "source.py"
    source.write_text("def fetch_fasta_files(directory): return ['z.fa', 'a.fa']")
    enumerate_files = native_enumerator({"path": str(source), "sha256": digest(source)})
    assert enumerate_files("unused") == ["z.fa", "a.fa"]
    source.write_text("changed")
    with pytest.raises(ValueError, match="source changed"):
        enumerate_files("unused")
    with pytest.raises(ValueError, match="digest differs"):
        read_pinned(source, "0" * 64)
