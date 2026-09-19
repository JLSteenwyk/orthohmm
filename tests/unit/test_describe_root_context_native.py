from copy import deepcopy
import gzip
import json
from pathlib import Path

import pytest

from benchmark_tools import describe_root_context_native as module


@pytest.fixture
def example(tmp_path):
    source = Path(module.__file__).parent / "results/root_context_audit_22020_20260919.json.gz"
    controls = json.loads(gzip.decompress(source.read_bytes()))
    native = controls["trials"][0]["replay"]["measurement"]
    measured = deepcopy(native["lineage"]["measured"])
    row = dict(index=0, method=module.METHODS[0], status="validated", job_id=measured["job_id"],
        root_context=native["context"], original_flagged_intervals=[0, 2], narrow_flagged_intervals=[2],
        native_wall_s=native["lineage"]["native_wall_s"], output_equivalent=True)
    path = tmp_path / "root_context_native_v1/run_00/measurement/lineage_report.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(measured))
    audit = dict(status="root_context_native_audited_not_scientific_admission", scientific_timings_admitted=False,
        runs=[row, dict(index=1, method=module.METHODS[1], status="retained_failure"),
              dict(index=2, method=module.METHODS[2], status="retained_unrun")],
        issues=[dict(stage="retained_fixture")], inventory=[module.record(path)])
    return audit, measured


def test_actual_raw_context_values_partition_without_losing_flags(example):
    audit, measured = example
    result = module.describe(audit)
    row = result["runs"][0]
    subsets = row["subsets"]
    count = len(measured["points"])-1
    assert subsets["all"]["intervals"] == count
    assert subsets["original_flagged"]["intervals"] == 2
    assert subsets["original_unflagged"]["intervals"] == count-2
    assert subsets["narrow_flagged"]["intervals"] == 1
    assert subsets["narrow_unflagged"]["intervals"] == count-1
    comparisons = audit["runs"][0]["root_context"]["intervals"]
    for original, described in zip(comparisons, row["intervals"]):
        assert described["values"]["root_minus_three_named_children_cpu_usec"] == original["root_minus_three_named_children_cpu_usec"]
        assert described["values"]["host_ticks:guest"] == original["enclosing_host_category_ticks"]["guest"]
        assert described["values"]["host_enclosing_window_s"] > 0
    assert all(r["subsets"] is None for r in result["runs"][1:])
    assert result["panel_issues"] == audit["issues"]
    assert result["scientific_timings_admitted"] is False


def test_empty_flagged_subset_and_mismatch_remain_visible(example):
    audit, _ = example
    audit["runs"][0].update(status="output_mismatch", output_equivalent=False,
                            original_flagged_intervals=[], narrow_flagged_intervals=[])
    row = module.describe(audit)["runs"][0]
    assert row["status"] == "output_mismatch" and row["output_equivalent"] is False
    assert row["subsets"]["narrow_flagged"]["distributions"] is None
    assert row["subsets"]["all"]["intervals"] > 0


@pytest.mark.parametrize("flags", [[True], [0, 0], [-1], [999999], [1, 0]])
def test_bad_flag_indices(example, flags):
    audit, measured = example
    audit["runs"][0]["narrow_flagged_intervals"] = flags
    with pytest.raises(ValueError):
        module.intervals(audit["runs"][0], measured)


@pytest.mark.parametrize("fault", ["order", "context", "source", "missing", "duplicate", "admission"])
def test_changed_evidence_rejected(example, fault):
    audit, _ = example
    if fault == "order":
        audit["runs"].reverse()
    elif fault == "context":
        audit["runs"][0]["root_context"]["intervals"][0]["root_minus_system_cpu_usec"] += 1
    elif fault == "source":
        Path(audit["inventory"][0]["path"]).write_text("{}")
    elif fault == "missing":
        audit["inventory"] = []
    elif fault == "duplicate":
        audit["inventory"].append(audit["inventory"][0])
    else:
        audit["scientific_timings_admitted"] = 0
    with pytest.raises(ValueError):
        module.describe(audit)


def test_report_pins_audit_bytes(example, tmp_path):
    audit, _ = example
    path = tmp_path / "audit.json.gz"
    path.write_bytes(gzip.compress(json.dumps(audit).encode()))
    assert module.report(path, module.record(path)["sha256"])["status"] == "native_root_context_described"
    with pytest.raises(ValueError):
        module.report(path, "0"*64)
