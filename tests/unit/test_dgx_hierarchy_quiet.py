from pathlib import Path

import pytest


@pytest.mark.parametrize("old,new", [
    ("run_dgx_hierarchy_native_smoke.py", "run_dgx_hierarchy_quiet_smoke.py"),
    ("run_dgx_hierarchy_native_smoke.sh", "run_dgx_hierarchy_quiet_smoke.sh"),
    ("audit_hierarchy_native_smokes.py", "audit_hierarchy_quiet_smokes.py"),
])
def test_only_output_recipe_and_reporting_names_differ(old, new):
    root = Path(__file__).resolve().parents[2] / "benchmark_tools"
    expected = (root / old).read_text().replace("hierarchy_native_smoke_v1", "hierarchy_quiet_smoke_v1")
    expected = expected.replace("hierarchy_native_recipe_v1", "hierarchy_quiet_recipe_v1")
    expected = expected.replace("run_dgx_hierarchy_native_smoke", "run_dgx_hierarchy_quiet_smoke")
    expected = expected.replace("three_hierarchy_native_smokes_validated", "three_hierarchy_quiet_smokes_validated")
    assert (root / new).read_text() == expected


def test_quiet_recipe_matches_transferred_source_inventory():
    import hashlib
    import json

    root = Path(__file__).resolve().parents[2] / "benchmark_tools"
    manifest = json.loads((root / "results/dgx_hierarchy_quiet_recipe_v1_20260918.json").read_text())
    files = [r for r in manifest["records"] if r["kind"] == "file"]
    assert len(files) == 27
    for row in files:
        path = root / Path(row["path"]).name
        if not path.exists():
            path = root / "results" / path.name
        assert len(path.read_bytes()) == row["bytes"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]


def test_actual_quiet_control_preserves_remaining_flag():
    import json
    from benchmark_tools.measure_native_hierarchy_step import evaluate

    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    result = json.loads((root / "dgx_hierarchy_quiet_smokes_21820.json").read_text())
    assert result["scientific_timings_admitted"] == 0
    assert result["publication_ready"] is False
    flags = []
    for row in result["runs"]:
        measurement = row["verification"]["measurement"]
        assert evaluate(measurement["points"], measurement["native"], measurement["job_id"]) == measurement["screening"]
        assert measurement["native"]["exit_code"] == 0
        assert row["output_validation"]["input_genes"] == 645
        flags.append(measurement["screening"]["original_threshold_screen"]["flagged_intervals"])
        for token in ("JobState=COMPLETED", "Restarts=0", "CPUs/Task=20", "MinMemoryNode=96G", "OverSubscribe=NO"):
            assert token in row["scheduler"].split()
    assert flags == [[], [3], []]
    flagged = result["runs"][1]["verification"]["measurement"]["screening"]["hierarchy_intervals"][3]
    assert flagged["host_minus_job_outer_cpu_s"] == pytest.approx(.367417)
    assert flagged["step_cpu_s"]["step_batch"] == pytest.approx(.002749)
