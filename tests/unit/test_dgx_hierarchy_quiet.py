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
