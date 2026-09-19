import json
from pathlib import Path
import shutil
import subprocess

import pytest


@pytest.mark.parametrize("ratio", [None, -.005, .02])
def test_summary_preserves_missing_signed_values_and_coverage(ratio):
    jq = shutil.which("jq")
    if not jq:
        pytest.skip("jq unavailable")
    script = Path(__file__).resolve().parents[2] / "benchmark_tools/summarize_lineage_overhead.jq"
    data = dict(
        runs=[dict(index=0, mode="boundary", flagged_intervals=None,
                   narrow_flagged_intervals=None, interval_screening_available=False),
              dict(index=1, mode="periodic", flagged_intervals=[1, 2],
                   narrow_flagged_intervals=[2], interval_screening_available=True)],
        paired=dict(methods=[dict(median=ratio)], pairs=[dict(wall_ratio_minus_one=ratio)]),
        scientific_timings_admitted=False, environmental_validity_established=False,
        publication_ready=False,
    )
    result = subprocess.run([jq, "--arg", "sha256", "a" * 64, "-f", str(script)],
                            input=json.dumps(data), text=True, capture_output=True, check=True)
    summary = json.loads(result.stdout)
    expected = None if ratio is None else ratio * 100
    assert summary["methods"][0]["median_percent"] == expected
    assert summary["pairs"][0]["percent"] == expected
    assert summary["runs"][0]["original_flags"] is None
    assert summary["runs"][0]["narrow_flags"] is None
    assert summary["runs"][0]["interval_screening_available"] is False
    assert summary["runs"][1]["original_flags"] == 2
    assert summary["runs"][1]["narrow_flags"] == 1
    assert summary["scientific_timings_admitted"] is False
    assert summary["environmental_validity_established"] is False
    assert summary["publication_ready"] is False
    assert summary["source_audit"]["sha256"] == "a" * 64
