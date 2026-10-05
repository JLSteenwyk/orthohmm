import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record


ROOT = Path(__file__).resolve().parents[2]
GUIDE = ROOT / "benchmark_tools/PUBLICATION_REPRODUCTION.md"


def current_entry():
    return GUIDE.read_text().split("The [native factorial adapter diagnostic]", 1)[0]


@pytest.mark.parametrize("module,flags", [
    ("prepare_native_factorial_qfo_pairs", ("--request", "--request-sha256", "--terminal-review", "--terminal-review-sha256", "--output-directory")),
    ("run_native_factorial_qfo_assessment", ("--root", "--pairs", "--pairs-sha256", "--conversion-job", "--check-only")),
    ("admit_native_factorial_qfo_assessment", ("--root", "--pairs", "--pairs-sha256", "--conversion-job", "--assessment-job", "--output-directory")),
    ("export_native_qfo_factorial_scores", ("--plan", "--plan-sha256", "--admission", "--output")),
])
def test_documented_clis_expose_required_bindings_without_execution(module, flags):
    name = "benchmark_tools." + module
    assert name in current_entry()
    process = subprocess.run([sys.executable, "-B", "-m", name, "--help"],
                             cwd=ROOT, text=True, capture_output=True, check=True)
    assert all(flag in process.stdout for flag in flags)
    assert not process.stderr


def test_guide_routes_to_committed_native_evidence():
    files = (
        "results/native_factorial_progress_20261005_v7/report.md",
        "results/NATIVE_FACTORIAL_SHARED_ATTEMPT_22437.md",
        "results/NATIVE_QFO_TERMINAL_REVIEW_22440.md",
        "results/NATIVE_QFO_CONVERSION_AND_ASSESSMENT_22439.md",
        "results/NATIVE_QFO_ADMISSION_AND_REPORTING_22442.md",
        "results/NATIVE_QFO_SCORE_RESULT_22442.md",
        "results/native_qfo_scores_20261005_v1/scores.md",
    )
    text = current_entry()
    for name in files:
        assert "(" + name + ")" in text
        assert (GUIDE.parent / name).is_file()
    snapshot = json.loads((GUIDE.parent / "results/native_factorial_progress_20261005_v7/report.json").read_text())
    assert len(snapshot["rows"]) == 13
    assert sum(r["job_id"] is not None for r in snapshot["rows"]) == 7
    assert "seven reviewed attempts" in text and "thirteen planned identities" in text


def test_current_status_does_not_override_live_gates_or_archives():
    text = current_entry()
    assert "actual scheduler state supersede dated live-status" in text
    assert "Do not duplicate the currently queued/running jobs" in text
    assert "actual completed export from successful22442" in text
    assert "Six\nothers remain unavailable, not zero" in text
    assert "They do not already contain these later native executions" in text
    assert "does not establish publication readiness" in text
    assert "not portable fresh-install" in text


def test_metric_coverage_and_uncertainty_limits_remain_distinct():
    text = current_entry()
    assert "VGNC/SwissTrees/TreeFam-A as F1, GO/EC similarity and FAS separately" in text
    assert "project-defined secondary summary" in text
    assert "not paired-family confidence intervals" in text
    assert "not\nreference coverage or accuracy" in text
    assert "R-on uses resolved native pairs instead" in text
    assert "reused P1C0R0 and older cached scores are not" in text


def test_binary_record_is_not_an_invocation_environment(tmp_path):
    binary = tmp_path / "binary"
    binary.write_bytes(b"synthetic executable bytes")
    invocation = tmp_path / "environment-python"
    invocation.symlink_to(binary)
    assert record(invocation) == record(binary)
    assert record(invocation)["path"] != str(invocation)
    text = current_entry()
    assert "Preserve the Python invocation path" in text
    assert 'record(python)["path"]' in text
    assert "not installed packages" in text
    assert "retain the original virtual-environment" in text


def test_diagnostic_pending_state_is_qualified_as_historical():
    text = GUIDE.read_text()
    assert "At that diagnostic checkpoint" in text
    assert "route above supersedes that pending state" in text
    assert "New full-cost executor/handoff\nremain required" not in text


def test_completed_native_table_not_confused_with_complete_factorial():
    table = json.loads((GUIDE.parent / "results/native_qfo_scores_20261005_v1/report.json").read_text())
    assert table["supplied_admissions"] == 1 and len(table["rows"]) == 7
    assert table["rows"][0]["index"] == 6 and table["rows"][0]["status"] == "supplied_native_admission"
    assert all(r["scores"]["VGNC"] is None for r in table["rows"][1:])
    assert "0.6900998166155209" in current_entry()
    assert table["new_scoring_or_admission"] is False and table["publication_ready"] is False


def test_manuscript_preserves_actual_native_metrics_and_partial_scope():
    results = GUIDE.parent / "results"
    manuscript = (results / "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md").read_text()
    section = manuscript.split("### First Fresh Native QfO Ablation", 1)[1].split(
        "The original-release QfO factorial", 1)[0]
    row = json.loads((results / "native_qfo_scores_20261005_v1/report.json").read_text())["rows"][0]
    for value in row["scores"].values():
        assert str(value) in section
    assert str(row["secondary_mean"]) in section
    assert "one of seven fresh native QfO identities" in section
    assert "not a non-HMM control" in section
    assert "No cached result or its family intervals" in section
    assert "not replace the retained" in section
    assert "unknown, tool-dependent" in section
    assert "not establish isolated performance" in section


def test_claim_checklist_does_not_promote_partial_native_result():
    text = (GUIDE.parent / "results/PUBLICATION_CLAIMS_20260916.md").read_text()
    row = next(line for line in text.splitlines()
               if line.startswith("| The first fresh native QfO ablation"))
    assert "Supported for P0C0R0 only" in row
    assert "One of seven fresh identities is scored; six are unavailable" in row
    assert "not official QfO F1" in row
    assert "no inherited cached intervals" in row
    assert "Earlier cached scores and comparator claims remain separate" in row
    for name in ("NATIVE_QFO_SCORE_RESULT_22442.md", "native_qfo_scores_20261005_v1/scores.md",
                 "native_qfo_assessment_admission_22442.json", "native_qfo_score_readback_22442.json"):
        assert "(" + name + ")" in row
        assert (GUIDE.parent / "results" / name).is_file()
    assert "| The package is publication-ready | All sections below | Not achieved |" in text
