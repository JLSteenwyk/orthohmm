"""Guard the current instructions, not the platform's continuation dispatcher."""

import hashlib
import os
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
CURRENT = ROOT / "benchmark_tools/PUBLICATION_GOAL_CURRENT.txt"
HISTORICAL = ROOT / "benchmark_tools/PUBLICATION_GOAL_20261003.txt"
MARKER = b"  > 1. Establish and freeze the publication baseline\n"
SCOPE_SHA256 = "5a4a46210cbeae65d33c5ebf827d24d9265fbd71122368c7323bba9c7b01656c"


def test_scientific_scope_and_completion_requirements_preserved_byte_for_byte():
    data = CURRENT.read_bytes()
    assert data.count(MARKER) == 1
    scope = MARKER + data.split(MARKER, 1)[1]
    assert hashlib.sha256(scope).hexdigest() == SCOPE_SHA256


def test_historical_evidence_is_not_replaced_by_current_instructions():
    assert hashlib.sha256(HISTORICAL.read_bytes()).hexdigest() == (
        "7d99ecb39a740b689101e885ca9a8e8d337aaa51d2aa78efb5295e4de27acde0"
    )
    text = CURRENT.read_text(encoding="utf-8")
    assert "not an execution-manifest input" in text
    assert "Do not substitute the editable current prompt for a frozen goal binding" in text


def test_current_contract_uses_live_ledger_not_stale_task_or_job_checkpoints():
    contract = CURRENT.read_bytes().split(MARKER, 1)[0].decode("ascii")
    assert "TOP of benchmark_tools/results/PUBLICATION_PROGRESS.md" in contract
    assert "actual scheduler/accounting and retained outcomes supply job state" in contract
    assert "Do not copy transient task or job-state assertions into this prompt" in contract
    for obsolete in ("Live checkpoint", "23902", "23910", "uncommitted and untested",
                     "finish its independent reader", "no review is queued"):
        assert obsolete not in contract


def test_continuation_checks_evidence_and_reassesses_wait_only_checkpoints():
    contract = CURRENT.read_bytes().split(MARKER, 1)[0].decode("ascii")
    for instruction in (
        "at the start of every automatic goal continuation",
        "use a concrete tool call",
        "Do not produce a no-tool acknowledgement or plan-only continuation",
        "not dummy activity to keep the loop alive",
        "no independent work is available is provisional, not a permanent gate",
        "check the outstanding requirements against the retained evidence",
        "Do not impose stronger requirements than the seven-part goal",
        "Record any newly identified action at the ledger TOP",
    ):
        assert instruction in contract


@pytest.mark.parametrize("instruction", [
    "No DGX, dedicated timing host, quiet window or renewed contention approval is required",
    "Do not request another routine authorization or another user resume at a milestone",
    "Do not end routine goal work with a status-only final response",
    "Running/dependency-pending jobs are waiting states, not blockers",
    "use an available wait/sleep tool",
    "Competing analyses alone NEVER block the goal",
    "only the affected launch",
    "keep the full goal ACTIVE",
    "only validated successful native output",
    "no-duplicate and no-automatic-retry rules",
    "unknown and potentially tool-dependent impact",
    "Prompt instructions cannot prevent platform interruptions",
])
def test_continuation_authority_and_scientific_safety_remain_explicit(instruction):
    assert instruction in CURRENT.read_text(encoding="utf-8")


def test_actual_goal_linked_attachment_matches_version_controlled_contract():
    attachment = os.environ.get("ORTHOHMM_CURRENT_GOAL_ATTACHMENT")
    if attachment is None:
        pytest.skip("Optional local goal attachment is not part of portable repository tests")
    assert Path(attachment).read_bytes() == CURRENT.read_bytes()
