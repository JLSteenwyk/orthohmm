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


def test_pending_dependency_does_not_gate_independent_scientific_work():
    contract = CURRENT.read_bytes().split(MARKER, 1)[0].decode("ascii")
    for instruction in (
        "gates ONLY its dependent conversion, scoring, admission",
        "does not gate independent error analyses, manuscript work",
        "distinguish executable independent work from dependency-waiting work",
        "identify the actual evidence needed to unblock a dependency",
        "into a whole-goal stop",
    ):
        assert instruction in contract


def test_waiting_is_bounded_without_repeated_administrative_work():
    contract = CURRENT.read_bytes().split(MARKER, 1)[0].decode("ascii")
    for instruction in (
        "execute a supported independent action if one exists",
        "bounded actual waits (normally 30-300 seconds)",
        "not by rereading the full ledger or reauditing unchanged evidence on every poll",
        "not a new commit for every unchanged wait",
        "Do not grow the scientific scope",
    ):
        assert instruction in contract


def test_wait_only_turn_has_an_action_and_not_a_user_resume_gate():
    contract = CURRENT.read_bytes().split(MARKER, 1)[0].decode("ascii")
    assert contract.index("Immediate continuation decision:") < contract.index(
        "Current execution contract"
    )
    for instruction in (
        "Waiting is an authorized next action, not an absence of a next action",
        "A completed monitoring window, unchanged poll or checkpoint commit",
        "every automatic wait-only turn must include a real wait and a post-wait state check",
        "unless a fresh user instruction or actual interruption preempts it",
        "If the job becomes terminal, replace waiting with the required review or failure diagnosis",
        "Full goal ACTIVE; next automatic continuation:",
        "not authorization to stop the goal",
        "Do not invent work, resubmit jobs or repeat completed analyses to avoid waiting",
    ):
        assert instruction in contract


def test_stale_checkpoint_is_reconciled_without_repeating_completed_execution():
    contract = CURRENT.read_bytes().split(MARKER, 1)[0].decode("ascii")
    for instruction in (
        "before executing a ledger next action, reconcile it",
        "existing output artifacts, execution receipts, source commits and live job handles",
        "If the action has already executed, do not rerun it",
        "Preserve any uncommitted result",
        "actual next unfinished validation or integration step",
        "A stale ledger is a recoverable bookkeeping issue",
    ):
        assert instruction in contract


def test_direct_user_answer_does_not_complete_or_pause_the_full_goal():
    contract = CURRENT.read_bytes().split(MARKER, 1)[0].decode("ascii")
    for instruction in (
        "answer a requested status check or bounded maintenance task promptly",
        "not completion, a pause request or a request for another resume",
        "Leave the full goal ACTIVE",
        "an actionable checkpoint for its next automatic continuation",
        "never promise that editing this attachment repairs the dispatcher",
    ):
        assert instruction in contract


def test_failed_validation_triggers_scoped_diagnosis_not_a_global_stop():
    contract = CURRENT.read_bytes().split(MARKER, 1)[0].decode("ascii")
    for instruction in (
        "Failure recovery is work, not an automatic whole-goal stop",
        "stop ONLY dependent conversion/scoring/admission",
        "immediately diagnose the failure using retained inputs, outputs",
        "Distinguish pipeline defects from parser/validator defects",
        "Do not leave a terminal job in a waiting checkpoint",
        "no fabricated score, weakened validation or automatic retry",
        "never edit a source or artifact bound to an existing attempt",
        "a separate prospective version",
        "check the existing history/authorization contract",
        "bounded checks using already permitted tools",
        "only if no meaningful permitted work remains",
    ):
        assert instruction in contract


def test_independent_task_switch_is_persisted_before_work_can_be_lost():
    contract = CURRENT.read_bytes().split(MARKER, 1)[0].decode("ascii")
    for instruction in (
        "replace the ledger TOP handoff before substantial implementation",
        "distinguish prepared/tested from committed/executed",
        "Do not leave a wait-only or already-committed action at the TOP",
        "On continuation, reconcile that work first",
        "not a new scientific gate",
    ):
        assert instruction in contract


def test_native_identity_state_is_observed_not_declared_statically():
    contract = CURRENT.read_bytes().split(MARKER, 1)[0].decode("ascii")
    assert "Determine attempted, running and genuinely unrun native identities from actual history" in contract
    assert "Genuinely unrun native identities 11 and 12" not in contract
    assert "Any genuinely unrun identity remains sequential behind the reviewed existing history" in contract
