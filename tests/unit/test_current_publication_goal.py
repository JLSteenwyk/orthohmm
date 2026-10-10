"""Guard the editable execution contract, not the platform dispatcher."""

import hashlib
import os
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
CURRENT = ROOT / "benchmark_tools/PUBLICATION_GOAL_CURRENT.txt"
HISTORICAL = ROOT / "benchmark_tools/PUBLICATION_GOAL_20261003.txt"
MARKER = b"  > 1. Establish and freeze the publication baseline\n"
SCOPE_SHA256 = "5a4a46210cbeae65d33c5ebf827d24d9265fbd71122368c7323bba9c7b01656c"


@pytest.fixture
def contract():
    return " ".join(CURRENT.read_bytes().split(MARKER, 1)[0].decode("ascii").split())


def test_scientific_scope_and_completion_requirements_preserved_byte_for_byte():
    data = CURRENT.read_bytes()
    assert data.count(MARKER) == 1
    scope = MARKER + data.split(MARKER, 1)[1]
    assert hashlib.sha256(scope).hexdigest() == SCOPE_SHA256


def test_historical_evidence_is_not_replaced_by_current_instructions(contract):
    assert hashlib.sha256(HISTORICAL.read_bytes()).hexdigest() == (
        "7d99ecb39a740b689101e885ca9a8e8d337aaa51d2aa78efb5295e4de27acde0"
    )
    assert "not an execution-manifest input" in contract
    assert "any historical frozen goal binding" in contract


def test_actual_goal_linked_attachment_matches_version_controlled_contract():
    attachment = os.environ.get("ORTHOHMM_CURRENT_GOAL_ATTACHMENT")
    if attachment is None:
        pytest.skip("Optional local attachment is not part of portable repository tests")
    assert Path(attachment).read_bytes() == CURRENT.read_bytes()


@pytest.mark.parametrize("instruction", [
    "At every automatic continuation",
    "use a concrete tool call for the next unfinished authorized task",
    "advance to unfinished validation or integration; never replay the producer",
    "wait and inspect that same job",
    "Do not ask the user to resume after a milestone",
    "full goal remains ACTIVE",
    "Older checkpoints are history, not commands",
    "Keep transient job/task state in the ledger, not in this prompt",
    "Every wait-only continuation includes a real wait and a post-wait observation",
    "prevent ONLY its dependent scoring or admission",
    "A recoverable parser, label, validation or reporting failure is not a whole-goal stop",
    "Deterministic postprocessing recovery is explicitly authorized after focused tests",
    "No duplicate submissions or automatic inference/timing retries",
    "Do not stop, suspend, renice or re-affinitize unrelated jobs",
    "do not commit large raw datasets",
    "Ordinary competing workloads never block the whole goal",
    "genuinely unsafe capacity defers only that launch",
    "CPU load, competing core counts and estimated timing distortion are NOT launch or evidence-admission thresholds",
    "safe memory and accounting requirements pass, even with competing analyses",
    "Scheduler-resource-PENDING is authorized waiting work",
    "unknown and potentially tool-dependent impact",
    "not estimates of isolated performance",
    "Pause only on an explicit user request",
    "mark complete only on verified completion",
    "use blocked only under the goal tool's sustained genuine impasse rule",
    "Distinguish goal status from a temporarily interrupted command runner",
    "Never claim that a prompt edit guarantees platform continuation",
])
def test_current_execution_safety_and_continuation_rules(contract, instruction):
    assert instruction in contract


def test_immediate_action_precedes_long_scientific_scope(contract):
    assert contract.index("Immediate Continuation Decision") < contract.index(
        "Authority And Current State"
    )
    assert "no DGX, dedicated host, quiet window, spare uncontended cores" in contract
    assert "Do not assume slight distortion" in contract
    assert "Do not modify private platform state to force continuation" in contract


def test_no_new_gates_or_placeholder_work_are_authorized(contract):
    assert "Do not add new all-tool hermeticity" in contract
    assert "or invent new analyses merely to keep the goal busy" in contract
    assert "Documenting a limitation does not falsely mark an unmet scientific requirement as fulfilled" in contract
    assert "Do not make more release candidates, copied-verifier layers, integrity receipts or prompt revisions without a specific remaining requirement" in contract
