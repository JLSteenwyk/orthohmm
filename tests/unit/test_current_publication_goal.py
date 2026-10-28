"""Guard the editable execution contract, not the platform dispatcher."""

import hashlib
import os
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
CURRENT = ROOT / "benchmark_tools/PUBLICATION_GOAL_CURRENT.txt"
HISTORICAL = ROOT / "benchmark_tools/PUBLICATION_GOAL_20261003.txt"
MARKER = b"  > 1. Establish and freeze the publication baseline\n"
SCOPE_SHA256 = "bd079a3be2651222dda73b1ad77eae1e118a3717ac17a5881c38a18b65d9e073"


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
    assert "historical checkpoints" in contract
    assert "Retain attempt-bound source/input/output bytes and frozen scientific protocols" in contract


def test_actual_goal_linked_attachment_matches_version_controlled_contract():
    attachment = os.environ.get("ORTHOHMM_CURRENT_GOAL_ATTACHMENT")
    if attachment is None:
        pytest.skip("Optional local attachment is not part of portable repository tests")
    assert Path(attachment).read_bytes() == CURRENT.read_bytes()


@pytest.mark.parametrize("instruction", [
    "Each automatic continuation with an executable next step must actually run it",
    "Finish pending integration before a new candidate",
    "never rerun its producer",
    "inspect the SAME handle afterward",
    "A commit, milestone, ledger update, completed wait, prompt repair or progress message is NOT a stop condition",
    "Keep the full goal ACTIVE while useful required work remains",
    "Older NEXT entries are history, not a queue to replay",
    "make the TOP handoff executable",
    "use an actual wait and inspect its handle",
    "block only its dependent scoring/admission",
    "diagnose it once with focused checks",
    "fresh paths is authorized after focused tests",
    "no replacement panel, duplicate submission or automatic failed inference/scoring/timing retry",
    "Do not stop, suspend, renice or re-affinitize unrelated processes",
    "do not commit large raw data",
    "No DGX, dedicated host, quiet window, free CPU cores or renewed approval for contention is required",
    "genuinely unsafe capacity defers only that launch",
    "CPU load, competing core counts and estimated timing distortion are NOT launch or evidence-admission thresholds",
    "Keep enforced CPU/RAM limits, safe memory and valid accounting",
    "If a necessary job is running/pending, advance independent unfinished work",
    "unknown and potentially tool-dependent impact",
    "not estimates of isolated performance",
    "the goal ACTIVE unless its actual lifecycle rules require another status",
    "an interruption is not evidence of completion or failure",
    "Never change private goal stores/application settings to force continuation",
])
def test_current_execution_safety_and_continuation_rules(contract, instruction):
    assert instruction in contract


def test_continuation_loop_precedes_finite_scientific_scope(contract):
    assert contract.index("One Continuation Loop") < contract.index("Finite Scope And Evidence")
    assert "No DGX, dedicated host, quiet window, free CPU cores" in contract
    assert "Do not assume slight distortion" in contract
    assert "Never change private goal stores/application settings to force continuation" in contract


def test_no_new_gates_or_placeholder_work_are_authorized(contract):
    assert "Do not invent substitute analyses, extra release candidates" in contract
    assert "all-tool hermeticity" in contract
    assert "An unavailable outcome can narrow claims but does not silently fulfill a requirement" in contract
    assert "Do not fabricate intervals, scores or readiness" in contract
