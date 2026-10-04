import hashlib
import json
import pytest

from benchmark_tools import audit_publication_requirements as audit


GOAL = "Goal: Preserve the full publication scope.\n\n> Accept shared timing.\n\n" + "".join(
    f"> {i}. Original section {i}\n\n> - Requirement {i}.\n\n" for i in range(1, 8)
) + "> Engineering and reporting requirements\n\n> - Preserve user changes.\n\n> Completion criteria\n\n> Complete only with evidence;\n> retain missing science.\n"


@pytest.fixture
def inputs(tmp_path):
    goal = tmp_path / "goal.txt"
    goal.write_text(GOAL)
    (tmp_path / "evidence.json").write_text('{"scope":"direct evidence only"}\n')
    review = tmp_path / "review.md"
    review.write_text("# Human review\n\n" + "\n\n".join(
        f"## {r['id']} | {'unmet' if r['id'].startswith('completion') else 'supported'}\n\n"
        "Bounded evidence, not transitive certification. [Receipt](evidence.json)."
        for r in audit.requirements(GOAL)
    ))
    return goal, review


def test_wrapped_completion_is_one_requirement():
    rows = audit.requirements(GOAL)
    assert len(rows) == 10
    assert rows[-1] == {"id": "completion.1", "requirement": "Complete only with evidence; retain missing science."}
    assert [r["id"] for r in rows] == ["scope.1", *[f"{i}.1" for i in range(1, 8)], "engineering.1", "completion.1"]


def test_headings_do_not_need_surrounding_blank_lines():
    compact = GOAL.replace("requirements\n\n> -", "requirements\n> -").replace(
        "criteria\n\n> Complete", "criteria\n> Complete")
    assert audit.requirements(compact) == audit.requirements(GOAL)


def test_adjacent_bullets_and_wrapped_bullet_are_distinct():
    expanded = GOAL.replace("> - Requirement 1.",
        "> - Requirement 1.\n> - Second requirement\n> continued here.\n> - Third requirement.")
    rows = audit.requirements(expanded)
    assert len(rows) == 12
    assert rows[2:4] == [
        {"id": "1.2", "requirement": "Second requirement continued here."},
        {"id": "1.3", "requirement": "Third requirement."},
    ]


@pytest.mark.parametrize("change", ["missing", "duplicate", "reordered", "unrecognized"])
def test_goal_section_integrity(change):
    altered = GOAL
    if change == "missing":
        altered = altered.replace("> 7. Original section 7\n\n> - Requirement 7.\n\n", "")
    elif change == "duplicate":
        altered = altered.replace("> 7. Original section 7", "> 6. Original section 7")
    elif change == "reordered":
        altered = altered.replace("> 6. Original section 6", "> 7. Original section 6").replace("> 7. Original section 7", "> 6. Original section 7")
    else:
        altered = altered.replace("> - Requirement 7.", "> Unexpected non-bullet instruction.")
    with pytest.raises(ValueError):
        audit.requirements(altered)


def test_exact_coverage_and_byte_pins(inputs, tmp_path):
    goal, review = inputs
    result = audit.build(goal, review, tmp_path / "new")
    assert result["requirement_count"] == 10
    assert result["status_counts"] == {"supported": 9, "unmet": 1}
    assert result["goal_completion_proven"] is False
    assert result["goal"]["sha256"] == hashlib.sha256(goal.read_bytes()).hexdigest()
    assert (tmp_path / "new" / "goal.txt").read_bytes() == goal.read_bytes()
    assert json.loads((tmp_path / "new" / "audit.json").read_text()) == result
    assert len(result["evidence"]) == 1


@pytest.mark.parametrize("change", ["omit", "extra", "duplicate", "complete", "no_evidence", "remote"])
def test_bad_review_cannot_emit_audit(inputs, tmp_path, change):
    goal, review = inputs
    text = review.read_text()
    if change == "omit":
        text = text[:text.index("## completion.1")]
    elif change == "extra":
        text += "\n\n## 1.2 | supported\n\n[Receipt](evidence.json)"
    elif change == "duplicate":
        text += "\n\n## 1.1 | supported\n\n[Receipt](evidence.json)"
    elif change == "complete":
        text = text.replace("| unmet", "| complete")
    elif change == "no_evidence":
        text = text.replace("[Receipt](evidence.json)", "No evidence")
    else:
        text = text.replace("(evidence.json)", "(https://example.org/unreviewed)")
    review.write_text(text)
    with pytest.raises(ValueError):
        audit.build(goal, review, tmp_path / "new")
    assert not (tmp_path / "new").exists()


def test_existing_output_is_preserved(inputs, tmp_path):
    goal, review = inputs
    output = tmp_path / "existing"
    output.mkdir()
    marker = output / "marker"
    marker.write_text("preserve")
    with pytest.raises(FileExistsError):
        audit.build(goal, review, output)
    assert marker.read_text() == "preserve"
    assert sorted(p.name for p in output.iterdir()) == ["marker"]


def test_evidence_identity_changes_without_certifying_science(inputs, tmp_path):
    goal, review = inputs
    first = audit.build(goal, review, tmp_path / "first")
    (tmp_path / "evidence.json").write_text('{"scope":"changed"}\n')
    second = audit.build(goal, review, tmp_path / "second")
    assert first["evidence"] != second["evidence"]
    assert second["goal_completion_proven"] is False


def test_missing_evidence_does_not_emit(inputs, tmp_path):
    goal, review = inputs
    (tmp_path / "evidence.json").unlink()
    with pytest.raises(FileNotFoundError):
        audit.build(goal, review, tmp_path / "new")
    assert not (tmp_path / "new").exists()


def test_all_supported_still_does_not_certify_goal(inputs, tmp_path):
    goal, review = inputs
    review.write_text(review.read_text().replace("| unmet", "| supported"))
    result = audit.build(goal, review, tmp_path / "new")
    assert result["goal_completion_proven"] is False
    assert result["status"] == "publication_requirement_audit_incomplete"
