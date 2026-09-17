import pytest

from benchmark_tools.score_simulation_tree_panel import check_artifact_gate, score
from benchmark_tools.summarize_simulation_tree_panel import CONDITIONS, SEEDS, METHODS, CONTRASTS


@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "wrong_admission", "status", "unknown", "different"])
def test_complete_artifact_gate_preserves_real_differences(problem):
    indexed = {(c, s, m, v): {} for c in CONDITIONS for s in SEEDS for m in METHODS for v in ("generating", "nni1", "nni2")}
    rows = [{"condition": c, "seed": s, "method": m, "target": a, "reference": b, "status": "retained_upstream_equivalent"}
            for c in CONDITIONS for s in SEEDS for m in METHODS for a, b in CONTRASTS]
    artifact = {"status": "tree_artifact_contrasts_checked", "accuracy_evaluated": False, "admission": {"sha256": "fixture"}, "contrasts": rows}
    if problem == "missing":
        rows.pop()
    elif problem == "duplicate":
        rows[-1] = rows[0]
    elif problem == "wrong_admission":
        artifact["admission"] = {}
    elif problem == "status":
        artifact["status"] = "partial"
    elif problem == "unknown":
        rows[0]["status"] = "unchecked"
    elif problem == "different":
        rows[0]["status"] = "retained_upstream_different"
    if problem not in {None, "different"}:
        with pytest.raises(ValueError):
            check_artifact_gate(artifact, {"sha256": "fixture"}, indexed)
    else:
        check_artifact_gate(artifact, {"sha256": "fixture"}, indexed)


def test_no_overwrite(tmp_path):
    with pytest.raises(FileExistsError):
        score(tmp_path, tmp_path, "unused", "unused")
