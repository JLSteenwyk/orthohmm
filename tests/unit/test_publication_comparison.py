import json
import hashlib

import pytest

from benchmark_tools.publication_comparison import METRICS, qfo_record, unique_index, validate_three_score
from benchmark_tools.publication_comparison import verify_final_group_workflow


def completed_workflow(directory):
    lines = []
    for name in ("pairs.tsv", "pairs.qfo.tsv"):
        path = directory / name
        path.write_text("a\tb\n")
        lines.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path}\n")
    (directory / "result.sha256").write_text("".join(lines))
    (directory / "run_metadata.tsv").write_text("exit_code\t0\nfinished\t2026-09-16T15:22:00Z\n")


def test_workflow_pending_without_manifest(tmp_path):
    assert verify_final_group_workflow(tmp_path) is None


def test_workflow_checks_actual_artifacts(tmp_path):
    completed_workflow(tmp_path)
    assert len(verify_final_group_workflow(tmp_path)["artifacts"]) == 2
    (tmp_path / "pairs.tsv").write_text("changed")
    with pytest.raises(ValueError, match="checksum mismatch"):
        verify_final_group_workflow(tmp_path)


@pytest.mark.parametrize("mutation", ["failure", "unfinished", "duplicate", "missing"])
def test_workflow_rejects_invalid_completion(tmp_path, mutation):
    completed_workflow(tmp_path)
    metadata = tmp_path / "run_metadata.tsv"
    manifest = tmp_path / "result.sha256"
    if mutation == "failure":
        metadata.write_text("exit_code\t1\nfinished\tnow\n")
    elif mutation == "unfinished":
        metadata.write_text("exit_code\t0\n")
    elif mutation == "duplicate":
        manifest.write_text(manifest.read_text() * 2)
    else:
        manifest.write_text(manifest.read_text().splitlines()[0] + "\n")
    with pytest.raises(ValueError):
        verify_final_group_workflow(tmp_path)


def test_duplicate_methods_rejected():
    with pytest.raises(ValueError, match="Duplicate"):
        unique_index([{"method": "a"}, {"method": "a"}], "method")


def test_missing_qfo_is_not_zero(tmp_path):
    result = qfo_record(tmp_path)
    assert result["status"] == "incomplete"
    assert len(result["missing_metrics"]) == 6
    assert "mean" not in result


def write_metrics(directory, participant):
    for relative, _ in METRICS.values():
        path = directory / "results" / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"datalink": {"inline_data": {
            "challenge_participants": [{"participant_id": participant, "metric_x": 0.8, "metric_y": 0.6}],
            "visualization": {"x_axis": "x", "y_axis": "y"}}}}))


def test_qfo_retains_axes_and_provenance(tmp_path):
    write_metrics(tmp_path, tmp_path.name)
    result = qfo_record(tmp_path)
    assert result["status"] == "metrics_available"
    assert result["scores"]["VGNC F"] == pytest.approx(2 * 0.8 * 0.6 / 1.4)
    assert result["metric_details"]["GO"]["participant"]["metric_x"] == 0.8
    assert len(result["metric_details"]["GO"]["source"]["sha256"]) == 64


def test_wrong_participant_rejected(tmp_path):
    write_metrics(tmp_path, "wrong")
    with pytest.raises(ValueError, match="Unexpected participant"):
        qfo_record(tmp_path)


def test_invalid_precision_cannot_be_hidden_by_harmonic_mean(tmp_path):
    write_metrics(tmp_path, tmp_path.name)
    path = tmp_path / "results/VGNC/VGNC.json"
    payload = json.loads(path.read_text())
    payload["datalink"]["inline_data"]["challenge_participants"][0]["metric_x"] = 1.1
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="Invalid QfO axis"):
        qfo_record(tmp_path)


def test_three_kingdoms_requires_exact_count_consistency():
    score = dict(true_positive_gene_pairs=2, false_positive_gene_pairs=1,
                 false_negative_gene_pairs=3, precision=2 / 3, recall=2 / 5, f_score=0.5)
    validate_three_score(score)
    score["f_score"] = 0.9
    with pytest.raises(ValueError, match="disagrees"):
        validate_three_score(score)
