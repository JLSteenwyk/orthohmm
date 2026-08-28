import json

import pytest

from benchmark_tools.qfo_summarize_scores import harmonic_mean, summarize


def _write_metric(path, x, y):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "datalink": {
                    "inline_data": {
                        "challenge_participants": [
                            {"participant_id": "test", "metric_x": x, "metric_y": y}
                        ]
                    }
                }
            }
        )
    )


def test_harmonic_mean_handles_zero_pair():
    assert harmonic_mean(0.0, 0.0) == 0.0
    assert harmonic_mean(0.8, 0.5) == pytest.approx(0.6153846154)


def test_summarize_uses_f_scores_and_scalar_metrics(tmp_path):
    scoring_dir = tmp_path / "method"
    _write_metric(scoring_dir / "results/VGNC/VGNC.json", 0.8, 0.5)
    _write_metric(scoring_dir / "results/SwissTrees/SwissTrees.json", 0.9, 0.6)
    _write_metric(scoring_dir / "results/TreeFam-A/TreeFam-A.json", 0.7, 0.7)
    _write_metric(scoring_dir / "results/EC/EC.json", 100, 0.4)
    _write_metric(scoring_dir / "results/GO/GO.json", 100, 0.5)
    _write_metric(scoring_dir / "results/FAS/FAS.json", 100, 0.6)

    summary = summarize(scoring_dir)

    assert summary["prediction"] == "method"
    assert summary["scores"]["VGNC F"] == pytest.approx(0.6153846154)
    assert summary["scores"]["SwissTrees F"] == pytest.approx(0.72)
    assert summary["scores"]["TreeFam-A F"] == pytest.approx(0.7)
    assert summary["scores"]["EC"] == pytest.approx(0.4)
    assert summary["scores"]["GO"] == pytest.approx(0.5)
    assert summary["scores"]["FAS"] == pytest.approx(0.6)
