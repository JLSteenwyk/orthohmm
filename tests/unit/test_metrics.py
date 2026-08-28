import json

from orthohmm.metrics import PipelineMetrics, process_tree_rss_bytes


def test_process_tree_rss_includes_current_process():
    assert process_tree_rss_bytes(__import__("os").getpid()) > 0


def test_pipeline_metrics_writes_stage_counts_and_metadata(tmp_path):
    output = tmp_path / "metrics.json"
    with PipelineMetrics(str(output), sample_interval=0.001) as metrics:
        metrics.add_metadata(dataset="tiny")
        with metrics.stage("search"):
            payload = bytearray(1024 * 1024)
            assert payload
        metrics.add_counts(search_candidates=12, significant_hits=3)

    data = json.loads(output.read_text())
    assert data["status"] == "complete"
    assert data["metadata"]["dataset"] == "tiny"
    assert data["counts"] == {
        "search_candidates": 12,
        "significant_hits": 3,
    }
    assert data["stages"]["search"]["wall_s"] >= 0
    assert data["peak_process_tree_rss_bytes"] > 0
