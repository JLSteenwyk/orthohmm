from copy import deepcopy
import json

import pytest

from benchmark_tools import compare_orthohmm_generations as runner
from benchmark_tools import plot_orthohmm_generations as plot


def data():
    return dict(runs=[dict(r,status="validated",elapsed_s=100+r["index"],search_s=10+r["index"])
                      for r in runner.schedule()])


def test_plot_independently_recomputes_exact_runner_summaries():
    snapshot = data()
    original = deepcopy(snapshot)
    assert plot.summarize(snapshot) == runner.summarize(snapshot)
    assert snapshot == original


@pytest.mark.parametrize("kind", ["failed", "missing", "duplicate", "negative", "nonfinite"])
def test_plot_refuses_cherry_picked_or_invalid_runs(kind):
    snapshot = data()
    if kind == "failed": snapshot["runs"][0]["status"] = "timed_out"
    elif kind == "missing": snapshot["runs"].pop()
    elif kind == "duplicate": snapshot["runs"][0] = snapshot["runs"][1]
    elif kind == "negative": snapshot["runs"][0]["search_s"] = -1
    else: snapshot["runs"][0]["elapsed_s"] = float("nan")
    with pytest.raises(ValueError): plot.summarize(snapshot)


def test_plot_never_overwrites_existing_attempt(tmp_path):
    with pytest.raises(FileExistsError): plot.export(tmp_path / "missing.json", tmp_path)


def test_plot_rejects_fabricated_summary_before_any_figure_is_written(tmp_path):
    snapshot = data()
    snapshot["summary"] = []
    path = tmp_path / "snapshot.json"
    path.write_text(json.dumps(snapshot))
    output = tmp_path / "figures"
    with pytest.raises(ValueError): plot.export(path, output)
    assert not output.exists()
