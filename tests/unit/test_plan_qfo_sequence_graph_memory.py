from copy import deepcopy
from pathlib import Path

import pytest

from benchmark_tools import plan_qfo_sequence_graph_memory as module


def setup(monkeypatch):
    conversion = {"variants": {}}
    calls = []
    for label, hits in (("all_hits", 20), ("top100", 10)):
        conversion["variants"][label] = dict(checkpoint="/" + label,
            manifest=dict(path="/" + label + "/manifest.json", sha256=label),
            audit=dict(summary=dict(genes=5, hits=hits, self_hits=5)))
    monkeypatch.setattr(module, "check", lambda item: None)

    def estimate(checkpoint, sha, core):
        calls.append((checkpoint, sha, core))
        summary = conversion["variants"][checkpoint.name]["audit"]["summary"]
        return dict(estimate=dict(summary, graph_feasibility_admitted=False))

    monkeypatch.setattr(module, "estimate", estimate)
    return conversion, calls


def test_both_variants_exact_hashes(monkeypatch):
    conversion, calls = setup(monkeypatch)
    result = module.estimate_variants(conversion, Path("/core"))
    assert list(result) == ["all_hits", "top100"]
    assert calls == [(Path("/all_hits"), "all_hits", Path("/core")),
                     (Path("/top100"), "top100", Path("/core"))]


def test_missing_variant(monkeypatch):
    conversion, _ = setup(monkeypatch)
    del conversion["variants"]["all_hits"]
    with pytest.raises(ValueError, match="both"):
        module.estimate_variants(conversion, Path("/core"))


def test_path_disagreement(monkeypatch):
    conversion, _ = setup(monkeypatch)
    conversion["variants"]["all_hits"]["manifest"]["path"] = "/wrong/manifest.json"
    with pytest.raises(ValueError, match="location"):
        module.estimate_variants(conversion, Path("/core"))


@pytest.mark.parametrize("field,value", [("hits", 1), ("genes", 1), ("self_hits", 1),
                                         ("graph_feasibility_admitted", True)])
def test_bad_estimator_result(monkeypatch, field, value):
    conversion, _ = setup(monkeypatch)
    estimate = module.estimate

    def corrupt(*args):
        result = deepcopy(estimate(*args))
        result["estimate"][field] = value
        return result

    monkeypatch.setattr(module, "estimate", corrupt)
    with pytest.raises(ValueError):
        module.estimate_variants(conversion, Path("/core"))


def test_failed_scheduler_prevents_output(tmp_path, monkeypatch):
    def fail(*args):
        raise ValueError("not completed")
    monkeypatch.setattr(module, "completed", fail)
    output = tmp_path / "result.json"
    with pytest.raises(ValueError, match="not completed"):
        module.run(tmp_path, "21792", output)
    assert not output.exists()
