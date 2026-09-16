import json

import pytest

from benchmark_tools import assemble_sequence_search_control as assembler


def family(tp, fp, fn):
    return {"refog_records": [{"refog": "RefOG001", "genes": 3,
                              "true_positive": tp, "false_positive": fp, "false_negative": fn}]}


def test_fixed_three_methods_six_endpoints_and_direction():
    scores = {assembler.BASELINE: family(3, 0, 0), "all_hits": family(2, 0, 1), "top100": family(1, 0, 2)}
    result = assembler.bootstrap(scores)
    assert result["replicates"] == 20000
    assert result["seed"] == 20260918
    assert "6 reported" in result["multiplicity"]
    assert len(result["comparisons"]) == 2
    assert result["comparisons"]["all_hits"]["metrics"]["f_score"]["difference_percentage_points"] == -20
    assert result["comparisons"]["top100"]["family_f1_losses"] == 1


def test_partial_comparison_rejected():
    with pytest.raises(ValueError):
        assembler.bootstrap({assembler.BASELINE: family(3, 0, 0), "all_hits": family(2, 0, 1)})


@pytest.mark.parametrize("which", ["native", "coverage"])
def test_failed_gate_prevents_reference_access(tmp_path, monkeypatch, which):
    def fail(*args):
        raise ValueError("Incomplete gate")

    monkeypatch.setattr(assembler, "validate", fail if which == "native" else lambda *a: {})
    monkeypatch.setattr(assembler, "admit_coverage", fail)
    monkeypatch.setattr(assembler, "read_frozen", lambda *a: pytest.fail("Loaded benchmark before gates"))
    output = tmp_path / "out"
    with pytest.raises(ValueError, match="Incomplete"):
        assembler.assemble(tmp_path, output)
    assert not output.exists()


def test_live_coverage_cannot_be_admitted(tmp_path, monkeypatch):
    monkeypatch.setattr(assembler.subprocess, "check_output", lambda *a, **k:
                        "JobIDRaw|State|ExitCode|Elapsed\n21294|RUNNING|0:0|00:01:00\n")
    with pytest.raises(ValueError):
        assembler.admit_coverage(tmp_path, {})


def test_no_overwrite(tmp_path, monkeypatch):
    monkeypatch.setattr(assembler, "validate", lambda *a: pytest.fail("Validated over existing output"))
    with pytest.raises(FileExistsError):
        assembler.assemble(tmp_path, tmp_path)


@pytest.mark.parametrize("bad", [False, True])
def test_resources_preserve_scope_and_unmeasured_baseline(tmp_path, bad):
    path = tmp_path / "metrics.json"
    path.write_text(json.dumps({"wall_s": -1 if bad else 10, "peak_process_rss_gib": 2,
                                "timings": {"load_and_index_hits_s": 1}}))
    native = {"variants": {label: {"native_records": {"replay.json": assembler.file_provenance(path)}}
                           for label in assembler.VARIANTS}}
    if bad:
        with pytest.raises(ValueError):
            assembler.graph_resources(native)
    else:
        result = assembler.graph_resources(native)
        assert result[assembler.BASELINE]["wall_s"] is None
        assert result["all_hits"]["wall_s"] == 10
        assert "not process tree" in result["all_hits"]["scope"]
