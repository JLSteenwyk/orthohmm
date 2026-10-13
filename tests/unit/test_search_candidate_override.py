import json

import numpy as np
import pytest

from orthohmm.accuracy import resolve_accuracy_profile, resolve_search_candidate_cap
from orthohmm.args_processing import process_args
from orthohmm.helpers import StartStep, SubstitutionMatrix
from orthohmm.parser import create_parser


@pytest.mark.parametrize("profile,expected", [("standard", 0), ("high_sensitivity", 100)])
def test_implicit_cap_keeps_named_profile(profile, expected):
    config = resolve_accuracy_profile(profile)
    assert resolve_search_candidate_cap(config) == expected
    assert resolve_search_candidate_cap(config, 400) == 400
    assert config.max_candidates_per_query == expected


@pytest.mark.parametrize("cap", [0, 100, 400, np.int64(400), 2**31 - 1])
def test_explicit_cap_is_int32_safe(cap):
    assert resolve_search_candidate_cap(resolve_accuracy_profile("standard"), cap) == cap


@pytest.mark.parametrize("cap", [-1, 2**31, True, False, np.bool_(True), 1.5, 400.0, "400", float("nan")])
def test_invalid_cap_is_rejected(cap):
    with pytest.raises(ValueError, match="nonnegative int32"):
        resolve_search_candidate_cap(resolve_accuracy_profile("standard"), cap)


def test_cli_override_is_opt_in(tmp_path):
    default = process_args(create_parser().parse_args([str(tmp_path)]))
    assert "search_max_candidates_per_query" not in default
    configured = process_args(create_parser().parse_args([
        str(tmp_path), "--accuracy_profile", "high_sensitivity",
        "--search_max_candidates_per_query", "400",
    ]))
    assert configured["search_max_candidates_per_query"] == 400
    adaptive = process_args(create_parser().parse_args([
        str(tmp_path), "--search_max_candidates_per_query", "0",
    ]))
    assert adaptive["search_max_candidates_per_query"] == 0


@pytest.mark.parametrize("extra", [
    ["--search_max_candidates_per_query", "-1"],
    ["--search_max_candidates_per_query", str(2**31)],
    ["--search_max_candidates_per_query", "400", "--start", "search_res"],
    ["--search_max_candidates_per_query", "400", "--search_mode", "phmmer"],
])
def test_cli_rejects_invalid_or_reused_search(tmp_path, monkeypatch, extra):
    import orthohmm.args_processing as processing
    monkeypatch.setattr(processing, "find_executable", lambda name: name)
    with pytest.raises(SystemExit) as error:
        process_args(create_parser().parse_args([str(tmp_path), *extra]))
    assert error.value.code == 1


@pytest.mark.parametrize("profile,override,expected", [
    ("standard", None, 0), ("high_sensitivity", None, 100),
    ("high_sensitivity", 400, 400), ("standard", 0, 0),
])
def test_pipeline_handoff_and_metadata_match(tmp_path, monkeypatch, profile, override, expected):
    from orthohmm.orthohmm import execute
    import orthohmm.search.engine as engine

    class SearchReached(Exception):
        pass

    captured = {}

    def search(*args, **kwargs):
        captured.update(kwargs)
        raise SearchReached()

    monkeypatch.setattr(engine, "execute_builtin_search", search)
    (tmp_path / "species.faa").write_text(">gene\nACDEFGHIK\n")
    metrics = tmp_path / "metrics.json"
    kwargs = {} if override is None else {"search_max_candidates_per_query": override}
    with pytest.raises(SearchReached):
        execute(str(tmp_path), str(tmp_path / "out"), None, 2, 0.5, "mcl",
                1.5, None, None, SubstitutionMatrix.blosum62, 1e-4,
                accuracy_profile=profile, metrics_json=str(metrics), **kwargs)
    assert captured["max_candidates_per_query"] == expected
    payload = json.loads(metrics.read_text())
    assert payload["metadata"]["search_max_candidates_per_query"] == expected
    assert payload["metadata"]["accuracy_profile"] == profile


@pytest.mark.parametrize("search,start", [("phmmer", None), ("builtin", StartStep.search_res)])
def test_programmatic_pipeline_cannot_mislabel_cached_or_phmmer_search(tmp_path, search, start):
    from orthohmm.orthohmm import execute
    with pytest.raises(ValueError, match="fresh built-in search"):
        execute(str(tmp_path), str(tmp_path / "out"), None, 2, 0.5, "mcl",
                1.5, start, None, SubstitutionMatrix.blosum62, 1e-4,
                search_mode=search, search_max_candidates_per_query=400)
