import numpy as np
import pytest

from benchmark_tools import compare_search_hit_coverage as coverage


def test_directional_self_reciprocal_and_missing_query_counts():
    # a->a, a->b, b->a, b->c; d has no hit and c only receives one.
    report, codes = coverage.summarize(np.array([0, 0, 1, 1]), np.array([0, 1, 0, 2]),
                                       np.array([4., 3., 2., 1.]), np.array([0, 0, 1, 1]))
    assert codes.tolist() == [0, 1, 4, 6]
    assert report["directed_hits"] == 4
    assert report["self_hits"] == 1
    assert report["cross_species_hits"] == 1
    assert report["queries_without_hits"] == 2
    assert report["queries_without_nonself_hits"] == 2
    assert report["queries_without_cross_species_hits"] == 3
    assert report["targets_without_hits"] == 1
    assert report["reciprocal_nonself_directed_hits"] == 2
    assert report["reciprocal_nonself_unordered_pairs"] == 1
    assert report["normalized_score_quantiles"][2] == 2.5
    rows = report["species_directions"]
    assert [r["directed_hits"] for r in rows] == [3, 1, 0, 0]
    assert [r["queries_with_hits"] for r in rows] == [2, 1, 0, 0]
    assert all(r["query_species_genes"] == 2 for r in rows)


def test_overlap_keeps_direction_and_excludes_self_separately():
    result = coverage.overlap(np.array([0, 1, 4, 6]), np.array([0, 1, 9]), 4)
    assert result["all"]["intersection"] == 2
    assert result["all"]["union"] == 5
    assert result["all"]["jaccard"] == .4
    assert result["nonself"]["intersection"] == 1
    assert result["nonself"]["first_only"] == 2
    assert result["nonself"]["second_only"] == 1


def test_empty_hits_are_not_perfect_overlap_or_coverage():
    empty = np.array([], dtype=np.int32)
    report, codes = coverage.summarize(empty, empty, np.array([]), np.array([0, 1]))
    assert report["queries_without_hits"] == 2
    assert report["normalized_score_quantiles"] is None
    assert coverage.overlap(codes, codes, 2)["all"]["jaccard"] is None


@pytest.mark.parametrize("q,t,s", [([0, 0], [1, 1], [1., 2.]), ([2], [1], [1.]),
                                  ([-1], [1], [1.]), ([0.5], [1], [1.]),
                                  ([0], [1], [float("nan")]), ([0], [1], [0.]),
                                  ([0], [1, 0], [1.])])
def test_invalid_or_duplicate_hits_rejected(q, t, s):
    with pytest.raises(ValueError):
        coverage.summarize(np.array(q), np.array(t), np.array(s), np.array([0, 1]))


def test_species_codes_may_differ_but_membership_must_match():
    coverage.same_species_partition([0, 0, 1], [9, 9, 8])
    with pytest.raises(ValueError):
        coverage.same_species_partition([0, 0, 1], [9, 8, 8])
    with pytest.raises(ValueError):
        coverage.same_species_partition([0, 1], [9])


def test_live_conversion_cannot_load_hits(tmp_path, monkeypatch):
    monkeypatch.setattr(coverage.subprocess, "check_output", lambda *a, **k:
                        "JobIDRaw|State|ExitCode|Elapsed\n21292|RUNNING|0:0|00:01:00\n")
    monkeypatch.setattr(coverage, "load_replay_input", lambda *a, **k: pytest.fail("Loaded partial hits"))
    monkeypatch.setattr(coverage, "read_frozen", lambda *a, **k: pytest.fail("Passed terminal gate"))
    output = tmp_path / "report.json"
    with pytest.raises(ValueError):
        coverage.assemble(tmp_path, output)
    assert not output.exists()
