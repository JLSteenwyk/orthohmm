import pytest

from benchmark_tools.summarize_dgx_observation_gaps import resource_errors


def rows():
    return [{"snapshot": {"sampling_errors": errors, "sampled_sum_process_rss_bytes": rss}}
            for errors, rss in (([], 10), ([{"error_type": "NoSuchProcess"}, {"error_type": "AccessDenied"}], 30),
                                ([{"error_type": "NoSuchProcess"}], 20))]


def test_events_are_distinct_from_affected_samples():
    result = resource_errors(iter(rows()), {"observations": 3, "samples_with_process_errors": 2,
                                          "maximum_sampled_sum_rss_bytes": 30})
    assert result["error_samples"] == 2
    assert result["error_events"] == 3
    assert result["error_types"] == {"NoSuchProcess": 2, "AccessDenied": 1}
    assert result["error_samples_at_recorded_rss_maximum"] == 1
    assert result["fraction_of_samples_with_errors"] == 2 / 3


@pytest.mark.parametrize("observations,affected", [(2, 2), (3, 1)])
def test_replayed_counts_must_match(observations, affected):
    with pytest.raises(ValueError):
        resource_errors(iter(rows()), {"observations": observations, "samples_with_process_errors": affected,
                                      "maximum_sampled_sum_rss_bytes": 30})


def test_no_error_does_not_create_errors():
    result = resource_errors([rows()[0]], {"observations": 1, "samples_with_process_errors": 0,
                                          "maximum_sampled_sum_rss_bytes": 10})
    assert result["error_events"] == 0
    assert result["error_samples_at_recorded_rss_maximum"] == 0
