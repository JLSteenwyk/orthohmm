import pytest

from benchmark_tools.summarize_matched_resources import verbose_time


def sample(elapsed="0:01.25", exit_code="0"):
    return f"""Elapsed (wall clock) time (h:mm:ss or m:ss): {elapsed}
User time (seconds): 2.1
System time (seconds): 0.3
Maximum resident set size (kbytes): 12345
Exit status: {exit_code}
"""


@pytest.mark.parametrize("elapsed,expected", [("0:01.25", 1.25), ("1:02:03", 3723), ("61:05", 3665)])
def test_valid_elapsed(elapsed, expected):
    result = verbose_time(sample(elapsed))
    assert result["elapsed_seconds"] == expected
    assert result["max_process_rss_kib"] == 12345
    assert "not simultaneous" in result["semantics"]["max_process_rss_kib"]


@pytest.mark.parametrize("elapsed", ["1", "0:60", "1:60:00", "0:nan", "-1:02", "1:2:3:4"])
def test_bad_elapsed(elapsed):
    with pytest.raises(ValueError):
        verbose_time(sample(elapsed))


def test_failed_missing_and_duplicate():
    for text in (sample(exit_code="1"), sample().replace("Exit status: 0\n", ""), sample() + "Exit status: 0\n"):
        with pytest.raises(ValueError):
            verbose_time(text)
